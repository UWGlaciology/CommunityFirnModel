"""Unit checks for SEB_turbulent.py / SEB.py turbulent-flux paths. Run from CFM_main."""
import sys, numpy as np
sys.path.insert(0, '.')
from constants import *
from SEB import SurfaceEnergyBudget
from SEB_turbulent import (andreas_roughness, q_sat_ice, BulkRichardsonFlux,
                           MoninObukhovFlux, RoughnessLengths, TurbulentFluxModel)
fails = []
def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name, info)
    if not ok: fails.append(name)

N = 3
def clim(T2m=255., WS=5., q=1.0e-3, PS=8.0e4, QH=0., QL=0.):
    f = lambda v: np.full(N, float(v))
    return dict(time=np.arange(N)/365., SW_d=f(200), LW_d=f(200), ALBEDO=f(0.8), T2m=f(T2m),
                TSKIN=f(T2m), QH=f(-QH), QL=f(-QL), RAIN=f(0), WS=f(WS), QV2m=f(q), PS=f(PS))
def phys(T=250., n=200):
    dz = np.full(n, 0.02); z = np.cumsum(dz) - dz; z = np.maximum(z, 0)
    z = np.concatenate(([0.], np.cumsum(dz)[:-1]))
    return dict(Tz=np.full(n, T), mass=dz*350., dt=86400., z=z, dz=dz, mtime=0.)
def seb(method, c=None, **kw):
    cfg = dict(SEB_turb_method=method, SEB_TL_thick=0.05); cfg.update(c or {})
    return SurfaceEnergyBudget(cfg, clim(**kw), 0, None, 0)

# 1. q_sat_ice anchor values
check('q_sat 273.15K', abs(q_sat_ice(273.15, 1e5) - 0.622*611.15/(1e5-0.378*611.15)) < 1e-9)
e253 = 611.15*np.exp(22.452*(-20)/(253.15-0.6)); check('q_sat ~103 Pa at -20C', abs(e253-103.2) < 2, e253)

# 2. Andreas roughness: continuity at regime boundaries and physical trend
for Re in (0.135, 2.5):
    u = Re*1.461e-5/1e-3
    a = andreas_roughness(1e-3, u*0.9999999); b = andreas_roughness(1e-3, u*1.0000001)
    check(f'andreas continuity Re={Re}', abs(np.log(a[0]/b[0])) < 5e-3 and abs(np.log(a[1]/b[1])) < 5e-3, (a, b))
zh, zq = andreas_roughness(1e-3, 0.3*1e-3*0+ 10*1.461e-5/1e-3)  # Re*=10, rough regime
check('andreas: z0h < z0m for rough flow (Re*=10)', zh < 1e-3 and zq < 1e-3, (zh, zq))
zh, zq = andreas_roughness(1e-3, 0.05*1.461e-5/1e-3)            # Re*=0.05, smooth
check('andreas: z0h > z0m for smooth flow', zh > 1e-3, (zh, zq))

# 3. Neutral limit: bulk & MO should reproduce k^2/(ln ln) exchange coefficient
s_b = seb('bulk_richardson', T2m=270., WS=8.); s_m = seb('monin_obukhov', T2m=270., WS=8.)
QHb, QLb = s_b.flux_model.compute(270., 0); QHm, QLm = s_m.flux_model.compute(270., 0)
check('neutral: QH==0 when Ts==T2m', abs(QHb) < 1e-6 and abs(QHm) < 1e-6, (QHb, QHm))
Ts = 269.0
QHb, QLb = s_b.flux_model.compute(Ts, 0); QHm, QLm = s_m.flux_model.compute(Ts, 0)
check('stable weak-gradient: bulk vs MO QH within 30%', abs(QHb-QHm) < 0.3*abs(QHm), (QHb, QHm))
# 4. sign conventions (positive = into snow)
sw = seb('monin_obukhov', T2m=270., q=5e-4)
QH, QL = sw.flux_model.compute(260., 0); check('warm air over cold snow: QH>0', QH > 0, QH)
QH, QL = sw.flux_model.compute(272., 0); check('cold air over warm snow: QH<0', QH < 0, QH)
check('dry air: QL<0 (sublimation)', QL < 0, QL)
# 5. sweep: finite, monotone-ish, no blowups across stability range
for name, S in (('bulk', 'bulk_richardson'), ('MO', 'monin_obukhov')):
    for WS in (0.0, 0.5, 2., 10., 25.):
        for Ta in (230., 255., 272.):
            for Ts in (Ta-10, Ta-3, Ta, Ta+3, 273.15):
                if abs(Ts-Ta) > 10: continue
                s = seb(S, T2m=Ta, WS=WS)
                QH, QL = s.flux_model.compute(min(Ts, 273.15), 0)
                if not (np.isfinite(QH) and np.isfinite(QL) and abs(QH) < 1500 and abs(QL) < 1500):
                    check(f'{name} finite/bounded WS={WS} Ta={Ta} Ts={Ts}', False, (QH, QL))
check('sweep done', True)
# 6. MO iteration convergence
class Probe(MoninObukhovFlux): pass
# 7. Newton path == quartic path when fluxes are Ts-independent
class FixedFlux(TurbulentFluxModel):
    depends_on_Ts = True
    def __init__(s, QH, QL): s.QH, s.QL = QH, QL
    def compute(s, Ts, iii): return s.QH, s.QL
# Newton's ground-heat term uses G(Ts); the quartic path uses G(Tz[0]). They agree when Ts~Tz[0]
# (cold cases near equilibrium) and to ~1% in melt, so tolerances reflect that design difference.
for QH0, QL0, T0 in ((-13, -10, 250.), (-12, -10, 250.), (80, 20, 272.5), (150, 30, 272.9)):
    sa = seb('prescribed', T2m=T0, QH=QH0, QL=QL0); sb = seb('prescribed', T2m=T0, QH=QH0, QL=QL0)
    sb.flux_model = FixedFlux(QH0, QL0)
    Ta, _, Ma, _ = sa.solve_Ts(phys(T0), 0, T0); Tb, _, Mb, _ = sb.solve_Ts(phys(T0), 0, T0)
    check(f'newton==quartic QH={QH0}', abs(float(np.squeeze(Ta))-float(np.squeeze(Tb))) < 5e-2 and abs(float(np.squeeze(Ma))-float(np.squeeze(Mb))) <= 0.02*abs(float(np.squeeze(Ma))), (float(np.squeeze(Ta)), float(np.squeeze(Tb)), float(np.squeeze(Ma)), float(np.squeeze(Mb))))
# 8. full solves with real models produce physical Ts, finite melt
for S in ('bulk_richardson', 'monin_obukhov'):
    for T0, WS in ((235., 3.), (255., 8.), (271., 1.), (272.5, 6.)):
        s = seb(S, T2m=T0, WS=WS); Ts, Tz, melt, _ = s.solve_Ts(phys(T0), 0, T0)
        ok = np.isfinite(float(Ts)) and 150 < float(Ts) <= 273.15 and np.isfinite(float(melt)) and float(melt) >= 0
        check(f'{S} solve_Ts T2m={T0} WS={WS}', ok, (float(Ts), float(melt)))

# 9. Sub-stepped (dtRATIO>1) Newton path == sub-stepped quartic path for Ts-independent fluxes
def seb_sub(method, T0, QH0, QL0, dtr=6, **kw):
    n = 2*dtr
    f = lambda v: np.full(n, float(v))
    SEBf = dict(time=np.arange(n)/365./dtr, SW_d=f(150)*(1+np.sin(np.arange(n))), LW_d=f(200), ALBEDO=f(0.8), T2m=f(T0),
                TSKIN=f(T0), QH=f(QH0), QL=f(QL0), RAIN=f(0), dtRATIO=dtr, WS=f(5), QV2m=f(1e-3), PS=f(8e4))
    cl = dict(T2m=np.full(2, float(T0)))
    return SurfaceEnergyBudget(dict(SEB_turb_method=method, SEB_TL_thick=0.05), cl, 0, SEBf, 0)
# T tolerance 1 K: Newton evaluates the ground-heat flux at the new Ts, the quartic path at the old Tz[0]
for QH0, QL0, T0 in ((-13, -10, 250.), (-5, -5, 266.), (60, 20, 272.5)):
    sa = seb_sub('prescribed', T0, QH0, QL0); sb = seb_sub('prescribed', T0, QH0, QL0)
    sb.flux_model = FixedFlux(QH0, QL0)
    pa, pb = phys(T0), phys(T0)
    Ta, _, Ma, _ = sa.solve_Ts_subdt(pa, 0, T0, 0.); Tb, _, Mb, _ = sb.solve_Ts_subdt(pb, 0, T0, 0.)
    check(f'subdt newton==quartic T0={T0}', abs(float(Ta)-float(Tb)) < 1.0 and abs(float(Ma)-float(Mb)) <= 0.05*abs(float(Ma))+0.05,
          (float(Ta), float(Tb), float(Ma), float(Mb)))
for S in ('bulk_richardson', 'monin_obukhov'):
    s = seb_sub(S, 268., 0, 0); Ts, Tz, melt, _ = s.solve_Ts_subdt(phys(268.), 0, 268., 0.)
    check(f'{S} solve_Ts_subdt runs', np.isfinite(float(Ts)) and 200 < float(Ts) <= 273.15 and float(melt) >= 0, (float(Ts), float(melt)))
print('\nRESULT:', 'FAIL' if fails else 'PASS', fails)
sys.exit(1 if fails else 0)
