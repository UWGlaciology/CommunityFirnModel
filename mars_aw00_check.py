'''
Mars check of the vapor transport against Arthern, Winebrenner, and Waddington (2000, AW00; Icarus 144, 367-386), Figs. 9 and 10.

AW00 Sec. 4.2: a thin layer of density rho_s over a half space of ice, with the surface temperature
Ts = Tbar + dTy cos(wy t) + dTd cos(wd t) (Eq. 20), has at z = 0 the temperature gradient of Eq. 21
dT/dz = (k_i/k_s) [dTy/zeta_y (sin wy t - cos wy t) + dTd/zeta_d (sin wd t - cos wd t)]    (z positive down),
and the annual mean of the vapor flux q = -D phi d(rho_v)/dz (Eqs. 5 and 6) there gives the critical accumulation rate
b_crit = qbar(0) / (1 - rho_s/rho_i)    (Eq. 22; the equations were checked on the rendered page): the accumulation
that is needed to keep a porous snowpack (Fig. 9: Tbar = 175 K, no diurnal cycle: about 0.5 kg m-2 Mars year-1 for rho_s = 360
and dTy = 30 K; Fig. 10: dTd = 10 K: about 5).

This script evaluates that with the functions of CFM_main/vapor.py (CO2 saturation vapor density, Eq. 6; diffusivity of Clow 1987 with
the Yosida factor of 4) and prints the critical accumulation rate for the porosity of AW00 (phi) and for the pore space of the CFM
(phi times the Johnsen 1/tortuosity). The atmospheric pressure of AW00 is not given in the paper: 600 Pa is assumed (the mean
pressure at the pole of Mars); the flux is proportional to 1/P (AW00 Sec. 5).

    python mars_aw00_check.py            prints the table
Not part of the CFM; tests/test_vapor.py has a test with the numbers.
'''
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CFM_main'))
from vapor import firn_tortuosity_inverse, saturation_vapor_density, vapor_diffusivity

SOL = 88775.2                 # s
MARS_YEAR = 668.6 * SOL       # s
RHO_I = 917.0


def critical_accumulation(Tbar, dTy, dTd, rho_s, P=600.0, jonsen_pore=False, yosida=None, n_per_sol=64, uniform=False):
    '''b_crit [kg m-2 Mars year-1] of AW00 Eq. 22 (flux from Eqs. 5, 6, 21). uniform: a half space of snow of density rho_s instead of a thin layer
    over ice (the skin depth of the snow, no conductivity ratio), the case that the model column below reproduces without the first-cell problem.'''
    t = np.arange(int(round(668.6 * n_per_sol))) * (SOL / n_per_sol)
    wy, wd = 2 * np.pi / MARS_YEAR, 2 * np.pi / SOL
    ki = 9.828 * np.exp(-0.0057 * Tbar)                     # ice at the mean temperature (Paterson 1994)
    c = 152.5 + 7.122 * Tbar
    ks = ki * (rho_s / RHO_I) ** 2                          # conductivity of the surface layer, AW00 Eq. B3 without the vapor term
    zeta = (lambda w: np.sqrt(2 * ks / (rho_s * c * w))) if uniform else (lambda w: np.sqrt(2 * ki / (RHO_I * c * w)))   # AW00 Eq. 16
    dTdz = (1.0 if uniform else ki / ks) * (dTy / zeta(wy) * (np.sin(wy * t) - np.cos(wy * t)) + dTd / zeta(wd) * (np.sin(wd * t) - np.cos(wd * t)))
    Ts = Tbar + dTy * np.cos(wy * t) + dTd * np.cos(wd * t)
    Ts = np.maximum(Ts, 148.0)                              # the clip of AW00 Eq. B4 (the CO2 frost point)
    dT = 0.01
    drho_dT = (saturation_vapor_density(Ts + dT, 'CO2') - saturation_vapor_density(Ts - dT, 'CO2')) / (2 * dT)
    D = vapor_diffusivity(Ts, P, 'CO2', yosida)
    phi = 1 - rho_s / RHO_I
    if jonsen_pore:
        phi = phi * float(firn_tortuosity_inverse(rho_s))
    q_mean = np.mean(-D * phi * drho_dT * dTdz)             # downward positive [kg m-2 s-1]
    return q_mean * MARS_YEAR / (1 - rho_s / RHO_I)


if __name__ == '__main__':
    print('b_crit [kg m-2 Mars year-1], Tbar = 175 K, rho_s = 360, dTy = 30 K, P = 600 Pa')
    for dTd, fig in ((0.0, 'Fig. 9, about 0.5'), (10.0, 'Fig. 10, about 5')):
        a = critical_accumulation(175.0, 30.0, dTd, 360.0)
        b = critical_accumulation(175.0, 30.0, dTd, 360.0, jonsen_pore=True)
        print(f'  dTd = {dTd:4.1f} K ({fig}):  porosity of AW00 {a:7.3f},  with Johnsen 1/tau {b:7.3f}')


def column_flux(Tbar, dTy, dTd, rho_s, P=600.0, H=0.01, years=2, n_per_sol=16, uniform=False):
    '''
    The same mean flux q(0) [kg m-2 Mars year-1], from the model: a layer of density rho_s and thickness H over pore-free ice
    (layers growing by 10% with depth to 60 m), the heat solver of the model (transient_solve_TR, as heatDiff) and the top
    flux of vapor_mass_flux_divergence with vapor_surface_exchange (Johnsen pore space, CO2, Yosida factor 4): the flux just beneath
    the surface, from the gradient between the top two layers (AW00's q(0)). The densities do not change (the mass flux is only
    recorded). Mean over the last year, divided by (1 - rho_s/rho_i) as Eq. 22.
    '''
    import types
    from solver import transient_solve_TR
    from vapor import vapor_mass_flux_divergence
    dz = [H]
    while sum(dz) < 60.0:
        dz.append(dz[-1] * 1.1 if len(dz) > 1 else 0.01)
    dz = np.array(dz)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    rho = np.full(len(dz), RHO_I)
    rho[0] = rho_s
    if uniform:
        rho[:] = rho_s
    zd = np.concatenate((z, [z[-1] + dz[-1]]))             # layer edges with the dummy point, as heatDiff
    zP = (zd[1:] + zd[:-1]) / 2
    ki = 9.828 * np.exp(-0.0057 * Tbar)
    c = 152.5 + 7.122 * Tbar
    K = ki * (rho / RHO_I) ** 2
    cvol = rho * c
    dt = SOL / n_per_sol
    wy, wd = 2 * np.pi / MARS_YEAR, 2 * np.pi / SOL
    config = {'vapor_atmosphere': 'CO2', 'site_pressure': P / 100.0, 'vapor_surface_exchange': True, 'vapor_min_thickness': 0.0}
    T = np.full(len(dz), float(Tbar))
    nyear = int(round(668.6 * n_per_sol))
    q = []
    for i in range(years * nyear):
        t = i * dt
        Ts = max(Tbar + dTy * np.cos(wy * t) + dTd * np.cos(wd * t), 148.0)
        T = transient_solve_TR(zd, zP, dt, K, T, len(dz), Ts, cvol)
        if i >= (years - 1) * nyear:
            col = types.SimpleNamespace(z=z, dz=dz, Tz=T, rho=rho, mass=rho * dz, c=config)
            q.append(vapor_mass_flux_divergence(col, dt, max_fill=np.inf, max_loss=np.inf)[0].sum())
    return np.mean(q) * MARS_YEAR / (1 - rho_s / RHO_I)


if __name__ == '__main__' and 'column' in sys.argv:
    for dTd in (0.0, 10.0):
        n = 16 if dTd else 4
        print(f'  dTd = {dTd:4.1f} K: thin layer over ice, from the model column {column_flux(175.0, 30.0, dTd, 360.0, n_per_sol=n):7.3f}; '
              f'snow half space: model column {column_flux(175.0, 30.0, dTd, 360.0, n_per_sol=n, uniform=True):7.3f}, '
              f'Eq. 21 with the snow {critical_accumulation(175.0, 30.0, dTd, 360.0, jonsen_pore=True, uniform=True):7.3f}')
