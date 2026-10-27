"""
Checks for the optional 'eps_zz' return form of densification models (physics.get_drho_dt).
Run from CFM_main:  python ../release_checks/test_eps_zz_interface.py
Convention: eps_zz is the vertical COMPACTION rate (positive = compaction, like drho_dt); eps_xx/eps_yy are
extension-positive horizontal strain rates (strain.py). drho_dt = rho*(eps_zz - (eps_xx+eps_yy)).
"""
import sys, json, shutil
sys.path.insert(0, '.')
import numpy as np
from physics import get_drho_dt, FirnPhysics
fails = []
def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name, info)
    if not ok: fails.append(name)

rng = np.random.default_rng(0)
rho = rng.uniform(300, 900, 500)
drho = rng.uniform(0, 1e-7, 500)
exx, eyy = rng.uniform(-1e-10, 1e-10, 500), rng.uniform(-1e-10, 1e-10, 500)

# 1. legacy form is passed through untouched (same object, horizontal strain NOT applied)
RD = {'drho_dt': drho, 'viscosity': 1.}
out = get_drho_dt(RD, rho, exx, eyy)
check('legacy drho_dt returned unchanged', out is drho)
# 2. eps_zz form with zero horizontal strain reproduces a legacy model
out = get_drho_dt({'eps_zz': drho / rho}, rho)
check('eps_zz (no horizontal strain) == legacy drho_dt', np.allclose(out, drho, rtol=1e-13, atol=0), np.max(np.abs(out-drho)/drho))
# 3. sign: positive eps_zz densifies
check('positive eps_zz -> positive drho_dt', np.all(get_drho_dt({'eps_zz': np.full(500, 1e-9)}, rho) > 0))
# 4. incompressible flow: vertical compaction exactly balances horizontal extension -> no density change
out = get_drho_dt({'eps_zz': exx + eyy}, rho, exx, eyy)
check('incompressible flow: drho_dt == 0', np.allclose(out, 0, atol=1e-25), np.max(np.abs(out)))
# 5. horizontal extension (eps_xx+eps_yy>0) at fixed eps_zz lowers the density rate; convergence raises it
a = get_drho_dt({'eps_zz': np.full(500, 1e-9)}, rho, 1e-10, 0.); b = get_drho_dt({'eps_zz': np.full(500, 1e-9)}, rho, -1e-10, 0.)
check('extension reduces / convergence raises drho_dt', np.all(a < get_drho_dt({'eps_zz': np.full(500, 1e-9)}, rho)) and np.all(b > a))
# 6. error handling
for name, RDx, kw in (('both keys', {'drho_dt': drho, 'eps_zz': drho}, {}), ('neither key', {'viscosity': 1.}, {}),
                      ('eps_zz + strain_softening', {'eps_zz': drho / rho}, {'strain_softening': True})):
    try: get_drho_dt(RDx, rho, **kw); check(f'raises: {name}', False)
    except (ValueError, KeyError): check(f'raises: {name}', True)
check('legacy + strain_softening still allowed', get_drho_dt({'drho_dt': drho}, rho, strain_softening=True) is drho)

# 7. end-to-end: wrap legacy models so they return eps_zz INCLUDING the horizontal term the central step removes
#    (eps_zz = drho_dt/rho + eps_xx + eps_yy, with eps_xx/eps_yy read from the FirnPhysics instance, i.e. exercising the
#    PhysParams wiring in both spin-up and transient loops). Results must match the plain legacy run.
from firn_density_nospin import FirnDensityNoSpin
def run(tag, physRho, wrap, divergence):
    c = json.load(open('example.json'))
    c.update(resultsFolder=f'/tmp/cfm_eps_{tag}', physRho=physRho, horizontal_divergence=divergence, strain_softening=False)
    shutil.rmtree(c['resultsFolder'], ignore_errors=True)
    orig = getattr(FirnPhysics, physRho)
    if wrap:
        def f(self, orig=orig):
            RD = orig(self); d = RD.pop('drho_dt')
            RD['eps_zz'] = d / self.rho + (self.eps_xx + self.eps_yy)
            return RD
        setattr(FirnPhysics, physRho, f)
    try:
        json.dump(c, open(f'/tmp/cfm_eps_{tag}.json', 'w'))
        m = FirnDensityNoSpin(f'/tmp/cfm_eps_{tag}.json', NewSpin=True); m.time_evolve()
    finally:
        setattr(FirnPhysics, physRho, orig)
    import h5py
    f = h5py.File(c['resultsFolder'] + '/CFMresults.hdf5', 'r')
    return f['density'][:], f['depth'][:]
import contextlib, io
for physRho in ('KuipersMunneke2015', 'HLdynamic'):
    for div in (False, True):
        with contextlib.redirect_stdout(io.StringIO()):
            base = run(f'{physRho}_{div}_base', physRho, False, div); alt = run(f'{physRho}_{div}_eps', physRho, True, div)
        ok = np.allclose(base[0], alt[0], rtol=5e-7, atol=0, equal_nan=True) and np.allclose(base[1], alt[1], rtol=5e-7, atol=0, equal_nan=True)  # outputs are float32: drho/rho*rho round-off can flip the last bit
        d = np.nanmax(np.abs(base[0] - alt[0])) if base[0].shape == alt[0].shape else 'shape'
        check(f'end-to-end {physRho} horizontal_divergence={div}', ok, f'max|d density|={d}')
print('\nRESULT:', 'FAIL' if fails else 'PASS', fails)
sys.exit(1 if fails else 0)
