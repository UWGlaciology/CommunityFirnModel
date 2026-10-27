"""
Checks for the horizontal-strain input sources (strain.get_strain_input / load_strain).
Run from CFM_main:  python ../release_checks/test_strain_inputs.py
"""
import sys, json, shutil, types, subprocess, importlib.util, contextlib, io
sys.path.insert(0, '.')
import numpy as np
import strain as S
fails = []
def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name, info)
    if not ok: fails.append(name)

def mock(spin=False, **c):
    m = types.SimpleNamespace()
    m.c = dict(int_type='linear', InputFileFolder='CFMinput_example', InputFileNameStrain='example_STRAIN.csv',
               spinup_climate_type='mean', strain_softening=True, horizontal_divergence=True)
    m.c.update(c); m.modeltime = np.linspace(1980, 1990, 121); m.stp = 50; m.gridLen = 300
    return m
def load(m, spin=False, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return S.load_strain(m, spin=spin, **kw)

# 1. csv path reproduces the original (git HEAD) implementation
src = subprocess.check_output(['git', 'show', 'HEAD:CFM_main/strain.py']).decode()
spec = importlib.util.spec_from_loader('strain_head', loader=None); old = importlib.util.module_from_spec(spec)
exec(compile(src, 'strain_head', 'exec'), old.__dict__)
for spin in (False, True):
    for ncols in ('example_STRAIN.csv',):
        a, b = mock(), mock()
        with contextlib.redirect_stdout(io.StringIO()):
            ea, da, _ = old.load_strain(a, spin=spin); eb, db, _ = S.load_strain(b, spin=spin)
        n = 121 if not spin else 50
        check(f'csv == original (spin={spin})', np.allclose(da[:n], db[:n], rtol=0, atol=0) and np.allclose(ea[:n], eb[:n], rtol=0, atol=0), (da[:2], db[:2]))
        shape_ok = len(db) == n
        check(f'array length = number of steps (spin={spin})', shape_ok, len(db))

# 1b. 2- and 3-row csv files (first row time, then components) also reproduce the original
tt = np.linspace(1975, 1995, 241)
for nrow, name in ((2, 'strain2.csv'), (3, 'strain3.csv')):
    rows = [tt] + [1e-4 * (i + 1) * (1 + 0.3 * np.sin(tt * (i + 1))) for i in range(nrow)]
    np.savetxt('/tmp/' + name, np.vstack(rows), delimiter=',')
    for spin in (False, True):
        for sc in ('mean', 'initial'):
            a_, b_ = mock(InputFileFolder='/tmp', InputFileNameStrain=name, spinup_climate_type=sc), mock(InputFileFolder='/tmp', InputFileNameStrain=name, spinup_climate_type=sc)
            with contextlib.redirect_stdout(io.StringIO()):
                ea, da, _ = old.load_strain(a_, spin=spin); eb, db, _ = S.load_strain(b_, spin=spin)
            n = 50 if spin else 121
            check(f'{nrow}-row csv == original (spin={spin}, {sc})', np.array_equal(da[:n], db[:n]) and np.array_equal(ea[:n], eb[:n]))

# 2. sources
m = mock(); load(m, strain_rates=1e-4)
check('scalar -> divergence', np.allclose(m.eps_divergence, 1e-4) and len(m.eps_divergence) == 121 and m.c['strain_softening'] is False)
m = mock(); load(m, strain_rates=0.0)
check('scalar 0.0 is used (not mistaken for "no input")', np.all(m.eps_divergence == 0.0) and len(m.eps_divergence) == 121)
m = mock(); load(m, strain_rates=[1e-4, 2e-4, 5e-5])
check('3-list -> xx,yy,xy', np.allclose(m.eps_xx, 1e-4) and np.allclose(m.eps_yy, 2e-4) and np.allclose(m.eps_xy, 5e-5) and np.allclose(m.eps_divergence, -3e-4))
m = mock(); load(m, strain_rates=[1e-4, 2e-4])
check('2-list -> principal', np.allclose(m.eps_divergence, -3e-4) and np.allclose(m.eps_eff_hor_2, (1e-8 + 4e-8)/2))
t = np.linspace(1970, 2000, 361)
ts = {'time': t, 'eps_xx': 1e-4 * np.ones_like(t), 'eps_yy': 3e-4 * np.sin(t), 'eps_xy': np.zeros_like(t)}
m = mock(); load(m, climateTS=ts)
check('climateTS named keys, time-interpolated', np.allclose(m.eps_yy, np.interp(m.modeltime, t, ts['eps_yy']), atol=2e-5) and np.allclose(m.eps_xx, 1e-4))
m = mock(); load(m, strain_rates={'time': t, 'eps': np.vstack([ts['eps_xx'], ts['eps_yy'], ts['eps_xy']])})
check('dict with eps (3xN) == climateTS keys', np.allclose(m.eps_xx, 1e-4) and np.allclose(m.eps_divergence, -(m.eps_xx + m.eps_yy)))
m = mock(); load(m, strain_rates={'time': t, 'eps': np.vstack([ts['eps_xx'], ts['eps_yy'], ts['eps_xy']]).T})
check('dict with eps (Nx3) accepted', np.allclose(m.eps_xx, 1e-4))
m = mock(strain_rate_const=2e-4); load(m)
check('json strain_rate_const (number)', np.allclose(m.eps_divergence, 2e-4))
m = mock(strain_rate_const=[1e-4, 2e-4, 0.0]); load(m)
check('json strain_rate_const (list of 3)', np.allclose(m.eps_divergence, -3e-4))
# spin
m = mock(); load(m, spin=True, strain_rates=1e-4)
check('spin: constant, length = stp', np.allclose(m.eps_divergence, 1e-4) and len(m.eps_divergence) == 50)
m = mock(spinup_climate_type='initial'); load(m, spin=True, climateTS={'time': t, 'eps_divergence': np.linspace(1e-4, 2e-4, len(t))})
check("spin 'initial' uses first value", np.allclose(m.eps_divergence, 1e-4))
m = mock(); load(m, spin=True, climateTS={'time': t, 'eps_divergence': np.linspace(1e-4, 2e-4, len(t))})
check("spin 'mean' uses mean", np.allclose(m.eps_divergence, 1.5e-4))

# 3. priority: argument > climateTS > json > csv
m = mock(strain_rate_const=7e-4); load(m, strain_rates=1e-4, climateTS={'time': t, 'eps_divergence': np.full(len(t), 3e-4)})
check('priority: argument first', np.allclose(m.eps_divergence, 1e-4))
m = mock(strain_rate_const=7e-4); load(m, climateTS={'time': t, 'eps_divergence': np.full(len(t), 3e-4)})
check('priority: climateTS over json', np.allclose(m.eps_divergence, 3e-4))
m = mock(strain_rate_const=7e-4); load(m)
check('priority: json over csv', np.allclose(m.eps_divergence, 7e-4))

# 4. errors
for name, kw in (('missing time', dict(strain_rates={'eps': np.ones(5)})), ('length mismatch', dict(strain_rates={'time': np.arange(4.), 'eps': np.ones(5)})),
                 ('4 components', dict(strain_rates={'time': np.arange(4.), 'eps': np.ones((4, 4))}))):
    try: load(mock(), **kw); check(f'raises: {name}', False)
    except (KeyError, ValueError): check(f'raises: {name}', True)
m = mock(); del m.c['InputFileNameStrain']
try: load(m); check('raises: no source at all', False)
except KeyError: check('raises: no source at all', True)

# 5. horizontal_strain_rates accessor
m = mock(); load(m, strain_rates=[1e-4, 2e-4, 0.])
exx, eyy, exy = S.horizontal_strain_rates(m, 3)
check('accessor in 1/s, on', np.isclose(exx, 1e-4 / S.S_PER_YEAR) and np.isclose(eyy, 2e-4 / S.S_PER_YEAR))
m.c['horizontal_divergence'] = False
check('accessor zero when horizontal_divergence off', S.horizontal_strain_rates(m, 3) == (0.0, 0.0, 0.0))

# 6. end to end (dataframe path): strain columns in the forcing dataframe reach the model via climateTS
import pandas as pd, h5py
import RCMpkl_to_spin as RCM
from firn_density_nospin import FirnDensityNoSpin
def run(tag, edit=None, df_extra=None, **kw):
    c = json.load(open('example_df.json')); c.update(resultsFolder=f'/tmp/cfm_strain_{tag}', horizontal_divergence=True, strain_softening=False)
    if edit: c.update(edit)
    shutil.rmtree(c['resultsFolder'], ignore_errors=True)
    df = pd.read_pickle('CFMinput_example/' + c['DFfile'])
    if df_extra is not None:
        for k, v in df_extra.items(): df[k] = v
    with contextlib.redirect_stdout(io.StringIO()):
        cts, spy, d1, d2, db, SEBf = RCM.makeSpinFiles(df, timeres=c['DFresample'], melt=c['MELT'], desired_depth=c['H'] - c['HbaseSpin'])
        json.dump(c, open(f'/tmp/cfm_strain_{tag}.json', 'w'))
        FirnDensityNoSpin(f'/tmp/cfm_strain_{tag}.json', climateTS=cts, NewSpin=True, **kw).time_evolve()
    f = h5py.File(c['resultsFolder'] + '/CFMresults.hdf5', 'r'); return f['density'][:], f['depth'][:]
base = run('const', edit={'strain_rate_const': 3e-4})
arg = run('arg', strain_rates=3e-4)
col = run('col', df_extra={'eps_divergence': 3e-4})
none = run('none', edit={'horizontal_divergence': False})
same = lambda a, b: np.allclose(a[0], b[0], rtol=1e-6, equal_nan=True) and np.allclose(a[1], b[1], rtol=1e-6, equal_nan=True)
check('e2e: json constant == strain_rates argument', same(base, arg))
check('e2e: json constant == dataframe column', same(base, col))
check('e2e: divergence changes the result (sanity)', not same(base, none))
print('\nRESULT:', 'FAIL' if fails else 'PASS', fails)
sys.exit(1 if fails else 0)
