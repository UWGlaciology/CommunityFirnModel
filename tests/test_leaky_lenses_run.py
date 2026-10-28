'''
End-to-end tests of the LeakyLenses option: short CFM runs (bucket scheme)
driven by a small synthetic forcing with enough melt to form ice lenses.
Each run takes a few seconds.
'''

import json
import subprocess
import sys

import h5py
import numpy as np
import pytest

from conftest import CFM_MAIN

RHO_W, RHO_I = 1000.0, 917.0
N_YEARS = 3
START = 2000.0


def write_forcing(folder):
    '''Daily forcing: seasonal cycle plus weather noise around -9 C, 0.4 m w.e./yr
    of snow and degree-day melt, so the firn gets ice lenses. CFM reads
    accumulation and melt as m ice eq./yr.'''
    n_days = int(N_YEARS * 365.25)
    t = START + np.arange(n_days) / 365.25
    rng = np.random.default_rng(0)
    T = -9.0 - 8.0 * np.cos(2 * np.pi * t) + rng.normal(0.0, 3.8, n_days)
    forcing = {
        'temperature.csv': np.minimum(T, 0.0) + 273.15,
        'accumulation.csv': np.full(n_days, 0.4 * RHO_W / RHO_I),
        'melt.csv': 0.004 * np.maximum(T, 0.0) * 365.25 * RHO_W / RHO_I,
        # SUBLIM is kept on with a zero sublimation file: with MELT on and SUBLIM off,
        # the run currently fails on an unrelated issue
        'sublimation.csv': np.zeros(n_days),
    }
    for name, values in forcing.items():
        np.savetxt(folder / name, np.vstack([t, values]), delimiter=',')


def run_cfm(folder, **options):
    '''Run CFM in `folder` with example.json modified for a short melt run plus
    `options`; a value of None removes the key. Returns the results file.'''
    folder.mkdir()
    write_forcing(folder)
    with open(CFM_MAIN / 'example.json') as f:
        c = json.load(f)
    c.update({
        'input_type': 'csv',
        'InputFileFolder': str(folder),
        'InputFileNameTemp': 'temperature.csv',
        'InputFileNamebdot': 'accumulation.csv',
        'InputFileNamemelt': 'melt.csv',
        'InputFileNameSublim': 'sublimation.csv',
        'resultsFolder': str(folder / 'results'),
        'physRho': 'Crocus',
        'MELT': True,
        'RAIN': False,
        'SUBLIM': True,
        'liquid': 'bucket',
        'isoDiff': False,
        'FirnAir': False,
        'AutoSpinUpTime': False,
        'yearSpin': 1.0,
        'H': 50.0,
        'HbaseSpin': 40.0,              # 10 m domain
        'doublegrid': False,            # the example grid2bottom (10 m) needs a deeper domain
        'TWriteStart': START,
        'outputs': ['density', 'depth', 'LWC', 'meltoutputs'],
    })
    for key, value in options.items():
        if value is None:
            c.pop(key, None)
        else:
            c[key] = value
    config = folder / 'config.json'
    with open(config, 'w') as f:
        json.dump(c, f, indent=2)
    result = subprocess.run([sys.executable, str(CFM_MAIN / 'main.py'), str(config), '-n'],
                            cwd=folder, capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    assert 'Mass conservation error' not in result.stdout
    return folder / 'results' / c['resultsFileName']


def read(results):
    with h5py.File(results, 'r') as f:
        return {k: f[k][:] for k in ('density', 'LWC', 'meltvol', 'refreeze', 'runoff')}


def assert_identical(a, b):
    for key in a:
        np.testing.assert_array_equal(a[key], b[key], err_msg=key)


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    base = tmp_path_factory.mktemp('leaky_lenses')
    return {
        'default': read(run_cfm(base / 'default', LeakyLenses=None, LensDecayLength=None)),
        'off': read(run_cfm(base / 'off', LeakyLenses=False, LensDecayLength=0.3)),
        'thickimp0': read(run_cfm(base / 'thickimp0', LeakyLenses=False, ThickImp=0.0)),
        'tiny_decay': read(run_cfm(base / 'tiny_decay', LeakyLenses=True, LensDecayLength=1e-9)),
        'leaky': read(run_cfm(base / 'leaky', LeakyLenses=True, LensDecayLength=0.1)),
    }


def test_off_is_identical_to_omitting_the_keys(runs):
    assert_identical(runs['off'], runs['default'])


def test_tiny_decay_length_is_identical_to_thickimp_zero(runs):
    assert_identical(runs['tiny_decay'], runs['thickimp0'])


def test_leaky_lenses_change_the_result(runs):
    # the forcing must actually form lenses, or the tests above prove little
    assert runs['leaky']['runoff'][:, 1].sum() < runs['thickimp0']['runoff'][:, 1].sum()
    assert runs['leaky']['refreeze'][:, 1].sum() > runs['thickimp0']['refreeze'][:, 1].sum()


def test_leaky_run_closes_the_water_budget(runs):
    r = runs['leaky']
    melt = r['meltvol'][1:, 1].sum()
    out = r['refreeze'][1:, 1].sum() + r['runoff'][1:, 1].sum()
    stored_change = np.nansum(r['LWC'][-1, 1:]) - np.nansum(r['LWC'][0, 1:])
    assert out + stored_change == pytest.approx(melt, rel=1e-3)
