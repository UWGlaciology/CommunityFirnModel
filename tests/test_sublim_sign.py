'''
Sign convention of SUBLIM: negative = sublimation (mass loss), positive = deposition. It goes with the
sign of QL (positive into the snow): sublimation cools the surface, so SUBLIM and QL are positively
correlated. Before 4.1.0 the example forcing held the raw MERRA-2 sign (EVAP positive for sublimation).
'''

import contextlib
import io

import numpy as np
import pandas as pd
import pytest

import RCMpkl_to_spin as RCM
from cfm_testing import CFM_MAIN

SITES = ['66.5_-46.25', '72.5_-38.75']   # DYE-2 and Summit, Greenland
WARNING = 'SUBLIM and QL are negatively correlated'


def _example(site, years=('1980', '1989')):
    return pd.read_pickle(CFM_MAIN / 'CFMinput_example' / f'CFM_example_{site}.pkl').loc[years[0]:years[1]].copy()


def _spin(df, **kw):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        res = RCM.makeSpinFiles(df, timeres='1d', spin_date_st=1980.0, spin_date_end=1982.0, melt=True,
                                desired_depth=20, num_reps=2, **kw)
    return res, out.getvalue()


@pytest.mark.unit
@pytest.mark.parametrize('site', SITES)
def test_example_dataframes_have_the_cfm_sign(site):
    df = _example(site)
    assert df['SUBLIM'].corr(df['QL']) > 0.99


@pytest.mark.unit
def test_example_csv_matches_the_dataframe_sign():
    # the csv-route sublimation (example_SUBLIM.csv, m ie/yr) has Summit-like values; the sign of its mean
    # must agree with that of the Summit dataframe (kg m-2 per day), which is in the CFM convention
    csv = np.loadtxt(CFM_MAIN / 'CFMinput_example' / 'example_SUBLIM.csv', delimiter=',')
    df = _example('72.5_-38.75', years=('1980', '2019'))
    assert np.sign(np.mean(csv[1])) == np.sign(df['SUBLIM'].mean())


@pytest.mark.unit
def test_legacy_evap_column_is_converted_with_the_sign_change():
    # the legacy example pickle holds the raw MERRA-2 name and sign
    legacy = pd.read_pickle(CFM_MAIN / 'CFMinput_example' / 'CFM_example_input_df.pkl').loc['1980':'1989'].copy()
    (cd, *_), _ = _spin(legacy, SEB=False)
    summit = _example('72.5_-38.75')
    (cd2, *_), _ = _spin(summit[['TSKIN', 'SUBLIM', 'BDOT', 'SMELT', 'RAIN']], SEB=False)
    assert np.mean(cd['SUBLIM']) * np.mean(cd2['SUBLIM']) > 0   # same sign as the CFM-convention dataframe
    n = min(len(cd['SUBLIM']), len(cd2['SUBLIM']))
    assert np.corrcoef(cd['SUBLIM'][-n:], cd2['SUBLIM'][-n:])[0, 1] > 0.99


@pytest.mark.unit
def test_sublim_column_is_not_changed_by_makeSpinFiles():
    df = _example('66.5_-46.25')
    (cd, *_), _ = _spin(df[['TSKIN', 'SUBLIM', 'BDOT', 'SMELT', 'RAIN']], SEB=False)
    raw = df['SUBLIM'].resample('1D').sum().values * 365.25 / 917
    assert np.corrcoef(cd['SUBLIM'][-len(raw):], raw)[0, 1] > 0.999
    assert np.mean(cd['SUBLIM'][-len(raw):]) * np.mean(raw) > 0


@pytest.mark.unit
@pytest.mark.parametrize('site', SITES)
def test_no_warning_with_the_cfm_sign(site):
    _, text = _spin(_example(site), SEB=True)
    assert WARNING not in text


@pytest.mark.unit
@pytest.mark.parametrize('site', SITES)
def test_warning_when_sublim_has_the_raw_merra_sign(site):
    df = _example(site)
    df['SUBLIM'] *= -1
    _, text = _spin(df, SEB=True)
    assert WARNING in text and 'multiplied by -1' in text


@pytest.mark.unit
def test_sublim_sign_check_edge_cases():
    rng = np.random.default_rng(0)
    q = rng.normal(size=200)
    assert RCM.sublim_sign_check(q, q)['consistent'] is True
    assert RCM.sublim_sign_check(-q, q)['consistent'] is False
    assert RCM.sublim_sign_check(q[:10], q[:10])['consistent'] is None          # too short
    assert RCM.sublim_sign_check(np.zeros(200), q)['consistent'] is None        # constant
    a = q.copy(); a[:5] = np.nan
    assert RCM.sublim_sign_check(a, q)['consistent'] is True                    # NaNs ignored
