'''
Sign convention of the turbulent fluxes in the SEB forcing: QH and QL are positive INTO the snow.
flux_sign_check warns when QH looks positive upward (the sign of MERRA-2 HFLUX and EFLUX), for
example when the multiplication by -1 was forgotten while making the forcing dataframe.
'''

import contextlib
import io

import numpy as np
import pandas as pd
import pytest

import RCMpkl_to_spin as RCM
from cfm_testing import CFM_MAIN
from SEB import SurfaceEnergyBudget, flux_sign_check

SITES = ['66.5_-46.25', '72.5_-38.75']   # DYE-2 and Summit, Greenland (daily MERRA-2, fluxes stored positive into the snow)
WARNING = 'WARNING (SEB)'


def _example(site, years=('1980', '1989')):
    return pd.read_pickle(CFM_MAIN / 'CFMinput_example' / f'CFM_example_{site}.pkl').loc[years[0]:years[1]].copy()


def _seb(site, sign, with_sebfluxes=True):
    '''SurfaceEnergyBudget made from the example forcing through makeSpinFiles, with QH and QL
    multiplied by `sign` (-1 reproduces forgetting the -1 on MERRA-2 fluxes). Returns (seb, printed text).'''
    df = _example(site)
    df['QH'] *= sign
    df['QL'] *= sign
    with contextlib.redirect_stdout(io.StringIO()):
        CD, _, _, _, _, SEBfluxes = RCM.makeSpinFiles(
            df, timeres='1d', spin_date_st=1980.0, spin_date_end=1982.0, melt=True,
            desired_depth=20, SEB=True, num_reps=2)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        seb = SurfaceEnergyBudget(dict(SEB_TL_thick=0.05), CD, 0, SEBfluxes if with_sebfluxes else None, 0)
    return seb, out.getvalue()


@pytest.mark.unit
@pytest.mark.parametrize('site', SITES)
def test_no_warning_when_fluxes_are_positive_into_the_snow(site):
    seb, text = _seb(site, +1)
    assert WARNING not in text
    check = flux_sign_check(seb.QH, seb.T2m, seb.TSKIN, seb.time_in)
    assert check['judged'] and not check['looks_upward']
    assert check['r'] > 0.2 and check['mean'] > 0


@pytest.mark.unit
@pytest.mark.parametrize('site', SITES)
def test_warning_when_the_minus_one_was_forgotten(site):
    # MERRA-2 HFLUX and EFLUX are positive upward: the dataframe has to hold -1 times them
    seb, text = _seb(site, -1)
    assert WARNING in text
    assert 'positive UPWARD' in text and 'multiplied by -1' in text


@pytest.mark.unit
def test_no_sign_change_when_climateTS_is_given_without_sebfluxes():
    # before 4.1.0 this path multiplied QH and QL by -1; all paths now expect fluxes positive into the snow
    seb, text = _seb(SITES[0], +1, with_sebfluxes=False)
    assert WARNING not in text
    assert np.mean(seb.QH) > 0


@pytest.mark.unit
def test_the_check_is_only_made_for_prescribed_fluxes():
    # with calculated turbulent fluxes the QH of the forcing is not used
    df = _example(SITES[0])
    df['QH'] *= -1
    with contextlib.redirect_stdout(io.StringIO()):
        CD, _, _, _, _, SEBfluxes = RCM.makeSpinFiles(df, timeres='1d', spin_date_st=1980.0, spin_date_end=1982.0,
                                                      melt=True, desired_depth=20, SEB=True, num_reps=2)
        for key, val in (('WS', 3.0), ('QV2m', 1e-3), ('PS', 70000.0)):
            SEBfluxes[key] = np.full_like(SEBfluxes['QH'], val)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        SurfaceEnergyBudget(dict(SEB_TL_thick=0.05, SEB_turb_method='bulk_richardson'), CD, 0, SEBfluxes, 0)
    assert WARNING not in out.getvalue()


@pytest.mark.unit
@pytest.mark.parametrize('site', SITES)
def test_short_records_are_not_judged(site):
    # the correlation is negative in winter even for correctly signed fluxes
    df = _example(site).loc['1980-12':'1981-02']
    t = np.arange(len(df)) / 365.25
    for sign in (+1, -1):
        check = flux_sign_check(sign * df['QH'], df['T2m'], df['TSKIN'], t)
        assert not check['judged'] and not check['looks_upward']


@pytest.mark.unit
def test_constant_or_missing_data_are_not_judged():
    t = np.arange(1000) / 365.25 * 3
    assert not flux_sign_check(np.ones(1000), np.arange(1000.0), np.zeros(1000), t)['judged']
    assert not flux_sign_check(np.full(1000, np.nan), np.arange(1000.0), np.zeros(1000), t)['judged']
    assert not flux_sign_check(np.arange(50.0), np.arange(50.0), np.zeros(50))['judged']


@pytest.mark.unit
def test_fluxes_that_agree_with_the_temperature_difference_are_accepted_and_the_opposite_flagged():
    rng = np.random.default_rng(0)
    n = 4000
    t = np.arange(n) / 365.25 / 2
    dT = rng.normal(0.5, 2.0, n)
    QH = 8.0 * dT + 10.0 + rng.normal(0, 15, n)    # positive into the snow
    assert not flux_sign_check(QH, 250 + dT, np.full(n, 250.0), t)['looks_upward']
    assert flux_sign_check(-QH, 250 + dT, np.full(n, 250.0), t)['looks_upward']
