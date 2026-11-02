'''
Tests of the bare-ice melt enhancement (IceMeltFactor, IceMeltRho): the melt forcing is
multiplied by IceMeltFactor at time steps where the density of the top model node is
above IceMeltRho. Short DYE-2 runs, bucket scheme, using the end-of-run totals that
time_evolve prints. A DYE-2 surface is snow (about 350 kg m-3) at nearly every step, so
an IceMeltRho below the snow density is used to make the condition true at every step.
'''

import contextlib
import io
import os

import pytest

from cfm_testing import CFM_MAIN, ShortRun
from test_water_budget import DYE2_FORCING, _printed_totals

FACTOR = 1.27
CASES = {
    'off': {},
    'factor_one': dict(IceMeltFactor=1.0, IceMeltRho=300.),
    'snow_surface': dict(IceMeltFactor=FACTOR),                       # default IceMeltRho = 800: snow on top does not trigger
    'high_threshold': dict(IceMeltFactor=FACTOR, IceMeltRho=1000.),   # nothing is denser than this
    'always': dict(IceMeltFactor=FACTOR, IceMeltRho=300.),            # every surface node counts as bare ice
    'always_ponding': dict(IceMeltFactor=FACTOR, IceMeltRho=300., Ponding=True),
}


@pytest.fixture(scope='module')
def totals(tmp_path_factory):
    '''Run each case once; returns {case: end-of-run totals [m w.e.]}.'''
    base = tmp_path_factory.mktemp('ice_melt_factor')
    results = {}
    cwd = os.getcwd()
    os.chdir(CFM_MAIN)  # some modules open files relative to CFM_main
    try:
        for name, options in CASES.items():
            run = ShortRun(base / name, years=2, start_year=2011, liquid='bucket',
                           DFfile=DYE2_FORCING, **options)
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                run.run()
            results[name] = _printed_totals(out.getvalue())
    finally:
        os.chdir(cwd)
    return results


def test_factor_one_is_identical_to_omitting_the_keys(totals):
    assert totals['factor_one'] == totals['off']


@pytest.mark.parametrize('case', ['snow_surface', 'high_threshold'])
def test_snow_surface_does_not_trigger_the_factor(totals, case):
    assert totals[case] == totals['off']


def test_melt_is_scaled_by_the_factor_when_the_surface_counts_as_ice(totals):
    base, scaled = totals['off'], totals['always']
    assert base['meltvol'] > 0
    assert scaled['meltvol'] == pytest.approx(FACTOR * base['meltvol'], rel=1e-6)
    # the rain is not scaled
    assert scaled['rainvol'] == pytest.approx(base['rainvol'], rel=1e-9)
    assert scaled['melt_rain'] == pytest.approx(FACTOR * base['meltvol'] + base['rainvol'], rel=1e-6)


@pytest.mark.conservation
@pytest.mark.parametrize('case', ['always', 'always_ponding'])
def test_water_budget_closes_with_the_factor(totals, case):
    t = totals[case]
    assert t['refreeze'] + t['runoff'] + t['lwc'] == pytest.approx(t['melt_rain'], rel=1e-3)


@pytest.mark.parametrize('options, expected', [
    ({}, None),
    ({'IceMeltFactor': 1.0, 'SEB': True}, None),
    ({'IceMeltFactor': 1.27}, None),
    ({'IceMeltFactor': 1.27, 'SEB': True}, "not used when 'SEB' is true"),
    ({'IceMeltFactor': 1.27, 'MELT': False}, "no effect when 'MELT' is false"),
])
def test_validator_warns_when_ice_melt_factor_has_no_effect(capsys, options, expected):
    from config_validation import validate_config
    from cfm_testing import example_config
    validate_config({**example_config(), **options})
    out = capsys.readouterr().out
    if expected is None:
        assert 'IceMeltFactor' not in out
    else:
        assert expected in out
