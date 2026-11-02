'''
End-to-end tests of the LeakyLenses option (bucket scheme) over a short DYE-2
run (2011-2012, which includes the large 2012 melt), using the end-of-run totals that time_evolve prints. The option and the
routing function are by Aslak Grinsted (pull request #25).
'''

import contextlib
import io
import json
import os

import pytest

from cfm_testing import CFM_MAIN, ShortRun
from test_water_budget import DYE2_FORCING, _printed_totals

CASES = {
    'default': {},                                                   # keys absent from the config (popped below)
    'off': dict(LeakyLenses=False, LensDecayLength=0.3),
    'thickimp0': dict(LeakyLenses=False, ThickImp=0.0),
    'tiny_decay': dict(LeakyLenses=True, LensDecayLength=1e-9),
    'leaky': dict(LeakyLenses=True, LensDecayLength=0.1),
    'leaky_ponding': dict(LeakyLenses=True, LensDecayLength=0.1, Ponding=True),
}


@pytest.fixture(scope='module')
def totals(tmp_path_factory):
    '''Run each case once; returns {case: end-of-run totals [m w.e.]}.'''
    base = tmp_path_factory.mktemp('leaky_lenses')
    results = {}
    cwd = os.getcwd()
    os.chdir(CFM_MAIN)  # some modules open files relative to CFM_main
    try:
        for name, options in CASES.items():
            run = ShortRun(base / name, years=2, start_year=2011, liquid='bucket', DFfile=DYE2_FORCING, **options)
            if name == 'default':
                for key in ('LeakyLenses', 'LensDecayLength'):
                    run.config.pop(key, None)
                run.config_path.write_text(json.dumps(run.config))
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                run.run()
            results[name] = _printed_totals(out.getvalue())
    finally:
        os.chdir(cwd)
    return results


def test_off_is_identical_to_omitting_the_keys(totals):
    assert totals['off'] == totals['default']


def test_tiny_decay_length_is_identical_to_thickimp_zero(totals):
    assert totals['tiny_decay'] == totals['thickimp0']


def test_leaky_lenses_change_the_result(totals):
    # the forcing must form lenses that matter, or the tests above prove little
    assert totals['leaky']['runoff'] < totals['thickimp0']['runoff']
    assert totals['leaky']['refreeze'] > totals['thickimp0']['refreeze']


@pytest.mark.conservation
@pytest.mark.parametrize('case', ['leaky', 'leaky_ponding'])
def test_leaky_run_closes_the_water_budget(totals, case):
    t = totals[case]
    assert t['refreeze'] + t['runoff'] + t['lwc'] == pytest.approx(t['melt_rain'], rel=1e-3)
