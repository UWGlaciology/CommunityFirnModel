'''
Shared fixtures for the CFM test suite. The helpers are in cfm_testing.py.
'''

import itertools

import pytest

from cfm_testing import CFM_MAIN, ShortRun, SyntheticColumn


@pytest.fixture(autouse=True)
def cfm_main_cwd(monkeypatch):
    '''
    Run each test in CFM_main/. Some modules open files with paths relative to
    that folder (for example AirConfig.json). This is a temporary workaround
    while the CFM is a set of flat scripts and not an installable package.
    '''
    monkeypatch.chdir(CFM_MAIN)


@pytest.fixture
def synthetic_column():
    '''
    Factory for SyntheticColumn. Usage:
        column = synthetic_column(n=40, rho_profile=lambda z: 350 + 10*z, melt=0.01)
    '''
    return SyntheticColumn


@pytest.fixture
def short_config(tmp_path):
    '''
    Factory for ShortRun. Each call gets its own folder in tmp_path. Usage:
        run = short_config(years=2, liquid='darcy')
        model = run.run()
    '''
    counter = itertools.count()

    def make(years=2, start_year=1980, num_reps=None, **overrides):
        folder = tmp_path / f'run{next(counter)}'
        return ShortRun(folder, years=years, start_year=start_year,
                        num_reps=num_reps, **overrides)

    return make
