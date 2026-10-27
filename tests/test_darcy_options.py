'''
Options of the Darcy scheme. They were read in one block, so one missing key made all of them take their
defaults, and eps_cvg_darcy was read and then overwritten.
'''

import numpy as np
import pytest

import melt

DARCY_KEYS = list(melt._DARCY_DEFAULTS)


@pytest.fixture
def solver_calls(monkeypatch):
    '''Record the convergence tolerance that darcyscheme gives to the two flux solvers.'''
    calls = []
    for name in ('flux_bisection', 'flux_newtonraphson'):
        original = getattr(melt, name)

        def wrapper(*args, _original=original, **kwargs):
            calls.append(args[-1])   # eps_cvg is the last argument
            return _original(*args, **kwargs)

        monkeypatch.setattr(melt, name, wrapper)
    return calls


def _wet_cold_column(synthetic_column, **config):
    # Melt on a cold column with water flowing down through it.
    return synthetic_column(n=80, rho_profile=lambda z: 350.0 + 60.0 * z, T_profile=268.15, melt=0.03, config=config)


@pytest.mark.scheme
def test_eps_cvg_darcy_is_used(synthetic_column, solver_calls):
    melt.darcyscheme(_wet_cold_column(synthetic_column, eps_cvg_darcy=7.5e-4), 0)
    assert len(solver_calls) > 0
    assert set(solver_calls) == {7.5e-4}


@pytest.mark.scheme
def test_a_missing_option_does_not_reset_the_others(synthetic_column, solver_calls, capsys):
    melt._darcy_missing_reported.clear()
    column = _wet_cold_column(synthetic_column, eps_cvg_darcy=7.5e-4)
    del column.c['dtmax_darcy']
    melt.darcyscheme(column, 0)
    assert set(solver_calls) == {7.5e-4}   # it used to fall back to 0.1e-3 with the missing key
    printed = capsys.readouterr().out
    assert '"dtmax_darcy" is not in the .json' in printed
    assert printed.count('is not in the .json') == 1   # only the missing one is reported


@pytest.mark.scheme
def test_a_missing_option_is_reported_once_per_run(synthetic_column, capsys):
    melt._darcy_missing_reported.clear()
    for step in range(2):
        column = _wet_cold_column(synthetic_column)
        del column.c['slope_proxy_darcy']
        melt.darcyscheme(column, 0)
    assert capsys.readouterr().out.count('"slope_proxy_darcy" is not in the .json') == 1


@pytest.mark.unit
def test_defaults_are_the_documented_ones():
    assert melt._DARCY_DEFAULTS['eps_cvg_darcy'] == pytest.approx(0.1e-3)
    assert melt._DARCY_DEFAULTS['RhoImp_darcy'] == 873.0
    assert melt._DARCY_DEFAULTS['rholim_dcy'] == 910.0
    assert melt._DARCY_DEFAULTS['runoff_method_darcy'] == 'Darcy'
