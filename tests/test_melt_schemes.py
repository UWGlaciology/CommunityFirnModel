'''
Single calls of the meltwater schemes on a synthetic column.
'''

import numpy as np
import pytest

import melt
from constants import RHO_I


@pytest.mark.scheme
@pytest.mark.parametrize('scheme', ['bucket', 'darcyscheme'])
def test_dh_melt_is_melted_thickness(synthetic_column, scheme):
    # With a uniform density, the melted thickness is melt mass / density.
    # darcyscheme melts in sub-steps, so this also checks its sum. The column
    # is at the melting point: refreeze between sub-steps would make the
    # surface denser and the later sub-steps would melt a thinner layer.
    density, melt_ice_equivalent = 400.0, 0.02
    column = synthetic_column(n=60, rho_profile=density, T_profile=273.15, melt=melt_ice_equivalent)
    scheme_outputs = getattr(melt, scheme)(column, 0)
    assert len(scheme_outputs) == 13
    dh_melt = scheme_outputs[-1]
    assert dh_melt == pytest.approx(-melt_ice_equivalent * RHO_I / density, rel=1e-9)


class _FluxSolverCalls:
    '''Wraps the two flux solvers that darcyscheme calls, and records each call.'''

    def __init__(self, monkeypatch):
        import darcy_funcs
        self.calls = []
        for name in ('flux_bisection', 'flux_newtonraphson'):
            monkeypatch.setattr(melt, name, self._wrap(name, getattr(darcy_funcs, name)))

    def _wrap(self, name, solver):
        import darcy_funcs

        def wrapper(guess, lwc_available, lwc_accommodated, th_i, th_s, lwc, dz, alpha, n, m, eps_cvg):
            saturation = darcy_funcs.thetae_update(guess, th_i, th_s, lwc, dz)
            head = darcy_funcs.phead_vG(alpha, n, m, saturation)
            flux = solver(guess, lwc_available, lwc_accommodated, th_i, th_s, lwc, dz, alpha, n, m, eps_cvg)
            self.calls.append(dict(
                solver=name, f_eq=head[0] - head[1] - dz.sum() / 2, flux=flux,
                n_nodes=len(lwc_available), upper_bound=min(lwc_available[0], lwc_accommodated[1])))
            return flux
        return wrapper


@pytest.fixture
def darcy_wetting_front(synthetic_column, monkeypatch):
    '''
    One day of melt into dry, temperate firn, with the flux solver calls recorded.
    The front moves into dry firn, so many interfaces start with f_eq < -1.
    '''
    calls = _FluxSolverCalls(monkeypatch)
    column = synthetic_column(n=60, rho_profile=lambda z: 350.0 + 60.0 * z, T_profile=273.15,
                              dz=0.05, r2=(5e-4) ** 2, melt=0.02, dt=86400.0)
    melt.darcyscheme(column, 0)
    assert len(calls.calls) > 100
    return column, calls.calls


@pytest.mark.scheme
def test_darcy_flux_solvers_get_the_two_interface_nodes(darcy_wetting_front):
    # The upper bound of the flux at interface j1 is min(LWCav[j1], glwcacm[j1+1]).
    # The solvers read it from elements 0 and 1, so they must get the two nodes
    # at the interface, not the full column.
    _, calls = darcy_wetting_front
    assert {call['n_nodes'] for call in calls} == {2}


@pytest.mark.scheme
def test_darcy_fluxes_are_within_the_available_water_and_pore_space(darcy_wetting_front):
    column, calls = darcy_wetting_front
    for call in calls:
        assert 0.0 <= call['flux'] <= call['upper_bound'] * (1 + 1e-12), call
    saturated_water = (917.0 - column.rho) / 917.0 * 0.917 * column.dz
    assert np.all(column.LWC <= saturated_water)
