'''
Tests for transient_solve_enthalpy when a layer holds liquid water below 0 C.

The prefsnowpack and resingledomain schemes can return such layers: water in
the preferential-flow domain, round-off residue of the prewetting, and the
temperature average in combineCFM. The solver must refreeze that water with
its latent heat, as if the layer was in equilibrium before the step.
'''

import numpy as np
import pytest

import solver
from constants import CP_I, CP_W, LF_I

N_LAYERS = 20
DZ = 0.1  # layer thickness [m]
RHO = 450.0  # layer density [kg m-3]
DT = 15 * 86400.0  # step length [s]
COLD_WET_LAYER = 5  # index of the layer with liquid water below 0 C


def _solve(T_C, th_liquid):
    z_edges = np.arange(N_LAYERS + 1) * DZ
    Z_P = z_edges[:-1] + DZ / 2
    Gamma_P = np.full(N_LAYERS, 0.5)
    th_solid = RHO - th_liquid
    return solver.transient_solve_enthalpy(z_edges, Z_P, DT, Gamma_P, T_C, th_liquid,
                                           th_solid, 0, bc_u=T_C[0])


def _cold_column(cold_liquid):
    T_C = np.linspace(-5.0, -12.0, N_LAYERS)
    th_liquid = np.zeros(N_LAYERS)
    th_liquid[COLD_WET_LAYER] = cold_liquid
    return T_C, th_liquid


def _refrozen_first(T_C, th_liquid):
    '''Refreeze all liquid at constant enthalpy (valid when the layer stays below 0 C).'''
    enthalpy = (RHO - th_liquid) * CP_I * T_C + th_liquid * (CP_W * T_C + LF_I)
    T_refrozen_C = enthalpy / (RHO * CP_I)
    assert np.all(T_refrozen_C <= 0.0)
    return T_refrozen_C, np.zeros(N_LAYERS)


@pytest.mark.unit
@pytest.mark.parametrize('cold_liquid', [1e-19, 2.0, 20.0])  # kg m-3
def test_cold_liquid_refreezes_with_its_latent_heat(cold_liquid):
    result = _solve(*_cold_column(cold_liquid))
    reference = _solve(*_refrozen_first(*_cold_column(cold_liquid)))
    np.testing.assert_allclose(result['TzC_return'], reference['TzC_return'], rtol=0, atol=1e-8)
    np.testing.assert_allclose(result['th_liquid'], reference['th_liquid'], rtol=0, atol=1e-10)
    np.testing.assert_allclose(result['th_solid'], reference['th_solid'], rtol=0, atol=1e-8)


@pytest.mark.unit
@pytest.mark.parametrize('cold_liquid', [1e-19, 2.0])  # kg m-3
def test_no_clamp_reported_for_cold_liquid(cold_liquid):
    result = _solve(*_cold_column(cold_liquid))
    assert result['claw_mushy'] == 0.0
