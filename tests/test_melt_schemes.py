'''
Single calls of the meltwater schemes on a synthetic column.
'''

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
