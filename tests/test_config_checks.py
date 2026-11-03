"""
Checks made by config_validation.validate_config that change or report a setting.
"""

import pytest

from config_validation import validate_config
from cfm_testing import example_config


@pytest.mark.unit
@pytest.mark.parametrize('melt', [False, None])
def test_seb_without_melt_switches_melt_on(capsys, melt):
    c = example_config()
    c['SEB'] = True
    if melt is None:
        c.pop('MELT')
    else:
        c['MELT'] = melt
    validate_config(c)
    assert c['MELT'] is True
    out = capsys.readouterr().out
    assert "'SEB' is true but 'MELT' is false" in out
    assert "Switching 'MELT' to true" in out


@pytest.mark.unit
@pytest.mark.parametrize('seb, melt', [(True, True), (False, False), (False, True)])
def test_other_seb_melt_combinations_are_left_alone(capsys, seb, melt):
    c = example_config()
    c['SEB'], c['MELT'] = seb, melt
    validate_config(c)
    assert c['MELT'] is melt
    assert "Switching 'MELT'" not in capsys.readouterr().out


@pytest.mark.unit
@pytest.mark.parametrize('method, sublim, warns', [
    ('bulk_richardson', False, True), ('monin_obukhov', False, True),
    ('bulk_richardson', True, False), ('prescribed', False, False)])
def test_seb_turbulent_flux_with_sublimation_off_warns(capsys, method, sublim, warns):
    # the latent heat flux is in the energy balance, but without SUBLIM the sublimation that goes with it is not applied
    c = example_config()
    c.update(SEB=True, MELT=True, SEB_turb_method=method, SUBLIM=sublim)
    validate_config(c)
    assert ("'SUBLIM' is false but 'SEB_turb_method'" in capsys.readouterr().out) is warns
