'''
Smoke tests for the shared fixtures. They show that the fixtures give a state
that the model code accepts.
'''

import numpy as np
import pytest

import melt


@pytest.mark.unit
def test_synthetic_column_has_state_attributes(synthetic_column):
    column = synthetic_column(n=20, rho_profile=lambda z: 350.0 + 50.0 * z, dz=0.05)
    for name in column.STATE_ATTRIBUTES:
        assert hasattr(column, name), name
    assert column.z[0] == 0.0
    np.testing.assert_allclose(np.diff(column.z), column.dz[:-1])
    np.testing.assert_allclose(column.mass, column.rho * column.dz)
    assert column.rho[0] == 350.0


@pytest.mark.scheme
def test_bucket_closes_water_budget_on_synthetic_column(synthetic_column):
    column = synthetic_column(n=60, rho_profile=lambda z: 350.0 + 40.0 * z,
                              T_profile=265.0, melt=0.02, config={'melt_strict': True})
    melt.bucket(column, 0)
    diagnostics = column.melt_diagnostics
    assert all(check['ok'] for check in diagnostics['checks'].values())
    assert diagnostics['refrozen'] > 0.0
    water_out_and_stored = diagnostics['lwc_end'] + diagnostics['refrozen'] + diagnostics['runoff']
    assert water_out_and_stored == pytest.approx(diagnostics['water_in'], abs=1e-12)


@pytest.mark.scheme
def test_darcyscheme_runs_on_synthetic_column(synthetic_column):
    column = synthetic_column(n=60, rho_profile=lambda z: 350.0 + 40.0 * z,
                              T_profile=265.0, melt=0.02, config={'melt_strict': True})
    melt.darcyscheme(column, 0)
    assert all(check['ok'] for check in column.melt_diagnostics['checks'].values())


@pytest.mark.conservation
def test_short_run_completes_in_strict_mode(short_config):
    run = short_config(years=2, melt_strict=True)
    model = run.run()
    assert run.results_path.exists()
    assert (run.folder / 'diagnostics_enthalpy.csv').exists()
    assert model.melt_diagnostics['scheme'] == 'bucket'
