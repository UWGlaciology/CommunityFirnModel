'''
Tests for the check recording in melt.py when a check fails. The fixture smoke
tests cover only checks that pass.
'''

from types import SimpleNamespace

import numpy as np
import pytest

import melt


def _stub_model(**config):
    '''A stub with the two attributes that the diagnostics helpers use.'''
    model = SimpleNamespace(c=config)
    melt._start_diagnostics(model, 'bucket', 7)
    return model


@pytest.mark.unit
def test_failed_check_raises_in_strict_mode():
    model = _stub_model(melt_strict=True)
    with pytest.raises(melt.MeltCheckError, match="bucket: check 'budget_1' failed at step 7"):
        melt._record_check(model, 'budget_1', False, residual=2e-3, threshold=1e-3)
    assert model.melt_diagnostics['checks']['budget_1']['ok'] is False


@pytest.mark.unit
def test_failed_check_is_recorded_without_strict_mode():
    model = _stub_model()  # melt_strict is not in the config: the default is False
    melt._record_check(model, 'cold_1', False, n_wet_cold_layers=np.int64(3), cold_lwc_removed=np.float64(0.0))
    entry = model.melt_diagnostics['checks']['cold_1']
    assert entry == {'ok': False, 'n_wet_cold_layers': 3, 'cold_lwc_removed': 0.0}
    assert type(entry['n_wet_cold_layers']) is int  # stored as a Python scalar, not a numpy one
