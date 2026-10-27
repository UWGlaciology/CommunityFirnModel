'''
Short model runs. Each test checks that the model completes a run with one
config combination.
'''

import pytest


@pytest.mark.conservation
def test_melt_with_sublimation_completes(short_config):
    # SUBLIM true runs the sublimation block, which sets subLWCvol at each step.
    run = short_config(years=2, MELT=True, SUBLIM=True, melt_strict=True)
    run.run()
    assert run.results_path.exists()


@pytest.mark.conservation
def test_melt_without_sublimation_completes(short_config):
    # SUBLIM false skips the sublimation block, so subLWCvol keeps the value
    # that __init__ gives it. See UWGlaciology/CommunityFirnModel#19.
    run = short_config(years=2, MELT=True, SUBLIM=False, melt_strict=True)
    run.run()
    assert run.results_path.exists()
