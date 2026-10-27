'''
Total compaction without the horizontal strain (UWGlaciology/CommunityFirnModel#9):
the thickness change of the column by densification alone is kept apart from the
thinning by horizontal strain. The surface elevation change needs both.
'''

import h5py
import numpy as np
import pytest

from constants import S_PER_YEAR

SUMMIT = 'CFM_example_input_df.pkl'


def _run(short_config, **overrides):
    run = short_config(years=2, liquid='bucket', DFfile=SUMMIT, TWriteInt=1,
                       outputs=['density', 'depth', 'temperature', 'age', 'DIP', 'compaction'], **overrides)
    model = run.run()
    return run, model


@pytest.mark.conservation
def test_without_horizontal_strain_all_the_thickness_change_is_compaction(short_config):
    run, model = _run(short_config, horizontal_divergence=False)
    assert model.dsdz_strain == 0.0
    assert model.dsdz_comp == model.dsdz_step < 0.0


@pytest.mark.conservation
@pytest.mark.parametrize('divergence', [-2.0e-2, 2.0e-2], ids=['thins', 'thickens'])
def test_horizontal_strain_is_kept_out_of_the_compaction(short_config, divergence):
    # A single number is eps_divergence [1/yr]; the layer mass is scaled by 1 + eps_divergence*dt.
    run, model = _run(short_config, horizontal_divergence=True, strain_rate_const=divergence)
    dt_yr = model.dt[-1] / S_PER_YEAR
    thickness_after_densification = model.sdz_old + model.dsdz_comp
    assert model.dsdz_strain == pytest.approx(divergence * dt_yr * thickness_after_densification, rel=1e-6)
    assert model.dsdz_step == pytest.approx(model.dsdz_comp + model.dsdz_strain, rel=1e-12)
    assert model.dsdz_comp < 0.0   # the strain is not in it, compaction alone thins the column


@pytest.mark.conservation
def test_compaction_layers_add_up_to_the_compaction_in_dip(short_config):
    # Written every step: the layers add up to minus DIP column 4, also with horizontal strain on.
    # (The first written row is zero by design; steps with melt are left out, see outputs.rst.)
    run, model = _run(short_config, horizontal_divergence=True, strain_rate_const=-2.0e-2)
    with h5py.File(run.results_path, 'r') as f:
        compaction = f['compaction'][:]
        dip = f['DIP'][:]
        melt = f['meltvol'][:, 1]
    layer_sum = np.nansum(compaction[:, 1:], axis=1)
    rows = (melt == 0.0)
    rows[:2] = False
    assert rows.sum() > 5
    np.testing.assert_allclose(layer_sum[rows], -dip[rows, 4], rtol=1e-3, atol=2e-5)
