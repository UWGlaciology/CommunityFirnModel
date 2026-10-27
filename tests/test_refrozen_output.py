'''
The per-layer 'refrozen' output (UWGlaciology/CommunityFirnModel#20). It was set
once in __init__ and never updated by the bucket and darcy schemes, so the
output was all zeros.
'''

import h5py
import numpy as np
import pytest

import melt


@pytest.mark.scheme
@pytest.mark.parametrize('scheme', ['bucket', 'darcyscheme'])
def test_scheme_sets_refrozen_per_layer(synthetic_column, scheme):
    # Melt on a cold column: the water refreezes in the layers below the surface.
    column = synthetic_column(n=60, rho_profile=400.0, T_profile=253.15, melt=0.02)
    outputs = getattr(melt, scheme)(column, 0)
    refrozen_total = outputs[10]  # refrozentot (bucket) or refr_tot (darcyscheme) [m we]
    assert refrozen_total > 0.0
    assert column.refrozen.shape == column.dz.shape
    assert np.all(column.refrozen >= 0.0)
    assert column.refrozen.sum() == pytest.approx(refrozen_total, rel=1e-9)


@pytest.mark.conservation
@pytest.mark.parametrize('liquid', ['bucket', 'darcy'])
@pytest.mark.parametrize('grid_outputs', [False, True], ids=['native_grid', 'gridded'])
def test_refrozen_output_matches_refreeze(short_config, liquid, grid_outputs):
    # Every step is written. The refrozen layers must add up to the refreeze
    # output of the same step: the melt scheme and the refreezing diffusion
    # both count in the two. Gridded output regrids the amounts and keeps the total.
    run = short_config(years=2, liquid=liquid, DFfile='CFM_example_66.5_-46.25.pkl', TWriteInt=1,
                       grid_outputs=grid_outputs,
                       outputs=['density', 'depth', 'temperature', 'age', 'LWC', 'refrozen'])
    run.run()
    with h5py.File(run.results_path, 'r') as f:
        refrozen = f['refrozen'][:]    # column 0 is the time
        refreeze = f['refreeze'][:]    # [time, total refreezing of the step in m we]
    layer_sum = np.nansum(refrozen[:, 1:], axis=1)
    assert refreeze[:, 1].sum() > 1e-4  # the run has melt
    assert layer_sum.sum() > 0.0
    if liquid == 'bucket':
        np.testing.assert_allclose(layer_sum, refreeze[:, 1], rtol=1e-3, atol=1e-7)
    else:
        # darcy melts in sub-steps, so water that refroze in a surface layer can melt again in a later
        # sub-step of the same step. The refreeze output counts it, the layers do not (it is no longer there).
        assert np.all(layer_sum <= refreeze[:, 1] * (1 + 1e-3) + 1e-7)
        assert layer_sum.sum() >= 0.95 * refreeze[:, 1].sum()
