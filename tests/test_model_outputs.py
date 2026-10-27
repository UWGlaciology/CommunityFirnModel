'''
Regridding of the outputs that are an amount per layer (UWGlaciology/CommunityFirnModel#16):
compaction [m], LWC and refrozen [m w.e.]. They were interpolated like a profile of
temperature, so the sum over the cells was not the sum over the layers.
'''

import h5py
import numpy as np
import pytest

from ModelOutputs import ModelOutputs, EXTENSIVE_OUTPUTS


def _layers(n=40):
    '''Layers that get thicker with depth, as in the model: top depths z and amounts.'''
    dz = np.linspace(0.05, 0.6, n)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    amount = np.random.default_rng(1).uniform(0.0, 1.0, n) * dz
    return z, dz, amount


@pytest.mark.unit
@pytest.mark.parametrize('spacing', [0.1, 0.25, 1.0])
def test_regrid_keeps_the_total(spacing):
    z, dz, amount = _layers()
    n_pts = int(np.floor((z[-1] - z[0]) / spacing)) + 1
    grid = z[0] + np.arange(n_pts) * spacing
    cells = ModelOutputs.RGfun(None, z, amount, grid)
    assert cells.sum() == pytest.approx(amount.sum(), rel=1e-12)
    assert np.all(cells >= 0.0)


@pytest.mark.unit
def test_regrid_puts_the_amount_in_the_depth_range_of_the_layer():
    # One layer holds everything: all the amount has to land in the cells that overlap that layer.
    z, dz, _ = _layers()
    amount = np.zeros_like(z)
    k = 20
    amount[k] = 1.0
    grid = np.arange(0.0, z[-1], 0.05)
    cells = ModelOutputs.RGfun(None, z, amount, grid)
    top, bottom = z[k], z[k] + dz[k]
    inside = (grid + 0.05 > top) & (grid < bottom)
    assert cells[inside].sum() == pytest.approx(1.0, rel=1e-12)
    assert cells[~inside].sum() == pytest.approx(0.0, abs=1e-12)


@pytest.mark.unit
def test_regrid_sum_between_two_depths_is_the_amount_between_them():
    # Layer thicknesses that are multiples of the cell size, so the layer edges are cell edges.
    dz = 0.05 * np.arange(1, 41)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    amount = np.random.default_rng(2).uniform(0.0, 1.0, 40) * dz
    grid = np.arange(0.0, z[-1], 0.05)
    cells = ModelOutputs.RGfun(None, z, amount, grid)
    # Between layer edges 10 and 25 the cells add up to those layers
    in_range = (grid >= z[10] - 1e-9) & (grid < z[25] - 1e-9)
    assert cells[in_range].sum() == pytest.approx(amount[10:25].sum(), rel=1e-9)


@pytest.mark.unit
def test_amounts_per_layer_are_regridded_as_amounts():
    assert {'LWC', 'refrozen', 'compaction'} <= set(EXTENSIVE_OUTPUTS)


@pytest.mark.conservation
def test_gridded_compaction_adds_up_to_the_native_compaction(short_config):
    # The same run on the model grid and on the output grid: the layers and the cells add up to the same.
    sums = {}
    for grid_outputs in (False, True):
        run = short_config(years=2, liquid='bucket', DFfile='CFM_example_66.5_-46.25.pkl', TWriteInt=1,
                           grid_outputs=grid_outputs, outputs=['density', 'depth', 'temperature', 'age', 'compaction', 'LWC'])
        run.run()
        with h5py.File(run.results_path, 'r') as f:
            for name in ('compaction', 'LWC'):
                sums[(grid_outputs, name)] = np.nansum(f[name][:, 1:], axis=1)
    assert np.abs(sums[(False, 'compaction')]).max() > 1e-3
    np.testing.assert_allclose(sums[(True, 'compaction')], sums[(False, 'compaction')], rtol=1e-4, atol=1e-7)
    np.testing.assert_allclose(sums[(True, 'LWC')], sums[(False, 'LWC')], rtol=1e-4, atol=1e-9)
