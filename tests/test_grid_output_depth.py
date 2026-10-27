'''
Depth of the gridded outputs (grid_output_max_depth). There were two defaults (200 m in ModelOutputs,
which could not be reached, and the depth of the initial domain in the model); a smaller value than the
domain was silently raised; and a domain deeper than the grid was cut off without a warning.
'''

import numpy as np
import pytest

from ModelOutputs import ModelOutputs


def _outputs(**config):
    c = dict(grid_outputs=True, grid_output_res=0.25, output_bits='float32')
    c.update(config)
    z = np.linspace(0.0, 30.0, 61)           # a 30 m domain, 0.5 m nodes
    MOd = dict(z=z, density=np.linspace(350.0, 900.0, 61))
    return ModelOutputs(c, MOd, 3, 1980.0, len(z)), MOd


@pytest.mark.unit
def test_default_depth_is_the_depth_of_the_domain(capsys):
    outputs, MOd = _outputs()
    assert outputs.c['grid_output_max_depth'] == 30.0
    assert outputs.max_grid_len == 121   # 0, 0.25, ..., 30 m
    assert 'Defaulting to the depth of the model domain, 30.0 m' in capsys.readouterr().out
    row = outputs.Mout_dict['density'][0]
    assert np.isfinite(row[1:]).sum() == 121   # the whole domain is on the grid


@pytest.mark.unit
def test_a_domain_deeper_than_the_grid_is_cut_off_with_one_warning(capsys):
    outputs, MOd = _outputs(grid_output_max_depth=10.0)
    outputs.updateMO(MOd, 1980.1, 1)
    outputs.updateMO(MOd, 1980.2, 2)
    printed = capsys.readouterr().out
    assert printed.count('is deeper than grid_output_max_depth') == 1
    assert np.isfinite(outputs.Mout_dict['density'][1, 1:]).sum() == 41   # 0, 0.25, ..., 10 m


@pytest.mark.conservation
def test_model_uses_a_smaller_depth_as_given(short_config, capsys):
    run = short_config(years=2, liquid='bucket', grid_outputs=True, grid_output_max_depth=5.0)
    model = run.build()
    assert model.c['grid_output_max_depth'] == 5.0
    assert 'is shallower than the model domain' in capsys.readouterr().out


@pytest.mark.conservation
def test_model_default_depth_is_the_initial_domain(short_config):
    run = short_config(years=2, liquid='bucket', grid_outputs=True)
    assert 'grid_output_max_depth' not in run.config
    model = run.build()
    assert model.c['grid_output_max_depth'] == np.ceil(model.z[-1])
