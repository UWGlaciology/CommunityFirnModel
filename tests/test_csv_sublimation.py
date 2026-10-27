'''
Sublimation input with input_type "csv" and MELT true (UWGlaciology/CommunityFirnModel#19).
'''

import json
import shutil

import numpy as np
import pytest

from cfm_testing import CFM_MAIN, example_config


def _run(tmp_path, input_folder=None, **overrides):
    '''Run example.json (csv input) in tmp_path with the overrides; return the model.'''
    from firn_density_nospin import FirnDensityNoSpin
    config = json.load(open(CFM_MAIN / 'example.json'))
    config.update(resultsFolder=str(tmp_path / 'results'), NewSpin=True)
    if input_folder is not None:
        config['InputFileFolder'] = str(input_folder)
    config.update(overrides)
    for key in [k for k, v in overrides.items() if v is None]:
        del config[key]
    config_path = tmp_path / 'config.json'
    config_path.write_text(json.dumps(config))
    model = FirnDensityNoSpin(str(config_path), NewSpin=True)
    model.time_evolve()
    return model


@pytest.mark.conservation
@pytest.mark.parametrize('overrides', [
    {'SUBLIM': False},
    {'SUBLIM': True},
    {'SUBLIM': True, 'InputFileNameSublim': None},
], ids=['sublim_off', 'sublim_from_file', 'sublim_from_negative_bdot'])
def test_csv_melt_runs_with_each_sublimation_source(tmp_path, overrides):
    # SUBLIM false used to fail at sublwc2check[iii] = self.subLWCvol (a 1-element
    # array with NumPy 2.5). SUBLIM true without InputFileNameSublim used to fail
    # on climateTS['time'], and climateTS is None for csv input.
    _run(tmp_path, MELT=True, **overrides)
    assert (tmp_path / 'results' / 'CFMresults.hdf5').exists()


@pytest.mark.conservation
def test_csv_sublimation_is_taken_from_negative_bdot(tmp_path):
    # Negative accumulation is sublimation (an outflux); bdot itself is clipped at zero.
    folder = tmp_path / 'input'
    shutil.copytree(CFM_MAIN / example_config()['InputFileFolder'], folder)
    bdot_file = folder / example_config()['InputFileNamebdot']
    data = np.loadtxt(bdot_file, delimiter=',')  # row 0: time, row 1: bdot
    n = data.shape[1]
    data[1, n - 40:n - 30] = -0.05
    np.savetxt(bdot_file, data, delimiter=',')
    model = _run(tmp_path, input_folder=folder, MELT=True, SUBLIM=True, InputFileNameSublim=None)
    sublim = model.forcing_dict['SUBLIM']
    bdot = model.forcing_dict['BDOT']
    assert sublim.min() == pytest.approx(-0.05)
    assert np.count_nonzero(sublim) == 10
    assert bdot.min() >= 0.0
