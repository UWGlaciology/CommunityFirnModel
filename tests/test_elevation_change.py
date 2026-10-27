'''
Surface elevation change in a column that is in equilibrium (UWGlaciology/CommunityFirnModel#10).

The forcing of a few years is repeated for 500 years. When the model is in equilibrium the
elevation change over a repeat has to be zero: what is added at the surface (new snow), lost
by melt and sublimation, and lost by densification has to equal what leaves at the bottom.
If melt were counted twice, the drift would be the thickness of the melt.
'''

import json

import h5py
import numpy as np
import pytest

YEARS = 3         # length of the repeated forcing
REPEATS = 170     # 170 * 3 = 510 years: the column is renewed about once at 0.2 m ie/yr to 100 m


def _elevation_change_per_repeat(short_config, forcing):
    run = short_config(years=YEARS, num_reps=REPEATS, liquid='bucket', DFfile=forcing, TWriteInt=1,
                       outputs=['density', 'depth', 'temperature', 'age', 'DIP', 'LWC'])
    # write every step from the start of the repeats, and do not restart from a spin-up update
    run.config['TWriteStart'] = float(run.climateTS['time'][0])
    run.config['spinUpdate'] = False
    with open(run.config_path, 'w') as f:
        json.dump(run.config, f)
    run.run()
    with h5py.File(run.results_path, 'r') as f:
        dip = f['DIP'][:]
    time, dH_uncorrected, dH_corrected = dip[:, 0], dip[:, 2], dip[:, 5]
    first_of_last_repeat = time[0] + YEARS * np.floor((time[-1] - time[0]) / YEARS)
    last = time >= first_of_last_repeat   # the last repeat of the forcing
    return dH_uncorrected[last].sum(), dH_corrected[last].sum(), dip[0, 0], time[-1]


@pytest.mark.conservation
@pytest.mark.parametrize('forcing, accumulation', [
    ('CFM_example_input_df.pkl', 0.20),       # Summit: little melt
    ('CFM_example_66.5_-46.25.pkl', 0.53),    # DYE-2: melt and refreezing every summer
], ids=['summit', 'dye2'])
def test_elevation_change_is_zero_in_equilibrium(short_config, forcing, accumulation):
    dH_uncorrected, dH_corrected, t0, t1 = _elevation_change_per_repeat(short_config, forcing)
    assert t1 - t0 > 500.0
    accumulated = accumulation * YEARS   # m ice eq. in a repeat
    # The corrected elevation change takes the density of the layer that leaves at the bottom into account.
    assert abs(dH_corrected) < 0.01 * accumulated
    # The uncorrected one assumes that the domain reaches ice density; it is biased by the bottom density.
    assert abs(dH_uncorrected) < 0.15 * accumulated
