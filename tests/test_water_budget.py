'''
Water budget over a short model run, from the totals that time_evolve prints
at the end of the run.
'''

import re

import pytest

DYE2_FORCING = 'CFM_example_66.5_-46.25.pkl'  # DYE-2, Greenland: more melt than the default Summit forcing


def _printed_totals(text):
    '''Return the end-of-run totals [m w.e.] that time_evolve prints.'''
    labels = {'melt_rain': 'Melt+Rain', 'meltvol': 'meltvol', 'rainvol': 'rainvol',
              'refreeze': 'Refreezing', 'runoff': 'Runoff', 'lwc': 'LWC (current)'}
    totals = {}
    for key, label in labels.items():
        match = re.search(rf'^{re.escape(label)}:\s+(\S+)$', text, re.MULTILINE)
        totals[key] = float(match.group(1))
    return totals


@pytest.mark.conservation
@pytest.mark.parametrize('liquid', ['bucket', 'darcy', 'prefsnowpack', 'resingledomain'])
def test_water_budget_closes(short_config, capsys, liquid):
    # The prefsnowpack and resingledomain branches did not set the step totals
    # refreeze, runoff, meltvol and rainvol, so the end-of-run totals were wrong.
    run = short_config(years=2, liquid=liquid, DFfile=DYE2_FORCING)
    run.run()
    totals = _printed_totals(capsys.readouterr().out)
    melt_rain = totals['melt_rain']
    assert totals['meltvol'] + totals['rainvol'] == pytest.approx(melt_rain, rel=1e-6)
    water_out_and_stored = totals['refreeze'] + totals['runoff'] + totals['lwc']
    assert water_out_and_stored == pytest.approx(melt_rain, rel=1e-3)
