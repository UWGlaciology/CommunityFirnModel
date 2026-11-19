'''
The forcing file (CFMforcing.hdf5) is enough, with the configuration file, to repeat a run. The
repeats of the reference climate interval (RCI) are not stored; reader.read_forcing puts them back.
'''

import types

import h5py
import numpy as np
import pandas as pd
import pytest

import RCMpkl_to_spin as RCM
from cfm_testing import CFM_MAIN
from reader import expand_rci, read_forcing
from writer import forcing_writer

DYE2 = 'CFM_example_66.5_-46.25.pkl'
OLD = 'CFM_example_input_df.pkl'


def _same(a, b):
    assert set(a) == set(b)
    for k in a:
        x, y = np.asarray(a[k]), np.asarray(b[k])
        assert x.dtype == y.dtype and x.shape == y.shape, k
        assert np.array_equal(x, y, equal_nan=True), k


def _forcing(pickle, years, **kw):
    df = pd.read_pickle(CFM_MAIN / 'CFMinput_example' / pickle).loc['1980':f'{1979 + years}'].copy()
    clim, _, _, _, _, seb = RCM.makeSpinFiles(df, melt=True, spin_date_st=1980.0,
                                              spin_date_end=1979.0 + years, num_reps=4, **kw)
    return clim, seb


def _round_trip(tmp_path, clim, seb, strain=None):
    writer_self = types.SimpleNamespace(c={'resultsFolder': str(tmp_path)})
    forcing_writer(writer_self, dict(clim), seb, strain)
    return read_forcing(tmp_path)


def _without_rci_keys(clim):
    return {k: v for k, v in clim.items() if k not in ('sds', 'sde', 'num_reps')}


@pytest.mark.unit
def test_expand_rci_puts_the_repeats_in_front():
    record = {'time': np.arange(1980.0, 1983.0, 0.5), 'x': np.arange(6.0), 'start': 1980.0}
    out = expand_rci(record, 1980.0, 1980.0, 2)   # the RCI is the year 1980 (2 steps), twice
    assert list(out['time'][:4]) == [1978.0, 1978.5, 1979.0, 1979.5]
    assert list(out['x']) == [0, 1, 0, 1] + list(range(6))
    assert out['start'] == 1980.0


@pytest.mark.unit
@pytest.mark.parametrize('kw, pickle, resolution', [
    (dict(SEB=True), DYE2, '1d'),
    (dict(SEB=False), OLD, '1d'),
    (dict(SEB=False), OLD, '1ME'),
])
def test_forcing_file_gives_back_the_forcing(tmp_path, kw, pickle, resolution):
    # SEB: the repeats are exactly the record. Without SEB the first step of the mass fluxes has the
    # mean step length (flux_to_rate), so those elements are stored as a patch.
    clim, seb = _forcing(pickle, 2, timeres=resolution, **kw)
    climateTS, SEBfluxes, strain = _round_trip(tmp_path, clim, seb)
    _same(_without_rci_keys(clim), climateTS)
    if seb is None:
        assert SEBfluxes is None
    else:
        _same(seb, SEBfluxes)
    assert strain is None
    with h5py.File(tmp_path / 'CFMforcing.hdf5', 'r') as f:
        assert f['main'].attrs['rci_trimmed']
        assert len(f['main/time']) < len(clim['time']) / 2   # the repeats are not in the file


@pytest.mark.unit
def test_forcing_file_stores_everything_when_the_repeats_cannot_be_removed(tmp_path):
    # sds is not the start of the record: the repeats are not copies of the record
    df = pd.read_pickle(CFM_MAIN / 'CFMinput_example' / OLD).loc['1980':'1983'].copy()
    clim, *_ = RCM.makeSpinFiles(df, melt=True, timeres='1d', spin_date_st=1982.0,
                                 spin_date_end=1982.0, num_reps=3)
    climateTS, _, _ = _round_trip(tmp_path, clim, None)
    _same(_without_rci_keys(clim), climateTS)
    with h5py.File(tmp_path / 'CFMforcing.hdf5', 'r') as f:
        assert not f['main'].attrs['rci_trimmed']
        assert len(f['main/time']) == len(clim['time'])


@pytest.mark.unit
def test_forcing_file_holds_the_strain_rates_argument(tmp_path):
    clim, seb = _forcing(OLD, 2, timeres='1d', SEB=False)
    strain = {'time': np.array([1900.0, 2000.0]), 'eps_divergence': np.array([1e-4, 2e-4])}
    _, _, got = _round_trip(tmp_path, clim, seb, strain)
    _same(strain, got)
    _, _, got = _round_trip(tmp_path, clim, seb, [1e-4, 2e-4])
    assert np.array_equal(got, [1e-4, 2e-4])


@pytest.mark.conservation
def test_run_from_the_forcing_file_repeats_the_run(short_config):
    # The model run from the forcing file and the same configuration gives the same results as the run
    # that wrote the file: every result dataset is identical.
    kw = dict(years=2, liquid='bucket', DFfile=DYE2, SEB=True, MELT=True)
    first = short_config(**kw)
    first.run()

    second = short_config(**kw)
    second.climateTS, second.SEBfluxes, _ = read_forcing(first.folder)
    second.run()

    with h5py.File(first.results_path, 'r') as a, h5py.File(second.results_path, 'r') as b:
        assert set(a.keys()) == set(b.keys())
        for k in a.keys():
            assert np.array_equal(a[k][()], b[k][()], equal_nan=True), k
