'''
A run is repeated exactly by its configuration file and forcing. The random numbers
(surface-density noise, synthetic isotope signal) come from the key random_seed.
'''

import h5py
import numpy as np
import pytest

from config_loader import get_rng

DYE2 = 'CFM_example_66.5_-46.25.pkl'


@pytest.mark.unit
def test_get_rng_repeats_with_the_same_seed():
    a = get_rng({'random_seed': 7}, 'srho').normal(size=5)
    b = get_rng({'random_seed': 7}, 'srho').normal(size=5)
    assert np.array_equal(a, b)


@pytest.mark.unit
def test_get_rng_differs_with_seed_and_with_stream():
    base = get_rng({'random_seed': 7}, 'srho').normal(size=5)
    assert not np.array_equal(base, get_rng({'random_seed': 8}, 'srho').normal(size=5))
    assert not np.array_equal(base, get_rng({'random_seed': 7}, 'dD').normal(size=5))


@pytest.mark.unit
def test_get_rng_default_seed_is_zero():
    assert np.array_equal(get_rng({}, 'srho').normal(size=5),
                          get_rng({'random_seed': 0}, 'srho').normal(size=5))


def _results(path):
    with h5py.File(path, 'r') as f:
        return {k: f[k][()] for k in f.keys()}


def _run_twice(short_config, **overrides):
    out = []
    for _ in range(2):
        run = short_config(years=2, liquid='bucket', DFfile=DYE2, SEB=True, MELT=True, **overrides)
        run.run()
        out.append(_results(run.results_path))
    return out


@pytest.mark.conservation
def test_seb_run_repeats_exactly(short_config):
    a, b = _run_twice(short_config)
    for k in a:
        assert np.array_equal(a[k], b[k], equal_nan=True), k


@pytest.mark.conservation
def test_noisy_surface_density_repeats_with_the_seed(short_config):
    noise = dict(variable_srho=True, srho_type='noise')
    a, b = _run_twice(short_config, random_seed=3, **noise)
    for k in a:
        assert np.array_equal(a[k], b[k], equal_nan=True), k

    # another seed gives another realisation
    run = short_config(years=2, liquid='bucket', DFfile=DYE2, SEB=True, MELT=True, random_seed=4, **noise)
    run.run()
    assert not np.array_equal(a['density'], _results(run.results_path)['density'], equal_nan=True)
