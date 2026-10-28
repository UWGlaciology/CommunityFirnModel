'''
Unit tests for melt.leaky_lens_routing (partially permeable ice lenses).
'''

import numpy as np
import pytest

from melt import leaky_lens_routing


def column(dense_nodes, nnd=10, dz=0.1, stcap=0.01):
    '''A column of nnd nodes with the given dense node indices; dense nodes and
    the bottom node have no storage capacity (as in bucket).'''
    dense = np.zeros(nnd, dtype=bool)
    dense[list(dense_nodes)] = True
    dz = np.full(nnd, dz)
    cap = np.where(dense, 0.0, stcap)
    cap[-1] = 0.0
    return dense, dz, cap


def blocked_above(node, amount, nnd=10):
    LWCblocked = np.zeros(nnd)
    LWCblocked[node] = amount
    return LWCblocked


def test_water_is_conserved():
    rng = np.random.default_rng(1)
    nnd = 50
    dense = rng.random(nnd) < 0.3
    dz = rng.uniform(0.01, 0.2, nnd)
    stcap = np.where(dense, 0.0, rng.uniform(0.0, 0.02, nnd))
    LWCblocked = np.where(rng.random(nnd) < 0.2, rng.uniform(0.0, 0.05, nnd), 0.0)
    for decay_length in [0.0, 0.05, 0.3, 10.0]:
        storage, blocked = leaky_lens_routing(LWCblocked, stcap, dense, dz, decay_length)
        assert storage.sum() + blocked.sum() == pytest.approx(LWCblocked.sum(), rel=1e-12)
        assert np.all(storage <= stcap + 1e-15)
        assert np.all(storage >= 0) and np.all(blocked >= 0)


def test_zero_decay_length_blocks_everything():
    # as ThickImp = 0: water blocked above a lens stays where it is
    dense, dz, stcap = column([4, 5])
    LWCblocked = blocked_above(3, 0.02)
    storage, blocked = leaky_lens_routing(LWCblocked, stcap, dense, dz, 0.0)
    assert np.all(storage == 0)
    np.testing.assert_array_equal(blocked, LWCblocked)


def test_transmitted_fraction_fills_storage_below():
    # one dense node of 0.1 m with decay length 0.1 m passes exp(-1) of the water
    dense, dz, stcap = column([4], stcap=1.0)
    LWCblocked = blocked_above(3, 0.02)
    storage, blocked = leaky_lens_routing(LWCblocked, stcap, dense, dz, 0.1)
    assert storage[5] == pytest.approx(0.02 * np.exp(-1.0))
    assert blocked[3] == pytest.approx(0.02 * (1 - np.exp(-1.0)))
    assert storage[:5].sum() == 0


def test_large_decay_length_passes_all_lenses():
    dense, dz, stcap = column([2, 5, 6], stcap=0.004)
    LWCblocked = blocked_above(1, 0.02)
    storage, blocked = leaky_lens_routing(LWCblocked, stcap, dense, dz, 1e9)
    # everything below the first lens fills up; only what does not fit is left at the bottom
    np.testing.assert_allclose(storage[[3, 4, 7, 8]], 0.004)
    assert blocked[8] == pytest.approx(0.02 - 4 * 0.004)
    assert blocked[:8].sum() == pytest.approx(0.0, abs=1e-9)


def test_lens_thickness_adds_up():
    # two adjacent 0.05 m dense nodes transmit the same as one 0.1 m node
    nnd = 10
    dense_two = np.zeros(nnd, dtype=bool)
    dense_two[[4, 5]] = True
    dz_two = np.full(nnd, 0.1)
    dz_two[[4, 5]] = 0.05
    dense_one, dz_one, _ = column([4])
    stcap_two = np.where(dense_two, 0.0, 1.0)
    stcap_one = np.where(dense_one, 0.0, 1.0)
    stcap_two[-1] = stcap_one[-1] = 0.0
    LWCblocked = blocked_above(3, 0.02)
    s_two, b_two = leaky_lens_routing(LWCblocked, stcap_two, dense_two, dz_two, 0.07)
    s_one, b_one = leaky_lens_routing(LWCblocked, stcap_one, dense_one, dz_one, 0.07)
    assert s_two.sum() == pytest.approx(s_one.sum())
    assert b_two[3] == pytest.approx(b_one[3])


def test_water_inside_a_lens_only_passes_the_rest_of_the_lens():
    dense, dz, stcap = column([4, 5], stcap=1.0)
    LWCblocked = blocked_above(5, 0.02) # excess water in the lower lens node
    storage, blocked = leaky_lens_routing(LWCblocked, stcap, dense, dz, 0.1)
    assert storage[6] == pytest.approx(0.02 * np.exp(-1.0))


def test_nothing_blocked_returns_zeros():
    dense, dz, stcap = column([4])
    storage, blocked = leaky_lens_routing(np.zeros(10), stcap, dense, dz, 0.1)
    assert np.all(storage == 0) and np.all(blocked == 0)
