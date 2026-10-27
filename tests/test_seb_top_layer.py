'''
Thickness of the top layer that takes the surface energy in the SEB module. If the melt it would produce is
more than its own mass, the layer is made thicker. The layer used to be the one that holds the melt of the
first, thin guess, which can be a node thicker than needed; it is now the thinnest layer that is at least as
heavy as its own melt. Melt that is more than the mass of the whole column raised IndexError.
'''

import numpy as np
import pytest

from SEB import SurfaceEnergyBudget


def _seb(SW, LW=300.0, T2m=275.0, method='prescribed', cls=SurfaceEnergyBudget):
    f = lambda v: np.full(2, float(v))
    climate = dict(time=np.arange(2) / 365.25, SW_d=f(SW), LW_d=f(LW), ALBEDO=f(0.7), T2m=f(T2m), TSKIN=f(T2m),
                   QH=f(30.0), QL=f(5.0), RAIN=f(0.0), WS=f(5.0), QV2m=f(3e-3), PS=f(8e4))
    return cls(dict(SEB_TL_thick=0.05, SEB_turb_method=method), climate, 0, None, 0)


def _column(T, n=300, dz0=0.02, gradient=3.0):
    dz = np.full(n, dz0)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    return dict(Tz=T - np.linspace(0.0, gradient, n), mass=dz * np.linspace(330.0, 450.0, n), dt=86400.0,
                z=z, dz=dz, mtime=0.0)


class Probe(SurfaceEnergyBudget):
    '''Keeps the function that gives the melt for each top layer, to check the layer that was chosen.'''
    def _thicken_top_layer(self, mass, iTL, melt_for):
        self.first = iTL
        self.cum_mass = np.cumsum(mass)
        self.melt_for = melt_for
        self.chosen = super()._thicken_top_layer(mass, iTL, melt_for)
        return self.chosen


@pytest.mark.unit
@pytest.mark.parametrize('seed', range(30))
def test_thickening_finds_the_thinnest_layer_that_can_supply_its_melt(seed):
    rng = np.random.default_rng(seed)
    n = 60
    mass = rng.uniform(5.0, 15.0, n)
    cum = np.cumsum(mass)
    cold = np.cumsum(rng.uniform(0.0, 4.0, n))                 # cold content of the layer, grows with thickness
    energy = rng.uniform(0.0, 300.0)
    melt = lambda i: max(energy - cold[i], 0.0)                 # melt gets smaller with a thicker layer
    first = int(rng.integers(0, 5))
    chosen = _seb(0.0)._thicken_top_layer(mass, first, melt)
    consistent = [i for i in range(first, n) if melt(i) <= cum[i]]
    expected = consistent[0] if consistent else n - 1
    assert chosen == expected


@pytest.mark.unit
def test_a_first_layer_that_is_heavy_enough_is_kept():
    mass = np.full(10, 10.0)
    assert _seb(0.0)._thicken_top_layer(mass, 2, lambda i: 5.0) == 2


@pytest.mark.unit
@pytest.mark.parametrize('method', ['prescribed', 'bulk_richardson'])
@pytest.mark.parametrize('SW', [600.0, 750.0, 900.0])
def test_solver_uses_the_thinnest_consistent_layer(method, SW):
    seb = _seb(SW, method=method, cls=Probe)
    column = _column(268.0)
    solve = seb._solve_Ts_quartic if method == 'prescribed' else seb._solve_Ts_newton
    _, Tz, melt, _ = solve(column, 0, 268.0)
    chosen = seb.chosen
    assert chosen > seb.first                                  # these forcings melt more than the first layer holds
    assert seb.melt_for(chosen) <= seb.cum_mass[chosen]
    assert seb.melt_for(chosen - 1) > seb.cum_mass[chosen - 1]  # a thinner layer could not supply its melt
    assert float(np.squeeze(melt)) == pytest.approx(seb.melt_for(chosen))
    assert np.all(Tz[:chosen + 1] == Tz[0]) and Tz[chosen + 1] != Tz[0]   # the layer takes the surface temperature


@pytest.mark.unit
def test_more_melt_than_the_mass_of_the_column_gives_the_whole_column():
    # The old iteration looked for the node that holds the melt and raised IndexError here.
    mass = np.full(10, 10.0)
    assert _seb(0.0)._thicken_top_layer(mass, 1, lambda i: 500.0) == 9
