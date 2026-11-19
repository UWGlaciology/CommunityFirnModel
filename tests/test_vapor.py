'''
Vapor diffusion in the pore space (vapor.py, configuration key vaporTransport).

Vapor transport changes the mass of a layer at fixed volume: the column mass is conserved
(zero flux at the surface and the bottom), and the thickness of a layer does not change.
'''

import types

import numpy as np
import pytest

import firn_density_nospin
import vapor
from vapor import (block_ids, vapor_latent_conductivity, firn_tortuosity_inverse, saturation_vapor_density, saturation_vapor_pressure, vapor_diffusivity,
                   vapor_mass_flux_divergence)


def _column(rho, T, dz=None, c=None):
    '''A column in the form that vapor_mass_flux_divergence reads (z = depth of the top of each node).'''
    n = len(rho)
    dz = np.full(n, 0.05) if dz is None else np.asarray(dz, dtype=float)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    rho = np.asarray(rho, dtype=float)
    return types.SimpleNamespace(z=z, dz=dz, Tz=np.asarray(T, dtype=float), rho=rho, mass=rho * dz, c=c or {})


DAY = 86400.0


def _flux(column, dt=DAY, **kw):
    return vapor_mass_flux_divergence(column, dt, **kw)[0]


@pytest.mark.unit
def test_saturation_vapor_pressure_of_air_is_the_one_of_johnsen():
    T = np.array([223.15, 243.15, 263.15, 273.15])
    assert saturation_vapor_pressure(T) == pytest.approx(3.454e12 * np.exp(-6133 / T), rel=1e-14)   # the formula of Johnsen et al. (2000)
    assert saturation_vapor_pressure(273.15) == pytest.approx(611.0, rel=0.005)
    assert np.all(np.diff(saturation_vapor_density(T)) > 0)
    assert saturation_vapor_density(253.15) == pytest.approx(saturation_vapor_pressure(253.15) * 0.018 / (8.314 * 253.15))


@pytest.mark.unit
def test_saturation_vapor_pressure_of_co2_is_the_one_of_arthern():
    assert saturation_vapor_pressure(273.16, 'CO2') == pytest.approx(611.0)                  # AW00 Eq. 6: 611 Pa at the triple point
    assert saturation_vapor_density(273.16, 'CO2') == pytest.approx(4.84e-3, rel=2e-3)
    T = np.linspace(223.15, 273.15, 11)
    assert np.allclose(saturation_vapor_pressure(T, 'CO2'), saturation_vapor_pressure(T, 'air'), rtol=0.0036)   # the two agree
    assert not np.array_equal(saturation_vapor_pressure(T, 'CO2'), saturation_vapor_pressure(T, 'air'))
    with pytest.raises(ValueError):
        saturation_vapor_pressure(250.0, 'argon')


@pytest.mark.unit
def test_firn_tortuosity():
    rho = np.array([350.0, 600.0, 700.0, 804.0, 805.0, 850.0, 917.0])
    inv = firn_tortuosity_inverse(rho)
    assert inv[1] == pytest.approx(1 - 1.3 * (600 / 917) ** 2)
    assert np.all(inv[:4] > 0) and np.all(np.diff(inv[:4]) < 0)
    assert np.all(inv[4:] == 0)                                                           # zero from 804.26 kg m-3
    assert np.all(firn_tortuosity_inverse(rho, b=2.0)[3:] == 0)                           # a larger b closes earlier (648 kg m-3)
    assert firn_tortuosity_inverse(700.0, b=1.0) > firn_tortuosity_inverse(700.0, b=1.3)


@pytest.mark.unit
def test_vapor_diffusivity_in_air():
    D = lambda T, P, **k: vapor_diffusivity(T, P, atmosphere='air', **k)
    assert D(273.15, 101325.0) == pytest.approx(2.11e-5)                                   # the free-air value
    assert D(253.15, 101325.0) == pytest.approx(2.11e-5 * (253.15 / 273.15) ** 1.94)
    assert D(273.15, 101325.0) / D(273.15, 50662.5) == pytest.approx(0.5)                   # 1/P
    assert D(260.0, 101325.0, yosida_factor=4.0) == pytest.approx(4 * D(260.0, 101325.0))
    assert vapor_diffusivity(260.0, 101325.0) == D(260.0, 101325.0)                          # air is the default


@pytest.mark.unit
def test_vapor_diffusivity_in_co2_is_that_of_clow():
    D = lambda T, P, **k: vapor_diffusivity(T, P, atmosphere='CO2', **k)
    assert D(273.0, 101300.0) == pytest.approx(4 * 1.37e-5)                                 # the Yosida factor of Clow (1987)
    assert D(273.0, 101300.0, yosida_factor=1.0) == pytest.approx(1.37e-5)
    assert D(2 * 273.0, 101300.0) / D(273.0, 101300.0) == pytest.approx(4.0)                # T^2
    with pytest.raises(ValueError):
        vapor_diffusivity(273.0, 101300.0, atmosphere='argon')


@pytest.mark.unit
def test_vapor_configuration_selects_the_diffusivity():
    T = np.linspace(250.0, 262.0, 20)
    rho = np.full(20, 600.0)
    air = _flux(_column(rho, T))                                                              # the default: air at 1013.25 hPa
    assert np.array_equal(air, _flux(_column(rho, T, c={'vapor_atmosphere': 'air', 'site_pressure': 1013.25})))
    D = lambda c: vapor_diffusivity(260.0, 101325.0, atmosphere=c.get('vapor_atmosphere', 'air'),
                                    yosida_factor=c.get('vapor_yosida_factor'))
    co2 = _flux(_column(rho, T, c={'vapor_atmosphere': 'CO2'}))
    assert np.allclose(co2 / air, D({'vapor_atmosphere': 'CO2'}) / D({}), rtol=0.03)          # the temperature dependences differ a little
    half = _flux(_column(rho, T, c={'site_pressure': 506.625}))
    assert np.allclose(half, 2 * air, rtol=1e-12, atol=0)                                     # half the pressure: twice the diffusivity
    four = _flux(_column(rho, T, c={'vapor_yosida_factor': 4.0}))
    assert np.allclose(four, 4 * air, rtol=1e-12, atol=0)


@pytest.mark.unit
def test_vapor_flux_conserves_mass():
    rng = np.random.default_rng(1)
    n = 80
    dz = rng.uniform(0.005, 0.3, n)                               # layers of different thickness
    rho = np.sort(rng.uniform(330.0, 900.0, n))
    T = 255.0 + 10.0 * rng.standard_normal(n)                     # a noisy temperature profile
    d = _flux(_column(rho, T, dz))
    assert abs(d.sum()) < 1e-12 * np.abs(d).sum()                 # zero flux at the top and the bottom: the sum is zero to rounding
    assert np.abs(d).max() > 0


@pytest.mark.unit
def test_vapor_moves_from_the_warm_to_the_cold_layers():
    n = 20
    T = np.concatenate((np.full(n // 2, 250.0), np.full(n // 2, 258.0)))    # cold above, warm below
    d = _flux(_column(np.full(n, 400.0), T))
    assert d[n // 2 - 1] > 0 and d[n // 2] < 0                    # the cold layer at the step gains, the warm one loses
    assert np.all(d[:n // 2 - 1] == 0) and np.all(d[n // 2 + 1:] == 0)   # uniform T away from the step: no flux
    Tw = T[::-1]                                                  # warm above, cold below: the other way
    dw = _flux(_column(np.full(n, 400.0), Tw))
    assert dw[n // 2 - 1] < 0 and dw[n // 2] > 0


@pytest.mark.unit
def test_vapor_flux_is_zero_in_an_isothermal_column():
    d = _flux(_column(np.linspace(350.0, 800.0, 40), np.full(40, 255.0)))
    assert np.abs(d).max() < 1e-18                                           # (the mean temperature of a block can differ in the last bit)


@pytest.mark.unit
def test_vapor_flux_stops_where_the_tortuosity_does():
    n = 20
    T = np.linspace(250.0, 262.0, n)
    open_col = _flux(_column(np.full(n, 600.0), T))
    shut = _flux(_column(np.full(n, 810.0), T))                   # above 804 kg m-3
    assert np.abs(open_col).max() > 0
    assert np.abs(shut).max() < 1e-6 * np.abs(open_col).max()


@pytest.mark.unit
def test_vapor_tortuosity_option():
    n = 20
    rho, T = np.full(n, 700.0), np.linspace(250.0, 262.0, n)
    default = _flux(_column(rho, T))
    more = _flux(_column(rho, T, c={'vapor_tortuosity_b': 1.0}))
    assert np.abs(more).sum() > np.abs(default).sum()
    ratio = (1 - 1.0 * (700 / 917) ** 2) / (1 - 1.3 * (700 / 917) ** 2)                   # only the pore space changes
    assert np.allclose(more, ratio * default, rtol=1e-12, atol=0)


### In the model ###
####################

def _run(short_config, **kw):
    run = short_config(years=2, num_reps=2, liquid='bucket', MELT=False, **kw)
    return run.run()


@pytest.mark.conservation
def test_vapor_conserves_the_mass_of_the_column(short_config):
    off = _run(short_config)
    on = _run(short_config, vaporTransport=True)
    assert on.mass.sum() == pytest.approx(off.mass.sum(), rel=1e-12)
    n = min(len(on.rho), len(off.rho))
    assert np.abs(on.rho[:n] - off.rho[:n]).max() > 0.1           # and vapor does something


@pytest.mark.conservation
def test_vapor_does_not_change_the_thickness_of_layers(short_config, monkeypatch):
    # Without densification, a change of the density by vapor (at fixed mass it would change dz) must leave the thickness alone.
    monkeypatch.setattr(firn_density_nospin, 'get_drho_dt', lambda RD, rho, *a, **k: np.zeros_like(rho))
    off = _run(short_config)
    on = _run(short_config, vaporTransport=True)
    assert np.abs(on.rho - off.rho).max() > 0.1
    assert np.allclose(on.dz, off.dz, rtol=1e-9, atol=0)
    assert on.dz.sum() == pytest.approx(off.dz.sum(), rel=1e-12)


@pytest.mark.conservation
def test_no_densification_also_freezes_the_vapor_mass(short_config):
    off = _run(short_config, no_densification=True)
    on = _run(short_config, no_densification=True, vaporTransport=True)
    assert np.array_equal(on.rho, off.rho) and np.array_equal(on.mass, off.mass)


@pytest.mark.conservation
@pytest.mark.parametrize('vapor, isotopes, warned', [(True, True, True), (True, False, False), (False, True, False)])
def test_warning_when_vapor_transport_and_isotope_diffusion_are_both_on(short_config, capsys, vapor, isotopes, warned):
    # the isotope diffusion does not include the isotopes that the vapor moves with the mass
    run = short_config(years=2, num_reps=2, liquid='bucket', MELT=False, vaporTransport=vapor, isoDiff=isotopes, iso=['dD'])
    run.build()
    out = capsys.readouterr().out
    assert out.count('ISOTOPE DIFFUSION DOES NOT CURRENTLY INCLUDE VAPOR TRANSPORT') == (1 if warned else 0)


### The limiter of the flux ###
###############################

def _thin_column(c=None):
    '''Layers as at the surface of a run with daily steps: some almost without mass (days without snowfall).'''
    n = 40
    dz = np.concatenate((np.full(10, 3e-5), np.full(n - 10, 0.005)))       # 0.03 mm and 5 mm
    rho = np.full(n, 150.0)
    T = np.linspace(262.0, 268.0, n) + 2.0 * np.sin(np.arange(n))           # a gradient between neighbours
    return _column(rho, T, dz, c=c)


@pytest.mark.unit
def test_vapor_flux_limiter_keeps_the_mass_positive_and_conserved():
    column = _thin_column({'vapor_min_thickness': 0})                         # one layer per block
    unlimited = _flux(column, max_loss=np.inf, max_fill=np.inf)
    assert np.any(DAY * unlimited < -column.mass)                             # without the limiter a layer would lose more than all its mass
    d, n_limited = vapor_mass_flux_divergence(column, DAY)
    assert n_limited > 0
    assert np.all(column.mass + DAY * d >= 0.5 * column.mass - 1e-15)         # no layer loses more than half its mass
    assert abs(d.sum()) < 1e-12 * np.abs(d).sum()                             # and the column mass is still conserved
    d2, _ = vapor_mass_flux_divergence(column, DAY, max_loss=0.1)
    assert np.all(column.mass + DAY * d2 >= 0.9 * column.mass - 1e-15)
    room = 917.0 * column.dz - column.mass                                   # the free pore space
    gain = DAY * d
    assert np.all(gain <= 0.5 * room + 1e-15)                                 # no layer fills more than half of its pores
    assert np.any(DAY * unlimited > 0.5 * room)                               # (and without the limiter one would)


@pytest.mark.unit
def test_vapor_flux_limiter_does_not_touch_normal_layers():
    n = 30
    column = _column(np.linspace(350.0, 700.0, n), np.linspace(250.0, 262.0, n), np.full(n, 0.05))   # layers of 5 cm
    d, n_limited = vapor_mass_flux_divergence(column, DAY)
    assert n_limited == 0
    assert np.array_equal(d, _flux(column, max_loss=np.inf, max_fill=np.inf))


@pytest.mark.unit
def test_vapor_flux_with_layers_without_mass():
    # a day without snowfall, or a layer that the sublimation emptied: mass and dz are zero. No flux through it, no NaN.
    n = 12
    column = _column(np.full(n, 350.0), np.linspace(255.0, 265.0, n), np.full(n, 0.05))
    column.mass[4] = 0.0
    column.dz[4] = 0.0
    d, n_limited = vapor_mass_flux_divergence(column, DAY)
    assert np.all(np.isfinite(d)) and d[4] == 0.0
    assert abs(d.sum()) < 1e-12 * np.abs(d).sum()


### The blocks ###
##################

@pytest.mark.unit
def test_block_ids_group_thin_layers_with_their_neighbours():
    dz = np.array([0.03, 3e-5, 4e-4, 2e-3, 0.02, 0.0, 1e-3, 0.5, 0.0004, 0.6])
    assert list(block_ids(dz, 0.01)) == [0, 1, 1, 1, 1, 2, 2, 2, 3, 3]       # a 3 cm layer is alone; the thin ones are grouped
    assert list(block_ids(dz, 0)) == list(range(len(dz)))                    # 0: one layer per block
    assert list(block_ids(np.full(5, 0.05), 0.01)) == [0, 1, 2, 3, 4]        # thick layers are their own blocks


@pytest.mark.unit
def test_block_ids_of_random_layers():
    rng = np.random.default_rng(3)
    for _ in range(20):
        dz = 10 ** rng.uniform(-5, 0, 200) * (rng.uniform(size=200) > 0.05)  # 5 micrometres to 1 m, and some of zero thickness
        ids = block_ids(dz, 0.01)
        assert ids[0] == 0 and np.all(np.diff(ids) >= 0) and np.all(np.diff(ids) <= 1)   # consecutive blocks from the top
        thickness = np.bincount(ids, weights=dz)
        assert np.all(thickness[:-1] >= 0.005 - 1e-12)                       # no thin block is left, except the last
        assert thickness.sum() == pytest.approx(dz.sum())


@pytest.mark.unit
def test_vapor_does_not_empty_thin_layers_when_they_are_in_blocks():
    # the column that empties the thin layers when each of them is a block (see above)
    alone = _thin_column({'vapor_min_thickness': 0})
    d0, limited0 = vapor_mass_flux_divergence(alone, DAY)
    assert limited0 > 0 and np.max(np.abs(DAY * d0 / alone.mass)) >= 0.5 - 1e-12   # the layers are at the limit
    blocks = _thin_column()                                                         # blocks of at least 1 cm (the default)
    d, limited = vapor_mass_flux_divergence(blocks, DAY)
    assert limited == 0
    change = DAY * d / blocks.mass
    assert np.max(np.abs(change)) < 0.2                                             # nothing is near the limit
    assert np.ptp(change[:10]) < 1e-9                                               # the ten thin layers change by the same fraction
    assert abs(d.sum()) < 1e-12 * np.abs(d).sum()


@pytest.mark.unit
def test_layers_of_a_block_lose_the_same_fraction_of_their_mass():
    dz = np.array([0.02] + [0.002] * 5 + [0.05, 0.05, 0.05])                         # layers 1 to 5 are one block of 1 cm
    T = np.array([262.0] + [263.0] * 5 + [255.0, 250.0, 248.0])
    column = _column(np.linspace(250.0, 400.0, 9), T, dz)
    assert list(block_ids(dz, 0.01)) == [0, 1, 1, 1, 1, 1, 2, 3, 4]
    d, _ = vapor_mass_flux_divergence(column, DAY)
    change = d / column.mass
    assert np.ptp(change[1:6]) < 1e-12 * np.abs(change).max() and change[1] != 0


@pytest.mark.unit
def test_blocks_do_not_change_the_flux_between_thick_layers():
    n = 30
    thick = _column(np.linspace(350.0, 700.0, n), np.linspace(250.0, 262.0, n), np.full(n, 0.05))      # all layers are 5 cm
    by_layer = _flux(_column(thick.rho, thick.Tz, np.full(n, 0.05), c={'vapor_min_thickness': 0}))
    assert np.allclose(_flux(thick), by_layer, rtol=1e-9, atol=0)                    # the default (blocks of 1 cm) is the same


@pytest.mark.unit
def test_vapor_min_thickness_is_read_from_the_configuration():
    thin = _thin_column({'vapor_min_thickness': 0})
    coarse = _thin_column({'vapor_min_thickness': 0.05})
    assert not np.allclose(vapor_mass_flux_divergence(thin, DAY)[0], vapor_mass_flux_divergence(coarse, DAY)[0], rtol=0.1, atol=0)
    assert abs(vapor_mass_flux_divergence(coarse, DAY)[0].sum()) < 1e-12


### The latent heat of the vapor flux (vapor_latent_heat): an added conductivity

def _kv_by_difference(T, rho, c, h=1e-3):
    atmosphere = c.get('vapor_atmosphere', 'air')
    D = vapor_diffusivity(T, 100.0 * c.get('site_pressure', 1013.25), atmosphere=atmosphere, yosida_factor=c.get('vapor_yosida_factor'))
    pore = (1.0 - rho / 917.0) * firn_tortuosity_inverse(rho, c.get('vapor_tortuosity_b', 1.3))
    L = 2.834e6 if atmosphere == 'air' else 2.83e6
    slope = (saturation_vapor_density(T + h, atmosphere) - saturation_vapor_density(T - h, atmosphere)) / (2 * h)
    return L * D * pore * slope


@pytest.mark.unit
@pytest.mark.parametrize('c', [{}, {'vapor_atmosphere': 'CO2'}, {'vapor_yosida_factor': 2.0, 'site_pressure': 700.0}])
def test_latent_conductivity_is_the_latent_heat_of_the_vapor_flux(c):
    T = np.array([200.0, 230.0, 255.0, 270.0])
    rho = np.array([150.0, 350.0, 600.0, 750.0])
    assert np.allclose(vapor_latent_conductivity(T, rho, c), _kv_by_difference(T, rho, c), rtol=1e-4)


@pytest.mark.unit
def test_latent_conductivity_size_and_cutoff():
    # a few percent of the conduction of firn in surface snow near the melting point, nothing above the tortuosity cutoff
    K = vapor_latent_conductivity(np.array([270.0, 270.0, 270.0, 270.0]), np.array([150.0, 450.0, 803.0, 850.0]), {})
    assert 0.01 < K[0] < 0.02 and 0.004 < K[1] < 0.01           # 0.014 and 0.006 W m-1 K-1 at 270 K
    assert K[2] < 1e-4 and K[3] == 0.0
    assert np.all(K >= 0)


@pytest.mark.unit
def test_latent_conductivity_gives_the_heat_of_the_mass_flux():
    # In a column with a smooth temperature gradient the vapor mass flux through an interior face times L is K_v times the gradient.
    n, dz = 40, 0.25
    T = 255.0 + 0.5 * np.arange(n) * dz                            # 0.5 K per metre, warm below
    rho = np.full(n, 350.0)
    col = _column(rho, T, np.full(n, dz), {'vapor_min_thickness': 0.0})
    d = vapor_mass_flux_divergence(col, DAY, max_loss=np.inf, max_fill=np.inf)[0]
    # the flux divergence is the difference of the face fluxes, so the flux through a face is the sum of the divergences above it
    q = -np.cumsum(d)[:-1] / 1.0                                    # [kg m-2 s-1], downward positive (the sum above is minus the flux)
    mid = slice(5, n - 6)
    K = vapor_latent_conductivity(0.5 * (T[1:] + T[:-1]), rho[1:], {})
    gradient = 0.5
    # heat flux upward (cold above) = K_v * gradient; mass flux is upward (negative) with the same ratio L
    assert np.allclose(2.834e6 * q[mid], -K[mid] * gradient * np.ones(len(K[mid])), rtol=0.05)
    assert np.abs(q[mid]).min() > 0


@pytest.mark.conservation
def test_vapor_latent_heat_is_off_by_default_and_changes_the_temperature(short_config):
    base = _run(short_config, vaporTransport=True)
    off = _run(short_config, vaporTransport=True, vapor_latent_heat=False)
    on = _run(short_config, vaporTransport=True, vapor_latent_heat=True)
    ignored = _run(short_config, vapor_latent_heat=True)           # without vaporTransport the key does nothing
    plain = _run(short_config)
    assert np.array_equal(base.Tz, off.Tz) and np.array_equal(ignored.Tz, plain.Tz)
    assert np.abs(on.Tz - off.Tz).max() > 1e-6                     # it does something
    assert np.abs(on.Tz - off.Tz).max() < 1.0                      # and it is a small correction


### The surface (vapor_surface_exchange; AW00 Eq. 4, q(0) = the flux just beneath the surface)

def _surface_column(c=None):
    n = 10
    T = np.linspace(262.0, 250.0, n)                               # warm surface, as a summer
    return _column(np.full(n, 400.0), T, np.full(n, 0.5), c)


@pytest.mark.unit
def test_surface_flux_is_off_by_default():
    closed = vapor_mass_flux_divergence(_surface_column(), DAY)[0]
    assert abs(closed.sum()) < 1e-12 * np.abs(closed).sum()        # zero flux through the top: the sum is zero
    on = vapor_mass_flux_divergence(_surface_column({'vapor_surface_exchange': True}), DAY)[0]
    assert on.sum() > 0 and on[0] > closed[0]                       # a flux enters the column through the top layer


@pytest.mark.unit
def test_surface_flux_is_the_flux_beneath_the_surface():
    ### warm above, cold below: vapor enters at the top with the gradient between the top two layers
    c = {'vapor_surface_exchange': True}
    col = _surface_column(c)
    d = vapor_mass_flux_divergence(col, DAY, max_fill=np.inf)[0]
    D = vapor_diffusivity(col.Tz[0], 101325.0)
    pore = (1 - 400.0 / 917.0) * firn_tortuosity_inverse(400.0)
    q0 = -D * pore * (saturation_vapor_density(col.Tz[1]) - saturation_vapor_density(col.Tz[0])) / 0.5
    assert q0 > 0 and d.sum() == pytest.approx(q0, rel=1e-9)
    ### a cold surface: the vapor leaves
    cold = _column(np.full(10, 400.0), np.linspace(240.0, 250.0, 10), np.full(10, 0.5), c)
    assert vapor_mass_flux_divergence(cold, DAY)[0].sum() < 0
    ### an isothermal column has none
    iso = _column(np.full(10, 400.0), np.full(10, 250.0), np.full(10, 0.5), c)
    assert np.abs(vapor_mass_flux_divergence(iso, DAY)[0]).max() == 0


@pytest.mark.unit
def test_surface_flux_is_limited_by_the_pore_space():
    ### a thin layer of snow over ice with a warm surface: the flux fills only half of the free pore space in a step
    c = {'vapor_surface_exchange': True, 'vapor_min_thickness': 0.0}
    col = _column(np.array([790.0, 917.0, 917.0]), np.array([270.0, 200.0, 200.0]), np.array([0.0005, 1.0, 1.0]), c)
    d, n = vapor_mass_flux_divergence(col, 1e8)
    room = (917.0 - 790.0) * 0.0005
    assert n >= 1 and d[0] * 1e8 == pytest.approx(0.5 * room, rel=1e-9)


@pytest.mark.conservation
def test_surface_exchange_in_a_run_changes_the_mass_of_the_column(short_config):
    off = _run(short_config, vaporTransport=True)
    on = _run(short_config, vaporTransport=True, vapor_surface_exchange=True)
    assert on.mass.sum() != pytest.approx(off.mass.sum(), rel=1e-12)    # the atmosphere gives or takes mass through the surface
    assert abs(on.mass.sum() / off.mass.sum() - 1) < 0.01               # a small amount
    assert np.all(np.isfinite(on.rho)) and np.all(on.mass >= 0)


### Mars: the critical accumulation rate of AW00 (Eqs. 21-22, Figs. 9 and 10); the script is mars_aw00_check.py

def _mars_check():
    import importlib.util
    import os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mars_aw00_check.py')
    spec = importlib.util.spec_from_file_location('mars_aw00_check', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.unit
def test_mars_critical_accumulation_is_that_of_aw00():
    # b_crit = qbar(0)/(1 - rho_s/rho_i) for Tbar = 175 K, rho_s = 360, dTy = 30 K: AW00 read about 0.5 (Fig. 9) and about 5 (Fig. 10, dTd = 10 K)
    # kg m-2 Mars year-1. Their pressure is not stated (600 Pa here; the flux goes as 1/P) and the pore space differs (their porosity, here Johnsen's 1/tau).
    m = _mars_check()
    b0 = m.critical_accumulation(175.0, 30.0, 0.0, 360.0, jonsen_pore=True)
    b10 = m.critical_accumulation(175.0, 30.0, 10.0, 360.0, jonsen_pore=True)
    assert 0.25 < b0 < 1.0 and 2.5 < b10 < 10.0
    assert 8.0 < b10 / b0 < 12.0                                   # "about a factor of 10 larger" with the diurnal cycle


@pytest.mark.unit
def test_mars_surface_flux_of_the_model_column_follows_eq_21():
    # the flux just beneath the surface of the model (heat solver + top flux of the vapor transport) in a snow half space is that of Eq. 21 (1.5%)
    m = _mars_check()
    model = m.column_flux(175.0, 30.0, 0.0, 360.0, n_per_sol=4, uniform=True)
    analytic = m.critical_accumulation(175.0, 30.0, 0.0, 360.0, jonsen_pore=True, uniform=True)
    assert model == pytest.approx(analytic, rel=0.05)
