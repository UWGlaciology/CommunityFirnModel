'''
Tests for hl_analytic, the steady-state density and age profiles of Herron and
Langway (1980, J. Glaciol. 25, 373-385).

This file calculates its reference values from the published equations below.
The reference calculation does not call hl_analytic and does not import
constants.py, so an error in the code cannot also be in the reference.

    k0 = 11 exp(-10160 / RT),  k1 = 575 exp(-21400 / RT)
    Zone 1 (rho < 0.55):  ln[rho/(rho_i - rho)] = rho_i k0 h + ln[rho_0/(rho_i - rho_0)]
    Zone 2 (rho > 0.55):  ln[rho/(rho_i - rho)] = rho_i k1 (h - h_0.55) / sqrt(A) + ln[0.55/(rho_i - 0.55)]
    Zone 1 age:  t = ln[(rho_i - rho_0)/(rho_i - rho)] / (k0 A)
    Zone 2 age:  t = ln[(rho_i - 0.55)/(rho_i - rho)] / (k1 sqrt(A)) + t_0.55

Densities are in Mg m-3, h in m, A in m w.e. per year and t in years.
hl_analytic takes the accumulation in m ice equivalent per year and returns
the age in seconds.
'''

import numpy as np
import pytest

from hl_analytic import hl_analytic

R_GAS = 8.314  # J mol-1 K-1
RHO_ICE = 0.917  # Mg m-3
RHO_CRITICAL = 0.550  # Mg m-3, boundary between zone 1 and zone 2
RHO_CLOSE_OFF = 0.830  # Mg m-3
SECONDS_PER_YEAR = 31557600.0

# (temperature [K], accumulation [m ice eq. per year]): from a cold, dry
# East Antarctic site to a warm, wet Greenland site.
SITES = [(218.0, 0.025), (242.0, 0.24), (253.0, 0.5), (263.0, 1.0)]
SURFACE_DENSITIES = [300.0, 350.0, 400.0]  # kg m-3


def _rate_constants(T):
    return 11.0 * np.exp(-10160.0 / (R_GAS * T)), 575.0 * np.exp(-21400.0 / (R_GAS * T))


def _log_ratio(rho):
    return np.log(rho / (RHO_ICE - rho))


def _depth_of(rho, rho_surface, T, A_ice):
    '''Depth [m] of the density horizon rho [Mg m-3].'''
    k0, k1 = _rate_constants(T)
    A = A_ice * RHO_ICE
    h_critical = (_log_ratio(RHO_CRITICAL) - _log_ratio(rho_surface)) / (RHO_ICE * k0)
    if rho <= RHO_CRITICAL:
        return (_log_ratio(rho) - _log_ratio(rho_surface)) / (RHO_ICE * k0)
    return h_critical + np.sqrt(A) * (_log_ratio(rho) - _log_ratio(RHO_CRITICAL)) / (RHO_ICE * k1)


def _age_at(rho, rho_surface, T, A_ice):
    '''Age [years] of the density horizon rho [Mg m-3].'''
    k0, k1 = _rate_constants(T)
    A = A_ice * RHO_ICE
    t_critical = np.log((RHO_ICE - rho_surface) / (RHO_ICE - RHO_CRITICAL)) / (k0 * A)
    if rho <= RHO_CRITICAL:
        return np.log((RHO_ICE - rho_surface) / (RHO_ICE - rho)) / (k0 * A)
    return t_critical + np.log((RHO_ICE - RHO_CRITICAL) / (RHO_ICE - rho)) / (k1 * np.sqrt(A))


def _profile(rhos0, T, A_ice, depth=300.0, dz=0.25):
    # Offset the grid so that no node falls exactly on the zone 1 to zone 2 boundary.
    h = np.arange(0.0, depth, dz) + np.pi * 1e-6
    h[0] = 0.0
    age, rho = hl_analytic(rhos0, h, T, A_ice)
    return h, age, rho


@pytest.mark.unit
@pytest.mark.parametrize('T, A_ice', SITES)
@pytest.mark.parametrize('rhos0', SURFACE_DENSITIES)
def test_surface_density_is_the_input(rhos0, T, A_ice):
    _, age, rho = _profile(rhos0, T, A_ice)
    assert rho[0] == pytest.approx(rhos0, rel=1e-12)
    assert age[0] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.unit
@pytest.mark.parametrize('T, A_ice', SITES)
@pytest.mark.parametrize('rhos0', SURFACE_DENSITIES)
def test_density_and_age_increase_with_depth(rhos0, T, A_ice):
    h, age, rho = _profile(rhos0, T, A_ice)
    assert len(age) == len(rho) == len(h)
    assert np.all(np.diff(rho) > 0)
    assert np.all(np.diff(age) >= 0)
    # hl_analytic holds the age constant where rho > rho_i - 0.01 Mg m-3.
    below_cap = rho[1:] < (RHO_ICE - 0.01) * 1000.0
    assert np.all(np.diff(age)[below_cap] > 0)


@pytest.mark.unit
@pytest.mark.parametrize('T, A_ice', SITES)
def test_density_approaches_ice_density(T, A_ice):
    z_close_off = _depth_of(RHO_CLOSE_OFF, 0.35, T, A_ice)
    _, _, rho = _profile(350.0, T, A_ice, depth=20.0 * z_close_off, dz=1.0)
    assert np.all(np.isfinite(rho))
    assert np.all(rho <= RHO_ICE * 1000.0 * (1 + 1e-12))  # round-off when Z / (1 + Z) is 1
    assert rho[-1] == pytest.approx(RHO_ICE * 1000.0, abs=0.5)


@pytest.mark.unit
@pytest.mark.parametrize('T, A_ice', SITES)
@pytest.mark.parametrize('rhos0', SURFACE_DENSITIES)
@pytest.mark.parametrize('rho_horizon', [0.45, 0.55 - 1e-9, 0.55 + 1e-9, 0.70, RHO_CLOSE_OFF, 0.90])
def test_horizon_depth_and_age_match_herron_langway(rhos0, T, A_ice, rho_horizon):
    rho_surface = rhos0 / 1000.0
    if rho_horizon <= rho_surface:
        pytest.skip('horizon is above the surface density')
    depth = _depth_of(rho_horizon, rho_surface, T, A_ice)
    age, rho = hl_analytic(rhos0, np.array([0.0, depth]), T, A_ice)
    assert rho[1] == pytest.approx(rho_horizon * 1000.0, rel=1e-9)
    assert age[1] / SECONDS_PER_YEAR == pytest.approx(_age_at(rho_horizon, rho_surface, T, A_ice), rel=1e-9)


@pytest.mark.unit
@pytest.mark.parametrize('T, A_ice', SITES)
def test_profile_matches_herron_langway(T, A_ice):
    '''The full profile, against the zone equations solved for density.'''
    rhos0 = 350.0
    h, _, rho = _profile(rhos0, T, A_ice)
    k0, k1 = _rate_constants(T)
    A = A_ice * RHO_ICE
    h_critical = _depth_of(RHO_CRITICAL, rhos0 / 1000.0, T, A_ice)
    log_ratio = np.where(h < h_critical,
                         RHO_ICE * k0 * h + _log_ratio(rhos0 / 1000.0),
                         RHO_ICE * k1 * (h - h_critical) / np.sqrt(A) + _log_ratio(RHO_CRITICAL))
    Z = np.exp(log_ratio)
    np.testing.assert_allclose(rho, 1000.0 * RHO_ICE * Z / (1 + Z), rtol=1e-12)


@pytest.mark.unit
@pytest.mark.xfail(strict=True, raises=IndexError,
                   reason='hl_analytic selects zone 1 with h < h0_55 and zone 2 with h > h0_55, '
                          'so a node at exactly h0_55 is dropped from the profile')
def test_node_on_the_zone_boundary_is_kept():
    rhos0, T, A_ice = 350.0, 242.0, 0.24
    h_critical = _depth_of(RHO_CRITICAL, rhos0 / 1000.0, T, A_ice)
    h = np.array([0.0, h_critical, h_critical + 1.0])
    age, rho = hl_analytic(rhos0, h, T, A_ice)
    assert len(rho) == len(age) == len(h)
    assert rho[1] == pytest.approx(RHO_CRITICAL * 1000.0, rel=1e-9)
