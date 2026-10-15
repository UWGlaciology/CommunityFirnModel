'''
Tests for the liquid-water functions of the bucket and Darcy schemes against
the equations of Verjans et al. (2019, The Cryosphere 13, 1819-1842).

The reference values are calculated in this file from the published equations,
with literal constants. They do not use constants.py.

    Water-holding capacity (Eqs. 9-10, in the form of Langen et al., 2017, Eqs. 3-4):
        W = 0.057 (rho_i - rho) / rho + 0.017
        theta_h = W / (1 - W) * rho rho_i / (rho_w (rho_i - rho))
    Porosity and saturated water content (Eqs. 16-17):
        P = 1 - rho / rho_i,  theta_sat = P rho_i / rho_w
    The irreducible water of a layer is theta_h theta_sat dz.
    van Genuchten and Yamaguchi et al. (2012) (Eqs. 12-15):
        S_e = (1 + (alpha h)^n)^-m
        alpha = 4.4e6 (rho / 2r)^-0.98,  n = 1 + 2.7e-3 (rho / 2r)^0.61,  m = 1 - 1/n
    Mualem relative conductivity (Eq. 20):  K / K_sat = S_e^1/2 (1 - (1 - S_e^1/m)^m)^2
    Calonne et al. (2012) (Eq. 21):  K_sat = 3 r^2 (g rho_w / mu) exp(-0.013 rho)
    Zuo and Oerlemans (1996) runoff (Eqs. 26-27):
        dR/dt = L_excess / tau,  tau = c1 + c2 exp(-c3 S),  c1 = 1.5 d, c2 = 25 d, c3 = 140
'''

import numpy as np
import pytest

import melt
from darcy_funcs import (flux_bisection, flux_newtonraphson, hydrconducsat_Calonne,
                         krel_vG, phead_vG, thetae_update, thetaeff_equaliser,
                         vG_Yama, vG_Yama_params)

RHO_ICE = 917.0  # kg m-3
RHO_WATER = 1000.0  # kg m-3
T_MELT = 273.15  # K
DAY = 86400.0  # s

DENSITIES = [350.0, 450.0, 550.0, 650.0]  # kg m-3
GRAIN_RADII = [2.0e-4, 5.0e-4, 1.0e-3, 1.5e-3]  # m


def _van_genuchten(rho, r):
    ratio = rho / (2 * r)
    alpha = 4.4e6 * ratio ** -0.98
    n = 1 + 2.7e-3 * ratio ** 0.61
    return alpha, n, 1 - 1 / n


def _saturated_water_content(rho):
    return (1 - rho / RHO_ICE) * RHO_ICE / RHO_WATER


def _water_holding_capacity(rho, coleou_lesaffre, irr_val=0.02):
    if not coleou_lesaffre:
        return np.full_like(rho, irr_val)
    W = 0.057 * (RHO_ICE - rho) / rho + 0.017
    return W / (1 - W) * rho * RHO_ICE / (RHO_WATER * (RHO_ICE - rho))


def _irreducible_water(rho, dz, coleou_lesaffre, rho_impermeable, irr_val=0.02):
    '''Irreducible water of each layer [m w.e.]. It is 0 where rho >= rho_impermeable.'''
    lwc = _water_holding_capacity(rho, coleou_lesaffre, irr_val) * _saturated_water_content(rho) * dz
    return np.where(rho >= rho_impermeable, 0.0, lwc)


# ---------------------------------------------------------------------------
# Hydraulic properties


@pytest.mark.unit
@pytest.mark.parametrize('r', GRAIN_RADII)
def test_saturated_conductivity_matches_calonne(r):
    rho = np.array(DENSITIES)
    g, mu = 9.81, 0.001792  # Verjans et al. round g to 9.8 m s-2; the code uses 9.81
    expected = 3.0 * r ** 2 * g * RHO_WATER / mu * np.exp(-0.013 * rho)
    np.testing.assert_allclose(hydrconducsat_Calonne(r, rho), expected, rtol=1e-12)


@pytest.mark.unit
@pytest.mark.parametrize('r', GRAIN_RADII)
def test_van_genuchten_parameters_match_yamaguchi(r):
    rho = np.array(DENSITIES)
    for got, expected in zip(vG_Yama_params(r, rho), _van_genuchten(rho, r)):
        np.testing.assert_allclose(got, expected, rtol=1e-12)


@pytest.mark.unit
@pytest.mark.parametrize('rho', DENSITIES)
@pytest.mark.parametrize('r', GRAIN_RADII)
def test_pressure_head_inverts_the_retention_curve(rho, r):
    alpha, n, m = vG_Yama_params(r, rho)
    saturation = np.logspace(-8, np.log10(1 - 1e-6), 200)
    head = phead_vG(alpha, n, m, saturation)
    assert np.all(head > 0)
    assert np.all(np.diff(head) < 0)  # a wetter layer has less suction
    np.testing.assert_allclose((1 + (alpha * head) ** n) ** -m, saturation, rtol=1e-9)


@pytest.mark.unit
@pytest.mark.parametrize('rho', DENSITIES)
@pytest.mark.parametrize('r', GRAIN_RADII)
def test_relative_conductivity_matches_mualem(rho, r):
    _, _, m = vG_Yama_params(r, rho)
    saturation = np.linspace(0.0, 1.0, 201)
    krel = krel_vG(m, saturation)
    expected = np.sqrt(saturation) * (1 - (1 - saturation ** (1 / m)) ** m) ** 2
    np.testing.assert_allclose(krel, expected, rtol=1e-12, atol=0)
    assert krel[0] == 0.0
    assert krel[-1] == pytest.approx(1.0, rel=1e-12)
    assert np.all(np.diff(krel) > 0)
    assert krel_vG(m, 1 - 1e-9) == pytest.approx(1.0, abs=1e-3)
    assert krel_vG(m, 1e-9) < 1e-6


@pytest.mark.unit
def test_vG_Yama_agrees_with_the_separate_functions():
    rho, r = np.array(DENSITIES), 5.0e-4
    saturation = np.array([0.01, 0.2, 0.5, 0.9])
    alpha, n, m = vG_Yama_params(r, rho)
    head, krel = vG_Yama(r, rho, saturation)
    np.testing.assert_array_equal(head, phead_vG(alpha, n, m, saturation))
    np.testing.assert_array_equal(krel, krel_vG(m, saturation))


# ---------------------------------------------------------------------------
# Flux between two nodes (Hirashima et al., 2010, Eq. 20)

EPS_CONVERGENCE = 0.1e-3  # m, the value in darcyscheme


def _random_two_node_states(seed, count):
    '''
    Two-node states as darcyscheme builds them. The node at index 0 is above the
    node at index 1. Water content is between the residual and the saturated value.
    '''
    rng = np.random.default_rng(seed)
    for _ in range(count):
        rho = rng.uniform(350.0, 650.0, 2)
        r = rng.uniform(2.0e-4, 1.5e-3, 2)
        dz = rng.uniform(0.03, 0.3, 2)
        th_s = _saturated_water_content(rho) - 1e-6
        th_i = np.minimum(0.02, th_s - 1e-9)
        lwc = dz * (th_i + rng.uniform(0.0, 1.0, 2) * (th_s - th_i))
        lwc_available = np.maximum(0.0, lwc - th_i * dz)
        lwc_accommodated = dz * th_s - lwc
        alpha, n, m = vG_Yama_params(r, rho)
        upper_bound = min(lwc_available[0], lwc_accommodated[1])
        guess = min(max(thetaeff_equaliser(th_i, th_s, lwc, dz), 0.0), upper_bound)
        args = (lwc_available, lwc_accommodated, th_i, th_s, lwc, dz, alpha, n, m, EPS_CONVERGENCE)
        yield guess, args, upper_bound


def _equilibrium(flux, args):
    '''f_eq: the head difference minus the node spacing after the flux. The flux is the root.'''
    _, _, th_i, th_s, lwc, dz, alpha, n, m, _ = args
    saturation = thetae_update(flux, th_i, th_s, lwc, dz)
    head = phead_vG(alpha, n, m, saturation)
    return head[0] - head[1] - dz.sum() / 2, saturation


@pytest.mark.unit
@pytest.mark.parametrize('flux_solver', [flux_bisection, flux_newtonraphson])
def test_flux_solver_stops_at_equilibrium_or_at_a_stop_condition(flux_solver):
    # Where the root is inside [0, upper bound], each result must be the root, or
    # satisfy one of the stop conditions of the solvers: the upper node is dry, the
    # lower node is saturated (S_e > 0.95), or the flux is zero.
    from scipy.optimize import brentq
    checked = 0
    for guess, args, upper_bound in _random_two_node_states(seed=1, count=1500):
        if not _equilibrium(0.0, args)[0] < 0 < _equilibrium(upper_bound, args)[0]:
            continue
        root = brentq(lambda g: _equilibrium(g, args)[0], 0.0, upper_bound, xtol=1e-14)
        flux = flux_solver(guess, *args)
        f_eq, saturation = _equilibrium(flux, args)
        at_root = abs(f_eq) < EPS_CONVERGENCE or abs(flux - root) < 2e-6
        stopped = (f_eq < 0 and (saturation[0] < 1e-8 or saturation[1] > 0.95)) or (flux == 0 and f_eq > 0)
        assert at_root or stopped, f'flux {flux}, root {root}, f_eq {f_eq}, S_e {saturation}'
        checked += 1
    assert checked > 500


@pytest.mark.unit
def test_bisection_flux_is_within_the_available_water_and_pore_space():
    for guess, args, upper_bound in _random_two_node_states(seed=2, count=1500):
        flux = flux_bisection(guess, *args)
        assert 0.0 <= flux <= upper_bound


@pytest.mark.unit
def test_newton_raphson_flux_is_within_the_available_water_and_pore_space():
    for guess, args, upper_bound in _random_two_node_states(seed=2, count=1500):
        flux = flux_newtonraphson(guess, *args)
        assert 0.0 <= flux <= upper_bound * (1 + 1e-12)


# ---------------------------------------------------------------------------
# Irreducible water and Zuo-Oerlemans runoff in bucket


def _bucket_column(synthetic_column, rho, dt=DAY, melt_ie=0.15, **config):
    '''A temperate column (no refreeze) with 0.1 m nodes. One bucket call.'''
    settings = dict(Ponding=False, DirectRunoff=0.0, RunoffZuoOerlemans=False, Slope=0.0,
                    ColeouLesaffre=True, IrrVal=0.02, RhoImp=830.0, DownToIce=False, ThickImp=0.15)
    settings.update(config)
    column = synthetic_column(n=len(rho), rho_profile=rho, T_profile=T_MELT, dz=0.1,
                              melt=melt_ie, dt=dt, config=settings)
    melt.bucket(column, 0)
    return column


def _column_with_lenses(density):
    '''Firn of one density, a thin lens at node 10 and a thick lens at nodes 20-21.'''
    rho = np.full(40, density)
    rho[10] = 850.0  # 0.1 m, thinner than ThickImp: permeable
    rho[20:22] = 900.0  # 0.2 m: impermeable
    return rho


@pytest.mark.unit
@pytest.mark.parametrize('coleou_lesaffre', [True, False])
@pytest.mark.parametrize('density', [400.0, 500.0, 600.0])
def test_bucket_leaves_the_irreducible_water_above_a_lens(synthetic_column, coleou_lesaffre, density):
    # The melt is more than the column above the lens can hold. With no refreeze
    # and no ponding, each node above the lens keeps its irreducible water and the
    # rest runs off.
    column = _bucket_column(synthetic_column, _column_with_lenses(density), ColeouLesaffre=coleou_lesaffre)
    assert column.melt_diagnostics['runoff'] > 0
    above_lens = slice(0, np.argmax(column.rho >= 880.0))
    expected = _irreducible_water(column.rho, column.dz, coleou_lesaffre, rho_impermeable=830.0)
    np.testing.assert_allclose(column.LWC[above_lens], expected[above_lens], rtol=1e-9, atol=1e-15)
    assert column.LWC[np.argmax(column.rho == 850.0)] == 0.0  # no irreducible water at rho >= RhoImp


@pytest.mark.unit
@pytest.mark.parametrize('coleou_lesaffre', [True, False])
@pytest.mark.parametrize('slope', [0.0, 0.01])
def test_bucket_runoff_matches_zuo_oerlemans(synthetic_column, coleou_lesaffre, slope):
    # Two calls that differ only in RunoffZuoOerlemans. The runoff is applied last,
    # so the difference in LWC is the runoff of Eqs. 26-27 over one step.
    tau = 1.5 * DAY + 25.0 * DAY * np.exp(-140.0 * slope)
    dt = DAY
    assert dt < tau
    rho = _column_with_lenses(500.0)
    common = dict(Ponding=True, ColeouLesaffre=coleou_lesaffre, Slope=slope, dt=dt)
    ponded = _bucket_column(synthetic_column, rho, RunoffZuoOerlemans=False, **common)
    drained = _bucket_column(synthetic_column, rho, RunoffZuoOerlemans=True, **common)

    irreducible = _irreducible_water(ponded.rho, ponded.dz, coleou_lesaffre, rho_impermeable=830.0)
    excess = np.maximum(ponded.LWC - irreducible, 0.0)
    assert excess.max() > 1e-3  # water ponds above the lens
    expected_runoff = dt * excess / tau
    np.testing.assert_allclose(ponded.LWC - drained.LWC, expected_runoff, rtol=1e-9, atol=1e-15)
    assert drained.melt_diagnostics['runoff'] - ponded.melt_diagnostics['runoff'] == pytest.approx(
        expected_runoff.sum(), rel=1e-9)


@pytest.mark.unit
def test_bucket_runoff_does_not_overshoot_when_the_step_is_longer_than_tau(synthetic_column):
    # Slope 0.1 gives tau = 1.5 d. With an explicit update, a 30-day step removed
    # about 20 times the excess water. The LWC became negative, the refreeze that
    # followed was negative, and the density became negative. The water budget
    # check did not see this. A step longer than tau now removes the excess only.
    rho = _column_with_lenses(500.0)
    common = dict(dt=30 * DAY, Ponding=True, Slope=0.1)
    ponded = _bucket_column(synthetic_column, rho, RunoffZuoOerlemans=False, **common)
    drained = _bucket_column(synthetic_column, rho, RunoffZuoOerlemans=True, **common)
    irreducible = _irreducible_water(drained.rho, drained.dz, True, rho_impermeable=830.0)
    excess = np.maximum(ponded.LWC - irreducible, 0.0)
    assert drained.melt_diagnostics['refrozen'] >= 0.0
    np.testing.assert_array_equal(drained.rho, ponded.rho)
    np.testing.assert_allclose(drained.LWC, np.minimum(ponded.LWC, irreducible), rtol=1e-9, atol=1e-15)
    assert drained.melt_diagnostics['runoff'] - ponded.melt_diagnostics['runoff'] == pytest.approx(
        excess.sum(), rel=1e-9)
    assert all(check['ok'] for check in drained.melt_diagnostics['checks'].values())
