'''
Heat carried by rain in the surface energy balance. It was computed in J m-2 (rain mass times heat
capacity times temperature difference) and added to the other terms in W m-2, which made it too
large by the length of the time step (a factor 86400 for a daily step), and the sub-daily solvers
left it out.
'''

import numpy as np
import pandas as pd
import pytest

from constants import CP_W, LF_I, RHO_I, S_PER_YEAR, T_MELT
from SEB import SurfaceEnergyBudget

DAY = 86400.0


def _rate(kg_per_m2, seconds):
    '''Rain of kg_per_m2 over a period of the given length, as a rate in m ie/yr (as in climateTS).'''
    return kg_per_m2 / RHO_I * S_PER_YEAR / seconds


def _seb(rain_rate, T2m=278.0, dtratio=1, n=None):
    '''A warm day with melt in every sub-step, less than the mass of the top layer (so the top layer is not
    thickened and the melt is linear in the flux). Rain rate and T2m are constant.'''
    n = n or dtratio
    f = lambda v: np.full(n, float(v))
    sub = dict(time=np.arange(n) / 365.25 / dtratio, SW_d=f(190.), LW_d=f(300.), ALBEDO=f(0.7), T2m=f(T2m),
               TSKIN=f(T2m), QH=f(0.), QL=f(0.), RAIN=f(rain_rate))
    if dtratio == 1:
        return SurfaceEnergyBudget(dict(SEB_TL_thick=0.05), sub, 0, None, 0)
    sub['dtRATIO'] = dtratio
    daily = dict(T2m=np.full(2, float(T2m)))
    return SurfaceEnergyBudget(dict(SEB_TL_thick=0.05), daily, 0, sub, 0)


def _column(T=270.0, n=200):
    dz = np.full(n, 0.02)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    return dict(Tz=np.full(n, T), mass=dz * 350.0, dt=DAY, z=z, dz=dz, mtime=0.0)


@pytest.mark.unit
def test_rain_heat_is_a_flux_in_watts_per_square_metre():
    seb = _seb(0.0)
    rate = _rate(1.0, DAY)   # 1 kg m-2 of rain in a day, at 5 K above melting
    expected = CP_W * (1.0 / DAY) * 5.0
    assert seb._rain_heat(rate, T_MELT + 5.0) == pytest.approx(expected, rel=1e-12)
    assert expected < 0.3    # a small flux, not the 7.7 kW m-2 that the missing time step gave


@pytest.mark.unit
@pytest.mark.parametrize('T_air', [T_MELT, T_MELT - 5.0])
def test_no_rain_heat_when_the_air_is_not_above_melting(T_air):
    assert _seb(0.0)._rain_heat(_rate(10.0, DAY), T_air) == 0.0


@pytest.mark.unit
@pytest.mark.parametrize('method', ['quartic', 'newton'])
def test_daily_melt_from_rain_is_the_rain_heat_over_the_latent_heat(method):
    # With melt, the melt energy is linear in the net flux: the extra melt is Q * dt / LF.
    rain_kg = 20.0
    melts = []
    for rain in (0.0, _rate(rain_kg, DAY)):
        seb = _seb(rain)
        solve = seb._solve_Ts_quartic if method == 'quartic' else seb._solve_Ts_newton
        _, _, melt, _ = solve(_column(), 0, 270.0)
        melts.append(float(np.squeeze(melt)))
    expected = CP_W * (rain_kg / DAY) * (278.0 - T_MELT) * DAY / LF_I
    assert melts[1] - melts[0] == pytest.approx(expected, rel=1e-6)
    assert expected < 5.0   # kg m-2: about a millimetre of melt for 20 mm of warm rain


@pytest.mark.unit
@pytest.mark.parametrize('method', ['quartic', 'newton'])
def test_sub_daily_melt_from_rain_is_the_rain_heat_over_the_latent_heat(method):
    # The same day in six 4-hour sub-steps: each sub-step carries a sixth of the rain.
    rain_kg = 20.0
    melts = []
    for rain in (0.0, _rate(rain_kg / 6, DAY / 6)):
        seb = _seb(rain, dtratio=6)
        if method == 'newton':
            seb.flux_model.depends_on_Ts = True
        _, _, melt, _ = seb.solve_Ts_subdt(_column(), 0, 270.0, 0.0)
        melts.append(float(np.squeeze(melt)))
    expected = CP_W * (rain_kg / DAY) * (278.0 - T_MELT) * DAY / LF_I
    assert melts[1] - melts[0] == pytest.approx(expected, rel=1e-6)


@pytest.mark.unit
def test_sub_daily_rain_in_the_seb_fields_is_a_rate_per_seb_step():
    # makeSpinFiles with 4-hourly SEB fields: the SEB rain rate has to be the same rate as the
    # daily climate rain rate (it used the number of daily steps per year for the 4-hourly mass).
    import RCMpkl_to_spin as RCM
    index = pd.date_range('1980-01-01', '1984-12-31 20:00', freq='4h')
    n = len(index)
    df = pd.DataFrame(dict(T2m=np.full(n, 270.), TSKIN=np.full(n, 268.), SUBLIM=np.zeros(n), QH=np.zeros(n),
                           QL=np.zeros(n), SW_d=np.full(n, 100.), LW_d=np.full(n, 250.), RAIN=np.full(n, 0.5),
                           BDOT=np.full(n, 0.8), ALBEDO=np.full(n, 0.8)), index=index)
    climateTS, _, _, _, _, SEBfluxes = RCM.makeSpinFiles(
        df, timeres='1D', spin_date_st=1980.0, spin_date_end=1981.0, melt=True, SEB=True,
        desired_depth=60., num_reps=2)
    assert SEBfluxes['dtRATIO'] == 6
    expected = 0.5 * 6 * 365.25 / RHO_I * 6 / 6   # 0.5 kg per 4 h step = 3 kg per day, as a rate
    assert np.mean(SEBfluxes['RAIN'][-100:]) == pytest.approx(3.0 * 365.25 / RHO_I, rel=2e-3)
    assert np.mean(SEBfluxes['RAIN'][-100:]) == pytest.approx(np.mean(climateTS['RAIN'][-100:]), rel=2e-3)
