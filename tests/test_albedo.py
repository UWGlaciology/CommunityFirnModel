'''
The albedo of the surface energy balance (SEB): the albedo.py methods, how SEB.py uses them, and the
output 'albedo'.
'''

import h5py
import numpy as np
import pytest

import albedo
from albedo import AlbedoModel, PrescribedAlbedo, make_albedo_model
from SEB import SurfaceEnergyBudget
from constants import T_MELT

DAY = 86400.0


def _seb(config=None, dtratio=1, albedo_values=None, SW=300.0):
    n = 3 * dtratio
    f = lambda v: np.full(n, float(v))
    climate = dict(time=np.arange(n) / 365.25 / dtratio, SW_d=f(SW), LW_d=f(250.0), T2m=f(268.0), TSKIN=f(268.0),
                   QH=f(0.0), QL=f(0.0), RAIN=f(0.0),
                   ALBEDO=np.asarray(albedo_values if albedo_values is not None else f(0.8), dtype=float))
    config = dict(SEB_TL_thick=0.05, SEB_turb_method='prescribed', **(config or {}))
    keys = ('SW_d', 'LW_d', 'ALBEDO', 'T2m', 'TSKIN', 'QH', 'QL', 'RAIN')
    sub = None if dtratio == 1 else dict(time=climate['time'], dtRATIO=dtratio, **{k: climate[k] for k in keys})
    seb = SurfaceEnergyBudget(config, climate, 0, sub, 0)
    seb.dtRATIO = dtratio
    return seb


def _column(T=268.0, n=300, dz0=0.02):
    dz = np.full(n, dz0)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    return dict(Tz=T - np.linspace(0.0, 3.0, n), mass=dz * np.linspace(330.0, 450.0, n), dt=DAY, z=z, dz=dz, mtime=0.0)


@pytest.mark.unit
def test_prescribed_albedo_fills_nan_and_applies_the_factor():
    forcing = {'ALBEDO': np.array([0.8, np.nan, 0.6], dtype='float32')}
    model = PrescribedAlbedo({'albedo_factor': 0.5}, forcing)
    assert np.allclose(model.step(1, 0, 3, {}), [0.4, 0.35, 0.3])
    assert np.isnan(forcing['ALBEDO'][1])   # the forcing is not changed


@pytest.mark.unit
def test_prescribed_albedo_gives_the_values_of_the_sub_steps():
    model = PrescribedAlbedo({}, {'ALBEDO': np.arange(6.0) / 10})
    assert np.array_equal(model.step(1, 3, 3, {}), [0.3, 0.4, 0.5])


@pytest.mark.unit
def test_albedo_method_is_checked():
    with pytest.raises(ValueError, match='albedo_method'):
        make_albedo_model({'albedo_method': 'not_a_method'}, {})
    with pytest.raises(KeyError, match='ALBEDO'):
        make_albedo_model({}, {'SW_d': np.zeros(3)})
    assert isinstance(make_albedo_model({}, {'ALBEDO': np.zeros(3)}), PrescribedAlbedo)   # the default


class ConstantAlbedo(AlbedoModel):
    '''An albedo that does not come from the forcing, and records what the SEB asks for.'''

    def __init__(self, config, forcing):
        super().__init__(config, forcing)
        self.calls = []

    def step(self, iii, i0, n, PhysParams):
        self.calls.append((iii, i0, n, sorted(PhysParams)))
        return np.linspace(0.2, 0.9, n) if n > 1 else np.array([0.5])


@pytest.fixture
def constant_albedo(monkeypatch):
    monkeypatch.setitem(albedo.ALBEDO_METHODS, 'constant', ConstantAlbedo)


@pytest.mark.unit
@pytest.mark.parametrize('dtratio', [1, 4])
@pytest.mark.parametrize('solver', ['quartic', 'newton'])
def test_seb_uses_the_albedo_of_the_model(constant_albedo, dtratio, solver):
    # Same forcing, but the 'constant' model has a lower albedo than the forcing: more absorbed shortwave radiation
    # gives a warmer surface, and (SW is large) melt.
    dark = _seb(dict(albedo_method='constant'), dtratio, albedo_values=np.full(3 * dtratio, 0.95))
    bright = _seb({}, dtratio, albedo_values=np.full(3 * dtratio, 0.95))
    if dtratio == 1:
        solve = dark._solve_Ts_quartic if solver == 'quartic' else dark._solve_Ts_newton
        solve_b = bright._solve_Ts_quartic if solver == 'quartic' else bright._solve_Ts_newton
        a, b = solve(_column(), 0, 268.0), solve_b(_column(), 0, 268.0)
    else:
        a, b = dark.solve_Ts_subdt(_column(), 0, 268.0, 0.0), bright.solve_Ts_subdt(_column(), 0, 268.0, 0.0)
    assert a[0] > b[0] or a[2] > b[2]   # warmer, or more melt
    assert dark.albedo_out == pytest.approx(0.5 if dtratio == 1 else np.mean(np.linspace(0.2, 0.9, dtratio)))
    iii, i0, n, state = dark.albedo_model.calls[0]
    assert (iii, i0, n) == (0, 0, dtratio)
    assert {'Tz', 'mass', 'z', 'dz', 'dt'} <= set(state)   # the model state is given to the albedo model


@pytest.mark.unit
def test_albedo_out_is_the_mean_albedo_of_the_step():
    seb = _seb({}, dtratio=3, albedo_values=[0.9, 0.8, 0.7, 0.5, 0.5, 0.5, 0.4, 0.4, 0.4])
    seb.solve_Ts_subdt(_column(), 0, 268.0, 0.0)
    assert seb.albedo_out == pytest.approx(0.8)
    seb.solve_Ts_subdt(_column(), 1, 268.0, 0.0)
    assert seb.albedo_out == pytest.approx(0.5)


@pytest.mark.conservation
@pytest.mark.parametrize('resolution', ['1D', '1ME'])   # the model step has the SEB forcing, or has 30 SEB steps
def test_albedo_is_an_output_of_a_seb_run(short_config, resolution):
    run = short_config(years=2, num_reps=2, liquid='bucket', DFfile='CFM_example_66.5_-46.25.pkl', SEB=True,
                       MELT=True, DFresample=resolution, albedo_factor=0.9)
    run.run()
    with h5py.File(run.results_path, 'r') as f:
        alb = f['albedo'][()]
    assert alb.shape[1] == 2
    assert np.isnan(alb[0, 1])   # the initial state
    assert np.all((alb[1:, 1] > 0) & (alb[1:, 1] < 1))
    if resolution == '1D':   # one SEB step per model step: the albedo of the forcing, times the factor
        t = run.climateTS['time']
        step = np.argmin(np.abs(t[None, :] - alb[1:, 0, None]), axis=1)   # the forcing step of each written time
        expect = run.climateTS['ALBEDO'][step] * 0.9
        assert np.allclose(alb[1:, 1], expect, rtol=1e-5)


@pytest.mark.conservation
def test_albedo_is_not_an_output_without_seb(short_config):
    run = short_config(years=2, num_reps=2, liquid='bucket', outputs=['density', 'depth', 'temperature', 'albedo'])
    run.run()
    with h5py.File(run.results_path, 'r') as f:
        assert 'albedo' not in f


### Gardner and Sharp (2010) and Greuell and Konzelmann (1994) ###
###################################################################

from albedo import (GardnerSharpAlbedo, GreuellKonzelmannAlbedo, cos_zenith, gardner_sharp_albedo,
                    specific_surface_area)


@pytest.mark.unit
def test_specific_surface_area():
    assert specific_surface_area(1e-4) == pytest.approx(329.67, rel=1e-4)   # 3 / (910 kg m-3 * 0.1 mm), in cm2 g-1


@pytest.mark.unit
def test_gardner_sharp_albedo_reproduces_the_values_of_the_paper():
    # Values read from the figures of Gardner and Sharp (2010) (Fig 9b, 10c); the paper's standard error is 0.01.
    S = specific_surface_area
    assert gardner_sharp_albedo(S(1e-4)) == pytest.approx(0.815, abs=0.003)                     # new snow, re = 0.1 mm
    assert gardner_sharp_albedo(S(1e-3), c=0.3) == pytest.approx(0.563, abs=0.006)              # old dirty snow
    assert gardner_sharp_albedo(S(1e-4), u=np.cos(np.radians(85))) == pytest.approx(0.885, abs=0.003)
    assert gardner_sharp_albedo(S(1e-3), tau=10) == pytest.approx(0.81, abs=0.012)              # Fig 10c, re = 1 mm


@pytest.mark.unit
def test_gardner_sharp_albedo_responds_in_the_physical_direction():
    S = specific_surface_area(np.array([5e-5, 1e-4, 5e-4, 2e-3]))      # finer to coarser snow
    a = gardner_sharp_albedo(S)
    assert np.all(np.diff(a) < 0)                                      # coarser snow is darker
    assert np.all(gardner_sharp_albedo(S, c=0.1) < a)                  # carbon lowers it
    u = np.cos(np.radians([0, 40, 70, 85]))
    assert np.all(np.diff(gardner_sharp_albedo(S[1], u=u)) > 0)        # a lower sun raises it
    tau = np.array([0, 1, 5, 30])
    assert np.all(np.diff(gardner_sharp_albedo(S[1], tau=tau)) > 0)    # cloud raises it
    assert gardner_sharp_albedo(S[1], c=0.0) == pytest.approx(gardner_sharp_albedo(S[1]))


@pytest.mark.unit
def test_cos_zenith():
    t = lambda doy, hour: 2020.0 + (doy - 1 + hour / 24) / 366     # 2020 is a leap year
    assert cos_zenith(t(173, 12), 66.5, 0.0) == pytest.approx(np.cos(np.radians(66.5 - 23.44)), abs=2e-3)   # June solstice, noon
    assert cos_zenith(t(173, 0), 80.0, 0.0) > 0                      # polar day: the sun is up at midnight
    assert cos_zenith(t(356, 12), 80.0, 0.0) < 0                     # polar night
    assert cos_zenith(t(80, 12), 0.0, 0.0) == pytest.approx(1.0, abs=0.01)   # equinox, equator, noon
    assert cos_zenith(t(173, 12), 66.5, 90.0) < cos_zenith(t(173, 12), 66.5, 0.0)   # at 90 E it is 6 h before noon (UTC)


def _gs_forcing(n=366, step=1 / 366, start=2020.0, **extra):
    return dict(time=start + np.arange(n) * step, **extra)


@pytest.mark.unit
def test_gardner_sharp_model_is_the_mean_over_the_day_weighted_by_the_cosine():
    lat = 66.5
    model = GardnerSharpAlbedo(dict(albedo_latitude=lat, r2s0=1e-8), _gs_forcing())
    i = 171                                                      # a day in June
    t0 = _gs_forcing()['time'][i]
    hours = t0 + (np.arange(24) + 0.5) / 24 / 366                # the day that starts at the time of the step
    u = cos_zenith(hours, lat, 0.0)
    a = gardner_sharp_albedo(specific_surface_area(1e-4), 0.0, np.maximum(u, np.cos(np.radians(85))), 0.0)
    w = np.maximum(u, 0)
    assert model.step(i, i, 1, dict(r2=None, r2s0=1e-8))[0] == pytest.approx((a * w).sum() / w.sum(), rel=1e-4)


@pytest.mark.unit
def test_gardner_sharp_model_inputs():
    forcing = _gs_forcing()
    state = dict(r2=None, r2s0=1e-8)
    summer, winter = 171, 0
    clear = GardnerSharpAlbedo(dict(albedo_latitude=66.5), forcing)
    a = clear.step(summer, summer, 1, state)[0], clear.step(winter, winter, 1, state)[0]
    assert 0.8 < a[0] < a[1] < 0.9                                # the low winter sun gives the higher albedo
    cloudy = GardnerSharpAlbedo(dict(albedo_latitude=66.5, albedo_cloud_tau=10.0), forcing)
    assert cloudy.step(summer, summer, 1, state)[0] > a[0] + 0.05
    from_forcing = GardnerSharpAlbedo(dict(albedo_latitude=66.5, albedo_cloud_tau=0.0), _gs_forcing(TAU=np.full(366, 10.0)))
    assert from_forcing.step(summer, summer, 1, state)[0] == pytest.approx(cloudy.step(summer, summer, 1, state)[0])
    dirty = GardnerSharpAlbedo(dict(albedo_latitude=66.5, albedo_carbon=0.5), forcing)
    assert dirty.step(summer, summer, 1, state)[0] < a[0] - 0.05
    # the grain size of the surface node, when the grain growth is on
    z, mass = np.arange(0.0, 1.0, 0.01), np.full(100, 4.0)       # 1 cm nodes
    coarse = dict(r2=np.full(100, 1e-6), r2s0=1e-8, z=z, mass=mass)   # r = 1 mm
    assert clear.step(summer, summer, 1, coarse)[0] < a[0] - 0.05
    # the mean is over the top layer (5 cm: nodes 0 to 5) and is weighted by mass: S, not the radius, is averaged
    r2 = np.full(100, 1e-8); r2[3:] = 1e-6
    mixed = dict(r2=r2, r2s0=1e-8, z=z, mass=mass)
    S = clear.surface_S(mixed)
    assert S == pytest.approx((3 * specific_surface_area(1e-4) + 3 * specific_surface_area(1e-3)) / 6)
    shallow = GardnerSharpAlbedo(dict(albedo_latitude=66.5, albedo_grain_depth=0.0), forcing)
    assert shallow.surface_S(mixed) == pytest.approx(specific_surface_area(1e-4))   # the top node only
    assert clear.step(summer, summer, 3, state).shape == (3,)


@pytest.mark.unit
def test_gardner_sharp_model_steps_longer_than_a_day():
    forcing = _gs_forcing(n=24, step=1 / 12)                      # monthly steps
    model = GardnerSharpAlbedo(dict(albedo_latitude=66.5, r2s0=1e-8), forcing)
    a = np.array([model.step(i, i, 1, dict(r2=None, r2s0=1e-8))[0] for i in range(12)])
    assert np.all((a > 0.8) & (a < 0.9)) and a.argmin() in (4, 5, 6, 7)   # lowest around the summer


@pytest.mark.unit
def test_gardner_sharp_model_needs_the_latitude():
    with pytest.raises(KeyError, match='latitude'):
        make_albedo_model(dict(albedo_method='gardner_sharp'), _gs_forcing())
    model = make_albedo_model(dict(albedo_method='gardner_sharp', lat_val=66.5), _gs_forcing())   # lat_val is used
    assert isinstance(model, GardnerSharpAlbedo)


@pytest.mark.unit
def test_greuell_konzelmann_is_linear_in_the_density_between_snow_and_ice():
    cfg = dict(albedo_method='greuell_konzelmann', albedo_snow=0.84, albedo_ice=0.4)
    model = make_albedo_model(cfg, {})
    state = lambda rho: dict(rho=np.array([rho, 500.0]))
    assert model.step(0, 0, 1, state(300.0))[0] == pytest.approx(0.84)
    assert model.step(0, 0, 1, state(910.0))[0] == pytest.approx(0.4)
    assert model.step(0, 0, 1, state(605.0))[0] == pytest.approx(0.62)
    assert model.step(0, 0, 1, state(100.0))[0] == pytest.approx(0.84)     # fresher than fresh snow: not above albedo_snow
    assert model.step(0, 0, 1, state(917.0))[0] == pytest.approx(0.4)      # ice
    assert model.step(0, 0, 4, state(450.0)).shape == (4,)
    other = make_albedo_model(dict(cfg, albedo_rho_snow=100.0, albedo_rho_ice=900.0), {})
    assert other.step(0, 0, 1, state(500.0))[0] == pytest.approx(0.62)


@pytest.mark.unit
def test_greuell_konzelmann_needs_its_constants():
    with pytest.raises(KeyError, match='albedo_snow'):
        make_albedo_model(dict(albedo_method='greuell_konzelmann'), {})
    with pytest.raises(ValueError, match='albedo_rho'):
        make_albedo_model(dict(albedo_method='greuell_konzelmann', albedo_snow=0.8, albedo_ice=0.4,
                               albedo_rho_snow=900.0, albedo_rho_ice=300.0), {})


@pytest.mark.unit
@pytest.mark.parametrize('method, extra, expect', [
    ('greuell_konzelmann', dict(albedo_snow=0.84, albedo_ice=0.4), 0.84 - (350.0 - 300.0) / (910.0 - 300.0) * 0.44),
    ('gardner_sharp', dict(albedo_latitude=66.5, r2s0=1e-8), None),
])
def test_seb_runs_with_the_albedo_methods_without_an_albedo_forcing(method, extra, expect):
    seb = _seb(dict(albedo_method=method, **extra))
    del seb.albedo_model   # rebuild without ALBEDO in the forcing: the methods do not need it
    n = 3
    f = lambda v: np.full(n, float(v))
    forcing = dict(time=np.arange(n) / 365.25, SW_d=f(300.0))
    seb.albedo_model = make_albedo_model(seb.c, forcing)
    column = _column()
    column.update(rho=column['mass'] / column['dz'], r2=None, r2s0=1e-8)
    column['rho'][0] = 350.0
    Ts, _, melt, _, _ = seb._solve_Ts_quartic(column, 0, 268.0)
    assert np.isfinite(Ts) and 0.7 < seb.albedo_out < 0.95
    if expect is not None:
        assert seb.albedo_out == pytest.approx(expect)


@pytest.fixture
def no_albedo_forcing(tmp_path):
    '''A folder with the DYE-2 dataframe without the column ALBEDO (the new methods do not need it).'''
    import pandas as pd
    from cfm_testing import CFM_MAIN
    df = pd.read_pickle(CFM_MAIN / 'CFMinput_example' / 'CFM_example_66.5_-46.25.pkl').drop(columns=['ALBEDO'])
    folder = tmp_path / 'forcing'
    folder.mkdir()
    df.to_pickle(folder / 'no_albedo.pkl')
    return dict(InputFileFolder=str(folder), DFfile='no_albedo.pkl')


@pytest.mark.conservation
@pytest.mark.parametrize('method, extra', [
    ('gardner_sharp', dict(albedo_latitude=66.5, albedo_longitude=-46.25, albedo_carbon=0.05)),
    ('greuell_konzelmann', dict(albedo_snow=0.84, albedo_ice=0.4)),
])
def test_albedo_methods_in_a_seb_run(short_config, no_albedo_forcing, method, extra):
    run = short_config(years=2, num_reps=2, liquid='bucket', SEB=True, MELT=True, DFresample='1D',
                       albedo_method=method, **extra, **no_albedo_forcing)
    run.run()
    with h5py.File(run.results_path, 'r') as f:
        alb = f['albedo'][()]
    a = alb[1:, 1]
    assert np.all(np.isfinite(a)) and a.std() > 0       # the albedo changes in time
    if method == 'gardner_sharp':
        assert 0.7 < a.min() and a.max() < 0.95
    else:
        assert 0.4 <= a.min() and a.max() <= 0.84
