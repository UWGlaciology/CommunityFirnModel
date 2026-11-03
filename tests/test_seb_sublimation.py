'''
Sublimation and deposition from the latent heat flux of the surface energy balance (SEB).

With the turbulent-flux methods (bulk_richardson, monin_obukhov) the latent heat flux QL depends on the
surface temperature calculated by the SEB. The solvers return the mass exchange that goes with it,
QL / L * dt (positive = deposition, negative = sublimation), so that the energy and the mass flux are
consistent. With prescribed fluxes the solvers return None and the SUBLIM forcing is used.
'''

import numpy as np
import pytest

from SEB import SurfaceEnergyBudget
from SEB_turbulent import TurbulentFluxModel
from constants import LS_SUBLIMATION, LV_LIQUID, T_MELT

DAY = 86400.0


class FixedQL(TurbulentFluxModel):
    '''Flux model that depends on Ts (so the Newton solvers are used) and has a known QL for each step.'''
    depends_on_Ts = True

    def __init__(self, QL_series, QH=0.0):
        self.QL_series = np.asarray(QL_series, dtype=float)
        self.QH = QH

    def compute(self, Ts, iii):
        return self.QH, self.QL_series[iii]


def _seb(SW, QL_series, dtratio=1, LW=250.0, T2m=260.0):
    n = 2 * dtratio
    f = lambda v: np.full(n, float(v))
    climate = dict(time=np.arange(n) / 365.25, SW_d=f(SW), LW_d=f(LW), ALBEDO=f(0.8), T2m=f(T2m), TSKIN=f(T2m),
                   QH=f(0.0), QL=f(0.0), RAIN=f(0.0))
    seb = SurfaceEnergyBudget(dict(SEB_TL_thick=0.05, SEB_turb_method='prescribed'), climate, 0,
                              None if dtratio == 1 else dict(time=climate['time'], dtRATIO=dtratio, **{k: climate[k] for k in
                                  ('SW_d', 'LW_d', 'ALBEDO', 'T2m', 'TSKIN', 'QH', 'QL', 'RAIN')}), 0)
    seb.flux_model = FixedQL(QL_series)
    seb.dtRATIO = dtratio
    return seb


def _column(T, n=300, dz0=0.02, gradient=3.0):
    dz = np.full(n, dz0)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    return dict(Tz=T - np.linspace(0.0, gradient, n), mass=dz * np.linspace(330.0, 450.0, n), dt=DAY,
                z=z, dz=dz, mtime=0.0)


@pytest.mark.unit
@pytest.mark.parametrize('QL', [-25.0, -3.0, 0.0, 4.0])
def test_cold_surface_mass_is_QL_over_the_latent_heat_of_sublimation(QL):
    seb = _seb(SW=0.0, QL_series=[QL, QL])
    Ts, _, melt, _, sublim = seb.solve_Ts(_column(255.0), 0, 255.0)
    assert Ts < T_MELT and melt == 0
    assert sublim == pytest.approx(QL / LS_SUBLIMATION * DAY, rel=1e-12)
    assert np.sign(sublim) == np.sign(QL)   # QL < 0 (heat lost) = sublimation = mass loss


@pytest.mark.unit
def test_melting_surface_uses_the_latent_heat_of_vaporization():
    seb = _seb(SW=900.0, QL_series=[-10.0, -10.0], T2m=275.0)
    Ts, _, melt, _, sublim = seb.solve_Ts(_column(272.0), 0, 272.0)
    assert Ts == T_MELT and melt > 0
    assert sublim == pytest.approx(-10.0 / LV_LIQUID * DAY, rel=1e-12)


@pytest.mark.unit
def test_sub_steps_are_summed_with_their_own_flux():
    QL = np.array([-30.0, -10.0, 0.0, 5.0, -20.0, -15.0])   # six 4-hour sub-steps of one model step
    seb = _seb(SW=0.0, QL_series=np.concatenate((QL, QL)), dtratio=6)
    _, _, _, _, sublim = seb.solve_Ts_subdt(_column(255.0), 0, 255.0, 0.0)
    assert sublim == pytest.approx(np.sum(QL) / LS_SUBLIMATION * DAY / 6, rel=1e-12)


@pytest.mark.unit
def test_prescribed_fluxes_return_none():
    f = lambda v: np.full(12, float(v))
    climate = dict(time=np.arange(12) / 365.25, SW_d=f(100.0), LW_d=f(250.0), ALBEDO=f(0.8), T2m=f(260.0),
                   TSKIN=f(260.0), QH=f(10.0), QL=f(-5.0), RAIN=f(0.0))
    seb = SurfaceEnergyBudget(dict(SEB_TL_thick=0.05, SEB_turb_method='prescribed'), climate, 0, None, 0)
    assert seb.solve_Ts(_column(260.0), 0, 260.0)[4] is None
    sub = dict(time=np.arange(12) / 365.25, dtRATIO=6, **{k: climate[k] for k in ('SW_d', 'LW_d', 'ALBEDO', 'T2m', 'TSKIN', 'QH', 'QL', 'RAIN')})
    seb = SurfaceEnergyBudget(dict(SEB_TL_thick=0.05, SEB_turb_method='prescribed'), climate, 0, sub, 0)
    seb.dtRATIO = 6
    assert seb.solve_Ts_subdt(_column(260.0), 0, 260.0, 0.0)[4] is None


### Run level: the sublimation of a turbulent-flux SEB run is the one that goes with QL.

def _run_with_fixed_QL(tmp_path, monkeypatch, QL, method='bulk_richardson'):
    '''Two years of DYE-2 forcing (daily, with the surface energy balance) with the flux model replaced by one
    that gives a constant QL. Returns the model after the run and the climateTS dictionary.'''
    import json
    import pandas as pd
    import RCMpkl_to_spin as RCM
    import SEB
    from cfm_testing import CFM_MAIN, example_config
    from firn_density_nospin import FirnDensityNoSpin

    c = example_config()
    c.update(InputFileFolder=str(CFM_MAIN / c['InputFileFolder']), resultsFolder=str(tmp_path), NewSpin=True,
             variable_srho=False, srho_type='userinput', TWriteStart=1980.0, spinUpdateDate=1980.0, SEB=True,
             MELT=True, SUBLIM=True, SEB_turb_method=method, DFresample='1d', rhos0=350.0,
             outputs=['density', 'depth', 'temperature', 'age', 'DIP', 'LWC', 'climate'])
    df = pd.read_pickle(CFM_MAIN / 'CFMinput_example' / 'CFM_example_66.5_-46.25.pkl').loc['1980':'1981'].copy()
    df['WS'], df['QV2m'], df['PS'] = 5.0, 1e-3, 8e4   # only needed by the turbulent-flux models
    (tmp_path / 'config.json').write_text(json.dumps(c))
    CT, _, _, _, _, SEBfluxes = RCM.makeSpinFiles(df, timeres='1d', spin_date_st=1980.0, spin_date_end=1981.0,
                                                  melt=True, desired_depth=c['H'] - c['HbaseSpin'], SEB=True, num_reps=2)

    class Fixed(TurbulentFluxModel):
        depends_on_Ts = True

        def compute(self, Ts, iii):
            return 0.0, QL

    original = SEB.SurfaceEnergyBudget.__init__

    def init(self, *a, **k):
        original(self, *a, **k)
        if method != 'prescribed':
            self.flux_model = Fixed()

    monkeypatch.setattr(SEB.SurfaceEnergyBudget, '__init__', init)
    model = FirnDensityNoSpin(str(tmp_path / 'config.json'), climateTS=CT, NewSpin=True, SEBfluxes=SEBfluxes)
    model.time_evolve()
    return model, CT


@pytest.mark.conservation
def test_run_sublimation_follows_the_latent_heat_flux(tmp_path, monkeypatch):
    from constants import RHO_I, S_PER_YEAR
    QL = -20.0
    model, _ = _run_with_fixed_QL(tmp_path, monkeypatch, QL)
    assert model.sublim_from_SEB
    rate = lambda L: QL / L / RHO_I * S_PER_YEAR   # m ice eq. per year
    allowed = np.isclose(model.sublim, rate(LS_SUBLIMATION), rtol=1e-4) | np.isclose(model.sublim, rate(LV_LIQUID), rtol=1e-4)
    assert np.all(allowed)                               # the latent heat of sublimation, or of vaporization at the melting point
    assert np.any(np.isclose(model.sublim, rate(LS_SUBLIMATION), rtol=1e-4))   # the cold steps use the latent heat of sublimation
    assert np.allclose(model.sublimSec * S_PER_YEAR * (S_PER_YEAR / model.dt), model.sublim, rtol=1e-6)


@pytest.mark.conservation
def test_run_deposition_becomes_accumulation(tmp_path, monkeypatch):
    from constants import RHO_I, S_PER_YEAR
    QL = 20.0
    model, CT = _run_with_fixed_QL(tmp_path, monkeypatch, QL)
    assert np.all(model.sublim == 0)                     # nothing sublimates
    extra = model.bdot - CT['BDOT'][:len(model.bdot)]    # the deposition has been added to the accumulation
    rate = lambda L: QL / L / RHO_I * S_PER_YEAR
    assert np.all(np.isclose(extra, rate(LS_SUBLIMATION), rtol=1e-4) | np.isclose(extra, rate(LV_LIQUID), rtol=1e-4))


@pytest.mark.conservation
def test_run_with_prescribed_fluxes_keeps_the_sublim_forcing(tmp_path, monkeypatch):
    model, CT = _run_with_fixed_QL(tmp_path, monkeypatch, 0.0, method='prescribed')
    assert not model.sublim_from_SEB
    assert np.allclose(model.sublim, np.minimum(CT['SUBLIM'][:len(model.sublim)], 0))   # deposition moves into bdot
