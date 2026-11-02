'''
Densification models (physics.FirnPhysics): every option of physRho on a plausible
dry-firn column, called directly (no model run).
'''

import numpy as np
import pytest

from config_validation import ALLOWED_VALUES
from constants import GRAVITY, S_PER_YEAR
from physics import FirnPhysics

N = 300
Z = np.linspace(0.0, 100.0, N)
RHO = np.clip(350.0 + (917.0 - 350.0) * (1 - np.exp(-Z / 35.0)), 350.0, 915.0)
DZ = np.gradient(Z)
SIGMA = np.cumsum(RHO * DZ * GRAVITY)   # overburden stress [Pa]


def params(rho=RHO, **over):
    p = dict(iii=0, steps=12, gridLen=N, rho=rho.copy(), sigma=SIGMA, Tz=np.full(N, 243.15),
             T_mean=np.array([243.15]), T10m=243.15, Ts=np.array([243.15]), z=Z, dz=DZ,
             mass=rho * DZ, age=np.linspace(0.0, 800.0, N) * S_PER_YEAR,
             bdotSec=np.array([0.2 / S_PER_YEAR / 12]), bdot_mean=np.full(N, 0.2), bdot_type='mean',
             dt=S_PER_YEAR / 12, r2=np.full(N, 1e-8), physGrain=True, FirnAir=False,
             LWC=np.zeros(N), Hx=np.ones(N), QMorris=110e3, ind1_old=0, Gamma_old2_Gou=0.0,
             Gamma_old_Gou=0.5 / S_PER_YEAR, Gamma_Gou=0.5 / S_PER_YEAR, bdot_av=np.full(N, 0.2),
             MELT=False, rhos0=350.0, r2s0=1e-8, calcGrainSize=False, GrGrowPhysics='Arthern')
    p.update(over)
    return p


def rate(name, **over):
    '''Densification rate [kg m-3 a-1] of the model `name` on the test column.'''
    return getattr(FirnPhysics(params(**over)), name)()['drho_dt'] * S_PER_YEAR


@pytest.mark.unit
@pytest.mark.parametrize('name', sorted(ALLOWED_VALUES['physRho']))
def test_every_model_gives_a_finite_nonnegative_rate(name):
    d = rate(name)
    assert d.shape == (N,)
    assert np.all(np.isfinite(d))
    assert np.all(d >= 0)
    assert d[(RHO > 400) & (RHO < 800)].max() > 0   # firn is densifying


@pytest.mark.unit
def test_breant2017_is_close_to_goujon2003():
    # Breant et al. (2017) change the temperature dependence of the Goujon creep parameter only
    # (equivalent activation energy close to the 60 kJ/mol of Goujon et al.). With the third term of
    # their Eq. 7 coded with a1 instead of a3 the rate was about 1e10 times that of Goujon.
    goujon, breant = rate('Goujon2003'), rate('Breant2017')
    for rho0 in (600.0, 700.0, 800.0, 880.0):
        i = np.argmin(np.abs(RHO - rho0))
        assert 0.3 < breant[i] / goujon[i] < 3.0, rho0


def viscosity_check(name, **over):
    '''Return (viscosity, expected) for the nodes below ice density, with expected from
    (1/rho) drho/dt = sigma / (2 eta), i.e. eta = rho sigma / (2 drho/dt).'''
    RD = getattr(FirnPhysics(params(**over)), name)()
    d = RD['drho_dt']
    ok = (RHO < 917.0) & (d != 0)
    assert ok.sum() > 100
    return RD['viscosity'][ok], (RHO * SIGMA / (2 * d))[ok]


@pytest.mark.unit
@pytest.mark.parametrize('name', sorted(ALLOWED_VALUES['physRho']))
def test_viscosity_is_newtonian(name):
    # The viscosity output is that of a Newtonian viscous fluid with the model's densification rate.
    v, expected = viscosity_check(name)
    np.testing.assert_allclose(v, expected, rtol=1e-10)


@pytest.mark.unit
@pytest.mark.parametrize('name', ['HLdynamic', 'Li2004', 'Li2011', 'Helsen2008', 'Arthern2010S',
                                  'Ligtenberg2011', 'Simonsen2013', 'KuipersMunneke2015', 'Brils2022'])
def test_viscosity_is_newtonian_with_instant_accumulation(name):
    v, expected = viscosity_check(name, bdot_type='instant')
    np.testing.assert_allclose(v, expected, rtol=1e-10)
