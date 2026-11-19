'''
Isotope diffusion (isotopeDiffusion.py, Johnsen et al. 2000): golden values of the current code.

The test drives isoDiff directly on a synthetic column (the isotope module makes its own surface values when it
has no input file). The expected values are from the code before the diffusivity of water vapor was shared
with vapor.py (see vapor_changes_notes.md): the shared code must give the same numbers.
'''

import numpy as np
import pytest

from isotopeDiffusion import isotopeDiffusion

INDEX = [0, 1, 10, 40, 80, 110, 130, 149]
GOLDEN = {   # isotope: (del_z at INDEX, sum(del_z), iso_sig2_z at INDEX, sum(iso_sig2_z))
    'dD': ([-40.5555555556, -41.6666666667, -39.6658248229, -34.6108557293, -37.7696865369, -46.4678770901,
            -32.7359114872, -42.2720606076], -5716.79971047,
           [0.0, 0.000149794163113, 0.00198217549694, 0.00156202900283, 0.00135481308211, 0.0013029580152,
            0.00128822540537, 0.00128750348158], 0.213380444427),
    'd18O': ([-40.5555555556, -41.6666666667, -39.6658997183, -34.6108752189, -37.7696856011, -46.4678770901,
              -32.7359114872, -42.2720606076], -5716.82753203,
             [0.0, 0.000171048752681, 0.0020811407831, 0.00160254280657, 0.00136515686035, 0.00130538447119,
              0.00128837998133, 0.00128754660681], 0.217271953549),
}


def _column():
    n = 150
    dz = np.linspace(0.02, 0.4, n)
    z = np.concatenate(([0.0], np.cumsum(dz)[:-1]))
    return n, dz, z, np.linspace(330.0, 880.0, n), 250.0 + 6.0 * np.exp(-z / 4.0)


def run_isotope_diffusion(isotope, config_extra=None, steps=5):
    '''Five monthly steps of isoDiff on a synthetic column (330 to 880 kg m-3); returns del_z and iso_sig2_z.'''
    n, dz, z, rho, T = _column()
    config = dict(site_pressure=1013.25, spinup_climate_type='mean')
    config.update(config_extra or {})
    iso = isotopeDiffusion(True, config, isotope, None, 10, z)      # no input file: synthetic surface values
    iso.del_z = -40.0 + 8.0 * np.sin(z / 3.0)
    iso.iso_sig2_z = 1e-3 * np.ones(n)
    iso.del_s = np.linspace(-45.0, -35.0, 10)
    iso.iso_sig2_s = np.zeros(10)
    params = dict(Tz=T, rho=rho, dt=86400.0 * 30, z=z, rhos0=350.0, dz=dz, drho_dt=-1e-8 * rho, bdot=1e-8)
    for i in range(steps):
        del_z, sig2 = iso.isoDiff(params, i)
    return del_z, sig2


@pytest.mark.unit
@pytest.mark.parametrize('isotope', ['dD', 'd18O'])
def test_isotope_diffusion_golden_values(isotope):
    del_z, sig2 = run_isotope_diffusion(isotope)
    d, d_sum, s, s_sum = GOLDEN[isotope]
    assert del_z[INDEX] == pytest.approx(d, rel=1e-9)
    assert del_z.sum() == pytest.approx(d_sum, rel=1e-9)
    assert sig2[INDEX] == pytest.approx(s, rel=1e-9, abs=1e-15)
    assert sig2.sum() == pytest.approx(s_sum, rel=1e-9)


@pytest.mark.unit
def test_isotope_diffusion_stops_where_the_tortuosity_does():
    # 1/tau = 1 - 1.3 (rho/rho_i)^2 is zero from 804 kg m-3. The difference with the module's own pure advection
    # ('NoDiffusion') is the smoothing: large in the upper firn, none from the cutoff down (1e-5 is a boundary effect of the solver).
    n, dz, z, rho, T = _column()
    diffusion = np.abs(run_isotope_diffusion('dD')[0] - run_isotope_diffusion('NoDiffusion')[0])
    assert diffusion[rho < 700.0].max() > 0.5
    assert diffusion[rho >= 804.0].max() < 1e-3


@pytest.mark.unit
def test_isotope_diffusion_uses_the_tortuosity_parameter_of_the_configuration():
    # vapor_tortuosity_b is shared with the vapor transport; a smaller b is a larger diffusivity and a later cutoff
    n, dz, z, rho, T = _column()
    default = run_isotope_diffusion('dD')[0]
    same = run_isotope_diffusion('dD', {'vapor_tortuosity_b': 1.3})[0]
    more = run_isotope_diffusion('dD', {'vapor_tortuosity_b': 1.0})[0]
    advection = run_isotope_diffusion('NoDiffusion')[0]
    assert np.array_equal(default, same)
    assert np.abs(more - advection).sum() > np.abs(default - advection).sum()
