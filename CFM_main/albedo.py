#!/usr/bin/env python

'''
albedo.py

Surface albedo for the surface energy balance (SEB.py).

The albedo is chosen with the configuration key 'albedo_method'. SurfaceEnergyBudget
builds the selected model with make_albedo_model() and, at each model time step, asks it for
the albedo of the shortwave radiation: SW_net = SW_d * (1 - albedo). The albedo that was used
is also an output of the model (output 'albedo').

Methods
-------
'prescribed'          the albedo is the forcing column ALBEDO (default)
'gardner_sharp'       Gardner and Sharp (2010), physically based, from the grain size of the surface, the
                      solar zenith angle, the cloud optical thickness and the carbon content
'greuell_konzelmann'  Greuell and Konzelmann (1994), from the density of the surface snow

HOW TO ADD A PARAMETERIZATION
-----------------------------
1. Write a class that inherits from AlbedoModel and implements step().
2. Add it to ALBEDO_METHODS (the key is the value of 'albedo_method'), and to the options of
   'albedo_method' in config_validation.py, in the albedo_method_options helper of the example
   configuration files, and in docs/running/json.rst and docs/extras/seb.rst.
3. Put the reference of the parameterization, and its coefficients, in the docstring of the class.
   Check the coefficients against the paper: mark them VERIFY until that is done.
4. Add a test (tests/test_albedo.py).

What step() gets: the model step iii, the index i0 of the first forcing element of the step
and the number n of SEB sub-steps in the step (n is 1 when the SEB forcing has the model time
step), and PhysParams, the dictionary of the model state at the start of the step (the same
dictionary the SEB solvers get: for example rho, dz, mass, Tz, z, LWC, age, r2, Ts, dt, mtime;
see FirnDensityNoSpin.time_evolve). It must return an array of n albedos, one per SEB
sub-step. The state is the one before the surface temperature of this step is solved, so
an albedo that depends on the melt or the temperature uses the last step's values. The
constructor gets the forcing (below), so a parameterization that needs, for example, the
snowfall or the air temperature takes it from there. A model that keeps a state (the age of the
snow, for example) is updated in step(), which the SEB calls once for each model step, in
order. (The sub-step solvers call it once per model step, not once per sub-step.)
'''

import numpy as np


class AlbedoModel:
    '''
    Base class of the albedo models.

    :param config: the configuration dictionary
    :param forcing: dictionary of the SEB forcing arrays (the keys of the forcing dataframe, for
        example SW_d, T2m, BDOT, RAIN, ALBEDO), starting at the first step of the run, on
        the time axis of the SEB forcing ('time' is in it)
    '''

    # Names of the forcing arrays and of the configuration keys that the model needs; make_albedo_model
    # checks that they are present.
    required_forcing = ()
    required_config = ()

    def __init__(self, config, forcing):
        self.c = config

    def step(self, iii, i0, n, PhysParams):
        '''
        Albedo of the n SEB steps of the model step iii (forcing indices i0 to i0 + n - 1).

        :return: array of n albedos (dimensionless, between 0 and 1)
        '''
        raise NotImplementedError
### end AlbedoModel ###
#######################


class PrescribedAlbedo(AlbedoModel):
    '''
    The albedo of the forcing (column ALBEDO). Missing values (NaN) are replaced by the mean of the
    series (MERRA-2 gives no albedo when there is no shortwave radiation), and the
    series is multiplied by the configuration key 'albedo_factor' (default 1), which is a simple way
    to tune the melt.
    '''

    required_forcing = ('ALBEDO',)

    def __init__(self, config, forcing):
        super().__init__(config, forcing)
        self.ALBEDO = np.array(forcing['ALBEDO'])  # a copy: the forcing is not changed
        self.ALBEDO[np.isnan(self.ALBEDO)] = np.nanmean(self.ALBEDO)
        if 'albedo_factor' in self.c:
            self.albedo_factor = self.c['albedo_factor']
            self.ALBEDO = self.ALBEDO * self.albedo_factor
            print(f'ALBEDO is multipled by {self.albedo_factor}')
        else:
            print('"albedo_factor" not defined in .json. Using 1')

    def step(self, iii, i0, n, PhysParams):
        return self.ALBEDO[i0:i0 + n]
### end PrescribedAlbedo ###
############################


RHO_ICE_GS = 910.0   # density of ice [kg m-3] in the definition of the effective radius of Gardner and Sharp (2010)
U_MIN_GS = np.cos(np.radians(85.0))   # the parameterization is fitted for solar zenith angles up to 85 degrees
DAY = 1.0 / 365.25   # a day in years (the model time is a decimal date)


def specific_surface_area(r_e):
    '''
    Specific surface area S [cm2 g-1] of snow with the optical effective radius r_e [m]: S = 3 / (rho_i r_e),
    with rho_i = 910 kg m-3 (Gardner and Sharp, 2010, section 2.1). A radius of 0.1 mm is 330 cm2 g-1.
    '''
    return 3.0 / (RHO_ICE_GS * np.asarray(r_e, dtype=float)) * 10.0   # m2 kg-1 to cm2 g-1
### end specific_surface_area ###
#################################


def gardner_sharp_albedo(S, c=0.0, u=1.0, tau=0.0):
    '''
    Broadband albedo of semi-infinite snow or ice, equations 7 to 12 of Gardner and Sharp (2010), J. Geophys.
    Res. 115, F01009, doi:10.1029/2009JF001444. The equations were checked against the typeset paper. The
    snow-depth correction (their equation 13, for a snow layer over ice) is not included.

    Valid for S of 0.07 to 1300 cm2 g-1, c of 0 to 2 ppmw, solar zenith angles up to 85 degrees (u >= 0.087)
    and tau of 0 to 30; the inputs are not limited here. All inputs can be arrays (they are broadcast).

    :param S: specific surface area [cm2 g-1] (see specific_surface_area)
    :param c: concentration of light-absorbing carbon [ppmw]
    :param u: cosine of the solar zenith angle
    :param tau: cloud optical thickness
    :return: broadband albedo
    '''
    S, c, u, tau = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (S, c, u, tau)))
    a_S = 1.48 - S ** -0.07                                                                       # (7)
    da_c = np.maximum(0.04 - a_S, -c ** 0.55 / (0.16 + 0.6 * S ** 0.5 + 1.8 * c ** 0.6 * S ** -0.25))   # (8)
    a_c = a_S + da_c
    x = np.minimum((tau / (3.0 * u)) ** 0.5, 1.0)                                                 # (10)
    u_eff = 0.64 * x + (1.0 - x) * u
    da_u = 0.53 * a_S * (1.0 - a_c) * (1.0 - u_eff) ** 1.2                                        # (9) with u'
    da_tau = 0.1 * tau * a_c ** 1.3 / (1.0 + 1.5 * tau) ** a_S                                    # (11)
    return a_S + da_c + da_u + da_tau                                                             # (12)
### end gardner_sharp_albedo ###
################################


def cos_zenith(time, lat, lon=0.0):
    '''
    Cosine of the solar zenith angle at the decimal date time (UTC) at latitude lat and longitude lon
    [degrees; east positive]. The solar declination is the Fourier series of Spencer (1971) as in Iqbal (1983).
    The equation of time (up to 16 minutes) is not included.

    :param time: decimal date (array)
    :return: cosine of the zenith angle (negative at night)
    '''
    time = np.asarray(time, dtype=float)
    year = np.floor(time)
    leap = ((year % 4 == 0) & (year % 100 != 0)) | (year % 400 == 0)
    ndays = np.where(leap, 366.0, 365.0)
    day = (time - year) * ndays                    # days since 1 January 00:00
    hour = (day % 1.0) * 24.0
    G = 2.0 * np.pi * (day - 0.5) / ndays          # day angle
    decl = (0.006918 - 0.399912 * np.cos(G) + 0.070257 * np.sin(G) - 0.006758 * np.cos(2 * G)
            + 0.000907 * np.sin(2 * G) - 0.002697 * np.cos(3 * G) + 0.00148 * np.sin(3 * G))
    h = np.radians((hour + lon / 15.0 - 12.0) * 15.0)   # hour angle
    phi = np.radians(lat)
    return np.sin(phi) * np.sin(decl) + np.cos(phi) * np.cos(decl) * np.cos(h)
### end cos_zenith ###
######################


class GardnerSharpAlbedo(AlbedoModel):
    '''
    Gardner and Sharp (2010) physically based broadband albedo of snow (gardner_sharp_albedo, equations 7 to 12).

    * S, the specific surface area, is from the grain size of the snow at the surface: the mean of S of the nodes of the
      top layer, weighted by their mass (S is per unit mass, so it adds up), for the nodes down to the depth
      albedo_grain_depth [m] (default 0.05, the top layer of the SEB). The grain radius of a node is sqrt(r2) when the grain
      growth is on (physGrain), else sqrt(r2s0), a constant. (The top node is new snow in a step with accumulation: only a layer
      deeper than the new snow shows the growth of the grains.) VERIFY: whether the grain radius of the CFM is the optical
      effective radius of the paper.
    * The solar zenith angle is calculated from the latitude (configuration key albedo_latitude, else lat_val)
      and longitude (albedo_longitude, else lon_val, else 0; it matters for steps shorter than a day), and the time of
      the forcing. The albedo of a step is the mean of the albedo at times in the step, weighted by the cosine of the
      zenith angle (the shortwave radiation of clear sky): hourly samples for a step up to a few days, daily
      samples (24 hourly values each, at most 10 days) for a longer one. A step is the interval that starts at its
      time (as for a daily mean of a dataframe). The sun is below the horizon: the albedo for the zenith angle of 85 degrees.
    * c, the carbon in the snow, is the configuration key albedo_carbon [ppmw] (default 0).
    * tau, the cloud optical thickness, is the forcing column TAU if there is one, else the configuration key
      albedo_cloud_tau (default 0: clear sky). It is not estimated from the radiation.
    '''

    required_forcing = ('time',)

    def __init__(self, config, forcing):
        super().__init__(config, forcing)
        lat = config.get('albedo_latitude', config.get('lat_val'))
        if lat is None:
            raise KeyError("albedo_method 'gardner_sharp' needs the latitude: the configuration key 'albedo_latitude' (or 'lat_val').")
        lon = config.get('albedo_longitude', config.get('lon_val', 0.0))
        self.carbon = float(config.get('albedo_carbon', 0.0))
        self.grain_depth = float(config.get('albedo_grain_depth', 0.05))
        t = np.asarray(forcing['time'], dtype=float)
        n = len(t)
        if 'TAU' in forcing:
            self.tau = np.asarray(forcing['TAU'], dtype=float)
        else:
            self.tau = np.full(n, float(config.get('albedo_cloud_tau', 0.0)))

        # the times in each step, for the mean over the zenith angles of the step
        dt = np.diff(t)
        step = np.median(dt[dt > 0]) if np.any(dt > 0) else DAY
        dt = np.append(dt, dt[-1] if len(dt) else step)
        dt = np.where(dt > 0, dt, step)                       # [yr]
        if step < DAY:
            m = max(1, int(round(step / DAY * 24)))
            off = (np.arange(m) + 0.5) / m * dt[:, None]
        else:
            nd = min(int(np.ceil(step / DAY)), 10)
            centres = (np.arange(nd) + 0.5) / nd * dt[:, None]                # centre of each sample day [yr]
            hours = (np.arange(24) + 0.5) / 24 * DAY - 0.5 * DAY               # a day around its centre [yr]
            off = (centres[:, :, None] + hours[None, None, :]).reshape(n, nd * 24)
        u = cos_zenith(t[:, None] + off, lat, lon)
        self._weight = np.maximum(u, 0.0)
        self._u = np.maximum(u, U_MIN_GS)
        if not config.get('physGrain', False):
            print("albedo_method 'gardner_sharp': physGrain is false, so the specific surface area is the constant of r2s0.")

    def surface_S(self, PhysParams):
        '''Specific surface area [cm2 g-1] of the snow of the top layer (mean of the nodes, weighted by mass).'''
        r2 = PhysParams.get('r2')
        if r2 is None:
            return specific_surface_area(np.sqrt(PhysParams['r2s0']))
        z, mass = PhysParams['z'], PhysParams['mass']
        i = min(np.searchsorted(z, self.grain_depth), len(z) - 1)   # the first node at or below the depth, as in SEB.py
        S = specific_surface_area(np.sqrt(r2[:i + 1]))
        return np.sum(mass[:i + 1] * S) / np.sum(mass[:i + 1])

    def step(self, iii, i0, n, PhysParams):
        S = self.surface_S(PhysParams)
        u, w = self._u[i0:i0 + n], self._weight[i0:i0 + n]
        a = gardner_sharp_albedo(S, self.carbon, u, self.tau[i0:i0 + n, None])
        wsum = w.sum(axis=1)
        return np.where(wsum > 0, (a * w).sum(axis=1) / np.where(wsum > 0, wsum, 1.0), a.mean(axis=1))
### end GardnerSharpAlbedo ###
##############################


class GreuellKonzelmannAlbedo(AlbedoModel):
    '''
    Greuell and Konzelmann (1994) albedo that depends on the density of the surface snow (their SOMARS
    model; equation 1 of Gardner and Sharp, 2010): the albedo changes linearly with the density from the
    albedo of fresh snow, alpha_s, at the density rho_s, to the albedo of ice, alpha_i, at the density rho_i:

        alpha = alpha_i + (rho_ss - rho_i) (alpha_s - alpha_i) / (rho_s - rho_i)

    with rho_ss the density of the surface node. The densities outside of the range from rho_s to rho_i are
    set to its ends, so that the albedo is between alpha_i and alpha_s.

    The configuration keys are albedo_snow (alpha_s) and albedo_ice (alpha_i), which have no default: they are
    constants "determined as part of the model tuning" (Gardner and Sharp, 2010), and the values depend on the
    application. albedo_rho_snow (rho_s, default 300 kg m-3) and albedo_rho_ice (rho_i, default 910 kg m-3) are
    the usual values of the densities of fresh snow and ice. The cloud dependence of Greuell and Konzelmann
    (1994) (alpha_cloud = alpha_clear + 0.05 (n_c - 0.5), equation 6 of Gardner and Sharp, 2010) is not included.
    '''

    required_config = ('albedo_snow', 'albedo_ice')

    def __init__(self, config, forcing):
        super().__init__(config, forcing)
        self.alpha_s = float(config['albedo_snow'])
        self.alpha_i = float(config['albedo_ice'])
        self.rho_s = float(config.get('albedo_rho_snow', 300.0))
        self.rho_i = float(config.get('albedo_rho_ice', 910.0))
        if not self.rho_s < self.rho_i:
            raise ValueError("albedo_rho_snow must be smaller than albedo_rho_ice.")

    def step(self, iii, i0, n, PhysParams):
        weight = np.clip((PhysParams['rho'][0] - self.rho_i) / (self.rho_s - self.rho_i), 0.0, 1.0)   # 1 at fresh snow, 0 at ice
        return np.full(n, self.alpha_i + weight * (self.alpha_s - self.alpha_i))
### end GreuellKonzelmannAlbedo ###
###################################


ALBEDO_METHODS = {
    'prescribed': PrescribedAlbedo,
    'gardner_sharp': GardnerSharpAlbedo,
    'greuell_konzelmann': GreuellKonzelmannAlbedo,
}


def make_albedo_model(config, forcing):
    '''
    Build the albedo model that the configuration key 'albedo_method' selects (default 'prescribed').

    :param config: the configuration dictionary
    :param forcing: dictionary of the SEB forcing arrays (see AlbedoModel)
    '''
    method = config.get('albedo_method', 'prescribed')
    if method not in ALBEDO_METHODS:
        raise ValueError(f"Unknown albedo_method: '{method}'. Options: {', '.join(repr(k) for k in ALBEDO_METHODS)}.")
    model = ALBEDO_METHODS[method]
    missing = [k for k in model.required_forcing if k not in forcing]
    if missing:
        raise KeyError(f"albedo_method '{method}' needs {missing} in the SEB forcing (the columns of the dataframe).")
    missing = [k for k in model.required_config if k not in config]
    if missing:
        raise KeyError(f"albedo_method '{method}' needs {missing} in the configuration file.")
    return model(config, forcing)
### end make_albedo_model ###
#############################
