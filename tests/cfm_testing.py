'''
cfm_testing.py
==============
Helpers for the CFM test suite. conftest.py exposes them as fixtures.

- SyntheticColumn: a small object with the state attributes that the meltwater
  schemes (melt.bucket, melt.darcyscheme) read. A scheme test calls, for example,
  bucket(column, 0) without a full model run.
- ShortRun: a short model run, built from example_df.json with the forcing cut
  to a few years.

The CFM modules are flat scripts in CFM_main/, not an installable package. As a
temporary workaround, this module puts CFM_main/ on sys.path so that they import.
'''

import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CFM_MAIN = Path(__file__).resolve().parent.parent / 'CFM_main'
if str(CFM_MAIN) not in sys.path:
    sys.path.insert(0, str(CFM_MAIN))

from constants import S_PER_YEAR  # noqa: E402

EXAMPLE_CONFIG = CFM_MAIN / 'example_df.json'


def example_config():
    '''Return the example_df.json config as a dict.'''
    with open(EXAMPLE_CONFIG) as f:
        return json.load(f)


def _values_per_node(value, z):
    '''Expand a scalar, an array, or a function of depth to one value per node.'''
    if callable(value):
        return np.asarray(value(z), dtype=float)
    return np.broadcast_to(np.asarray(value, dtype=float), z.shape).copy()


class SyntheticColumn:
    '''
    A firn column with the state attributes that the meltwater schemes read.

    STATE_ATTRIBUTES lists them. Units follow FirnDensityNoSpin:

    Per node (top node = surface):
        rho [kg m-3], Tz [K], dz [m], z [m, top edge of node], mass [kg m-2],
        LWC [m w.e.], age [s], r2 [m2], bdot_mean [m ie yr-1], Dcon, dzn,
        mass_sum [kg m-2]
    Per time step:
        dt [s], snowmeltSec and rainSec [m ie per step / S_PER_YEAR],
        modeltime [decimal year]
    Other:
        c (config dict), compboxes, doublegrid, gridtrack
    '''

    STATE_ATTRIBUTES = ('rho', 'Tz', 'dz', 'z', 'mass', 'LWC', 'age', 'r2',
                        'bdot_mean', 'Dcon', 'dzn', 'mass_sum', 'dt', 'snowmeltSec',
                        'rainSec', 'modeltime', 'c', 'compboxes', 'doublegrid',
                        'gridtrack')

    def __init__(self, n=50, rho_profile=400.0, T_profile=263.15, dz=0.1, LWC=0.0,
                 r2=1.0e-8, melt=0.0, rain=0.0, dt=S_PER_YEAR / 12, steps=1,
                 config=None):
        '''
        :param n: number of nodes
        :param rho_profile, T_profile, LWC, r2: a scalar, an array of length n, or a
            function of node depth z (top edge) that returns an array of length n
        :param dz: node thickness [m], a scalar or an array of length n
        :param melt, rain: surface melt and rain for each step [m ie], a scalar or
            an array of length steps
        :param dt: step length [s]
        :param steps: number of time steps (the length of the per-step arrays)
        :param config: overrides for self.c, which starts from example_df.json
        '''
        self.dz = np.broadcast_to(np.asarray(dz, dtype=float), (n,)).copy()
        self.z = np.concatenate(([0.0], np.cumsum(self.dz)[:-1]))

        self.rho = _values_per_node(rho_profile, self.z)
        self.Tz = _values_per_node(T_profile, self.z)
        self.LWC = _values_per_node(LWC, self.z)
        self.r2 = _values_per_node(r2, self.z)
        self.mass = self.rho * self.dz
        self.mass_sum = np.cumsum(self.mass)
        self.age = np.arange(n) * dt
        self.bdot_mean = np.full(n, 0.3)
        self.Dcon = np.flipud(np.arange(-n, 0))

        self.compboxes = n
        self.dzn = self.dz[0:self.compboxes]
        self.doublegrid = False
        self.gridtrack = None

        self.dt = np.full(steps, float(dt))
        self.snowmeltSec = np.broadcast_to(np.asarray(melt, dtype=float), (steps,)) / S_PER_YEAR
        self.rainSec = np.broadcast_to(np.asarray(rain, dtype=float), (steps,)) / S_PER_YEAR
        self.modeltime = 2000.0 + np.arange(steps) * dt / S_PER_YEAR

        self.c = example_config()
        self.c.update(config or {})

    def water_mass(self):
        '''Total liquid water in the column [kg m-2].'''
        return float(np.sum(self.LWC) * 1000.0)


class ShortRun:
    '''
    A short transient run built from example_df.json.

    The constructor:
      1. Reads example_df.json.
      2. Cuts the DataFrame forcing to the years [start_year, start_year + years - 1].
      3. Sets a constant surface density, NewSpin, and resultsFolder = folder.
      4. Applies the overrides.
      5. Writes the config to folder/config.json and builds climateTS with
         RCMpkl_to_spin.makeSpinFiles. The spin-up repeats the same years.

    Call build() to get a FirnDensityNoSpin (the spin-up runs in its __init__),
    or run() to also call time_evolve().
    '''

    def __init__(self, folder, years=2, start_year=1980, num_reps=None, **overrides):
        import RCMpkl_to_spin as RCM

        self.folder = Path(folder)
        self.folder.mkdir(parents=True, exist_ok=True)

        c = example_config()
        c.update({
            'InputFileFolder': str(CFM_MAIN / c['InputFileFolder']),
            'resultsFolder': str(self.folder),
            'NewSpin': True,
            'variable_srho': False,
            'srho_type': 'userinput',
            'TWriteStart': float(start_year),
            'spinUpdateDate': float(start_year),
        })
        c.update(overrides)
        self.config = c
        self.config_path = self.folder / 'config.json'
        with open(self.config_path, 'w') as f:
            json.dump(c, f, indent=1)

        if c['input_type'] != 'dataframe':
            raise ValueError('ShortRun supports input_type "dataframe" only')
        end_year = start_year + years - 1
        df = pd.read_pickle(os.path.join(c['InputFileFolder'], c['DFfile']))
        df = df.loc[f'{start_year}':f'{end_year}'].copy()  # makeSpinFiles changes df in place
        (self.climateTS, self.stepsperyear, _, _, _,
         self.SEBfluxes) = RCM.makeSpinFiles(
            df, timeres=c['DFresample'], melt=c['MELT'],
            desired_depth=c['H'] - c['HbaseSpin'],
            spin_date_st=float(start_year), spin_date_end=float(end_year),
            num_reps=num_reps)

    @property
    def results_path(self):
        return self.folder / self.config['resultsFileName']

    def build(self):
        '''Return a new FirnDensityNoSpin. Its __init__ runs the spin-up.'''
        import solver
        from firn_density_nospin import FirnDensityNoSpin
        solver._diag_reset()  # module-level list; it is not reset between runs
        return FirnDensityNoSpin(str(self.config_path), climateTS=self.climateTS,
                                 NewSpin=True, SEBfluxes=self.SEBfluxes)

    def run(self):
        '''Build the model, run time_evolve(), and return the model.'''
        model = self.build()
        model.time_evolve()
        # time_evolve writes the solver diagnostics to the working folder.
        diagnostics_file = Path(f'diagnostics_{self.config["meltwater_solver"]}.csv')
        if diagnostics_file.exists():
            shutil.move(diagnostics_file, self.folder / diagnostics_file.name)
        return model
