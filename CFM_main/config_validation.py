#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
config_validation.py
====================

Validation of the .json run-configuration dictionary for the Community Firn
Model. The goal is to fail fast with a single, clear message when a
configuration is missing a required key or uses an invalid option value,
rather than surfacing a cryptic KeyError deep inside a model run.

The validator is called once from ``FirnDensityNoSpin.__init__`` immediately
after the config is loaded; because that dictionary is passed on to
``FirnDensitySpin``, a single call covers the whole run.

Behavior:
  * ERROR (raises ConfigValidationError) on missing required keys, unmet
    conditional requirements, and invalid option values.
  * WARNING (printed, run continues) on unrecognized keys, which are usually
    typos.

Allowed option values are held here (authoritative, matching the densification
schemes in physics.py and docs/running/json.rst). They are deliberately NOT
read from the ``*_options`` helper arrays in example.json, because those have
drifted from the code (e.g. physRho_options omits GSFC2020, which is the scheme
used by run_CFM_example.py).

This module uses only the standard library so it can be imported anywhere
without creating circular imports.
'''

import json
import sys


class ConfigValidationError(Exception):
    '''Raised when the run configuration has one or more fatal problems.'''
    pass


# Keys that are read without a default guard anywhere in the model and so must
# always be present.
REQUIRED_KEYS = [
    'physRho',
    'isoDiff',
    'stpsPerYear',
    'H',
    'HbaseSpin',
    'resultsFolder',
    'InputFileFolder',
    'physGrain',
    'heatDiff',
    'outputs',
    'bdot_type',
    'spinup_climate_type',
]


# Authoritative allowed values for keys that select among a fixed set of
# options. For scalar keys the value must be in the set; for the list-valued
# keys (outputs, iso) each element must be in the set.
ALLOWED_VALUES = {
    'physRho': {
        'HLdynamic', 'HLSigfus', 'Li2004', 'Li2011', 'Li2015',
        'Helsen2008', 'Arthern2010S', 'Arthern2010T', 'Simonsen2013',
        'Ligtenberg2011', 'Barnola1991', 'Morris2014', 'KuipersMunneke2015',
        'Brils2022', 'Veldhuijsen2023', 'Goujon2003', 'Breant2017',
        'Crocus', 'GSFC2020',
    },
    'bdot_type': {'instant', 'mean', 'stress'},
    'conductivity': {
        'Schwander', 'Yen_fixed', 'Yen_var', 'Anderson', 'Yen_b', 'Sturm',
        'VanDusen', 'Schwerdtfeger', 'Riche', 'Jiawen', 'mix',
        'Calonne2011', 'Calonne2019',
    },
    'GrGrowPhysics': {'Arthern', 'Katsushima'},
    'input_type': {'csv', 'dataframe'},
    'int_type': {'nearest', 'linear'},
    'liquid': {
        'bucket', 'darcy', 'resingledomain', 'prefsnowpack',
        # accepted aliases for the Verjans bucket variant
        'bucketVV', 'percolation_bucket',
    },
    'spinup_climate_type': {'mean', 'initial'},
    'srho_type': {'userinput', 'param', 'noise', 'Brils22'},
    'timesetup': {'exact', 'interp', 'retmip'},
    'SeasonalThemi': {'north', 'south'},
    # list-valued keys (checked per element)
    'outputs': {
        'density', 'depth', 'temperature', 'age', 'Dcon', 'bdot_mean',
        'climate', 'compaction', 'grainsize', 'temp_Hx', 'isotopes', 'BCO',
        'DIPc', 'DIP', 'LWC', 'gasses', 'PLWC_mem', 'viscosity', 'runoff',
        'refrozen', 'meltoutputs',
    },
    'iso': {'d18O', 'dD', 'NoDiffusion'},
    'meltwater_solver': {'enthalpy', 'ahc', 'decp', 'ncz', 'legacy'},
}

# Keys that hold lists whose elements are validated against ALLOWED_VALUES.
LIST_VALUED_KEYS = {'outputs', 'iso'}


# Keys the model no longer reads, mapped to the message shown when one is
# found. These are recognized (so they are not reported as typos) but flagged,
# because silently ignoring them can change the physics of a run. Handling of
# the value itself lives at the point of use -- e.g. FirnDensityNoSpin.__init__
# migrates a lone 'LWC_heat' onto meltwater_solver='legacy'.
# claude, 26/09/11
DEPRECATED_KEYS = {
    'LWC_heat': ("replaced by 'meltwater_solver'; if 'meltwater_solver' is "
                 "absent the run falls back to the 'legacy' solver to "
                 "reproduce the old numerics, otherwise 'LWC_heat' is ignored"),
}


# Every configuration key the model recognizes. Used only to warn about
# unknown keys (likely typos), so it does not need to be exhaustive to be
# safe -- a missing entry only costs a spurious warning. Built from the keys
# documented in docs/running/json.rst, the example .json files, and keys that
# are injected at runtime by scripts such as run_CFM_example.py.
KNOWN_KEYS = {
    # core / IO
    'InputFileFolder', 'InputFileNameTemp', 'InputFileNamebdot',
    'InputFileNamemelt', 'InputFileNameRain', 'InputFileNameSublim',
    'InputFileNameIso', 'InputFileNamerho', 'InputFileNameStrain',
    'ManualTFilename', 'forcingFileName', 'resultsFolder', 'resultsFileName',
    'spinFileName', 'initfirnFile', 'initprofile', 'input_type',
    'DFresample', 'DFfile',
    # physics selection
    'physRho', 'MELT', 'RAIN', 'ReehCorrectedT', 'FirnAir', 'AirConfigName',
    'SUBLIM', 'bdm_sublim',
    # time / writing
    'TWriteInt', 'TWriteStart', 'int_type', 'timesetup',
    # seasonal cycle
    'SeasonalTcycle', 'SeasonalThemi', 'coreless', 'TAmp',
    # grain size
    'physGrain', 'calcGrainSize', 'GrGrowPhysics', 'r2s0',
    # heat
    'heatDiff', 'conductivity', 'LWCheat', 'LWCcorr_subdt',
    'correct_therm_prop',
    # surface density
    'variable_srho', 'srho_type', 'rhos0',
    # spin up / domain
    'AutoSpinUpTime', 'yearSpin', 'stpsPerYear', 'stpsPerYearSpin',
    'H', 'HbaseSpin', 'D_surf', 'bdot_type', 'iceblock', 'iceblock_rho',
    # outputs / grid
    'grid_outputs', 'grid_output_res', 'grid_output_max_depth', 'isoDiff', 'iso', 'spacewriteint',
    'outputs', 'output_bits', 'truncate_outputs',
    'doublegrid', 'nodestocombine', 'multnodestocombine', 'grid1bottom',
    'grid2bottom',
    # strain / divergence
    'strain', 'horizontal_divergence', 'strain_softening', 'residual_strain',
    'tuning_bias_correction', 'InputFileNamedudx',
    # climate
    'spinup_climate_type', 'manual_climate', 'deepT', 'bdot_long',
    'manual_iceout', 'iceout', 'QMorris', 'MQ', 'THist',
    # melt / liquid water
    'liquid', 'merging', 'merge_min', 'LWCcorrect', 'ColeouLesaffre',
    'IrrVal', 'RhoImp', 'ThickImp', 'DownToIce', 'Ponding', 'DirectRunoff',
    'RunoffZuoOerlemans', 'Slope', 'keep_firnthickness', 'meltwater_solver', 'meltwater_solver_options',
    # SEB
    'SEB', 'SEB_TL_thick', 'albedo_factor',
    # stage zero snow
    'stage_zero', 'snow_model', 's_zero_rho',
    # misc / measurement
    'manualT', 'no_densification', 'rad_pen', 'site_pressure',
    'spinUpdate', 'spinUpdateDate', 'DIPhorizon', 'NewSpin',
    # run metadata / site (injected by run scripts)
    'runID', 'lat_int', 'lon_int', 'lat_val', 'lon_val', 'runid', 'x_val', 'y_val', 'runloc', 'quad', 'x_int', 'y_int', 'resultspath'
}


# Conditional requirements: (predicate, [keys required when predicate True],
# human-readable reason). Predicates read from the config dict with .get so a
# missing controlling key is treated as its default.
CONDITIONAL_REQUIREMENTS = [
    (lambda c: c.get('input_type', 'csv') == 'dataframe',
     ['DFfile', 'DFresample'],
     "input_type is 'dataframe'"),
    (lambda c: c.get('MELT'),
     ['liquid'],
     "MELT is enabled"),
    (lambda c: (c.get('variable_srho') and c.get('srho_type') == 'userinput'
                and c.get('input_type', 'csv') == 'csv'),
     ['InputFileNamerho'],
     "variable_srho is True with srho_type 'userinput' (csv input)"),
    (lambda c: c.get('manualT'),
     ['ManualTFilename'],
     "manualT is enabled"),
    (lambda c: c.get('isoDiff'),
     ['iso'],
     "isoDiff is enabled"),
    (lambda c: c.get('isoDiff') and c.get('input_type', 'csv') == 'csv',
     ['InputFileNameIso'],
     "isoDiff is enabled (csv input)"),
    (lambda c: c.get('FirnAir'),
     ['AirConfigName'],
     "FirnAir is enabled"),
]


def _is_helper_key(key):
    '''Helper arrays such as ``physRho_options`` are not real config keys.'''
    return key.endswith('_options') or key == 'isoOptions'


def validate_config(c, config_path=None):
    '''
    Validate a CFM run-configuration dictionary.

    :param c: the configuration dictionary (parsed from the .json file).
    :param config_path: optional path to the source .json, used only to make
        error/warning messages more informative.

    Prints a warning for each unrecognized key, then raises
    ConfigValidationError (once, listing every problem found) if any required
    key is missing, a conditional requirement is unmet, or an option value is
    invalid. Returns None when the configuration is valid.
    '''
    errors = []
    warnings = []

    # 1. Required keys.
    for key in REQUIRED_KEYS:
        if key not in c:
            errors.append("missing required key '{}'".format(key))

    # 2. Conditional requirements.
    for predicate, required, reason in CONDITIONAL_REQUIREMENTS:
        if predicate(c):
            for key in required:
                if key not in c:
                    errors.append(
                        "key '{}' is required when {}".format(key, reason))

    # 3 & 4. Allowed option values (scalar and list-valued).
    for key, allowed in ALLOWED_VALUES.items():
        if key not in c:
            continue
        value = c[key]
        allowed_str = ', '.join(repr(v) for v in sorted(allowed))
        if key in LIST_VALUED_KEYS:
            if not isinstance(value, (list, tuple)):
                errors.append(
                    "key '{}' should be a list; got {!r}".format(key, value))
                continue
            for element in value:
                if element not in allowed:
                    errors.append(
                        "invalid value {!r} in '{}'; allowed values are: {}"
                        .format(element, key, allowed_str))
        else:
            if value not in allowed:
                errors.append(
                    "invalid value {!r} for '{}'; allowed values are: {}"
                    .format(value, key, allowed_str))

    # 5. Deprecated and unknown keys -> warnings only.
    # claude, 26/09/11: deprecated keys are checked first so they are not
    # mislabeled as typos.
    for key in c:
        if _is_helper_key(key):
            continue
        if key in DEPRECATED_KEYS:
            warnings.append(
                "deprecated key '{}': {}".format(key, DEPRECATED_KEYS[key]))
        elif key not in KNOWN_KEYS:
            warnings.append(
                "unrecognized key '{}' (possible typo?); it will be ignored"
                .format(key))

    src = " in {}".format(config_path) if config_path else ""

    for w in warnings:
        print("CFM config WARNING{}: {}".format(src, w))

    if errors:
        lines = ["Invalid CFM configuration{}:".format(src)]
        for i, e in enumerate(errors, 1):
            lines.append("  {}. {}".format(i, e))
        lines.append("Please fix the configuration and re-run. See the "
                     "documentation for the .json configuration file.")
        raise ConfigValidationError("\n".join(lines))


if __name__ == '__main__':
    # Standalone linter: `python config_validation.py <config.json>`
    if len(sys.argv) < 2:
        print("usage: python config_validation.py <config.json>")
        sys.exit(2)
    path = sys.argv[1]
    with open(path, "r") as _f:
        _c = json.loads(_f.read())
    try:
        validate_config(_c, config_path=path)
    except ConfigValidationError as err:
        print(err)
        sys.exit(1)
    print("{}: configuration is valid.".format(path))
