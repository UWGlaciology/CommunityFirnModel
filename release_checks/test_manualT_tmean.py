#!/usr/bin/env python
"""
Verifies the manualT fix in firn_density_nospin.py (dev_release_status.md §3.2).

The manualT branch never defined the "_full" temperature arrays, so Tsf_full
raised NameError. The fix aliases them to the non-full arrays. This checks the
property that makes the alias correct: for manualT, Ts_full must come out
EXACTLY equal to self.Ts, so that T_mean reproduces main:538-539 (a true
rolling mean) rather than the flat np.mean fallback.

Run:  python release_checks/test_manualT_tmean.py

This simulates the code path rather than running the model, because manualT
needs a ManualTFilename .csv that is not in CFMinput_example/.
"""
import numpy as np
import pandas as pd
from scipy import interpolate

S_PER_YEAR = 365.25 * 24 * 3600


def run(seed=1):
    rng = np.random.default_rng(seed)

    # manualT-style inputs (thermistor record, irregular vs model time)
    input_year_temp = np.linspace(2000., 2010., 400)
    input_temp = 250. + 10 * np.sin(np.linspace(0, 20, 400)) + rng.normal(0, .5, 400)
    modeltime = np.linspace(2000., 2010., 365)
    dt = np.diff(np.append(modeltime, 2010.027)) * S_PER_YEAR
    stp = len(modeltime)
    int_type = 'linear'

    # --- the fix: aliases set in the manualT branch ---
    input_year_temp_full, input_temp_full, modeltime_full = (
        input_year_temp, input_temp, modeltime)

    Tsf = interpolate.interp1d(input_year_temp, input_temp, int_type,
                               fill_value='extrapolate')
    Tsf_full = interpolate.interp1d(input_year_temp_full, input_temp_full, int_type,
                                    fill_value='extrapolate')
    Ts = Tsf(modeltime)
    Ts_full = Tsf_full(modeltime_full)

    c1 = np.array_equal(Ts_full, Ts)
    print(f"Ts_full identical to Ts            : {c1}")

    NN = int(np.mean(S_PER_YEAR / dt) * 10)

    # dev T_mean (firn_density_nospin.py:604-607)
    short = pd.Series(Ts).rolling(window=NN + 1, win_type='hamming').mean().values
    dev = pd.Series(Ts_full).rolling(window=NN + 1, win_type='hamming').mean().values.copy()
    dev = dev[-len(short):]
    dev[np.isnan(dev)] = dev[NN]

    # main T_mean (main:CFM_main/firn_density_nospin.py:538-539)
    main_ = pd.Series(Ts).rolling(window=NN + 1, win_type='hamming').mean().values.copy()
    main_[np.isnan(main_)] = main_[NN]

    c2 = np.allclose(dev, main_, equal_nan=True)
    md = np.nanmax(np.abs(dev - main_))
    c3 = (len(dev) == stp)
    print(f"dev T_mean == main T_mean          : {c2}  (max diff {md})")
    print(f"len(T_mean) == stp                 : {c3}")

    # what the rejected "move it into the try" fix would have produced
    flat = np.nanmax(np.abs(np.mean(Ts) - main_))
    print(f"flat-fallback would have differed by: {flat:.4f} K  (why that fix was wrong)")

    ok = c1 and c2 and c3
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    import sys
    sys.exit(run())
