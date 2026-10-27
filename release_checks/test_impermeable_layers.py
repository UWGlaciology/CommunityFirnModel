#!/usr/bin/env python
"""
Differential test: melt.py::_find_impermeable_layers vs. the pre-3.3 (main
branch) inline implementations it replaced.

Written 2026-09-15 while fixing the three release blockers in
dev_release_status.md §3.1 and §3.3. Kept so the verification is reproducible
rather than existing only in a chat transcript.

Run from CFM_main/:   python ../release_checks/test_impermeable_layers.py

Expected output: UNEXPLAINED == 0 for both schemes, and
"bottom node always in imp: True" for bucket. Any nonzero UNEXPLAINED means a
real behavior divergence from main that needs investigating.

The only tolerated difference is a lens whose thickness is *exactly* ThickImp:
main's darcy used a strict `>` while main's bucket used `>=`, and the unified
helper uses `>=`. That makes the two schemes agree, and is an intentional
change recorded in dev_release_status.md §3.6 item 6.
"""
import numpy as np
import sys, os

sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'CFM_main'))
from melt import _find_impermeable_layers as H


def main_bucket(rho, dz, RhoImp, ThickImp, DownToIce):
    """main:CFM_main/melt.py lines 243-272, transcribed verbatim."""
    nnd = len(rho)
    if DownToIce:
        if np.any(rho < RhoImp):
            imp = np.arange(np.where(rho < RhoImp)[0][-1] + 1, nnd, 1).astype(int)
        else:
            imp = np.arange(0, nnd, 1).astype(int)
    else:
        if ThickImp > 0:
            lens0 = np.array([ii for ii in range(1, nnd)
                              if (rho[ii] >= RhoImp and rho[ii - 1] < RhoImp)])
            lens1 = np.array([ii for ii in range(0, nnd - 1)
                              if (rho[ii] >= RhoImp and rho[ii + 1] < RhoImp)])
            if rho[0] >= RhoImp:
                lens0 = np.append(0, lens0).astype(int)
            lens1 = np.append(lens1, nnd - 1).astype(int)
            imp = np.array([])
            for ii in range(len(lens0)):
                lensdz = (sum(dz[lens0[ii]:lens1[ii] + 1]))
                if lensdz >= ThickImp or (ii == len(lens0) - 1):   # bucket used >=
                    imp = np.append(imp, np.arange(lens0[ii], lens1[ii] + 1, 1)).astype(int)
        else:
            imp = np.where(rho >= RhoImp)[0]
    imp = np.asarray(imp).astype(int)
    if len(imp) == 0:
        imp = np.array([nnd - 1])
    else:
        imp = np.append(imp, nnd - 1)          # unconditional, the key behavior
    return np.unique(imp)   # set-equivalent: stcap[imp]=0 is order/dup insensitive


def main_darcy(rho, dz, RhoImp, ThickImp):
    """main:CFM_main/melt.py lines 687-702, transcribed verbatim."""
    ncv = len(rho)
    l0 = np.array([i for i in range(1, ncv) if rho[i] >= RhoImp and rho[i - 1] < RhoImp])
    l1 = np.array([i for i in range(0, ncv - 1) if rho[i] >= RhoImp and rho[i + 1] < RhoImp])
    if rho[0] >= RhoImp:
        l0 = np.append(0, l0).astype(int)
    l1 = np.append(l1, ncv - 1).astype(int)
    imp, imptop = np.array([]), np.array([])
    for ii in range(len(l0)):
        if sum(dz[l0[ii]:l1[ii] + 1]) > ThickImp or (ii == len(l0) - 1):  # darcy used >
            imp = np.append(imp, np.arange(l0[ii], l1[ii] + 1, 1)).astype(int)
            imptop = np.append(imptop, l0[ii]).astype(int)
    if imp.size == 0:
        imp = np.array([ncv - 1])
        imptop = imp.copy()
    return np.unique(imp.astype(int)), np.sort(imptop.astype(int))


def _exact_thickness_case(rho, dz, RhoImp, ThickImp):
    """True if any lens is exactly ThickImp thick (the tolerated >=/> boundary)."""
    n = len(rho)
    l0 = np.array([i for i in range(1, n) if rho[i] >= RhoImp and rho[i - 1] < RhoImp])
    l1 = np.array([i for i in range(0, n - 1) if rho[i] >= RhoImp and rho[i + 1] < RhoImp])
    if rho[0] >= RhoImp:
        l0 = np.append(0, l0).astype(int)
    l1 = np.append(l1, n - 1).astype(int)
    return any(abs(sum(dz[l0[i]:l1[i] + 1]) - ThickImp) < 1e-12 for i in range(len(l0)))


def run(trials=6000, seed=7):
    rng = np.random.default_rng(seed)
    bmis = bexact = dmis = dexact = 0
    bottom_ok = True
    for _ in range(trials):
        n = rng.integers(5, 40)
        rho = rng.choice([400., 700., 850., 880., 900., 917.], size=n)
        dz = rng.choice([0.1, 0.25, 0.5, 0.6], size=n)

        # --- bucket: all three modes, force_bottom_impermeable=True ---
        for ThickImp, DownToIce in ((0.5, False), (0.0, False), (0.5, True)):
            RhoImp = 830.
            got, _ = H(rho, dz, RhoImp, ThickImp=ThickImp, DownToIce=DownToIce,
                       force_bottom_impermeable=True)
            exp = main_bucket(rho, dz, RhoImp, ThickImp, DownToIce)
            if not np.array_equal(got, exp):
                bmis += 1
                if ThickImp > 0 and not DownToIce and _exact_thickness_case(rho, dz, RhoImp, ThickImp):
                    bexact += 1
                elif bmis - bexact <= 3:
                    print("BUCKET unexplained:", ThickImp, DownToIce,
                          "\n  rho", rho, "\n  dz ", dz, "\n  got", got, "\n  exp", exp)
            if n - 1 not in got:
                bottom_ok = False

        # --- darcy: must remain byte-identical to main's darcy ---
        RhoImp = 873.
        g = H(rho, dz, RhoImp, ThickImp=0.5, DownToIce=False)
        e = main_darcy(rho, dz, RhoImp, 0.5)
        if not (np.array_equal(g[0], e[0]) and np.array_equal(g[1], e[1])):
            dmis += 1
            if _exact_thickness_case(rho, dz, RhoImp, 0.5):
                dexact += 1
            elif dmis - dexact <= 3:
                print("DARCY unexplained:\n  rho", rho, "\n  dz ", dz,
                      "\n  got", g, "\n  exp", e)

    print(f"trials={trials}")
    print(f"BUCKET: mismatches={bmis:5d}  explained={bexact:5d}  UNEXPLAINED={bmis-bexact}")
    print(f"        bottom node always in imp: {bottom_ok}")
    print(f"DARCY : mismatches={dmis:5d}  explained={dexact:5d}  UNEXPLAINED={dmis-dexact}")
    ok = (bmis - bexact == 0) and (dmis - dexact == 0) and bottom_ok
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(run())
