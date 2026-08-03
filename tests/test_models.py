"""Sanity and benchmark tests for src/models.py.

Benchmarks marked [S52] compare against Sondheimer, Adv. Phys. 1, 1
(1952), Table 1 (film rho/rho0 for p = 0 and p = 1/2), read from the
paper itself.
"""
import sys, os
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import models as M

# ------------------------- FS film ------------------------------------

# [S52] Table 1: kappa, rho/rho0 (p=0), rho/rho0 (p=1/2)
SOND_TABLE1 = [
    (0.001, 182.0, 73.5),
    (0.002, 100.4, 41.5),
    (0.005, 46.6, 20.0),
    (0.01, 26.5, 11.8),
    (0.02, 15.3, 7.1),
    (0.05, 7.69, 3.87),
    (0.1, 4.72, 2.62),
    (0.2, 3.00, 1.91),
    (0.5, 1.90, 1.402),
    (1.0, 1.462, 1.206),
    (2.0, 1.221, 1.102),
    (5.0, 1.081, 1.039),
    (10.0, 1.0390, 1.0191),
    (20.0, 1.0191, 1.0095),
    (50.0, 1.0076, 1.0038),
    (100.0, 1.0038, 1.0019),
]


@pytest.mark.parametrize("kappa,r0,rhalf", SOND_TABLE1)
def test_fs_film_vs_sondheimer_table1(kappa, r0, rhalf):
    """Reproduce Sondheimer (1952) Table 1.

    Note: in the crossover regime kappa ~ 0.05-0.2 the 1952 hand-computed
    table is itself inaccurate by up to ~3% (three independent numerical
    routes -- Gauss-Legendre in 1/t, adaptive quadrature in t, and the
    specularity series identity -- agree with each other to 6+ digits but
    deviate from the printed values there; e.g. kappa=0.1, p=0 gives
    4.7817 vs the printed 4.72).  The kappa=0.5, p=0 entry (printed 1.90,
    exact 1.9161) is off by 0.85% for the same reason.  Outside that
    window agreement is <0.3%.
    """
    tol = 0.04 if 0.02 <= kappa <= 0.5 else 7e-3
    got0 = M.fs_film(kappa, 1.0, p=0.0)
    goth = M.fs_film(kappa, 1.0, p=0.5)
    assert abs(got0 - r0) / r0 < tol
    assert abs(goth - rhalf) / rhalf < tol


def test_fs_film_bulk_recovery():
    assert abs(M.fs_film(1e5, 40.0, 0.0) - 1.0) < 1e-3
    assert abs(M.fs_film(50.0, 40.0, 1.0) - 1.0) < 1e-12


def test_fs_film_thick_asymptote():
    # rho/rho0 -> 1 + 3(1-p)/(8 kappa)
    for p in (0.0, 0.3, 0.6):
        kappa = 300.0
        expect = 1.0 + 3.0 * (1 - p) / (8 * kappa)
        got = M.fs_film(kappa, 1.0, p)
        assert abs(got - expect) < 4e-4


def test_fs_film_p_series_identity():
    for kappa in (0.3, 1.0, 5.0):
        for p in (0.3, 0.6):
            direct = M.fs_film(kappa, 1.0, p)
            series = M.film_p_series(kappa, 1.0, p)
            assert abs(direct - series) / direct < 2e-4


# ------------------------- Chambers wire ------------------------------

def test_wire_bulk_recovery():
    assert abs(M.chambers_wire_p0(4000.0, 4000.0, 1.0) - 1.0) < 2e-3


def test_wire_thick_square_asymptote():
    # Sondheimer Eq. (32): square wire, sigma0/sigma = 1 + (3/4)(1-p)/kappa
    kappa = 150.0
    got = M.chambers_wire_p0(kappa, kappa, 1.0)
    expect = 1.0 + 0.75 / kappa
    assert abs(got - expect) / (expect - 1.0) < 0.03


def test_wire_to_film_convergence():
    # w >> h: wire result must approach the film of thickness h
    lam = 40.0
    h = 230.0
    film = M.fs_film(h, lam, 0.0)
    wide = M.chambers_wire_p0(80 * h, h, lam)
    assert abs(wide - film) / (film - 1.0) < 0.03


def test_wire_p_series_monotone_in_p():
    vals = [M.wire_p_series(100.0, 230.0, 40.0, p) for p in (0.0, 0.3, 0.6, 0.9)]
    assert all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
    assert vals[-1] > 1.0


# ------------------------- Mayadas-Shatzkes ---------------------------

def test_ms_limits():
    assert abs(M.ms_rho_ratio(40.0, 1e9, 0.3) - 1.0) < 1e-6      # d -> inf
    assert abs(float(M.ms_sigma_ratio(0.0)) - 1.0) < 1e-12       # R -> 0
    big = M.ms_rho_ratio(40.0, 40.0, 0.999)                      # R -> 1
    assert big > 1e2


def test_ms_small_alpha_expansion():
    a = 1e-3
    exact = 3.0 * (1.0 / 3 - a / 2 + a * a - a ** 3 * np.log1p(1 / a))
    assert abs(float(M.ms_sigma_ratio(a)) - exact) < 1e-9


def test_ms_large_alpha():
    # f(alpha) -> 3/(4 alpha) for alpha >> 1
    a = 2000.0
    assert abs(float(M.ms_sigma_ratio(a)) - 0.75 / a) / (0.75 / a) < 2e-3


# ------------------------- combined & wrapper -------------------------

def test_combined_reduces_to_parts():
    lam, w, h = 40.0, 100.0, 230.0
    fs_only = M.wire_p_series(w, h, lam, 0.6)
    comb_no_gb = M.combined_wire(w, h, lam, 0.6, R=1e-9, d=1e9)
    assert abs(comb_no_gb - fs_only) < 1e-6
    ms_only = float(M.ms_rho_ratio(lam, w, 0.5))
    comb_no_fs = M.combined_wire(w, h, lam, 1.0, R=0.5, d=w)
    assert abs(comb_no_fs - ms_only) < 1e-9


def test_wrapper_matches_direct():
    lam, h = 40.0, 230.0
    for w in (50.0, 100.0, 400.0):
        ev = M.WireSizeEffect(h / w)
        direct = M.wire_p_series(w, h, lam, 0.6)
        fast = ev.rho_ratio(w, lam, 0.6)
        assert abs(fast - direct) / direct < 3e-3


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "--tb=short"]))


def test_rect_wire_interp_matches_direct():
    wi = M.RectWireInterp()  # no cache: small fresh table
    for (w, h) in ((60.0, 230.0), (300.0, 230.0)):
        for p in (0.0, 0.6):
            direct = M.wire_p_series(w, h, 40.0, p)
            fast = wi.rho_ratio(w, h, 40.0, p)
            assert abs(fast - direct) / direct < 3e-3
