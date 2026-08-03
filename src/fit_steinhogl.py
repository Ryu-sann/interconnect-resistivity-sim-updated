"""Fit the combined FS+MS model to the digitized Steinhoegl (2002) Fig. 3
data and compare FS-only / MS-only / combined variants.

Fixed by the paper (literature-constrained, not free):
    rho0 = 1.90 uOhm.cm (their fit of the bulk term at 295 K)
    lambda = 40 nm
    h = 230 nm, d = min(w, h)   (grain size tracks width up to the height)
Free within literature ranges:
    p in [0, 1]      (specularity; Cu literature ~0-0.6)
    R in [0.1, 0.9]  (GB reflection; Cu literature 0.24-0.8)

Outputs figures/fit_steinhogl.png and prints a fit report.
Run from the repository root:  python src/fit_steinhogl.py
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import models as M

HERE = os.path.dirname(__file__)
ROOT = os.path.abspath(os.path.join(HERE, ".."))

RHO0, LAM, H = 1.90, 40.0, 230.0
WIRE = M.RectWireInterp(cache_path=os.path.join(ROOT, "data", "cache",
                                                "chambers_p0_table.npz"))


def load_points():
    f = os.path.join(ROOT, "data", "steinhogl2002_fig3_points.csv")
    rows = []
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("w_nm"):
                continue
            rows.append([float(x) for x in line.split(",")])
    arr = np.array(rows)
    w, rho = arr[:, 0], arr[:, 1]
    sig = 0.5 * (arr[:, 2] + arr[:, 3])
    return w, rho, sig


def model_combined(w_arr, p, R):
    return np.array([RHO0 * WIRE.combined(w, H, LAM, p, R, min(w, H))
                     for w in w_arr])


def model_fs_only(w_arr, p):
    return np.array([RHO0 * WIRE.rho_ratio(w, H, LAM, p) for w in w_arr])


def model_ms_only(w_arr, R):
    return np.array([RHO0 * float(M.ms_rho_ratio(LAM, min(w, H), R))
                     for w in w_arr])


def chi2(res, sig):
    return float(np.sum((res / sig) ** 2))


def rmse(res):
    return float(np.sqrt(np.mean(res ** 2)))


def fit_grid(w, rho, sig):
    ps = np.linspace(0.0, 1.0, 51)
    Rs = np.linspace(0.05, 0.9, 69)
    C = np.empty((len(ps), len(Rs)))
    for i, p in enumerate(ps):
        for j, R in enumerate(Rs):
            C[i, j] = chi2(rho - model_combined(w, p, R), sig)
    i0, j0 = np.unravel_index(np.argmin(C), C.shape)
    # local quadratic refine
    from scipy.optimize import minimize
    r = minimize(lambda x: chi2(rho - model_combined(w, x[0], x[1]), sig),
                 [ps[i0], Rs[j0]], method="Nelder-Mead",
                 options=dict(xatol=1e-4, fatol=1e-8))
    return ps, Rs, C, r.x, r.fun


def main():
    w, rho, sig = load_points()
    ps, Rs, C, (p_fit, R_fit), c_min = fit_grid(w, rho, sig)
    best = model_combined(w, p_fit, R_fit)
    res = rho - best
    ndof = len(w) - 2
    print("\n=== combined FS+MS (additive, Steinhoegl Eq. 5) ===")
    print(f"  best fit: p = {p_fit:.3f}, R = {R_fit:.3f}")
    print(f"  chi2 = {c_min:.2f}  (ndof = {ndof}), chi2/ndof = {c_min/ndof:.2f}")
    print(f"  RMSE = {rmse(res):.4f} uOhm.cm")
    print(f"  paper values: p = 0.6, R = 0.50")
    paper = model_combined(w, 0.6, 0.50)
    print(f"  chi2 at paper values = {chi2(rho - paper, sig):.2f},"
          f" RMSE = {rmse(rho - paper):.4f}")

    # 1-sigma region: chi2 <= chi2_min + 2.30 (2 params)
    mask = C <= c_min + 2.30
    p_lo, p_hi = ps[mask.any(axis=1)].min(), ps[mask.any(axis=1)].max()
    R_lo, R_hi = Rs[mask.any(axis=0)].min(), Rs[mask.any(axis=0)].max()
    print(f"  1-sigma ranges: p in [{p_lo:.2f}, {p_hi:.2f}],"
          f" R in [{R_lo:.2f}, {R_hi:.2f}]")

    # multiplicative combination (model-form check):
    # rho = rho_MS(w) * (rho/rho0)_FS  instead of additive increments
    def model_mult(w_arr, p, R):
        out = np.empty_like(w_arr)
        for i, wi in enumerate(w_arr):
            d = min(wi, H)
            ms = float(M.ms_rho_ratio(LAM, d, R))
            out[i] = RHO0 * ms * WIRE.rho_ratio(wi, H, LAM, p)
        return out

    from scipy.optimize import minimize
    rm = minimize(lambda x: chi2(rho - model_mult(w, x[0], x[1]), sig),
                  [p_fit if p_fit > 0.05 else 0.3, R_fit],
                  method="Nelder-Mead")
    print("\n=== combined, multiplicative form (model-form check) ===")
    print(f"  best fit: p = {rm.x[0]:.3f}, R = {rm.x[1]:.3f}, "
          f"chi2 = {rm.fun:.2f}, "
          f"RMSE = {rmse(rho - model_mult(w, *rm.x)):.4f}")

    # single-mechanism fits
    from scipy.optimize import minimize_scalar
    rp = minimize_scalar(lambda p: chi2(rho - model_fs_only(w, p), sig),
                         bounds=(0.0, 0.999), method="bounded")
    rr = minimize_scalar(lambda R: chi2(rho - model_ms_only(w, R), sig),
                         bounds=(0.01, 0.99), method="bounded")
    fs_best = model_fs_only(w, rp.x)
    ms_best = model_ms_only(w, rr.x)
    print("\n=== FS only ===")
    print(f"  best p = {rp.x:.3f}, chi2 = {rp.fun:.1f}, "
          f"RMSE = {rmse(rho - fs_best):.4f}")
    print("=== MS only ===")
    print(f"  best R = {rr.x:.3f}, chi2 = {rr.fun:.1f}, "
          f"RMSE = {rmse(rho - ms_best):.4f}")

    print("\nresiduals of combined best fit (uOhm.cm):")
    for wi, ri, si in zip(w, res, sig):
        print(f"  w = {wi:6.1f} nm: {ri:+.3f}  ({ri/si:+.2f} sigma)")

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    wg = np.logspace(np.log10(32), 3, 120)
    # build curves on a common evaluator basis (use direct series; slower but exact)
    comb_fit = []
    comb_paper = []
    fs_c, ms_c = [], []
    for wi in wg:
        d = min(wi, H)
        comb_fit.append(RHO0 * WIRE.combined(wi, H, LAM, p_fit, R_fit, d))
        comb_paper.append(RHO0 * WIRE.combined(wi, H, LAM, 0.6, 0.50, d))
        fs_c.append(RHO0 * WIRE.rho_ratio(wi, H, LAM, p_fit))
        ms_c.append(RHO0 * float(M.ms_rho_ratio(LAM, d, R_fit)))

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(6.6, 6.8), sharex=True,
        gridspec_kw=dict(height_ratios=[3, 1], hspace=0.06))
    ax.errorbar(w, rho, yerr=sig, fmt="o", ms=5, mfc="white", mec="k",
                ecolor="k", elinewidth=1, capsize=2.5,
                label="Steinhoegl 2002, Fig. 3 (digitized)")
    ax.plot(wg, comb_fit, "-", c="tab:red", lw=1.8,
            label=f"combined fit: p={p_fit:.2f}, R={R_fit:.2f}")
    ax.plot(wg, comb_paper, "--", c="tab:gray", lw=1.3,
            label="combined, paper params (p=0.6, R=0.50)")
    ax.plot(wg, ms_c, "-.", c="tab:blue", lw=1.2,
            label=f"MS only (R={R_fit:.2f})")
    ax.plot(wg, fs_c, ":", c="tab:green", lw=1.4,
            label=f"FS only (p={p_fit:.2f})")
    ax.axhline(1.68, color="0.6", ls=":", lw=1)
    ax.text(650, 1.71, r"bulk $\rho_0$(295 K)", fontsize=8, color="0.4")
    ax.set_xscale("log")
    ax.set_ylabel(r"$\rho$  [$\mu\Omega\,$cm]")
    ax.set_ylim(1.5, 5.2)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("Cu wire resistivity vs. linewidth (h = 230 nm, 295 K)",
                 fontsize=10)

    axr.errorbar(w, res, yerr=sig, fmt="o", ms=4, mfc="white", mec="k",
                 ecolor="k", elinewidth=1, capsize=2)
    axr.axhline(0, color="tab:red", lw=1)
    axr.set_xlabel("linewidth w  [nm]")
    axr.set_ylabel("residual")
    fig.savefig(os.path.join(ROOT, "figures", "fit_steinhogl.png"),
                dpi=200, bbox_inches="tight")
    print("\nwrote figures/fit_steinhogl.png")


if __name__ == "__main__":
    main()
