"""Benchmark: reproduce the model curves of Steinhoegl et al. (2002),
Fig. 3 with their published parameter set (rho0 = 1.9 uOhm.cm,
lambda = 40 nm, p = 0.6, R = 0.50, h = 230 nm, d = min(w, h)) and
compare against curve samples digitized from the printed figure.

Also quantifies the combined-model form question at small w:
additive resistivity increments (their Eq. (5)) versus a
multiplicative combination.

Run from the repository root:  python src/benchmark_fig3.py
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import models as M

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RHO0, LAM, H, P, R = 1.90, 40.0, 230.0, 0.6, 0.50
WIRE = M.RectWireInterp(cache_path=os.path.join(ROOT, "data", "cache",
                                                "chambers_p0_table.npz"))


def load_curves():
    out = []
    with open(os.path.join(ROOT, "data",
                           "steinhogl2002_fig3_curves.csv")) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("curve"):
                continue
            c, w, r, bad = line.split(",")
            if w == "any":
                continue
            out.append((c, float(w), float(r), int(bad)))
    return out


def main():
    print(f"{'curve':8s} {'w[nm]':>7} {'digitized':>9} {'model':>8} "
          f"{'diff':>7}  flag")
    rows = load_curves()
    for c, w, r, bad in rows:
        d = min(w, H)
        if c == "dashed":       # FS model
            m = RHO0 * WIRE.rho_ratio(w, H, LAM, P)
        elif c == "dashdot":    # MS model
            m = RHO0 * float(M.ms_rho_ratio(LAM, d, R))
        elif c == "solid":      # combined (additive)
            m = RHO0 * WIRE.combined(w, H, LAM, P, R, d)
        else:
            continue
        flag = "unreliable-sample" if bad else ""
        print(f"{c:8s} {w:7.1f} {r:9.3f} {m:8.3f} {r - m:+7.3f}  {flag}")

    print("\ncombined-model form at small w (paper parameter set):")
    print(f"{'w[nm]':>7} {'additive':>9} {'multiplicative':>15} {'diff':>7}")
    for w in (40.0, 50.0, 65.0, 80.0, 100.0, 150.0, 300.0, 600.0):
        d = min(w, H)
        ms = float(M.ms_rho_ratio(LAM, d, R))
        fs = WIRE.rho_ratio(w, H, LAM, P)
        add = RHO0 * (ms + fs - 1.0)
        mult = RHO0 * ms * fs
        print(f"{w:7.0f} {add:9.3f} {mult:15.3f} {mult - add:+7.3f}")

    # ------- figure -------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    wg = np.logspace(np.log10(32), 3, 140)
    fs_c, ms_c, add_c, mult_c = [], [], [], []
    for w in wg:
        d = min(w, H)
        ms = float(M.ms_rho_ratio(LAM, d, R))
        fs = WIRE.rho_ratio(w, H, LAM, P)
        fs_c.append(RHO0 * fs)
        ms_c.append(RHO0 * ms)
        add_c.append(RHO0 * (ms + fs - 1.0))
        mult_c.append(RHO0 * ms * fs)

    # digitized data points
    pts = []
    with open(os.path.join(ROOT, "data",
                           "steinhogl2002_fig3_points.csv")) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("w_nm"):
                continue
            v = [float(x) for x in line.strip().split(",")]
            pts.append(v)
    pts = np.array(pts)

    fig, ax = plt.subplots(figsize=(6.6, 4.9))
    ax.plot(wg, add_c, "-", c="tab:red", lw=1.8,
            label="combined, additive (Eq. 5)")
    ax.plot(wg, mult_c, "-", c="tab:orange", lw=1.2, alpha=0.9,
            label="combined, multiplicative")
    ax.plot(wg, ms_c, "-.", c="tab:blue", lw=1.2, label="MS model")
    ax.plot(wg, fs_c, "--", c="tab:green", lw=1.2, label="FS model")
    ax.axhline(1.68, color="0.6", ls=":", lw=1)

    for c, mkr, col in (("solid", "s", "tab:red"),
                        ("dashdot", "D", "tab:blue"),
                        ("dashed", "^", "tab:green")):
        pw = [w for cc, w, r, bad in load_curves() if cc == c and not bad]
        pr = [r for cc, w, r, bad in load_curves() if cc == c and not bad]
        ax.plot(pw, pr, mkr, ms=6, mfc="none", mec=col, mew=1.4,
                ls="none",
                label=f"digitized {c} curve" if c == "solid" else None)
    ax.errorbar(pts[:, 0], pts[:, 1],
                yerr=0.5 * (pts[:, 2] + pts[:, 3]), fmt="o", ms=5,
                mfc="white", mec="k", ecolor="k", elinewidth=1, capsize=2,
                label="data (digitized)")
    ax.set_xscale("log")
    ax.set_xlabel("linewidth w  [nm]")
    ax.set_ylabel(r"$\rho$  [$\mu\Omega\,$cm]")
    ax.set_ylim(1.5, 5.3)
    ax.set_title("Reproduction of Steinhoegl (2002) Fig. 3, "
                 "paper parameter set", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(os.path.join(ROOT, "figures", "benchmark_fig3.png"),
                dpi=200, bbox_inches="tight")
    print("\nwrote figures/benchmark_fig3.png")


if __name__ == "__main__":
    main()
