"""Validation on an independent dataset: Yarimbiyik et al.,
Microelectron. Reliab. 49, 127 (2009) -- evaporated Cu films,
9-166 nm thick, with grain sizes measured per thickness (their
Tables 1 and 2).

Test protocol (per the project outline): parameters are constrained by
the literature, NOT refit to this dataset.
    p = 0        (their statement and Refs. therein: PVD Cu surfaces
                  scatter diffusely; also Kuan et al.)
    R = 0.32     (their g = 0.69 converted to MS reflection coefficient)
    lambda(301.75 K) = 38 nm   from rho*lambda = 6.6e-16 Ohm.m^2
                               (Gall, J. Appl. Phys. 119, 085101 (2016))
    rho0(301.75 K) = 1.736 uOhm.cm
                  = 1.678 (20 C, CRC) + 0.0067 uOhm.cm/C * 8.6 C
                    (temperature coefficient as used by the authors)
    d(t): linear interpolation/extrapolation of their measured
          in-plane grain size GS_xy versus thickness (Table 2)

Outputs figures/validate_yarimbiyik.png and a fit report.
Run from the repository root:  python src/validate_yarimbiyik.py
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import models as M

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

RHO0 = 1.736      # uOhm.cm at 28.6 C
LAM = 38.0        # nm at 28.6 C
P_LIT = 0.0
R_LIT = 0.32


def load(fname):
    rows = []
    with open(os.path.join(ROOT, "data", fname)) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line[0].isalpha():
                continue
            rows.append([float(x) for x in line.split(",")])
    return np.array(rows)


def main():
    films = load("yarimbiyik2009_table1_films.csv")     # t, rho
    grains = load("yarimbiyik2009_table2_grains.csv")   # t, GSxy
    t, rho = films[:, 0], films[:, 1]

    # grain size model: linear interpolation with linear extrapolation
    def d_of_t(tt):
        return np.interp(tt, grains[:, 0], grains[:, 1]) if \
            (grains[0, 0] <= tt <= grains[-1, 0]) else \
            _lin_extrap(tt, grains)

    def _lin_extrap(tt, g):
        if tt < g[0, 0]:
            s = (g[1, 1] - g[0, 1]) / (g[1, 0] - g[0, 0])
            return max(g[0, 1] + s * (tt - g[0, 0]), 1.0)
        s = (g[-1, 1] - g[-2, 1]) / (g[-1, 0] - g[-2, 0])
        return g[-1, 1] + s * (tt - g[-1, 0])

    d = np.array([d_of_t(x) for x in t])

    model = np.array([RHO0 * M.combined_film(ti, LAM, P_LIT, R_LIT, di)
                      for ti, di in zip(t, d)])
    fs_part = np.array([RHO0 * float(M.fs_film(ti, LAM, P_LIT)) for ti in t])
    ms_part = np.array([RHO0 * float(M.ms_rho_ratio(LAM, di, R_LIT))
                        for di in d])
    res = rho - model

    print("t [nm]   d [nm]   rho_exp   rho_model   resid")
    for a, b, c, e in zip(t, d, rho, model):
        print(f"{a:6.1f} {b:8.1f} {c:9.2f} {e:11.2f} {c - e:+8.2f}")
    print(f"\nRMSE (no refit) = {np.sqrt(np.mean(res**2)):.3f} uOhm.cm")
    print(f"mean |res|/rho  = {np.mean(np.abs(res)/rho)*100:.1f} %")

    # one-parameter sanity check (NOT used for conclusions): best-fit R
    from scipy.optimize import minimize_scalar
    def c2(R):
        m = np.array([RHO0 * M.combined_film(ti, LAM, P_LIT, R, di)
                      for ti, di in zip(t, d)])
        return float(np.sum((rho - m) ** 2))
    rb = minimize_scalar(c2, bounds=(0.05, 0.9), method="bounded")
    print(f"(diagnostic) best-fit R if freed: {rb.x:.3f}, "
          f"RMSE = {np.sqrt(rb.fun/len(t)):.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tg = np.logspace(np.log10(8), np.log10(200), 100)
    dg = np.array([d_of_t(x) for x in tg])
    mg = np.array([RHO0 * M.combined_film(ti, LAM, P_LIT, R_LIT, di)
                   for ti, di in zip(tg, dg)])
    m_lo = np.array([RHO0 * M.combined_film(ti, LAM, P_LIT, 0.25, di)
                     for ti, di in zip(tg, dg)])
    m_hi = np.array([RHO0 * M.combined_film(ti, LAM, P_LIT, 0.40, di)
                     for ti, di in zip(tg, dg)])
    fsg = np.array([RHO0 * float(M.fs_film(ti, LAM, P_LIT)) for ti in tg])
    msg = np.array([RHO0 * float(M.ms_rho_ratio(LAM, di, R_LIT))
                    for ti, di in zip(tg, dg)])

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.fill_between(tg, m_lo, m_hi, color="tab:red", alpha=0.15,
                    label="combined, R in [0.25, 0.40]")
    ax.plot(tg, mg, "-", c="tab:red", lw=1.8,
            label="combined (p=0, R=0.32, no refit)")
    ax.plot(tg, msg, "-.", c="tab:blue", lw=1.2, label="MS only (R=0.32)")
    ax.plot(tg, fsg, ":", c="tab:green", lw=1.4, label="FS only (p=0)")
    ax.plot(t, rho, "o", ms=6, mfc="white", mec="k",
            label="Yarimbiyik 2009 films")
    ax.axhline(RHO0, color="0.6", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("film thickness t  [nm]")
    ax.set_ylabel(r"$\rho$  [$\mu\Omega\,$cm]")
    ax.set_title("Independent validation: evaporated Cu films (302 K)",
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(os.path.join(ROOT, "figures", "validate_yarimbiyik.png"),
                dpi=200, bbox_inches="tight")
    print("wrote figures/validate_yarimbiyik.png")


if __name__ == "__main__":
    main()
