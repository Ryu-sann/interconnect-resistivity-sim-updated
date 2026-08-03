"""Exploratory case study: barrierless Ru versus Cu-with-liner for
narrow interconnects (outline item 5; qualitative, assumption-driven).

ASSUMPTIONS (all explicit; results are indicative only):
  * Square wire, w x w; effective resistivity referenced to the full
    w^2 cross-section.
  * Cu requires a diffusion barrier/liner: b = 1.5 nm per sidewall
    (2 x 1.5 nm total per dimension), treated as non-conducting; the Cu
    conductor is (w - 2b)^2.  Ru is assumed barrierless (b = 0).
  * Bamboo grain structure: d = conductor width for both metals.
  * Fully diffuse surfaces (p = 0) for both metals.
  * Same GB reflection coefficient R = 0.43 for both metals (our Cu fit;
    Ru literature values are similar in magnitude but uncertain).
  * Isotropic transport (real Ru is hcp and anisotropic; textured films
    differ).
  * Room-temperature parameters from Gall, J. Appl. Phys. 119, 085101
    (2016): Cu lambda = 39.9 nm, rho0 = 1.712 uOhm.cm;
            Ru lambda =  6.6 nm, rho0 = 7.1  uOhm.cm.

Run from the repository root:  python src/ru_case_study.py
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import models as M

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
WIRE = M.RectWireInterp(cache_path=os.path.join(ROOT, "data", "cache",
                                                "chambers_p0_table.npz"))
R_GB, P = 0.43, 0.0
B_LINER = 1.5     # nm per sidewall for Cu

CU = dict(lam=39.9, rho0=1.712)
RU = dict(lam=6.6, rho0=7.1)


def rho_eff(w, metal, liner):
    wc = w - 2.0 * liner
    if wc <= 2.0:
        return np.inf
    r = metal["rho0"] * WIRE.combined(wc, wc, metal["lam"], P, R_GB, d=wc)
    return r * (w / wc) ** 2


def main():
    print(f"{'w[nm]':>6} {'Cu+liner':>9} {'Ru bare':>8}")
    ws = np.array([8, 10, 12, 15, 18, 22, 28, 36, 48, 64, 90, 130, 200])
    cu = np.array([rho_eff(w, CU, B_LINER) for w in ws])
    ru = np.array([rho_eff(w, RU, 0.0) for w in ws])
    for w, a, b in zip(ws, cu, ru):
        print(f"{w:6.0f} {a:9.2f} {b:8.2f}")

    wg = np.logspace(np.log10(7), np.log10(220), 160)
    cug = np.array([rho_eff(w, CU, B_LINER) for w in wg])
    rug = np.array([rho_eff(w, RU, 0.0) for w in wg])
    diff = cug - rug
    idx = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]
    if len(idx):
        i = idx[-1]
        wx = np.interp(0.0, [diff[i], diff[i + 1]], [wg[i], wg[i + 1]])
        print(f"\ncrossover (Ru better below): w = {wx:.1f} nm "
              f"under the stated assumptions")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(wg, cug, "-", c="tab:orange", lw=1.8,
            label="Cu, 1.5 nm liner per wall")
    ax.plot(wg, rug, "-", c="tab:purple", lw=1.8, label="Ru, barrierless")
    if len(idx):
        ax.axvline(wx, color="0.7", lw=0.9, ls="--")
        ax.text(wx * 1.05, 25, f"crossover\n~{wx:.0f} nm", fontsize=8,
                color="0.35")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("drawn linewidth w  [nm]")
    ax.set_ylabel(r"effective $\rho$  [$\mu\Omega\,$cm]")
    ax.set_title("Exploratory: Cu+liner vs. barrierless Ru "
                 "(square wire, p=0, R=0.43, d=w)", fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(os.path.join(ROOT, "figures", "ru_case_study.png"),
                dpi=200, bbox_inches="tight")
    print("wrote figures/ru_case_study.png")


if __name__ == "__main__":
    main()
