"""One-at-a-time sensitivity analysis of the combined FS+MS model for
the Steinhoegl wire geometry (h = 230 nm, rho0 = 1.9 uOhm.cm,
lambda = 40 nm, d = s * min(w, h)).

Baseline: paper parameter set p = 0.6, R = 0.50, s = 1.
Literature ranges swept one at a time:
    p in [0, 0.6]     (Cu: mostly diffuse to moderately specular)
    R in [0.20, 0.65] (Mayadas 0.24; Kuan ~0.3; Rossnagel/Wu up to ~0.65)
    s in [0.5, 2]     (grain size between half and twice the linewidth)

Metric: change in rho at w = 50 nm (and full curves for the figure).
Run from the repository root:  python src/sensitivity.py
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import models as M

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RHO0, LAM, H = 1.90, 40.0, 230.0
WIRE = M.RectWireInterp(cache_path=os.path.join(ROOT, "data", "cache",
                                                "chambers_p0_table.npz"))
P0, R0, S0 = 0.6, 0.50, 1.0
W_REF = 50.0


def rho(w, p=P0, R=R0, s=S0):
    d = s * min(w, H)
    return RHO0 * WIRE.combined(w, H, LAM, p, R, d)


def main():
    base = rho(W_REF)
    print(f"baseline rho(w = {W_REF:.0f} nm) = {base:.3f} uOhm.cm "
          f"(p={P0}, R={R0}, s={S0})\n")

    sweeps = {
        "p (specularity)": [("p=0.0", dict(p=0.0)), ("p=0.3", dict(p=0.3)),
                            ("p=0.6", dict(p=0.6))],
        "R (GB reflection)": [("R=0.20", dict(R=0.20)), ("R=0.35", dict(R=0.35)),
                              ("R=0.50", dict(R=0.50)), ("R=0.65", dict(R=0.65))],
        "d scale s": [("s=0.5", dict(s=0.5)), ("s=1.0", dict(s=1.0)),
                      ("s=2.0", dict(s=2.0))],
    }
    spans = {}
    for name, entries in sweeps.items():
        print(f"--- {name} ---")
        vals = []
        for label, kw in entries:
            v = rho(W_REF, **kw)
            vals.append(v)
            print(f"  {label:8s} rho(50nm) = {v:6.3f}  (delta = {v-base:+.3f})")
        spans[name] = max(vals) - min(vals)
        print()

    print("total span of rho(50 nm) over each literature range:")
    for name, s in sorted(spans.items(), key=lambda kv: -kv[1]):
        print(f"  {name:22s} {s:6.3f} uOhm.cm")
    dom = max(spans, key=spans.get)
    print(f"\ndominant parameter at w = {W_REF:.0f} nm: {dom}")

    # ------- figure: curve families -------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    wg = np.logspace(np.log10(32), 3, 90)

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9), sharey=True)
    panels = [
        (axes[0], "p", [(0.0, "0"), (0.3, "0.3"), (0.6, "0.6")],
         lambda v: dict(p=v)),
        (axes[1], "R", [(0.20, "0.20"), (0.35, "0.35"), (0.50, "0.50"),
                        (0.65, "0.65")], lambda v: dict(R=v)),
        (axes[2], "s = d/min(w,h)", [(0.5, "0.5"), (1.0, "1"), (2.0, "2")],
         lambda v: dict(s=v)),
    ]
    cmap = plt.cm.viridis
    for ax, pname, vals, mk in panels:
        for i, (v, lab) in enumerate(vals):
            c = cmap(i / max(len(vals) - 1, 1) * 0.85)
            y = [rho(w, **mk(v)) for w in wg]
            ax.plot(wg, y, lw=1.6, color=c, label=f"{pname.split()[0]}={lab}")
        ax.axvline(W_REF, color="0.8", lw=0.8, zorder=0)
        ax.set_xscale("log")
        ax.set_xlabel("linewidth w [nm]")
        ax.set_title(f"vary {pname}", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
    axes[0].set_ylabel(r"$\rho$ [$\mu\Omega\,$cm]")
    fig.suptitle("Sensitivity, one-at-a-time around the Steinhoegl parameter set",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "sensitivity.png"),
                dpi=200, bbox_inches="tight")
    print("wrote figures/sensitivity.png")


if __name__ == "__main__":
    main()
