
import os
import sys
import csv
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import models as M

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

RHO0 = 1.736
LAM = 38.0
P_LIT = 0.0

R_LIT = 0.32
R_TRANSFER = 0.425


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
    films = load("yarimbiyik2009_table1_films.csv")
    grains = load("yarimbiyik2009_table2_grains.csv")

    t, rho = films[:, 0], films[:, 1]

    def _lin_extrap(tt, g):
        if tt < g[0, 0]:
            slope = (
                (g[1, 1] - g[0, 1]) /
                (g[1, 0] - g[0, 0])
            )

            return max(
                g[0, 1] + slope * (tt - g[0, 0]),
                1.0,
            )

        slope = (
            (g[-1, 1] - g[-2, 1]) /
            (g[-1, 0] - g[-2, 0])
        )

        return g[-1, 1] + slope * (tt - g[-1, 0])

    def d_of_t(tt):
        if grains[0, 0] <= tt <= grains[-1, 0]:
            return np.interp(
                tt,
                grains[:, 0],
                grains[:, 1],
            )

        return _lin_extrap(tt, grains)

    d = np.array([d_of_t(x) for x in t])

    # Experiment 4: Dataset 2 literature parameter
    model_lit = np.array([
        RHO0 * M.combined_film(
            ti,
            LAM,
            P_LIT,
            R_LIT,
            di,
        )
        for ti, di in zip(t, d)
    ])

    # Experiment 5: transfer R=0.425 from Steinhogl dataset
    model_transfer = np.array([
        RHO0 * M.combined_film(
            ti,
            LAM,
            P_LIT,
            R_TRANSFER,
            di,
        )
        for ti, di in zip(t, d)
    ])

    residual_lit = rho - model_lit
    residual_transfer = rho - model_transfer

    rmse_lit = np.sqrt(np.mean(residual_lit ** 2))
    rmse_transfer = np.sqrt(np.mean(residual_transfer ** 2))

    error_lit = (
        np.mean(np.abs(residual_lit) / rho) * 100
    )

    error_transfer = (
        np.mean(np.abs(residual_transfer) / rho) * 100
    )

    print("Yarimbiyik literature parameter")
    print(f"R = {R_LIT:.3f}")
    print(f"RMSE = {rmse_lit:.3f} uOhm.cm")
    print(f"mean relative error = {error_lit:.1f} %")

    print("\nSteinhogl parameter transfer")
    print(f"R = {R_TRANSFER:.3f}")
    print(f"RMSE = {rmse_transfer:.3f} uOhm.cm")
    print(f"mean relative error = {error_transfer:.1f} %")

    print("\nt [nm]   d [nm]   rho_exp   R=0.32   R=0.425   resid_032   resid_0425")

    for values in zip(
        t,
        d,
        rho,
        model_lit,
        model_transfer,
        residual_lit,
        residual_transfer,
    ):
        ti, di, exp, m_lit, m_transfer, r_lit, r_transfer = values

        print(
            f"{ti:6.1f} "
            f"{di:8.1f} "
            f"{exp:9.2f} "
            f"{m_lit:8.2f} "
            f"{m_transfer:9.2f} "
            f"{r_lit:+10.2f} "
            f"{r_transfer:+12.2f}"
        )

    figures_dir = os.path.join(ROOT, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    csv_path = os.path.join(
        figures_dir,
        "validate_yarimbiyik_transfer.csv",
    )

    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)

        writer.writerow([
            "t_nm",
            "d_nm",
            "rho_exp_uohm_cm",
            "rho_model_R032_uohm_cm",
            "rho_model_R0425_uohm_cm",
            "residual_R032_exp_minus_model",
            "residual_R0425_exp_minus_model",
        ])

        for row in zip(
            t,
            d,
            rho,
            model_lit,
            model_transfer,
            residual_lit,
            residual_transfer,
        ):
            writer.writerow(row)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tg = np.logspace(
        np.log10(8),
        np.log10(200),
        100,
    )

    dg = np.array([d_of_t(x) for x in tg])

    curve_lit = np.array([
        RHO0 * M.combined_film(
            ti,
            LAM,
            P_LIT,
            R_LIT,
            di,
        )
        for ti, di in zip(tg, dg)
    ])

    curve_transfer = np.array([
        RHO0 * M.combined_film(
            ti,
            LAM,
            P_LIT,
            R_TRANSFER,
            di,
        )
        for ti, di in zip(tg, dg)
    ])

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(7.2, 7.5),
        sharex=True,
    )

    axes[0].plot(
        tg,
        curve_lit,
        "-",
        linewidth=1.8,
        label="Yarimbiyik literature R=0.32",
    )

    axes[0].plot(
        tg,
        curve_transfer,
        "--",
        linewidth=1.8,
        label="Transferred Steinhogl R=0.425",
    )

    axes[0].plot(
        t,
        rho,
        "o",
        markersize=6,
        markerfacecolor="white",
        markeredgecolor="black",
        label="Yarimbiyik experiment",
    )

    axes[0].set_xscale("log")
    axes[0].set_ylabel(r"Resistivity $\rho$ [$\mu\Omega\cdot$cm]")
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].grid(alpha=0.25)

    axes[1].axhline(
        0,
        linestyle=":",
        linewidth=1,
    )

    axes[1].plot(
        t,
        residual_lit,
        "o-",
        label="R=0.32",
    )

    axes[1].plot(
        t,
        residual_transfer,
        "s--",
        label="R=0.425",
    )

    axes[1].set_xscale("log")
    axes[1].set_xlabel("Film thickness t [nm]")
    axes[1].set_ylabel(
        r"Residual: experiment $-$ model [$\mu\Omega\cdot$cm]"
    )
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].grid(alpha=0.25)

    fig.tight_layout()

    figure_path = os.path.join(
        figures_dir,
        "validate_yarimbiyik_transfer.png",
    )

    fig.savefig(
        figure_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    output_path = os.path.join(
        figures_dir,
        "validate_yarimbiyik_transfer_output.txt",
    )

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("Yarimbiyik literature parameter\n")
        fh.write(f"R = {R_LIT:.3f}\n")
        fh.write(f"RMSE = {rmse_lit:.3f} uOhm.cm\n")
        fh.write(
            f"mean relative error = {error_lit:.1f} %\n\n"
        )

        fh.write("Steinhogl parameter transfer\n")
        fh.write(f"R = {R_TRANSFER:.3f}\n")
        fh.write(
            f"RMSE = {rmse_transfer:.3f} uOhm.cm\n"
        )
        fh.write(
            f"mean relative error = {error_transfer:.1f} %\n"
        )

    print("\nwrote", os.path.relpath(figure_path, ROOT))
    print("wrote", os.path.relpath(csv_path, ROOT))
    print("wrote", os.path.relpath(output_path, ROOT))


if __name__ == "__main__":
    main()
