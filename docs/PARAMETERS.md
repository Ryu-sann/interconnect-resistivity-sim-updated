# Parameter table with sources

All parameters used anywhere in this repository, with provenance. "Fixed"
means taken from the cited source and never fitted; "fitted" means adjusted
here within the stated literature range.

## Copper, Dataset 1 (Steinhögl 2002 wires, 295 K)

| symbol | value | status | source / justification |
|---|---|---|---|
| ρ₀ | 1.90 µΩ·cm | fixed | Steinhögl et al. (2002) fit; includes Ta-liner dilution (their Ta-corrected value is 1.80; pure bulk Cu 1.75 at RT per their text, 1.678 at 20 °C per CRC Handbook) |
| λ | 40 nm | fixed | Steinhögl et al. (2002); consistent with ρλ = 6.7×10⁻¹⁶ Ωm² (Gall 2016) within 15 % |
| h | 230 nm | fixed | Steinhögl et al. (2002), wire height |
| d | min(w, h) | fixed (assumption) | Steinhögl et al. (2002): grains span the width, height-limited vertically |
| p | fitted in [0, 1] | fitted | result: unconstrained ([0,1] at 1σ); literature Cu: mostly diffuse to ~0.6 |
| R | fitted in [0.05, 0.9] | fitted | result: 0.43 [0.39, 0.55]; literature: Mayadas–Shatzkes 0.24, Kuan ~0.3, Steinhögl 0.50, up to ~0.65 (electroplated, Wu et al.) |

Temperature dependence (their Table I, for reference): (ρ₀, λ) =
(0.2, 330) at 77 K, (1.9, 40) at 300 K, (2.83, 26.5) at 423 K,
(3.97, 19) at 573 K. Caveat: the 77 K product ρ₀λ = 66 µΩ·cm·nm is
inconsistent with the 300–573 K product (75–76 µΩ·cm·nm); use with care.

## Copper, Dataset 2 (Yarimbiyik 2009 films, 301.75 K)

| symbol | value | status | source / justification |
|---|---|---|---|
| ρ₀ | 1.736 µΩ·cm | fixed | 1.678 µΩ·cm (20 °C, CRC) + 0.0067 µΩ·cm/°C × 8.6 °C (temperature coefficient as used by the authors, their Ref. [4]) |
| λ | 38 nm | fixed | ρλ = 6.6×10⁻¹⁶ Ωm² (Gall 2016) divided by ρ₀ above |
| t | Table 1, col. 2 | measured | Dektak, oxidation-corrected (−6.25 nm); σ_t = 3–4 nm |
| d(t) | Table 2 interp. | measured | XRD GS_z × 0.58 (EBSD-calibrated); linear inter-/extrapolation |
| p | 0 | fixed | authors' best fit; Kuan et al. (PVD Cu); Sondheimer refs. |
| R | 0.32 | fixed | authors' g = 0.69 via g = 0.2617·ln(R)+0.9913; Kuan R ≈ 0.3 |

## Ruthenium (exploratory case study, RT)

| symbol | value | status | source |
|---|---|---|---|
| ρ₀(Ru) | 7.1 µΩ·cm | fixed | commonly cited RT bulk value (see Gall 2016 and interconnect literature) |
| λ(Ru) | 6.6 nm | fixed | Gall (2016) |
| ρ₀(Cu) | 1.712 µΩ·cm | fixed | Gall (2016), 293 K |
| λ(Cu) | 39.9 nm | fixed | Gall (2016) |
| liner | 1.5 nm/wall (Cu), 0 (Ru) | assumption | representative TaN/Ta scaling; non-conducting |
| p, R, d | 0, 0.43, w | assumption | same for both metals; see src/ru_case_study.py header |

## Digitization (Dataset 1)

| quantity | value |
|---|---|
| x calibration | w[nm] = 100·10^((px−207)/169.5), verified on 16 minor ticks to <1 px |
| y calibration | ρ = 5 − (py−51.5)/51.667 µΩ·cm, verified on 4 major ticks to <1 px |
| pixel resolution | 0.019 µΩ·cm/px; 1.4 %/px in w |
| independent anchor | MS plateau: digitized 2.377 vs analytic 2.3758 µΩ·cm |
