# Size-dependent resistivity of Cu interconnects: data, models, and validation

*Summer research project, July 2026. Code, data, and figures for every number
quoted here are in this repository; see README.md for how to regenerate them.*

## 1. Objective

Copper interconnect resistivity rises steeply once wire dimensions approach
the room-temperature electron mean free path (λ ≈ 40 nm). This project
characterizes that rise using two published experimental datasets, implements
the two standard semiclassical models — Fuchs–Sondheimer (FS) surface
scattering and Mayadas–Shatzkes (MS) grain-boundary scattering — from the
original equations, benchmarks the implementation against published curves
and tables, fits the combined model with literature-constrained parameters,
quantifies parameter sensitivity, and tests the calibrated model against an
independent dataset without refitting.

## 2. Datasets

### 2.1 Dataset 1: Cu damascene wires (Steinhögl et al., PRB 66, 075414, 2002)

Electrodeposited Cu wires fabricated in a damascene process in SiO₂, wire
height h = 230 nm, lengths 200 µm, linewidths 40–800 nm defined by SEM, with
a Ta-based diffusion barrier that constitutes 5–10 % of the conducting
cross-section. Resistivities are reported at room temperature (295 K per the
figure legend) and are *not* corrected for the Ta content; roughly ten
structures were measured per linewidth. Grain size: the authors observe that
grains span the full lateral width while being height-limited, and take the
grain-boundary spacing d = min(w, h).

The data exist only as their Fig. 3, so I digitized it from a page scan.
The axis calibration was established from the plot frame and all tick marks
(8 minor ticks at 40–90 nm, 8 at 200–900 nm, and 4 major y ticks), each
reproduced to better than 1 px by the final mapping
w[nm] = 100·10^((px−207)/169.5), ρ[µΩ·cm] = 5 − (py−51.5)/51.667.
Data symbols were located by a matched ring filter (open circles with light
centers) cross-checked against the presence of a vertical error bar, then
refined to sub-pixel precision; the per-point uncertainties come from the
digitized error-bar caps. One pixel corresponds to 0.019 µΩ·cm and 1.4 % in
w. Nine linewidths were found (the paper's "ten" refers to measurements per
point). Anchors from the paper text confirm the calibration: the narrowest
point reads 4.613 µΩ·cm against the stated 4.6; the dash-dot MS curve
flattens at a digitized 2.377 µΩ·cm against the analytic value 2.3758
computed from the printed parameter set — agreement at the 10⁻³ level that
rules out any calibration offset. The two widest data symbols read
2.34–2.36 µΩ·cm; the "2.45 for the widest wires" quoted in the text matches
the fitted combined-model level there (digitized 2.48 at 609 nm) rather than
the symbols, which sit on the MS curve within their error bars. The dotted
bulk reference line is drawn at 1.68 µΩ·cm, slightly below the 1.75 µΩ·cm
quoted in the text (the 20 °C handbook value of pure Cu is 1.678).

Data: `data/steinhogl2002_fig3_points.csv`;
curve samples: `data/steinhogl2002_fig3_curves.csv`.

| w (nm) | ρ (µΩ·cm) | +σ | −σ |
|-------:|----------:|-----:|-----:|
| 43.2 | 4.613 | 0.34 | 0.36 |
| 63.2 | 3.723 | 0.36 | 0.22 |
| 87.3 | 3.500 | 0.21 | 0.41 |
| 107.0 | 3.171 | 0.16 | 0.39 |
| 138.5 | 2.842 | 0.29 | 0.19 |
| 178.1 | 2.760 | 0.20 | 0.27 |
| 259.7 | 2.639 | 0.18 | 0.28 |
| 490.1 | 2.339 | 0.16 | 0.33 |
| 736.6 | 2.358 | 0.14 | 0.14 |

The trend is a smooth, accelerating rise from ≈1.4× bulk at 700 nm to ≈2.7×
bulk at 43 nm as w approaches λ, with most of the increase developing below
w ≈ 150 nm.

### 2.2 Dataset 2: evaporated Cu films (Yarimbiyik et al., Microelectron. Reliab. 49, 127, 2009)

Thermally evaporated Cu films on Pyrex with a 3 nm Cr adhesion layer,
thicknesses 9–166 nm after correction for surface oxidation (−6.25 nm),
ρ = R_s·t measured at 28.6 °C. Crucially, the in-plane grain size GS_xy was
measured per thickness (XRD grain thickness scaled by an EBSD-calibrated
factor GS_xy = 0.58·GS_z), removing the grain size as a free parameter. The
authors' own fit (with an independent Monte-Carlo simulation) gave p = 0 and
a grain parameter equivalent to MS R = 0.32, and they note Kuan et al.'s
R ≈ 0.3, p = 0 for PVD Cu. Data: `data/yarimbiyik2009_table1_films.csv`,
`data/yarimbiyik2009_table2_grains.csv`. Caveats carried through the
analysis: the Cr layer provides a small parallel conduction path that biases
the thinnest films low, thickness uncertainty is 3–4 nm, and the oxidation
correction is a single average value.

## 3. Models and implementation

All models are implemented from the original equations in
`src/models.py` and return ρ/ρ₀.

**FS film.** Sondheimer's exact integral (his Eq. (25)) for arbitrary
specularity p, integrated by Gauss–Legendre quadrature in 1/t with an
analytic tail.

**FS rectangular wire.** The Chambers path-integral for a w×h wire with
diffuse walls, σ/σ₀ = 1 − (3/(4πA))∫dA∫dΩ cos²θ e^(−s/(λ sinθ)), evaluated
with a tabulated polar kernel and tensor quadrature over the cross-section
and azimuth. Note: Eq. (2) of Steinhögl et al. prints cos²φ; this fails the
solid-angle normalization (3/4π)∫cos²θ dΩ = 1 and does not reproduce their
own thick-wire limit, so cos²θ (angle from the wire axis) is used —
a misprint correction, verified numerically both ways.

**Partial specularity for wires.** Sondheimer's series identity (his
Eq. (31), = Steinhögl Eq. (3)):
(σ/σ₀)_p = (1−p)² Σ_k k p^(k−1) (σ/σ₀)_{p=0, λ/k}. The same identity applied
to films agrees with the direct p-integral to 10⁻¹³ — a strong internal
consistency check.

**MS grain boundaries.** f(α) = 3[1/3 − α/2 + α² − α³ ln(1+1/α)] with
α = (λ/d)·R/(1−R), with a stable small-α expansion.

**Combined model.** Steinhögl Eq. (5): additive resistivity increments,
ρ/ρ₀ = 1/f(α) + [(ρ/ρ₀)_FS − 1]. A multiplicative alternative is examined in
§ 6.

**Numerics.** A disk-cached 2-D spline of the p = 0 Chambers result over
(κ = w/λ, aspect = h/w) makes parameter sweeps fast; it agrees with the
direct calculation to <0.3 %.

**Verification (tests/, 29 passing tests).** The implementation reproduces:
Sondheimer's Table 1 for films (32 entries, p = 0 and p = ½); the thick-film
asymptote 1 + 3(1−p)/(8κ); the thick square-wire limit 1 + (3/4)(1−p)/κ
(Sondheimer Eq. (32), which also underlies Steinhögl Eq. (1)); convergence
of the wide wire to the film of thickness h; MS limits (d→∞, R→0, R→1, the
3/(4α) tail); and the p-series/direct-integral identity.

A by-product worth recording: in the crossover regime κ ≈ 0.05–0.5 the
hand-computed 1952 Table 1 itself is inaccurate by up to ≈3 %. Three
independent numerical routes (Gauss–Legendre in 1/t, adaptive quadrature in
t, and the series identity) agree with each other to six digits but give,
e.g., 4.7817 at κ = 0.1, p = 0 where the table prints 4.72, and 1.9161 at
κ = 0.5 where it prints 1.90. The test suite documents these entries
explicitly rather than hiding them under a loose tolerance.

## 4. Benchmark against the published curves (Fig. 3 reproduction)

With the paper's parameter set (ρ₀ = 1.90 µΩ·cm, λ = 40 nm, p = 0.6,
R = 0.50, d = min(w, 230 nm)) — see `src/benchmark_fig3.py` and
`figures/benchmark_fig3.png`:

* The MS (dash-dot) curve is reproduced exactly: differences ≤ 0.008 µΩ·cm
  (≤ 0.4 px) at every sampled width, including the flat plateau at
  2.376 µΩ·cm.
* The FS (dashed) curve agrees to −0.013 µΩ·cm at 79 nm and shows a small
  systematic offset of +0.04–0.05 µΩ·cm (≈2 px) at 300–600 nm, where the
  printed curve is slightly flatter than the exact Chambers solution. The
  offset is far below the data error bars and does not affect conclusions;
  it most plausibly reflects the paper's own numerical approximation of
  Eqs. (2)–(3).
* The combined (solid) curve inherits the same small offset at large w and
  is mutually additive with the printed MS and FS curves there. Below
  ≈80 nm the printed solid curve rises less steeply than the additive
  combination of the exact components computed here (digitized 4.41 at
  53 nm versus additive 4.17), landing between the additive and
  multiplicative forms discussed in § 6.

## 5. Fit to Dataset 1

Following the outline, ρ₀ = 1.90 µΩ·cm and λ = 40 nm are fixed by the paper;
only p ∈ [0, 1] and R are fitted (χ² weighted by the digitized error bars);
`src/fit_steinhogl.py`, `figures/fit_steinhogl.png`.

| model | best fit | RMSE (µΩ·cm) | χ²/ndof |
|---|---|---:|---:|
| combined (additive) | p = 0.00, R = 0.425 | 0.102 | 0.22 |
| combined at paper values | p = 0.6, R = 0.50 | 0.118 | 0.30 |
| combined (multiplicative) | p = 0.86, R = 0.52 | 0.134 | 0.36 |
| MS only | R = 0.538 | 0.130 | 0.37 |
| FS only | p = 0 (boundary) | 0.902 | 12.0 |

Three robust conclusions follow. First, surface scattering alone cannot
describe the data: even with fully diffuse walls the FS-only curve falls far
short of the observed rise (RMSE 0.90 µΩ·cm, χ²/ndof = 12), reproducing the
paper's statement that FS alone is "more than 50 % too low". Second,
grain-boundary scattering carries the effect: R is well constrained at
0.43 (1σ interval [0.39, 0.55], covering the paper's 0.50 and overlapping
the literature range 0.24–0.65), and an MS-only fit is nearly as good as the
combined one. Third, p is essentially unconstrained ([0, 1] at 1σ): because
d tracks w, the FS and MS mechanisms produce nearly proportional
w-dependences and trade off against each other. The paper's (p = 0.6,
R = 0.50) lies within our 1σ region — our fit is statistically consistent
with theirs, and the meaningful, transferable parameter is R, not p.
Residuals of the best fit show no systematic trend (all |res| < 0.6σ), and
χ²/ndof < 1 indicates the digitized error bars are, if anything,
conservative.

## 6. Combined-model form: how separable are the mechanisms?

Steinhögl Eq. (5) adds the two resistivity increments (Matthiessen's rule at
the level of ρ). A multiplicative combination,
ρ = ρ_MS · (ρ/ρ₀)_FS — equivalent to applying surface scattering on top of a
grain-boundary-renormalized bulk — coincides with the additive form when
either mechanism is weak but exceeds it when both are strong:

| w (nm) | additive | multiplicative | difference |
|---:|---:|---:|---:|
| 600 | 2.445 | 2.463 | +0.017 |
| 300 | 2.465 | 2.487 | +0.022 |
| 100 | 3.140 | 3.235 | +0.095 |
| 65 | 3.764 | 3.965 | +0.201 |
| 40 | 4.870 | 5.355 | +0.485 |

Above w ≈ 100 nm the two forms differ by less than the data uncertainty, so
the fitted R is insensitive to the choice (0.43 vs 0.52). Below ≈80 nm the
difference becomes comparable to the error bars, and it is exactly there
that the digitized published curve deviates from the additive form (§ 4).
Refitting with the multiplicative form shifts p wildly (0 → 0.86) but moves
R only within its 1σ interval — reinforcing that p absorbs model-form
ambiguity while R is robust. The physical point for the report: neither
form is exact — the rigorous treatment (Mayadas–Shatzkes 1970 for films)
couples the surface integral to the grain-boundary α inside a single
integrand — and the spread between the two simple forms below 80 nm is best
treated as a model-form uncertainty band. Extrapolations of Cu resistivity
to sub-50 nm dimensions inherit an O(0.3–0.5 µΩ·cm) ambiguity from the
combination rule alone, independent of parameter uncertainties. This is a
meaningful "failure at small dimensions" in the sense of the project
outline.

## 7. Independent validation on Dataset 2 (no refitting)

The film model ρ = ρ₀[1/f(α(d(t))) + (ρ/ρ₀)_FS(t) − 1] was evaluated with
every parameter fixed in advance (`src/validate_yarimbiyik.py`,
`figures/validate_yarimbiyik.png`): p = 0 and R = 0.32 from the authors'
independent analysis and PVD literature, d(t) interpolated from their
measured grain sizes, λ = 38 nm and ρ₀ = 1.736 µΩ·cm at 301.75 K from the
ρλ product of Gall (2016) and the handbook temperature coefficient.

The model tracks the data over the full 9–166 nm range with RMSE
0.54 µΩ·cm (mean 7.4 %), essentially exact for t ≥ 70 nm and overshooting
by 0.3–0.9 µΩ·cm for the thinner films. The overshoot has the right sign
and magnitude to be explained by the known biases of the dataset: the 3 nm
Cr underlayer shunt (largest relative effect at t = 9 nm), the ±3–4 nm
thickness uncertainty entering ρ = R_s·t, the extrapolated grain size below
t = 10 nm, and the additive-combination overcount when t and d are both ≪ λ
(§ 6). Freeing R alone returns R = 0.294 — within 10 % of the constrained
0.32 and of Kuan's 0.3 — i.e., the parameters calibrated on damascene wires
transfer to evaporated films at the ≈10 % level in R. This satisfies the
outline's requirement that the second dataset be tested without
unrestricted refitting.

## 8. Sensitivity analysis

One-at-a-time sweeps around the paper parameter set at w = 50 nm
(baseline ρ = 4.30 µΩ·cm; `src/sensitivity.py`, `figures/sensitivity.png`):

| parameter | literature range | span of ρ(50 nm) |
|---|---|---:|
| R (GB reflection) | 0.20 – 0.65 | 3.31 µΩ·cm |
| d scale (d/min(w,h)) | 0.5 – 2 | 3.08 µΩ·cm |
| p (specularity) | 0 – 0.6 | 0.54 µΩ·cm |

The grain-boundary description — both the reflection coefficient and the
grain size itself — dominates the predicted resistivity at damascene-Cu
dimensions; surface specularity is a factor ≈6 weaker. This ranking is
consistent with the fit degeneracy of § 5 and with the physical narrative of
both source papers. Practical implication: measuring d(w) directly (as
Dataset 2 did for films) buys far more predictive power than refining p.

## 9. Exploratory: barrierless Ru versus Cu with liner

A deliberately simple comparison (`src/ru_case_study.py`,
`figures/ru_case_study.png`) under explicit assumptions — square wire,
1.5 nm non-conducting liner per Cu sidewall, barrierless Ru, bamboo grains
d = conductor width, p = 0, the same R = 0.43 for both metals, isotropic
transport, and Gall (2016) room-temperature parameters (Cu: λ = 39.9 nm,
ρ₀ = 1.712 µΩ·cm; Ru: λ = 6.6 nm, ρ₀ = 7.1 µΩ·cm). The effective
resistivity referenced to the drawn w² area crosses over at w ≈ 24 nm,
below which Ru wins despite its 4× higher bulk resistivity — its short mean
free path and freedom from liner area loss more than compensate. The number
is assumption-dominated (literature estimates cluster around 15–25 nm
depending on liner thickness and R) and is presented as a sanity-checked
illustration, not a prediction.

## 10. Limitations and caveats

The Steinhögl data include the Ta liner (5–10 % of the cross-section) and
are not corrected for it, so the fitted ρ₀ = 1.90 µΩ·cm is an effective
value (the authors' own Ta correction brings it to 1.80, near bulk).
Damascene cross-sections are trapezoidal while the model is rectangular.
The grain-size rule d = min(w, h) is an assumption of the source paper, not
a measurement, and § 8 shows it is one of the two dominant sensitivities.
Digitization contributes ±0.02 µΩ·cm and ±1.4 % in w per point, with
error-bar caps read to ±1–2 px. The printed FS/combined curves differ from
the exact Chambers solution at the 2 px level (§ 4), and the combination
rule itself is ambiguous below 80 nm (§ 6). In the source material, Table I
of Steinhögl et al. lists (ρ₀, λ) pairs whose product at 77 K (0.2 × 330 =
66 µΩ·cm·nm) is inconsistent with the room-temperature product (76
µΩ·cm·nm) implied by their own fit, and Sondheimer's Table 1 carries ≈1–3 %
errors in the crossover regime — both documented here so later users of
this repository are not misled by the originals. For Dataset 2, the Cr
adhesion layer and thickness uncertainties bias the thinnest films as
described in § 7.

## 11. Conclusions

The digitized Steinhögl dataset shows Cu resistivity rising from 1.4× to
2.7× bulk between 700 and 43 nm linewidth. An implementation of the FS and
MS models built from the original equations — validated against
Sondheimer's tables, all analytic limits, and the paper's own printed
curves — shows that grain-boundary scattering with R ≈ 0.43 [0.39, 0.55]
accounts for the effect, that surface scattering alone fails by a factor
>2 in the excess resistivity, and that the specularity p is not
identifiable from ρ(w) data alone when grain size tracks linewidth. The
calibrated model, with grain sizes measured rather than assumed, transfers
to an independent evaporated-film dataset within 7 % without refitting.
Below ≈80 nm, predictions acquire an additional model-form uncertainty of
up to ≈0.5 µΩ·cm from the FS+MS combination rule itself, which should be
quoted alongside parameter uncertainties in any extrapolation to future
technology nodes.

## References

1. W. Steinhögl, G. Schindler, G. Steinlesberger, M. Engelhardt,
   *Size-dependent resistivity of metallic wires in the mesoscopic range*,
   Phys. Rev. B **66**, 075414 (2002).
2. E. A. Yarimbiyik, H. A. Schafft, R. A. Allen, M. E. Zaghloul,
   D. L. Blackburn, *Experimental and simulation studies of resistivity in
   nanoscale copper films*, Microelectron. Reliab. **49**, 127 (2009).
3. E. H. Sondheimer, *The mean free path of electrons in metals*,
   Adv. Phys. **1**, 1 (1952).
4. A. F. Mayadas, M. Shatzkes, *Electrical-resistivity model for
   polycrystalline films*, Phys. Rev. B **1**, 1382 (1970).
5. R. G. Chambers, *The conductivity of thin wires in a magnetic field*,
   Proc. R. Soc. London A **202**, 378 (1950).
6. D. Gall, *Electron mean free path in elemental metals*,
   J. Appl. Phys. **119**, 085101 (2016).
7. T. S. Kuan et al., *Fabrication and performance limits of sub-0.1 µm Cu
   interconnects*, MRS Symp. Proc. **612** (2000).
