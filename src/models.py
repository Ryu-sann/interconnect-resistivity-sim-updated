"""Size-effect resistivity models for thin metal films and wires.

Implements:
  * Fuchs-Sondheimer (FS) surface scattering for a thin film
    (exact integral, arbitrary specularity p)          -> fs_film()
  * Chambers / MacDonald-Sarginson surface scattering for a
    rectangular wire, diffuse walls (p = 0)            -> chambers_wire_p0()
  * Specularity extension of the p = 0 wire result via
    Sondheimer's series identity (his Eq. (31);
    Steinhoegl et al. 2002, Eq. (3))                   -> wire_p_series()
  * Mayadas-Shatzkes (MS) grain-boundary scattering    -> ms_rho_ratio()
  * Combined FS + MS model, additive in resistivity
    increments (Steinhoegl et al. 2002, Eq. (5))       -> combined_wire()

Conventions
-----------
All routines return the resistivity ratio rho/rho0 (dimensionless, >= 1),
where rho0 is the bulk ("infinite size") resistivity that enters through
the bulk mean free path lambda.  Lengths may be given in any consistent
unit (we use nm throughout).

References
----------
E. H. Sondheimer, Adv. Phys. 1, 1 (1952).           [film integral, Table 1]
R. G. Chambers, Proc. R. Soc. A 202, 378 (1950).    [path-integral formalism]
K. E. MacDonald & K. Sarginson, Proc. R. Soc. A 203, 223 (1950). [square wire]
A. F. Mayadas & M. Shatzkes, Phys. Rev. B 1, 1382 (1970).  [grain boundaries]
W. Steinhoegl et al., Phys. Rev. B 66, 075414 (2002).      [combined model]
"""

from __future__ import annotations

import os

import numpy as np
from functools import lru_cache

__all__ = [
    "RectWireInterp",
    "fs_film",
    "chambers_wire_p0",
    "wire_p_series",
    "ms_rho_ratio",
    "combined_wire",
    "combined_film",
    "WireSizeEffect",
]

@lru_cache(maxsize=8)
def _leggauss_cached(n: int):
    return np.polynomial.legendre.leggauss(n)


# ----------------------------------------------------------------------
# Fuchs-Sondheimer thin film, exact integral, arbitrary p
# ----------------------------------------------------------------------

def _fs_film_sigma_ratio(kappa: float, p: float) -> float:
    """sigma/sigma0 for a film of reduced thickness kappa = t/lambda.

    Sondheimer (1952), Eq. (25):
        sigma/sigma0 = 1 - (3(1-p)/(2 kappa)) *
                       Int_1^inf (1/t^3 - 1/t^5) (1-e^(-kappa t))/(1-p e^(-kappa t)) dt
    """
    if kappa <= 0:
        return 0.0
    if p >= 1.0:
        return 1.0
    # Split the integral: numerical part on t in [1, T], analytic tail.
    # For t > T with kappa*T >> 1 the exponentials vanish and the
    # integrand -> (1/t^3 - 1/t^5), whose tail integral is
    # 1/(2 T^2) - 1/(4 T^4).
    T = max(2.0, 60.0 / kappa)
    n = 1000
    # integrate in u = 1/t for better conditioning near t -> inf
    # t from 1..T  <=>  u from 1..1/T ; dt = -du/u^2
    u0, u1 = 1.0 / T, 1.0
    x, wgt = _leggauss_cached(n)
    u = 0.5 * (u1 - u0) * x + 0.5 * (u1 + u0)
    du = 0.5 * (u1 - u0) * wgt
    t = 1.0 / u
    e = np.exp(-kappa * t)
    integrand = (u - u ** 3) * (1.0 - e) / (1.0 - p * e)  # (1/t^3-1/t^5)*t^2 -> (u - u^3)
    val = float(np.sum(integrand * du))
    tail = 0.5 / T ** 2 - 0.25 / T ** 4  # exact, exp terms negligible by construction
    val += tail
    return 1.0 - 1.5 * (1.0 - p) / kappa * val


def fs_film(thickness, lam, p=0.0):
    """rho/rho0 of a thin film with thickness `thickness`, bulk mean free
    path `lam`, and specularity `p` (fraction of specular surface events).

    Exact Fuchs-Sondheimer result.  Vectorised over `thickness`.
    """
    th = np.atleast_1d(np.asarray(thickness, dtype=float))
    out = np.empty_like(th)
    for i, t in enumerate(th):
        out[i] = 1.0 / _fs_film_sigma_ratio(t / lam, p)
    return out if np.ndim(thickness) else float(out[0])


# ----------------------------------------------------------------------
# Chambers path-integral for a rectangular wire, diffuse walls (p = 0)
# ----------------------------------------------------------------------
#
# sigma/sigma0 = 1 - (3 / (4 pi A)) Int_A dA Int_0^{2pi} dphi  G(s(x,y,phi)/lambda)
# with the polar kernel
#     G(u) = 2 Int_0^{pi/2} sin(th) cos^2(th) exp(-u / sin(th)) d(th),
# G(0) = 2/3, and s(x,y,phi) the in-plane distance from (x,y) to the wire
# boundary along azimuth phi.  theta is measured from the wire axis; the
# printed Eq. (2) of Steinhoegl et al. contains cos^2(phi), which fails the
# solid-angle normalisation (3/(4pi)) Int cos^2(theta) dOmega = 1 and does
# not reproduce their own thick-wire limit; cos^2(theta) does both.

_GRID_U = np.concatenate(([0.0], np.logspace(-4, 3.2, 1200)))


@lru_cache(maxsize=None)
def _g_table():
    nth = 400
    x, wgt = _leggauss_cached(nth)
    th = 0.25 * np.pi * (x + 1.0)
    wth = 0.25 * np.pi * wgt
    s, c = np.sin(th), np.cos(th)
    base = 2.0 * s * c ** 2 * wth
    vals = np.empty(len(_GRID_U))
    for i, u in enumerate(_GRID_U):
        vals[i] = np.sum(base * np.exp(-u / s))
    return vals


def _G(u):
    """Kernel G(u), linear interpolation on a dense log grid."""
    tab = _g_table()
    return np.interp(u, _GRID_U, tab, right=0.0)


@lru_cache(maxsize=None)
def _shat_grid(aspect: float, nx: int = 72, ny: int = 72, nphi: int = 256):
    """Normalised boundary distances for a 1 x aspect rectangle.

    Returns (shat, weights) where shat[i,j,k] is the distance (in units of
    the width w) from interior Gauss node (x_i, y_j) to the boundary along
    azimuth phi_k, and weights are the correspondingly normalised
    quadrature weights with the phi average folded in.
    """
    gx, gwx = _leggauss_cached(nx)
    gy, gwy = _leggauss_cached(ny)
    X = 0.5 * (gx + 1.0)                      # in (0,1): units of w
    Y = 0.5 * (gy + 1.0) * aspect             # in (0,aspect)
    WX = 0.5 * gwx
    WY = 0.5 * aspect * gwy
    phi = (np.arange(nphi) + 0.5) * (2.0 * np.pi / nphi)
    cph, sph = np.cos(phi), np.sin(phi)

    with np.errstate(divide="ignore"):
        dx = np.where(cph > 0, (1.0 - X[:, None]) / cph,
                      np.where(cph < 0, X[:, None] / (-cph), np.inf))
        dy = np.where(sph > 0, (aspect - Y[:, None]) / sph,
                      np.where(sph < 0, Y[:, None] / (-sph), np.inf))
    shat = np.minimum(dx[:, None, :], dy[None, :, :])          # (nx, ny, nphi)
    warea = (WX[:, None] * WY[None, :]) / aspect               # integrates to 1
    wall = np.broadcast_to(warea[:, :, None] / nphi, shat.shape)  # phi average
    return shat.reshape(-1), wall.reshape(-1).copy()


def chambers_wire_p0(w: float, h: float, lam: float) -> float:
    """rho/rho0 for a rectangular wire w x h, diffuse surfaces (p = 0).

    Exact Chambers path integral (equivalent to Steinhoegl Eq. (2) after
    correcting the cos^2 phi misprint to cos^2 theta).
    """
    aspect = h / w
    shat, wt = _shat_grid(round(aspect, 9))
    kappa = w / lam
    deficit = np.sum(wt * _G(kappa * shat))
    sigma_ratio = 1.0 - (3.0 / 2.0) * deficit
    # note: (3/(4 pi A)) Int dA Int dphi -> with wt normalised so that
    # Sum wt = 1 (area x phi average), the prefactor becomes
    # (3/(4 pi)) * 2 pi = 3/2.
    return 1.0 / sigma_ratio


# ----------------------------------------------------------------------
# Specularity series (Sondheimer Eq. (31); Steinhoegl Eq. (3))
# ----------------------------------------------------------------------

def wire_p_series(w: float, h: float, lam: float, p: float,
                  tol: float = 1e-9, kmax: int = 500) -> float:
    """rho/rho0 for a rectangular wire with partially specular walls.

    (sigma/sigma0)_{p} = (1-p)^2 Sum_{k>=1} k p^{k-1} (sigma/sigma0)_{p=0, lambda/k}
    """
    if p <= 0.0:
        return chambers_wire_p0(w, h, lam)
    if p >= 1.0:
        return 1.0
    total = 0.0
    for k in range(1, kmax + 1):
        term = k * p ** (k - 1) / chambers_wire_p0(w, h, lam / k)
        total += term
        # remaining tail bound: terms decay ~ k p^(k-1) * 1
        if k > 10 and k * p ** (k - 1) < tol:
            break
    sigma_ratio = (1.0 - p) ** 2 * total
    return 1.0 / sigma_ratio


def film_p_series(t: float, lam: float, p: float, kmax: int = 500,
                  tol: float = 1e-9) -> float:
    """Same series identity applied to the film; cross-check of fs_film."""
    if p <= 0.0:
        return float(fs_film(t, lam, 0.0))
    if p >= 1.0:
        return 1.0
    total = 0.0
    for k in range(1, kmax + 1):
        total += k * p ** (k - 1) * _fs_film_sigma_ratio(k * t / lam, 0.0)
        if k > 10 and k * p ** (k - 1) < tol:
            break
    return 1.0 / ((1.0 - p) ** 2 * total)


# ----------------------------------------------------------------------
# Mayadas-Shatzkes grain-boundary scattering
# ----------------------------------------------------------------------

def ms_alpha(lam: float, d: float, R: float) -> float:
    """MS parameter alpha = (lambda/d) * R / (1 - R)."""
    return (lam / d) * R / (1.0 - R)


def ms_sigma_ratio(alpha) -> np.ndarray:
    """sigma/sigma0 = f(alpha) = 3 [1/3 - alpha/2 + alpha^2
                                   - alpha^3 ln(1 + 1/alpha)]."""
    a = np.atleast_1d(np.asarray(alpha, dtype=float))
    out = np.empty_like(a)
    small = a < 1e-4
    if np.any(small):
        s = a[small]
        # stable expansion: f = 1 - 3s/2 + 3 s^2 - 3 s^3 ln(1/s) - 3 s^4 ...
        with np.errstate(divide="ignore", invalid="ignore"):
            corr = np.where(s > 0.0, 3.0 * s ** 3 * np.log(1.0 / np.where(s > 0, s, 1.0)), 0.0)
        out[small] = 1.0 - 1.5 * s + 3.0 * s ** 2 - corr
    big = ~small
    ab = a[big]
    out[big] = 3.0 * (1.0 / 3.0 - ab / 2.0 + ab ** 2
                      - ab ** 3 * np.log1p(1.0 / ab))
    return out if np.ndim(alpha) else float(out[0])


def ms_rho_ratio(lam: float, d, R: float):
    """rho/rho0 from grain-boundary scattering alone."""
    d = np.asarray(d, dtype=float)
    return 1.0 / ms_sigma_ratio(ms_alpha(lam, d, R))


# ----------------------------------------------------------------------
# Combined models (Matthiessen-additive resistivity increments)
# ----------------------------------------------------------------------

def combined_wire(w: float, h: float, lam: float, p: float,
                  R: float, d: float) -> float:
    """rho/rho0, Steinhoegl Eq. (5): GB term + FS wire excess."""
    return float(ms_rho_ratio(lam, d, R)) + (wire_p_series(w, h, lam, p) - 1.0)


def combined_film(t: float, lam: float, p: float, R: float, d: float) -> float:
    """rho/rho0 for a film: MS term + FS film excess."""
    return float(ms_rho_ratio(lam, d, R)) + (float(fs_film(t, lam, p)) - 1.0)


# ----------------------------------------------------------------------
# Convenience wrapper with cached kappa-splines per aspect ratio
# ----------------------------------------------------------------------

class WireSizeEffect:
    """Fast evaluator: pre-tabulates the p = 0 Chambers result versus
    kappa = w/lambda for a fixed aspect ratio, then serves the p-series
    and combined model from the spline."""

    def __init__(self, aspect: float, kappa_min=5e-3, kappa_max=2e4, n=90):
        from scipy.interpolate import CubicSpline
        self.aspect = aspect
        kg = np.logspace(np.log10(kappa_min), np.log10(kappa_max), n)
        vals = np.array([1.0 / chambers_wire_p0(1.0, aspect, 1.0 / k)
                         for k in kg])       # sigma/sigma0 at kappa
        self._spl = CubicSpline(np.log(kg), vals)
        self._kmin, self._kmax = kappa_min, kappa_max

    def sigma_ratio_p0(self, kappa: float) -> float:
        if kappa >= self._kmax:
            # thick asymptote: 1 - c/kappa, anchor at kmax
            edge = float(self._spl(np.log(self._kmax)))
            return 1.0 - (1.0 - edge) * self._kmax / kappa
        if kappa <= self._kmin:
            kappa = self._kmin
        return float(self._spl(np.log(kappa)))

    def rho_ratio(self, w: float, lam: float, p: float,
                  tol: float = 1e-9, kmax: int = 500) -> float:
        kap = w / lam
        if p >= 1.0:
            return 1.0
        if p <= 0.0:
            return 1.0 / self.sigma_ratio_p0(kap)
        total = 0.0
        for k in range(1, kmax + 1):
            total += k * p ** (k - 1) * self.sigma_ratio_p0(k * kap)
            if k > 10 and k * p ** (k - 1) < tol:
                break
        return 1.0 / ((1.0 - p) ** 2 * total)

    def combined(self, w: float, lam: float, p: float, R: float,
                 d: float) -> float:
        return float(ms_rho_ratio(lam, d, R)) + (self.rho_ratio(w, lam, p) - 1.0)


# ----------------------------------------------------------------------
# Disk-cached 2-D interpolation of the p = 0 Chambers result over
# (kappa = w/lambda, aspect = h/w) -- makes curve sweeps fast.
# ----------------------------------------------------------------------

class RectWireInterp:
    """sigma/sigma0 (p = 0) on a (log kappa, log aspect) grid, with the
    Sondheimer series and the combined model built on top."""

    def __init__(self, cache_path=None, aspects=None, kappas=None):
        from scipy.interpolate import RectBivariateSpline
        if aspects is None:
            aspects = np.logspace(np.log10(0.18), np.log10(9.0), 28)
        if kappas is None:
            kappas = np.logspace(np.log10(5e-3), np.log10(2e4), 90)
        tab = None
        if cache_path and os.path.exists(cache_path):
            z = np.load(cache_path)
            if (len(z["aspects"]) == len(aspects)
                    and np.allclose(z["aspects"], aspects)
                    and np.allclose(z["kappas"], kappas)):
                tab = z["table"]
        if tab is None:
            tab = np.empty((len(kappas), len(aspects)))
            for j, a in enumerate(aspects):
                for i, k in enumerate(kappas):
                    tab[i, j] = 1.0 / chambers_wire_p0(1.0, a, 1.0 / k)
            if cache_path:
                os.makedirs(os.path.dirname(cache_path), exist_ok=True)
                np.savez(cache_path, aspects=aspects, kappas=kappas,
                         table=tab)
        from scipy.interpolate import CubicSpline
        self._kmax = kappas[-1]
        self._edge = CubicSpline(np.log(aspects), tab[-1, :])
        self._spl = RectBivariateSpline(np.log(kappas), np.log(aspects),
                                        tab, kx=3, ky=3)

    def sigma_p0(self, kappa, aspect):
        la = np.log(aspect)
        if kappa >= self._kmax:
            edge = float(self._edge(la))
            return 1.0 - (1.0 - edge) * self._kmax / kappa
        return float(self._spl(np.log(max(kappa, 5e-3)), la)[0, 0])

    def rho_ratio(self, w, h, lam, p, tol=1e-9, kmax=500):
        kap, a = w / lam, h / w
        if p >= 1.0:
            return 1.0
        if p <= 0.0:
            return 1.0 / self.sigma_p0(kap, a)
        total = 0.0
        for k in range(1, kmax + 1):
            total += k * p ** (k - 1) * self.sigma_p0(k * kap, a)
            if k > 10 and k * p ** (k - 1) < tol:
                break
        return 1.0 / ((1.0 - p) ** 2 * total)

    def combined(self, w, h, lam, p, R, d):
        return float(ms_rho_ratio(lam, d, R)) + (
            self.rho_ratio(w, h, lam, p) - 1.0)
