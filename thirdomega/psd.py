"""
Noise baseline characterisation.

The first measurement on any new cell is the noise floor with the drive OFF and
then ON at several current densities. Separating the instrument floor from the
bubble spectrum before either can be blamed for the other is the whole point.

Bubble noise is MULTIPLICATIVE, not additive. Nucleation and detachment modulate
the convective coefficient h; h sets dT; dT *is* the signal. No detrending column
touches it, because it scales the signal rather than adding to it. Detachment
events land at ~Hz -- in band. This is the terminal failure mode of the
experiment and it must be measured, not assumed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, signal


@dataclass
class NoiseBaseline:
    freq: np.ndarray
    psd: np.ndarray               # V^2/Hz
    fs: float
    label: str = ""

    def band_rms(self, f_lo: float, f_hi: float) -> float:
        """Integrated noise in a band, V rms."""
        m = (self.freq >= f_lo) & (self.freq <= f_hi)
        return float(np.sqrt(np.trapezoid(self.psd[m], self.freq[m])))

    def at(self, f: float) -> float:
        return float(np.interp(f, self.freq, self.psd))


def noise_psd(v: np.ndarray, fs: float, nperseg: int | None = None,
              label: str = "") -> NoiseBaseline:
    """Welch PSD of a voltage record."""
    v = np.asarray(v, dtype=float)
    if nperseg is None:
        nperseg = min(len(v), max(256, len(v) // 8))
    f, p = signal.welch(v, fs=fs, nperseg=nperseg, detrend="linear")
    return NoiseBaseline(f, p, fs, label)


def lorentzian_plus_flicker(f, s0, f_c, a_flicker, white):
    """
    Bubble departure spectrum model.

    A Lorentzian of corner f_c (the characteristic departure rate) over a 1/f
    flicker background and a white floor. Fitting f_c tells you directly whether
    a chosen drive frequency clears the bubble band.
    """
    f = np.asarray(f, dtype=float)
    fs = np.where(f <= 0, np.nan, f)
    return s0 / (1.0 + (fs / f_c) ** 2) + a_flicker / fs + white


def _log_bin(f, p, n_bins: int = 40):
    """Median-average a PSD into logarithmically spaced frequency bins."""
    m = f > 0
    f, p = f[m], p[m]
    edges = np.logspace(np.log10(f[0]), np.log10(f[-1]), n_bins + 1)
    idx = np.digitize(f, edges) - 1
    fb, pb = [], []
    for b in range(n_bins):
        sel = idx == b
        if sel.sum() >= 1:
            fb.append(np.sqrt(edges[b] * edges[b + 1]))
            pb.append(np.median(p[sel]))
    return np.asarray(fb), np.asarray(pb)


def fit_bubble_corner(baseline: NoiseBaseline, f_lo: float = 0.05,
                      f_hi: float | None = None) -> dict:
    """
    Extract the Lorentzian corner from a measured PSD.

    Identifiability caveat: at low frequency the Lorentzian plateau and the 1/f
    flicker term are partly degenerate, and if flicker is comparable to the
    plateau at the low end the corner is poorly constrained. Check
    f_corner_sigma before trusting the result, and extend the record to lower
    frequency if it is large.

    Returns the fitted parameters and, most usefully, f_corner -- the frequency
    above which bubble-driven h fluctuations fall away. Choose a drive frequency
    at least a decade above it, subject to staying below the wire's own thermal
    corner. If those two constraints do not leave a window, the geometry is
    wrong and no amount of signal processing will rescue it.
    """
    f_hi = f_hi or baseline.freq[-1]
    m = (baseline.freq >= f_lo) & (baseline.freq <= f_hi) & (baseline.freq > 0)
    f, p = baseline.freq[m], baseline.psd[m]

    # Seed f_c from the data: the frequency where the PSD falls to half its
    # low-frequency plateau. Estimated on a LOG-BINNED, median-averaged copy --
    # a raw Welch periodogram has chi-squared scatter of order 100 percent per
    # bin, and picking the first bin below half-power off unsmoothed data
    # triggers on noise and reports a corner far too low.
    fb, pb = _log_bin(f, p)
    plateau = float(np.median(pb[: max(2, len(pb) // 5)]))
    # Require the crossing to PERSIST for two consecutive bins. With only a
    # handful of Welch segments the low-frequency bins stay noisy even after
    # log-binning, and a single downward fluctuation otherwise triggers a
    # spuriously low corner.
    low = pb <= plateau / 2.0
    persist = low[:-1] & low[1:]
    hits = np.where(persist)[0]
    fc_seed = float(fb[hits[0]]) if len(hits) else float(np.sqrt(f[0] * f[-1]))
    p0 = [plateau, fc_seed, float(np.median(p) * np.median(f)), float(p.min())]
    try:
        # PSDs span many decades; fitting residuals in linear space lets the
        # low-frequency points dominate and pins f_c at its initial guess.
        def log_model(ff, *pars):
            return np.log10(lorentzian_plus_flicker(ff, *pars))

        popt, pcov = optimize.curve_fit(
            log_model, f, np.log10(p), p0=p0, maxfev=40000,
            # Bound f_c within a decade of the model-free half-power estimate.
            # Unbounded, the optimiser trades the Lorentzian against flicker+white
            # and wanders to a spurious minimum even on clean synthetic data.
            bounds=([0, fc_seed / 10.0, 0, 0],
                    [np.inf, fc_seed * 10.0, np.inf, np.inf]),
        )
        perr = np.sqrt(np.diag(pcov))
        ok = True
    except Exception:
        popt = np.array(p0)
        perr = np.full(4, np.nan)
        ok = False

    # The four-parameter fit is unstable even on clean synthetic data: the
    # Lorentzian plateau trades against flicker + white and the optimiser finds
    # spurious minima. The model-free half-power point is both more robust and
    # a direct estimate of the same quantity, so it is the primary answer. The
    # fitted value is reported alongside and used only when it is genuinely
    # constrained.
    constrained = bool(ok and np.isfinite(perr[1]) and perr[1] < 0.5 * popt[1])
    f_corner = float(popt[1]) if constrained else float(fc_seed)

    return {
        "f_corner": f_corner,
        "f_corner_method": "lorentzian_fit" if constrained else "half_power",
        "f_corner_halfpower": float(fc_seed),
        "fit_converged": ok,
        "fit_well_constrained": constrained,
        "fit_f_corner": float(popt[1]),
        "fit_f_corner_sigma": float(perr[1]),
        "s0": float(popt[0]), "flicker": float(popt[2]), "white": float(popt[3]),
        "recommended_drive_hz": 10.0 * f_corner,
    }


def frequency_window(f_bubble_corner: float, f_thermal_corner: float,
                     margin: float = 10.0) -> dict:
    """
    Reconcile the two competing frequency constraints.

    Lower bound: comfortably above the bubble corner.
    Upper bound: at or below the wire's thermal corner, or signal rolls off.

    If lo > hi the geometry cannot support the measurement -- go thinner (raises
    the thermal corner as 1/r) or reduce current density (lowers the bubble
    corner and the bubble amplitude together).
    """
    lo = margin * f_bubble_corner
    hi = f_thermal_corner
    return {
        "f_lo": lo, "f_hi": hi, "viable": lo < hi,
        "f_recommended": float(np.sqrt(lo * hi)) if lo < hi else float("nan"),
    }
