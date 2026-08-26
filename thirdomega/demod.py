"""
Coherent harmonic demodulation by simultaneous least squares.

Replaces sequential analog lock-in detection with a single linear solve over a
design matrix spanning the harmonic basis plus low-order drift polynomials.

Why this rather than a lock-in:
  - all harmonics recovered in one pass, with a shared noise estimate
  - proper covariance-derived uncertainties on every amplitude
  - baseline drift absorbed by explicit polynomial columns rather than filtered
  - the fundamental is fitted, not fought

Drift and the harmonic basis -- READ THIS
----------------------------------------
A tempting argument runs: over an integer number of drive cycles,

    int_0^T t * cos(w t) * cos(3w t) dt = 0     at  T = 2*pi*N/w

so linear baseline drift is orthogonal to the 3w basis and cannot contaminate
it. The integral is correct. The conclusion is WRONG for sampled data.

For N uniform samples spanning integer cycles, the discrete projections are

    sum_n n*cos(2*pi*k*n/N) = -N/2
    sum_n n*sin(2*pi*k*n/N) = -(N/2) * cot(pi*k/N)

The SINE projection is the problem. cot(pi*k/N) ~ N/(pi*k) grows linearly with
N, which exactly cancels the 1/N in the per-sample slope of a fixed ramp. The
leak is therefore N-INDEPENDENT: oversampling does not help at all.

Magnitude: a 20 mV ramp across a 20-cycle record leaks ~106 uV into the third
harmonic -- 58 percent of a 183 uV signal. Not a rounding error.

So the drift columns are MANDATORY, not a cheap precaution. Fitting the linear
term removes the leak to machine precision. Curvature is separately suppressed
by w^2 and is the smaller worry. Both behaviours are pinned in
tests/test_demod.py; do not set drift_order=0 on real data.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class HarmonicFit:
    """Result of a coherent harmonic fit."""

    harmonics: tuple[int, ...]
    amplitude: dict[int, float]          # volts, zero-to-peak
    amplitude_sigma: dict[int, float]    # 1-sigma, volts
    phase: dict[int, float]              # radians
    coeffs: np.ndarray                   # raw least-squares solution
    residual_rms: float                  # volts
    condition_number: float
    dof: int
    drift_coeffs: np.ndarray = field(default_factory=lambda: np.empty(0))

    def amplitude_rms(self, k: int) -> float:
        return self.amplitude[k] / np.sqrt(2.0)

    def dbc(self, k: int, reference: int = 1) -> float:
        """Amplitude of harmonic k relative to the fundamental, in dBc."""
        return 20.0 * np.log10(self.amplitude[k] / self.amplitude[reference])

    def snr(self, k: int) -> float:
        return self.amplitude[k] / self.amplitude_sigma[k]

    def summary(self) -> str:
        lines = [
            f"residual rms   {self.residual_rms:.3e} V",
            f"condition no.  {self.condition_number:.1f}",
            f"dof            {self.dof}",
        ]
        for k in self.harmonics:
            lines.append(
                f"  {k}w  {self.amplitude[k]:.6e} V "
                f"+/- {self.amplitude_sigma[k]:.2e}  "
                f"({self.dbc(k):+.1f} dBc, SNR {self.snr(k):.0f})"
            )
        return "\n".join(lines)


def design_matrix(
    t: np.ndarray,
    f0: float,
    harmonics=(1, 2, 3),
    drift_order: int = 2,
) -> tuple[np.ndarray, list[str]]:
    """
    Build the least-squares design matrix.

    Columns: [cos(k w t), sin(k w t) for each k] + [P_p(tau) for p = 0..drift_order]

    The drift columns use the time variable rescaled to tau in [-1, 1] and are
    built as Legendre polynomials. Raw powers of t form a Vandermonde block whose
    conditioning degrades sharply with order; Legendre on a symmetric interval is
    near-orthogonal and keeps the normal equations well conditioned.

    drift_order is capped at 3 deliberately. Higher orders begin to compete with
    the harmonics for the same variance and can absorb real signal.
    """
    if drift_order > 3:
        raise ValueError(
            "drift_order > 3 risks absorbing signal; cap at 3. "
            "If the baseline genuinely needs more, the window is too long."
        )

    t = np.asarray(t, dtype=float)
    w = 2.0 * np.pi * f0

    cols, names = [], []
    for k in harmonics:
        cols.append(np.cos(k * w * t))
        names.append(f"cos{k}w")
        cols.append(np.sin(k * w * t))
        names.append(f"sin{k}w")

    span = t[-1] - t[0]
    tau = 2.0 * (t - t[0]) / span - 1.0 if span > 0 else np.zeros_like(t)
    for p in range(drift_order + 1):
        c = np.zeros(p + 1)
        c[p] = 1.0
        cols.append(np.polynomial.legendre.legval(tau, c))
        names.append(f"P{p}")

    return np.column_stack(cols), names


def coherent_window(
    t: np.ndarray, v: np.ndarray, f0: float, min_cycles: int = 20
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Truncate a record to a whole number of drive cycles.

    Coherent windowing is what makes the linear-drift orthogonality above exact
    rather than approximate, and it is free. Always use it.
    """
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    dt = float(np.mean(np.diff(t)))
    duration = t[-1] - t[0] + dt          # record spans one dt past the last sample
    n_cycles = int(np.floor(duration * f0 + 1e-9))

    if n_cycles < min_cycles:
        raise ValueError(
            f"record holds {n_cycles} cycles, need >= {min_cycles}. "
            f"Short windows let the drift polynomials compete with the harmonics."
        )

    keep = t < t[0] + n_cycles / f0 - 1e-12
    return t[keep], v[keep], n_cycles


def fit(
    t: np.ndarray,
    v: np.ndarray,
    f0: float,
    harmonics=(1, 2, 3),
    drift_order: int = 2,
    enforce_coherent: bool = True,
    min_cycles: int = 20,
) -> HarmonicFit:
    """
    Fit harmonic amplitudes and phases to a voltage record.

    Parameters
    ----------
    t, v : sample times (s) and voltages (V)
    f0   : drive fundamental (Hz)

    Returns a HarmonicFit with covariance-derived 1-sigma uncertainties.

    Sign convention: v_k(t) = A_k * cos(k*w*t + phi_k)
    """
    if enforce_coherent:
        t, v, _ = coherent_window(t, v, f0, min_cycles)
    else:
        t = np.asarray(t, dtype=float)
        v = np.asarray(v, dtype=float)

    X, _ = design_matrix(t, f0, harmonics, drift_order)
    n, p = X.shape
    if n <= p:
        raise ValueError(f"{n} samples cannot constrain {p} parameters")

    coeffs, *_ = np.linalg.lstsq(X, v, rcond=None)
    resid = v - X @ coeffs
    dof = n - p
    sigma2 = float(resid @ resid) / dof

    XtX_inv = np.linalg.inv(X.T @ X)
    cov = sigma2 * XtX_inv

    amp, amp_sig, phase = {}, {}, {}
    for i, k in enumerate(harmonics):
        a, b = coeffs[2 * i], coeffs[2 * i + 1]
        A = float(np.hypot(a, b))
        amp[k] = A
        # v = a cos + b sin = A cos(kwt + phi) with phi = atan2(-b, a)
        phase[k] = float(np.arctan2(-b, a))

        if A > 0:
            g = np.array([a / A, b / A])
            sub = cov[np.ix_([2 * i, 2 * i + 1], [2 * i, 2 * i + 1])]
            amp_sig[k] = float(np.sqrt(g @ sub @ g))
        else:
            amp_sig[k] = float(np.sqrt(cov[2 * i, 2 * i]))

    return HarmonicFit(
        harmonics=tuple(harmonics),
        amplitude=amp,
        amplitude_sigma=amp_sig,
        phase=phase,
        coeffs=coeffs,
        residual_rms=float(np.sqrt(sigma2)),
        condition_number=float(np.linalg.cond(X)),
        dof=dof,
        drift_coeffs=coeffs[2 * len(harmonics):],
    )
