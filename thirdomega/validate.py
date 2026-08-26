"""
Systematic-error tests. Run these before believing any alpha value.

The 3w method has one classic way of returning a confident wrong answer:
third-harmonic content in the DRIVE current is in-band and indistinguishable
from the measurand. The discriminator is the power law.

    true 3w signal          scales as I^3   (dT ~ I^2, times I from V = IR)
    drive-harmonic leakage  scales as I^1
    sense-resistor artifact scales as I^3   -- NOT separable by exponent

That last line matters: the sense resistor artifact has the SAME exponent as
signal, so the current sweep will not catch it. It is caught only by the Z-foil
dummy substitution (dummy_substitution_expectation below).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ExponentTest:
    currents: np.ndarray
    amplitudes: np.ndarray
    exponent: float
    exponent_sigma: float
    r_squared: float

    @property
    def passes(self) -> bool:
        """Consistent with a cubic law within 3 sigma, and a good fit."""
        return abs(self.exponent - 3.0) < max(3.0 * self.exponent_sigma, 0.15) \
            and self.r_squared > 0.98

    def diagnosis(self) -> str:
        if self.passes:
            return (f"PASS  exponent {self.exponent:.2f} +/- {self.exponent_sigma:.2f}. "
                    f"Consistent with genuine 3w. Note this does NOT clear the sense "
                    f"resistor, which also scales as I^3 -- run the dummy substitution.")
        if self.exponent < 2.0:
            return (f"FAIL  exponent {self.exponent:.2f}. Drive harmonic distortion is "
                    f"dominating; an exponent near 1 means the current source IS the "
                    f"signal. Fix the source before proceeding.")
        return (f"MARGINAL  exponent {self.exponent:.2f} +/- {self.exponent_sigma:.2f}. "
                f"Mixed contributions. Widen the current range and repeat.")


def current_exponent_test(currents, amplitudes) -> ExponentTest:
    """
    Sweep drive current over ~a decade and fit log(V_3w) vs log(I).

    This is the same sweep that validates loading stationarity in an
    electrochemical cell -- one experiment, two systematics. Run it early.
    """
    i = np.asarray(currents, dtype=float)
    a = np.asarray(amplitudes, dtype=float)
    if np.any(i <= 0) or np.any(a <= 0):
        raise ValueError("currents and amplitudes must be positive")
    if len(i) < 4:
        raise ValueError("need at least 4 current points for a meaningful exponent")

    x, y = np.log(i), np.log(a)
    n = len(x)
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    s_slope = float(np.sqrt(ss_res / (n - 2) / ((x - x.mean()) ** 2).sum())) if n > 2 else np.nan

    return ExponentTest(i, a, float(slope), s_slope, r2)


def dummy_substitution_expectation(tcr_dummy: float = 0.2e-6,
                                   tcr_wire: float = 3.77e-3) -> float:
    """
    Expected signal suppression when the wire is replaced by a Z-foil dummy of
    equal resistance.

    At 0.2 ppm/K against palladium's 3770 ppm/K the true signal drops ~19000x.
    Anything still measured is the instrument: drive harmonics, ground loops,
    ADC nonlinearity, or the sense element.

    This is the only test that catches the sense-resistor artifact, because that
    artifact shares the I^3 exponent with real signal.
    """
    return tcr_wire / tcr_dummy


def platinum_control_expectation(alpha_pt: float = 3.92e-3,
                                 alpha_pd: float = 3.77e-3) -> float:
    """
    Expected V_3w ratio, Pt control against Pd, at matched geometry and drive.

    Within ~4 percent -- close enough that the control exercises the full chain
    at essentially the operating point, and Pt does not hydride. Any loading-
    dependent structure on Pt is instrumental.

    Match the diameter. rho_m*c_p differs between the metals, so tau differs
    (Pt is ~19 percent faster at equal radius); account for that in the roll-off
    when comparing, or drive well below both corners.
    """
    return alpha_pt / alpha_pd


def check_orthogonality(t, f0, harmonics=(1, 2, 3), drift_order: int = 2) -> dict:
    """
    Quantify leakage between drift columns and harmonic columns for a given
    window. Verifies in situ that the fit is well posed before it is trusted.
    """
    from .demod import design_matrix

    X, names = design_matrix(t, f0, harmonics, drift_order)
    Xn = X / np.linalg.norm(X, axis=0, keepdims=True)
    G = Xn.T @ Xn

    n_harm = 2 * len(harmonics)
    cross = np.abs(G[:n_harm, n_harm:])

    return {
        "condition_number": float(np.linalg.cond(X)),
        "max_drift_harmonic_overlap": float(cross.max()),
        "worst_pair": (names[int(np.argmax(cross) // cross.shape[1])],
                       names[n_harm + int(np.argmax(cross) % cross.shape[1])]),
        "well_posed": bool(np.linalg.cond(X) < 100 and cross.max() < 0.1),
    }
