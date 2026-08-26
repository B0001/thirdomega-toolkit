"""
Thermal and electrical model of a self-heated wire in a convective bath.

Establishes the forward model V_3w(geometry, drive, material) so that a measured
third harmonic can be turned into a temperature coefficient, and so that a rig
can be sized before anything is purchased.

SCALING CORRECTION
------------------
The thermal time constant of a wire in a convective bath scales as tau ~ r,
NOT as r^2.

For lumped capacitance with convective loss, tau = rho*c*V / (h*A_s), and for a
cylinder V/A_s = r/2, giving

    tau = rho * c * r / (2h)                          [linear in r]

The r^2 scaling belongs to *internal diffusion*, which is irrelevant here: the
Biot number Bi = h*r/(2k) is ~1e-4 for palladium at these radii, so the wire is
isothermal across its section and the lumped model governs.

Consequence: reducing 100 um -> 25 um raises the 2w thermal corner by 4x, not
16x (roughly 1.1 Hz -> 4.4 Hz). Signal scaling is unaffected and remains
V_3w ~ 1/r^5 at fixed current. Anyone planning a frequency budget from an r^2
assumption will overestimate the available band by 4x.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Material:
    name: str
    rho_e: float      # electrical resistivity, ohm*m, at reference T
    alpha: float      # temperature coefficient of resistance, 1/K
    rho_m: float      # mass density, kg/m^3
    c_p: float        # specific heat, J/(kg*K)
    k_th: float       # thermal conductivity, W/(m*K)


PALLADIUM = Material("Pd", 10.8e-8, 3.77e-3, 12023.0, 244.0, 71.8)
PLATINUM = Material("Pt", 10.6e-8, 3.92e-3, 21450.0, 133.0, 71.6)
"""
Platinum is the control specimen. Its resistivity and TCR are within a few
percent of palladium's, so it produces a comparable third harmonic through the
identical measurement chain -- and it does not form a hydride. Any structure
seen on Pt is the instrument.

Match the DIAMETER when using it. tau and dT both depend on geometry, so a
50 um Pt wire against a 25 um Pd wire is not a control, it is a different
experiment.
"""


@dataclass(frozen=True)
class WireGeometry:
    diameter: float   # m
    length: float     # m, the measured span

    @property
    def radius(self) -> float:
        return self.diameter / 2.0

    @property
    def area_cross(self) -> float:
        return np.pi * self.radius**2

    @property
    def area_surface(self) -> float:
        return np.pi * self.diameter * self.length


def resistance(geom: WireGeometry, mat: Material = PALLADIUM) -> float:
    """DC resistance of the span, ohms."""
    return mat.rho_e * geom.length / geom.area_cross


def tau_thermal(geom: WireGeometry, mat: Material = PALLADIUM, h: float = 1000.0) -> float:
    """
    Lumped thermal time constant, seconds.  tau = rho_m * c_p * r / (2h)

    h is the convective coefficient, W/(m^2*K). ~1e3 is typical for a thin wire
    under natural convection in aqueous electrolyte; measure it rather than
    trusting it, via the roll-off (see corner_frequency).
    """
    return mat.rho_m * mat.c_p * geom.radius / (2.0 * h)


def biot_number(geom: WireGeometry, mat: Material = PALLADIUM, h: float = 1000.0) -> float:
    """Bi = h*r/(2k). Must be << 1 for the lumped model above to hold."""
    return h * geom.radius / (2.0 * mat.k_th)


def corner_frequency(geom: WireGeometry, mat: Material = PALLADIUM, h: float = 1000.0) -> float:
    """
    Drive frequency at which the 2w temperature oscillation is 3 dB down.

    The roll-off is in 2w, so the corner sits at 2*w*tau = 1, i.e.
    f_c = 1 / (4*pi*tau).
    """
    return 1.0 / (4.0 * np.pi * tau_thermal(geom, mat, h))


def delta_t_mean(
    i_rms: float, geom: WireGeometry, mat: Material = PALLADIUM, h: float = 1000.0
) -> float:
    """
    Mean (DC) temperature rise of the wire above bath, kelvin.

    For I(t) = I0*cos(wt) with I_rms = I0/sqrt(2), instantaneous power is
    P = I_rms^2 * R * (1 + cos 2wt): the mean rise and the 2w amplitude are
    both I_rms^2 * R / (h*A_s) before roll-off.

    Keep this under ~2 K. Beyond that the alpha*dT linearity behind the whole
    method degrades, and the electrochemistry is no longer at bath temperature.
    """
    return i_rms**2 * resistance(geom, mat) / (h * geom.area_surface)


def delta_t_2w(
    i_rms: float, f0: float, geom: WireGeometry,
    mat: Material = PALLADIUM, h: float = 1000.0,
) -> float:
    """Amplitude of the 2w temperature oscillation, kelvin, including roll-off."""
    tau = tau_thermal(geom, mat, h)
    return delta_t_mean(i_rms, geom, mat, h) / np.sqrt(1.0 + (4.0 * np.pi * f0 * tau) ** 2)


def v_3w(
    i_rms: float, f0: float, geom: WireGeometry,
    mat: Material = PALLADIUM, h: float = 1000.0, alpha: float | None = None,
) -> float:
    """
    Third-harmonic voltage, volts RMS.

        V_3w,rms = I_rms * R0 * alpha * dT_2w / 2

    RMS throughout, matching what instruments report. demod.HarmonicFit returns
    zero-to-peak amplitudes; use .amplitude_rms(k) to compare against this.

    Pass alpha explicitly to model a hydride, whose TCR differs from the pure
    metal and is in fact the measurand.
    """
    a = mat.alpha if alpha is None else alpha
    R0 = resistance(geom, mat)
    return i_rms * R0 * a * delta_t_2w(i_rms, f0, geom, mat, h) / 2.0


def v_1w(i_rms: float, geom: WireGeometry, mat: Material = PALLADIUM) -> float:
    """Fundamental, volts RMS."""
    return i_rms * resistance(geom, mat)


def alpha_from_v3w(
    v3w_rms: float, i_rms: float, f0: float, geom: WireGeometry,
    mat: Material = PALLADIUM, h: float = 1000.0,
) -> float:
    """Invert the forward model for the temperature coefficient."""
    R0 = resistance(geom, mat)
    dT = delta_t_2w(i_rms, f0, geom, mat, h)
    return 2.0 * v3w_rms / (i_rms * R0 * dT)


def sense_resistor_artifact(
    r_sense: float, tcr_sense: float, i_rms: float,
    theta_ja: float = 100.0, tau_sense: float = 10.0, f0: float = 1.0,
) -> float:
    """
    Parasitic 3w fraction generated by the current-sense resistor.

    The sense element carries I^2*R, self-heats at 2w, and if it has any TCR its
    resistance modulates at 2w -- producing a third harmonic in series with, and
    indistinguishable from, the measurand.

    Two defences, both multiplicative:
      1. bulk metal foil (0.2 ppm/K vs 25 ppm/K thin film)  -> ~100x
      2. large thermal mass, so tau_sense >> tau_wire        -> ~2*w*tau_sense

    Defence 2 is the elegant one: it exploits a time-constant separation you
    already have. A potted foil resistor at tau ~ 10 s against a wire at
    tau ~ 20 ms is averaged by more than 2 orders of magnitude at 1 Hz, filtering
    the artifact while leaving the signal untouched.

    Returns the artifact as a fraction of the fundamental.
    """
    dT_sense = i_rms**2 * r_sense * theta_ja
    dT_2w = dT_sense / np.sqrt(1.0 + (4.0 * np.pi * f0 * tau_sense) ** 2)
    return 0.5 * tcr_sense * dT_2w


def rig_report(
    geom: WireGeometry, i_rms: float, f0: float,
    mat: Material = PALLADIUM, h: float = 1000.0,
) -> str:
    """Human-readable sizing check. Run this before ordering anything."""
    R0 = resistance(geom, mat)
    tau = tau_thermal(geom, mat, h)
    fc = corner_frequency(geom, mat, h)
    dT = delta_t_mean(i_rms, geom, mat, h)
    dT2 = delta_t_2w(i_rms, f0, geom, mat, h)
    v3 = v_3w(i_rms, f0, geom, mat, h)
    v1 = v_1w(i_rms, geom, mat)
    bi = biot_number(geom, mat, h)

    warn = []
    if dT > 2.0:
        warn.append(f"  !! mean rise {dT:.1f} K exceeds 2 K -- reduce drive current")
    if f0 > fc:
        warn.append(f"  !! drive {f0:.1f} Hz is above the {fc:.1f} Hz corner -- signal rolling off")
    if bi > 0.1:
        warn.append(f"  !! Biot {bi:.2e} too large -- lumped model invalid")

    out = [
        f"{mat.name} wire  d={geom.diameter*1e6:.0f} um  L={geom.length*1e3:.0f} mm",
        f"  R0             {R0:.3f} ohm",
        f"  tau            {tau*1e3:.1f} ms   (Bi = {bi:.1e}, lumped OK)",
        f"  2w corner      {fc:.2f} Hz",
        f"  drive          {i_rms*1e3:.1f} mA rms at {f0:.2f} Hz",
        f"  mean dT        {dT:.2f} K",
        f"  2w dT          {dT2:.2f} K",
        f"  V_1w           {v1*1e3:.2f} mV rms",
        f"  V_3w           {v3*1e6:.1f} uV rms   ({20*np.log10(v3/v1):+.1f} dBc)",
    ]
    return "\n".join(out + warn)
