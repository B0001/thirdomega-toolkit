"""
The first measurement on any new cell: characterise the noise floor before
trusting a single alpha value.

Take records with the drive OFF, then ON at several electrolysis current
densities. Fit the corner, then reconcile it against the wire's own thermal
corner to choose a drive frequency.

    python examples/noise_baseline.py
"""
import numpy as np

from thirdomega import psd, thermal as th

GEOM = th.WireGeometry(diameter=25e-6, length=10e-3)
FS = 200.0


def synth_bubble_noise(f_corner, amplitude, seconds=600.0, fs=FS, seed=0):
    """Lorentzian-shaped multiplicative noise from bubble departure."""
    rng = np.random.default_rng(seed)
    n = int(seconds * fs)
    white = rng.normal(0, 1, n)
    tau = 1.0 / (2 * np.pi * f_corner)
    a = np.exp(-1.0 / (tau * fs))
    out = np.zeros(n)
    for i in range(1, n):                      # one-pole filter -> Lorentzian
        out[i] = a * out[i - 1] + (1 - a) * white[i]
    out *= amplitude / (out.std() or 1.0)
    return out + rng.normal(0, amplitude / 20, n)


f_thermal = th.corner_frequency(GEOM)
print(f"wire thermal corner: {f_thermal:.2f} Hz  (25 um Pd)\n")

for j, (fc, amp) in enumerate([(0.15, 3e-6), (0.6, 1e-5), (2.5, 4e-5)]):
    v = synth_bubble_noise(fc, amp, seed=j)
    base = psd.noise_psd(v, FS, label=f"j{j}")
    r = psd.fit_bubble_corner(base, f_lo=0.02, f_hi=20.0)
    w = psd.frequency_window(r["f_corner"], f_thermal)

    print(f"current density step {j}   (true corner {fc} Hz)")
    print(f"  measured corner   {r['f_corner']:.3f} Hz  via {r['f_corner_method']}")
    print(f"  band rms 0.5-5 Hz {base.band_rms(0.5, 5.0)*1e6:.2f} uV")
    if w["viable"]:
        print(f"  usable window     {w['f_lo']:.2f} - {w['f_hi']:.2f} Hz"
              f"  -> drive at {w['f_recommended']:.2f} Hz\n")
    else:
        print(f"  NO VIABLE WINDOW  bubble corner x10 = {w['f_lo']:.2f} Hz "
              f"exceeds thermal corner {w['f_hi']:.2f} Hz")
        print("  -> go thinner (corner ~ 1/r) or reduce current density\n")
