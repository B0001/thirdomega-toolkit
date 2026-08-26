"""
End-to-end walkthrough on synthetic data: size a rig, generate a record,
demodulate it, and run the systematic tests.

    python examples/synthetic_demo.py
"""
import numpy as np

from thirdomega import demod, thermal as th, validate

GEOM = th.WireGeometry(diameter=25e-6, length=10e-3)
I_RMS, F0, FS = 0.025, 1.0, 1000.0


def synthesise(i_rms, f0, fs, cycles=20, drift_mv=20.0, noise_v=2e-6, seed=1):
    """A record with signal, a realistic baseline ramp, and white noise."""
    t = np.arange(0, cycles / f0, 1.0 / fs)
    w = 2 * np.pi * f0
    v1 = np.sqrt(2) * th.v_1w(i_rms, GEOM)
    v3 = np.sqrt(2) * th.v_3w(i_rms, f0, GEOM)
    v = v1 * np.cos(w * t) + v3 * np.cos(3 * w * t + 0.4)
    v += drift_mv * 1e-3 * (t - t.mean()) / (t.max() - t.mean())
    v += np.random.default_rng(seed).normal(0, noise_v, len(t))
    return t, v


print("=" * 66)
print(th.rig_report(GEOM, I_RMS, F0))
print("=" * 66)

t, v = synthesise(I_RMS, F0, FS)
print("\nexpected V_3w:", f"{th.v_3w(I_RMS, F0, GEOM)*1e6:.1f} uV rms\n")

print("--- fitted WITHOUT drift columns (drift_order=0) ---")
bad = demod.fit(t, v, F0, drift_order=0)
print(f"  3w = {bad.amplitude_rms(3)*1e6:.1f} uV rms   <-- contaminated\n")

print("--- fitted WITH drift columns (drift_order=2) ---")
good = demod.fit(t, v, F0, drift_order=2)
print(good.summary())
print(f"\n  3w = {good.amplitude_rms(3)*1e6:.1f} uV rms")
print(f"  alpha = {th.alpha_from_v3w(good.amplitude_rms(3), I_RMS, F0, GEOM)*1e3:.3f} e-3 /K")
print(f"  (true alpha = {th.PALLADIUM.alpha*1e3:.3f} e-3 /K)")

print("\n--- basis conditioning ---")
for k, val in validate.check_orthogonality(t, F0).items():
    print(f"  {k}: {val}")

print("\n--- current-exponent test ---")
currents = np.linspace(8e-3, 40e-3, 8)
amps = [demod.fit(*synthesise(i, F0, FS, seed=7), F0).amplitude_rms(3) for i in currents]
res = validate.current_exponent_test(currents, amps)
print(" ", res.diagnosis())
