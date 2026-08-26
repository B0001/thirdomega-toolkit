"""Tests for coherent harmonic demodulation."""
import numpy as np
import pytest

from thirdomega import demod


def synth(f0=1.0, fs=1000.0, cycles=20, amps=(55e-3, 0.0, 183e-6),
          drift_lin=0.0, drift_quad=0.0, noise=0.0, seed=0):
    t = np.arange(0, cycles / f0, 1 / fs)
    w = 2 * np.pi * f0
    v = np.zeros_like(t)
    for i, A in enumerate(amps, start=1):
        v += A * np.cos(i * w * t)
    span = t[-1] - t[0]
    tau = 2 * (t - t[0]) / span - 1
    v += drift_lin * tau + drift_quad * tau ** 2
    if noise:
        v += np.random.default_rng(seed).normal(0, noise, len(t))
    return t, v


def test_recovers_clean_amplitudes():
    t, v = synth()
    f = demod.fit(t, v, 1.0)
    assert f.amplitude[1] == pytest.approx(55e-3, rel=1e-5)
    assert f.amplitude[3] == pytest.approx(183e-6, rel=1e-3)


def test_dbc_matches_hand_calculation():
    t, v = synth()
    f = demod.fit(t, v, 1.0)
    assert f.dbc(3) == pytest.approx(-49.6, abs=0.3)


def test_linear_drift_leaks_badly_when_not_fitted():
    """
    Continuous orthogonality does NOT carry over to sampled data. The sine
    projection carries a cot(pi k/N) factor and the leak is large.
    """
    t, v = synth(drift_lin=10e-3)
    bare = demod.fit(t, v, 1.0, drift_order=0)
    leak = abs(bare.amplitude[3] - 183e-6) / 183e-6
    assert leak > 0.05


def test_fitting_the_linear_column_removes_the_residue():
    """This, not orthogonality, is why the drift columns are there."""
    t, v = synth(drift_lin=10e-3)
    f = demod.fit(t, v, 1.0, drift_order=1)
    assert f.amplitude[3] == pytest.approx(183e-6, rel=1e-9)


def test_leak_is_independent_of_sample_count():
    """Oversampling does not help: cot(pi k/N) cancels the 1/N in the slope."""
    leaks = []
    for fs in (500.0, 1000.0, 4000.0):
        t, v = synth(fs=fs, drift_lin=10e-3)
        f = demod.fit(t, v, 1.0, drift_order=0)
        leaks.append(abs(f.amplitude[3] - 183e-6))
    assert max(leaks) / min(leaks) < 1.1


def test_quadratic_drift_leaks_but_is_small_and_absorbed():
    """Curvature does project, suppressed by w^2, and vanishes once fitted."""
    t, v = synth(drift_quad=10e-3)
    bare = demod.fit(t, v, 1.0, drift_order=0)
    fitted = demod.fit(t, v, 1.0, drift_order=2)
    leak = abs(bare.amplitude[3] - 183e-6) / 183e-6
    assert leak < 0.05
    assert fitted.amplitude[3] == pytest.approx(183e-6, rel=2e-3)


def test_uncertainty_scales_with_noise():
    t, v = synth(noise=1e-5)
    f = demod.fit(t, v, 1.0)
    assert f.amplitude_sigma[3] > 0
    assert f.amplitude[3] == pytest.approx(183e-6, rel=0.2)
    assert f.snr(3) > 3


def test_phase_recovered():
    t = np.arange(0, 20.0, 1e-3)
    w = 2 * np.pi
    v = 1e-3 * np.cos(3 * w * t + 0.7)
    f = demod.fit(t, v, 1.0, harmonics=(3,))
    assert f.phase[3] == pytest.approx(0.7, abs=1e-3)


def test_design_matrix_stays_conditioned():
    t, _ = synth(cycles=20)
    f = demod.fit(t, np.ones_like(t) * 1e-9, 1.0, drift_order=2)
    assert f.condition_number < 10


def test_rejects_short_windows():
    t, v = synth(cycles=5)
    with pytest.raises(ValueError, match="cycles"):
        demod.fit(t, v, 1.0, min_cycles=20)


def test_rejects_high_drift_order():
    t, _ = synth()
    with pytest.raises(ValueError, match="drift_order"):
        demod.design_matrix(t, 1.0, drift_order=5)


def test_coherent_window_truncates_to_whole_cycles():
    t = np.arange(0, 20.7, 1e-3)
    v = np.zeros_like(t)
    tw, vw, n = demod.coherent_window(t, v, 1.0)
    assert n == 20
    assert tw[-1] <= 20.0 + 1e-9
