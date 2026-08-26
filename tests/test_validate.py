"""Tests for systematic-error discriminators and noise tooling."""
import numpy as np
import pytest

from thirdomega import psd, validate


def test_cubic_law_passes():
    i = np.linspace(5e-3, 50e-3, 8)
    r = validate.current_exponent_test(i, 2.0 * i ** 3)
    assert r.exponent == pytest.approx(3.0, abs=0.01)
    assert r.passes
    assert "sense resistor" in r.diagnosis()


def test_linear_leakage_fails():
    i = np.linspace(5e-3, 50e-3, 8)
    r = validate.current_exponent_test(i, 1e-3 * i)
    assert not r.passes
    assert "Drive harmonic" in r.diagnosis()


def test_mixed_contributions_land_between():
    i = np.linspace(5e-3, 50e-3, 8)
    r = validate.current_exponent_test(i, 2.0 * i ** 3 + 5e-4 * i)
    assert 1.0 < r.exponent < 3.0


def test_dummy_substitution_ratio():
    assert validate.dummy_substitution_expectation() == pytest.approx(18850, rel=0.01)


def test_platinum_control_close_to_unity():
    assert validate.platinum_control_expectation() == pytest.approx(1.04, abs=0.01)


def test_orthogonality_good_for_long_window():
    t = np.arange(0, 20.0, 1e-3)
    r = validate.check_orthogonality(t, 1.0, drift_order=2)
    assert r["well_posed"]
    assert r["condition_number"] < 10


def test_rejects_too_few_points():
    with pytest.raises(ValueError, match="at least 4"):
        validate.current_exponent_test([1e-3, 2e-3], [1e-9, 8e-9])


def test_psd_recovers_white_floor():
    fs = 1000.0
    v = np.random.default_rng(0).normal(0, 1e-6, 200000)
    b = psd.noise_psd(v, fs)
    expected = (1e-6 ** 2) / (fs / 2)
    assert b.at(100.0) == pytest.approx(expected, rel=0.2)


def test_bubble_corner_recovered_from_synthetic_lorentzian():
    f = np.logspace(-2, 2, 400)
    true_fc = 2.0
    p = psd.lorentzian_plus_flicker(f, 1e-10, true_fc, 1e-14, 1e-14)
    b = psd.NoiseBaseline(f, p, fs=1000.0)
    r = psd.fit_bubble_corner(b)
    assert r["f_corner"] == pytest.approx(true_fc, rel=0.1)
    assert r["recommended_drive_hz"] == pytest.approx(20.0, rel=0.1)


def test_half_power_estimator_is_the_robust_one():
    """The 4-parameter fit is unstable; the model-free estimator is not."""
    f = np.logspace(-2, 2, 400)
    p = psd.lorentzian_plus_flicker(f, 1e-10, 2.0, 1e-14, 1e-14)
    r = psd.fit_bubble_corner(psd.NoiseBaseline(f, p, fs=1000.0))
    assert r["f_corner_halfpower"] == pytest.approx(2.0, rel=0.05)
    assert r["f_corner_method"] in ("half_power", "lorentzian_fit")


def test_frequency_window_viable_case():
    w = psd.frequency_window(f_bubble_corner=0.2, f_thermal_corner=4.35)
    assert w["viable"]
    assert 2.0 < w["f_recommended"] < 4.35


def test_frequency_window_infeasible_case():
    """High bubble corner against a low thermal corner leaves no band."""
    w = psd.frequency_window(f_bubble_corner=1.0, f_thermal_corner=1.09)
    assert not w["viable"]
