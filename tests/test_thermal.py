"""Tests for the wire thermal model, pinned to hand calculations."""
import numpy as np
import pytest

from thirdomega import thermal as th

G100 = th.WireGeometry(diameter=100e-6, length=10e-3)
G25 = th.WireGeometry(diameter=25e-6, length=10e-3)


def test_resistance_hand_values():
    assert th.resistance(G100) == pytest.approx(0.1376, rel=0.01)
    assert th.resistance(G25) == pytest.approx(2.20, rel=0.01)


def test_tau_is_linear_in_radius_not_quadratic():
    """The correction: convective lumped tau ~ r, not r^2."""
    ratio = th.tau_thermal(G100) / th.tau_thermal(G25)
    assert ratio == pytest.approx(4.0, rel=0.01)      # 4x, not 16x


def test_tau_hand_values():
    assert th.tau_thermal(G100) == pytest.approx(0.0733, rel=0.01)
    assert th.tau_thermal(G25) == pytest.approx(0.0183, rel=0.01)


def test_biot_number_validates_lumped_model():
    assert th.biot_number(G25) < 1e-3
    assert th.biot_number(G100) < 1e-3


def test_corner_frequency_hand_values():
    assert th.corner_frequency(G25) == pytest.approx(4.35, rel=0.02)
    assert th.corner_frequency(G100) == pytest.approx(1.09, rel=0.02)


def test_delta_t_at_operating_point():
    assert th.delta_t_mean(0.025, G25) == pytest.approx(1.75, rel=0.02)


def test_v3w_at_operating_point():
    assert th.v_3w(0.025, 1.0, G25) == pytest.approx(177e-6, rel=0.05)


def test_v1w_at_operating_point():
    assert th.v_1w(0.025, G25) == pytest.approx(55e-3, rel=0.02)


def test_signal_scales_as_inverse_r_fifth():
    """V_3w ~ I^3 R^2/(hA) -> 1/r^5 at fixed current, well below both corners."""
    lo = th.v_3w(0.025, 0.02, G100)
    hi = th.v_3w(0.025, 0.02, G25)
    assert hi / lo == pytest.approx(4 ** 5, rel=0.05)


def test_100ma_on_thin_wire_would_cook_it():
    """Porting the thick-wire drive current to 25 um violates linearity."""
    assert th.delta_t_mean(0.100, G25) > 25.0


def test_rig_report_warns_on_excess_heating():
    assert "exceeds 2 K" in th.rig_report(G25, 0.100, 1.0)


def test_rig_report_warns_above_corner():
    assert "rolling off" in th.rig_report(G25, 0.025, 50.0)


def test_alpha_roundtrip():
    v = th.v_3w(0.025, 1.0, G25, alpha=2.1e-3)
    assert th.alpha_from_v3w(v, 0.025, 1.0, G25) == pytest.approx(2.1e-3, rel=1e-9)


def test_thermal_mass_suppresses_sense_artifact():
    fast = th.sense_resistor_artifact(10.0, 25e-6, 0.025, tau_sense=0.01)
    potted = th.sense_resistor_artifact(10.0, 0.2e-6, 0.025, tau_sense=10.0)
    assert fast / potted > 1e4


def test_platinum_matched_diameter_gives_comparable_signal():
    pd = th.v_3w(0.025, 0.5, G25, mat=th.PALLADIUM)
    pt = th.v_3w(0.025, 0.5, G25, mat=th.PLATINUM)
    assert 0.8 < pt / pd < 1.3
