"""
thirdomega -- coherent 3-omega thermal measurement toolkit.

A general-purpose instrument stack for third-harmonic thermal measurement on
self-heated conductors. Nothing here is specific to any one material system.

    from thirdomega import demod, thermal, psd, validate

Start with thermal.rig_report() to size a rig, psd for the noise baseline, then
demod.fit() on real records and validate for the systematics.
"""

from . import demod, psd, thermal, validate

__version__ = "0.1.0"
__all__ = ["demod", "thermal", "psd", "validate"]
