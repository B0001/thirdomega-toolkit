# Scope of this repository

This repository contains **general-purpose instrumentation only**: the 3-omega
forward model, harmonic demodulation, noise characterisation, and systematic
tests. All of it is applicable to any self-heated conductor and none of it is
specific to a material system or application.

## Deliberately not included

Analysis that inverts a measured `(R, alpha)` pair to a material state variable
is **not** in this repository. That inversion is the subject of a pending
disclosure and is kept private pending a filing decision.

The split is at a real seam, not an arbitrary one: harmonic extraction is a
generic signal-processing problem with a wide audience, while the inversion is
application-specific. Keep it that way when contributing.

## Why this matters

Publication is a statutory bar. Absolute novelty applies at the EPO and in most
jurisdictions outside the US; the US grace period is twelve months and not
renewable. A public commit carries an unambiguous timestamp and is prior art
against your own later filing.

If you are extending this toolkit toward an application, decide the IP question
*before* the first push, because that decision is made once.
