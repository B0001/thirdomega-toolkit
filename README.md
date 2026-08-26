# thirdomega

Coherent 3-omega thermal measurement toolkit: rig sizing, harmonic demodulation,
noise-floor characterisation, and systematic-error discriminators for
third-harmonic measurement on self-heated conductors.

Nothing here is specific to one material system. If you drive a wire, film, or
microbridge with an AC current and read the third harmonic, this applies.

```python
from thirdomega import thermal as th, demod, psd, validate

geom = th.WireGeometry(diameter=25e-6, length=10e-3)
print(th.rig_report(geom, i_rms=0.025, f0=1.0))   # size it before buying

fit = demod.fit(t, v, f0=1.0)                     # measure it
print(fit.summary())
```

## What it does

**`thermal`** — forward model. Resistance, lumped time constant, thermal corner,
temperature rise, and `V_3w` for a given geometry and drive. `rig_report()`
warns when the drive overheats the specimen or sits above the thermal corner.

**`demod`** — simultaneous least-squares fit over a harmonic basis plus low-order
drift polynomials. All harmonics in one solve, covariance-derived uncertainties,
condition-number reporting, enforced coherent windowing.

**`psd`** — noise-floor characterisation and corner extraction, for deciding a
drive frequency from measurement rather than assumption.

**`validate`** — the systematic-error discriminators. Current-exponent test,
dummy-substitution expectation, control-specimen ratio, basis conditioning.

## Three results worth knowing before you build

**Drift detrending is mandatory, not precautionary.** Over integer cycles the
*continuous* integral of a linear ramp against `cos(3wt)` vanishes, which invites
the conclusion that drift cannot contaminate the third harmonic. For sampled
data that conclusion is false. The discrete sine projection carries a factor
`cot(pi k/N)` which grows as `N/(pi k)`, exactly cancelling the `1/N` in a fixed
ramp's per-sample slope. **The leak is independent of sample count** — a 20 mV
ramp over a 20-cycle record puts ~106 uV into the third harmonic, 58% of a
183 uV signal. Oversampling does not help. Fitting the linear column does, to
machine precision. Pinned in `tests/test_demod.py`.

**The thermal time constant scales as `r`, not `r^2`.** For lumped capacitance
with convective loss, `tau = rho*c*r/(2h)`. The `r^2` law belongs to internal
diffusion, which is irrelevant when `Bi = h*r/(2k) ~ 1e-4`. Consequence: halving
the diameter twice raises the thermal corner 4x, not 16x. Signal scaling is
unaffected and remains `V_3w ~ 1/r^5` at fixed current — which is why the drive
current must come *down* when the specimen gets thinner, or the small-perturbation
premise behind `alpha*dT` fails.

**The current-exponent test does not clear the sense resistor.** True signal
scales as `I^3`; drive-harmonic leakage as `I^1`, so the sweep catches that. But
the sense resistor self-heats at `2w` and generates its own third harmonic with
the *same* `I^3` exponent. Only dummy substitution catches it. Use bulk metal
foil and pot it in thermal mass — a `tau_sense` of ~10 s against a wire at
~20 ms filters the artifact while leaving the signal untouched.

## Install

```bash
pip install -e .
pytest -q          # 39 tests
```

Requires numpy and scipy.

## Acquisition

Deliberately hardware-agnostic — no vendor bindings. Generate drive and sample
from **one clock**; coherent windowing is what makes the harmonic basis well
posed, and it is free. `validate.check_orthogonality()` verifies the fit is
well conditioned for your actual window before you trust it.

## Licence

MIT.
