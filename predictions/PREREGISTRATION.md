# Pre-registration

Committed **before** the corresponding measurement. The point is the timestamp:
a prediction registered in advance is worth considerably more than the same
claim made afterwards, and it costs one commit.

Do not amend entries. Add outcomes in a separate `OUTCOMES.md`, including the
ones that go against the prediction.

---

## P1 — Displacement of the TCR extremum from the resistivity maximum

**Registered:** _(fill in on commit)_
**Status:** open

### Claim

For a hydrided metal wire whose resistivity peaks at composition `x_p` (disorder
scattering maximum), the temperature coefficient `alpha = (1/R)(dR/dT)` has its
own extremum at a composition **displaced from** `x_p`, not coincident with it.

### Reasoning

Write `alpha = N/D` with `N = d(rho_ph)/dT` and `D = rho_total`. At the
resistivity maximum `D' = 0`, so

    d(alpha)/dx |_{x_p} = N'/D

`alpha` is therefore still moving where `R` is stationary. Expanding about the
peak with `n = N'/N` and `k` the relative curvature of the resistivity peak,

    delta ~ -n / (2k)

Coincidence requires `N' = 0` at exactly `x_p`, i.e. exact cancellation between
two competing contributions to the phonon-scattering temperature slope: optical
modes introduced by the interstitial sublattice (grows with content) against
band filling that depresses the density of states at the Fermi level
(suppresses all scattering).

### Predicted

| quantity | prediction | confidence |
|---|---|---|
| extrema coincide? | **no** | high — coincidence is a knife-edge |
| `\|delta\|` | ~0.04 in `x` | low — rests on an estimated `n` |
| sign of `delta` | **unknown** | none — registered as undetermined |

### Registered explicitly as NOT predicted

The sign. Which of the two competing terms dominates `N'` at the peak is not
determined by the argument above, and is not something the author has found in
the literature. Recording it as unknown in advance is the point; a sign
recovered after the fact would carry no weight.

### What would falsify this

The two extrema coinciding within measurement resolution, or `|delta|` so small
that the `(R, alpha)` locus cannot be resolved into separate branches in the
composition range of interest.

### Method

Single specimen, monotonic loading under steady current, logging `R` and
`alpha` continuously. The locus is plotted parametrically in the `(R, alpha)`
plane; `delta` is read as the separation between the two extrema along the
loading trajectory. No absolute composition calibration is required to
establish *whether* the extrema are separated — only to express `delta` in
composition units.

### Controls required before the result counts

- Platinum specimen, matched diameter, identical chain — no hydride, so any
  structure is instrumental
- Bulk-foil dummy resistor substitution — the only test that catches the
  sense-resistor artifact, which shares the `I^3` exponent with real signal
- Current-exponent sweep, fitted exponent consistent with 3
- Noise PSD at each current density, with the drive off and on

---

## P2 — Frequency window viability

**Registered:** _(fill in on commit)_
**Status:** open

**Claim:** at 25 um diameter the wire's thermal corner (~4.3 Hz for the
convective coefficient assumed) sits above ten times the bubble-departure
corner, leaving a usable drive band.

**Registered as uncertain:** the convective coefficient `h` is assumed, not
measured. It is extracted from the roll-off itself, so the corner is a
prediction only until the frequency sweep is run.

**Falsified if:** `10 x f_bubble > f_thermal`, i.e. no window exists. In that
case the geometry is wrong and no signal processing rescues it.
