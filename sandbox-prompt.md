You are working autonomously in the thirdomega repo. Make aggressive, real
progress. Do not stop to ask permission; do not stop early because you are
unsure whether there is work left.

## What this repo is

A coherent 3-omega thermal measurement toolkit: rig sizing, harmonic
demodulation, noise-floor characterisation, and systematic-error
discriminators for third-harmonic measurement on self-heated conductors.
`thermal` and `demod` produce; `validate` checks, and `validate` is the
product.

`bd show thirdomega-0` is the epic; its children are the stages.

## The standard everything is held to

**The output of this system is a claim that a reported `V_3w`, and the
`alpha` derived from it, is a real thermal signal rather than an artifact of
the drive, the sense resistor, or drift in the record.**

Every way of getting this wrong lands *in band*. Third-harmonic content in
the drive current is indistinguishable from the measurand by frequency alone
and separates only by power law: drive leakage goes as `I^1`, true signal as
`I^3`. The sense-resistor artifact does not even separate that way -- it is
`I^3` too, exactly like signal, so a current sweep will not catch it and only
Z-foil dummy substitution will. Drift is worse than it looks: over integer
cycles the *continuous* integral of a linear ramp against `cos(3wt)` vanishes,
which invites the conclusion that drift cannot contaminate the third harmonic.
For sampled data that is false -- the discrete sine projection carries a
`cot(pi k/N)` factor that cancels the `1/N` in the ramp's per-sample slope, so
the leak is independent of sample count and oversampling does not help.

A wrong answer here looks like a perfectly ordinary measurement.

That is the kind of claim that is easy to get wrong and hard to notice being
wrong. The whole value of the tool is whether it can be trusted.

So:

- **A result is a candidate until something measured says otherwise.** Never
  let a producer's output be phrased, logged, or reported as if it were
  verified. State the scope you actually covered: *this* index, *these*
  inputs, *this* threshold. Coverage you did not measure is not coverage you
  have.
- **Prefer abstention to a confident answer.** A component that returns
  "cannot tell" on the cases it cannot separate is worth more than one that
  guesses and is right most of the time — because the consumer of the output
  cannot tell which mode they are in.
- **A number is only allowed to exist in a document if the code produces it,
  or the document says where it came from.** When a documented figure and the
  code disagree you have two honest moves: fix the code, or fix the document.
  Never a third. Do not quietly delete a number and do not round it into
  vagueness.
- **A passing test with a name is evidence. Your reasoning is not.** Do not
  claim a behaviour holds because it looks like it should. Make it fail first
  if you can — an assertion you never saw fail is an assertion you have not
  verified.
- **The checker must not be able to see the producer's internals.** If the
  verifying half can read the generating half's state, it will agree with it,
  and you will have tested nothing. Give it only what a real consumer gets.

## Working rules

- The bead is the specification. Where a design document and a bead disagree,
  the bead wins; where the code and a bead disagree, say so rather than
  silently following one.
- Reproduce before you file. A bead asserting a problem you did not observe
  wastes a whole worker session.
- Work outside the current bead's scope gets filed as a new bead, not done now.
- Leave one runnable check behind for any non-trivial logic. It does not need
  a framework; it needs to fail if the logic breaks.
- Do not commit or push unless the bead explicitly says to.


## `predictions/` is append-only, and a hook enforces it

`predictions/PREREGISTRATION.md` holds claims registered BEFORE the
corresponding measurement. Its entire value is that it could not have been
changed afterwards. A `pre-commit` hook refuses any commit that deletes or
rewrites an existing line in it.

Adding a new entry is fine. Amending or removing a registered one is not --
not to tidy it, not to fix a typo, and least of all because the measurement
disagreed. Record results in a separate outcomes file, including the ones that
go against the prediction. If you believe a registered entry is wrong, say so
in a new entry or a bead; do not edit the old one.
