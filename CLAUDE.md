# Project Instructions for AI Agents

## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` for the
full workflow reference and session-close protocol.

- Use `bd` for ALL task tracking. Do not use TodoWrite or markdown TODO lists.
- Use `bd remember` for durable knowledge, not MEMORY.md files.
- Conservative git policy: do not commit, push, or sync without being asked.

## What this repo is

A coherent 3-omega thermal measurement toolkit: rig sizing, harmonic
demodulation, noise-floor characterisation, and systematic-error
discriminators for third-harmonic measurement on self-heated conductors.
Nothing here is specific to one material system.

## The claim this repo makes

That a reported `V_3w` — and the `alpha` derived from it — is a real thermal
signal rather than an artifact of the drive, the sense resistor, or a slow
drift in the record.

That claim is hard to verify because the ways of getting it wrong land
*in band*. Third-harmonic content in the drive current is indistinguishable
from the measurand by frequency alone, and separates only by its power law
(`I^1` for drive leakage, `I^3` for true signal). The sense-resistor artifact
does not even separate that way — it scales as `I^3`, exactly like signal, and
is caught only by dummy substitution. A wrong answer here looks like a
perfectly ordinary measurement.

## Producer / checker split

`thermal` and `demod` are the producer: a forward model and a least-squares
harmonic fit. `validate` is the checker, and it is the actual product — the
current-exponent test, dummy-substitution expectation, control-specimen ratio,
and basis conditioning are what turn a number into a claim someone can rely
on. `validate` works from measured quantities, not from the fit's internals;
a discriminator that could see the producer's state would agree with it and
discriminate nothing.

Where a discriminator genuinely cannot separate two causes, it says so. See
the sense-resistor case in `validate.py` — documenting the blind spot is the
correct outcome, not a gap to paper over.

## `predictions/` is append-only, and this is enforced

`predictions/PREREGISTRATION.md` records claims committed *before* the
corresponding measurement. Its value is entirely in the fact that it could not
have been changed afterwards, so a `pre-commit` hook refuses any commit that
deletes or rewrites an existing line in it. Adding a new entry is fine.

Record results in a separate outcomes file, **including the ones that go
against the prediction**. Never amend a registered entry, and never remove
one — not to tidy it, not to fix a typo in it, and least of all because the
measurement disagreed.

## Build & Test

```bash
uv sync --extra dev
uv run pytest
```


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:6cd5cc61 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Agent Context Profiles

The managed Beads block is task-tracking guidance, not permission to override repository, user, or orchestrator instructions.

- **Conservative (default)**: Use `bd` for task tracking. Do not run git commits, git pushes, or Dolt remote sync unless explicitly asked. At handoff, report changed files, validation, and suggested next commands.
- **Minimal**: Keep tool instruction files as pointers to `bd prime`; use the same conservative git policy unless active instructions say otherwise.
- **Team-maintainer**: Only when the repository explicitly opts in, agents may close beads, run quality gates, commit, and push as part of session close. A current "do not commit" or "do not push" instruction still wins.

## Session Completion

This protocol applies when ending a Beads implementation workflow. It is subordinate to explicit user, repository, and orchestrator instructions.

1. **File issues for remaining work** - Create beads for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Handle git/sync by active profile**:
   ```bash
   # Conservative/minimal/default: report status and proposed commands; wait for approval.
   git status

   # Team-maintainer opt-in only, unless current instructions forbid it:
   git pull --rebase
   git push
   git status
   ```
5. **Hand off** - Summarize changes, validation, issue status, and any blocked sync/commit/push step

**Critical rules:**
- Explicit user or orchestrator instructions override this Beads block.
- Do not commit or push without clear authority from the active profile or the current user request.
- If a required sync or push is blocked, stop and report the exact command and error.
<!-- END BEADS INTEGRATION -->
