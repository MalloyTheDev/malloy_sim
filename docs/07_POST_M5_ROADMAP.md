# 07 - Post-M5 Roadmap

M1-M8 are complete, so this is the live roadmap for what comes next.

Treat each entry as its own milestone: start it only when explicitly asked, and
build one milestone at a time. The all-in-one goal (`CLAUDE.md`,
`docs/decisions/0006-multi-domain-dispatch.md`) does not license building ahead.
Breadth is earned by finishing domains, not by scaffolding for them.

## The bar for "done"

A domain is finished when it has all four (`CLAUDE.md` rule 16):

1. validation of its own settings and state, returning status rather than throwing;
2. an invariant or conserved quantity checked by tests;
3. malformed and boundary input tests;
4. at least one scenario template in `scenarios/`.

## Complete

```text
M6  - stabilize N-body diagnostics      [done]
M7  - minimal scenario/config loading   [done]
M8  - simple 2D debug visualization     [done]
```

## Track 1: classical mechanics depth (active)

This is the committed next stretch. It turns one working domain into four and
produces the first real template library.

```text
M9  - collision primitives
M10 - 2D rigid body basics
M11 - ballistics/projectiles
M12 - springs and oscillators
```

## Track 2: the multi-domain shell

Only meaningful once Track 1 has produced a genuine second world type. This is
where the scenario `type` key and the dispatch switch are introduced, per ADR
0006. Introducing it earlier would be speculative infrastructure.

```text
M13 - multi-domain scenario dispatch
M14 - template library hardening
```

## Track 3: wider domains (direction, not commitment)

Candidates, roughly in order of how well they fit the existing foundation.
None of these is scheduled; the order will be revisited after Track 2.

```text
M15 - fluids (SPH, 2D)
M16 - thermodynamics / ideal gas
M17 - electromagnetism (charged particles)
M18 - vehicles, as an application of rigid bodies
```

## Track 4: dimension and presentation

```text
M19 - 3D math
M20 - 3D simulation experiments
M21 - rendering
```

Terminal-first still holds until M21 (ADR 0002). Raylib remains the likely first
visualization choice because it gets pixels on screen quickly without turning
this into a graphics project.

## Track 5: separate

```text
M22 - quantum
```

Quantum should start as a separate repo or isolated prototype. Only merge it as
a module if it develops clean boundaries and does not pollute classical
simulation types.
