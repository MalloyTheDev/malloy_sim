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
M9  - collision primitives              [done]
M10 - colliding particles               [done]
M11 - 2D rigid body basics              [done]
M12 - ballistics/projectiles            [done]
```

## Track 1: classical mechanics depth (active)

This is the committed next stretch. It turns one working domain into four and
produces the first real template library.

M9 shipped as a support library rather than a domain: `malloy_collide` works on
shapes, not bodies, so it has no world and no scenario template, and rule 16
does not apply to it any more than it does to `malloy_ascii`. M10 is what makes
collision demonstrable.

```text
M13 - springs and oscillators
M?? - rigid-body contact response       (promised, never numbered)
```

M12 shipped as uniform gravity in `malloy_particles` rather than a separate
library. A 2D projectile is a colliding particle under gravity, and a dedicated
domain would have duplicated almost all of `malloy_particles` for one extra
setting. Gravity is world configuration, not body state, so no body type
changed and ADR 0006 is unaffected. It defaults to zero, so every earlier
scenario runs bit-identically.

M11 promised that rigid-body contact response would be "its own milestone", but
no such milestone was ever numbered. It is the gap between M9 (collision
geometry) and M11 (rigid bodies): the two exist and nothing connects them. It
also needs static, infinite-mass bodies, which M11 deliberately left
unrepresented.

M11 shipped as free motion plus impulses: no persistent forces, no force or
torque accumulators, and no rigid-body contact response. Forces belong with
ballistics, where gravity is the point, and contact response is its own
milestone.

M10 was re-scoped from rigid bodies to colliding particles, so that collision
became demonstrable one milestone sooner and the project gained a second real
domain. Rigid bodies (orientation, angular velocity, torque) moved to M11.

M11's ownership boundary is settled in advance and recorded in
`docs/decisions/0007-rigid-bodies-own-their-state.md`: a dedicated
`RigidBody2D`, not a widened `Body2D`. 2D only, so scalar angle and scalar
inertia; quaternions and inertia tensors stay with 3D at M19.

## Track 2: the multi-domain shell

The scenario `type` key and the dispatch switch landed in M10, because that is
when a genuine second world type first existed. What remains is hardening.

```text
M14 - template library hardening        (dispatch itself shipped in M10)
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
