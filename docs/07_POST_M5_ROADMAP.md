# 07 - Post-M5 Roadmap

M1-M18 are complete, so this is the live roadmap for what comes next.

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
M13 - springs and oscillators           [done]
M14 - rigid-body contact response       [done]
M15 - gravity for rigid bodies          [done]
M16 - halfplanes, for true flat ground  [done]
M17 - Coulomb friction                  [done]
M18 - charged particles (E and B)       [done]
```

## Track 1: classical mechanics depth (active)

This is the committed next stretch. It turns one working domain into four and
produces the first real template library.

M9 shipped as a support library rather than a domain: `malloy_collide` works on
shapes, not bodies, so it has no world and no scenario template, and rule 16
does not apply to it any more than it does to `malloy_ascii`. M10 is what makes
collision demonstrable.

Track 1 is complete, including the rigid-body contact response that M11
promised and left unnumbered. It became M14.

M14 gave `RigidBody2D` a collision radius and made infinite mass and inertia
mean immovable, which is the representation M11 deliberately deferred on the
grounds that nothing needed statics until contact response did (ADR 0007).
Contacts are disc against disc: oriented boxes and SAT would be their own
milestone.

M15 added a uniform gravity field to the rigid domain, applied before the
position update so the integration stays semi-implicit Euler. It is an
acceleration rather than a force, so it does not scale with mass and static
bodies are skipped explicitly. Rotation under gravity is emergent: the field
acts through the centre of mass and generates no torque on its own, but a
contact away from the centre of mass does, so bodies rock and tumble without
any new machinery. `malloy_rigid` still has no force or torque accumulators.

M16 added `collide::Halfplane` and ground planes in `RigidSettings`. Before it,
a floor had to be built from overlapping discs, whose contact normal swings by
up to 14.48 degrees as a body moves along it, because it points at whichever disc
centre is nearest. A plane's normal never turns. It is also the only pair in
`malloy_collide` with no degenerate case, since the plane supplies the direction
instead of it being inferred from two centres.

M17 added Coulomb friction, which is what makes the flat normal worth having.
A tangential impulse clamped to the friction coefficient times the normal
impulse gives rolling, spin-down, backspin reversal and static holding, none of
which existed before, and it needs no new state on a body and no force
accumulator: friction is computed and consumed inside a contact exactly as
restitution is.

What M17 makes visible is the next pressure point. Friction makes stacking look
achievable, and the solver is a single pass with no iteration, so a stack will
not stand. That is where sequential-impulse iteration and then a constraint
solver start pulling, and a constraint solver is where a `Constraint` base class
grows. Rule 12 exists to stop exactly that, so it is named here and deferred
rather than left to drift in.

M13 shipped as a separate `malloy_springs` domain rather than as springs inside
`ParticleWorld`, because a spring network is interaction topology rather than
environmental configuration. It also introduced the first many-to-one force
pipeline here. `SpringWorld` is deliberately a local composition boundary and
not the engine-wide force architecture
(`docs/decisions/0008-spring-world-is-a-local-composition-boundary.md`).

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
inertia; quaternions and inertia tensors stay with the 3D milestone.

## Track 2: the multi-domain shell

The scenario `type` key and the dispatch switch landed in M10, because that is
when a genuine second world type first existed. What remains is hardening.

Template library hardening is what remains, and it is deliberately not
numbered: nothing here is scheduled, and a number written down before the
work is assigned is a number that collides with whatever actually ships.
This block used to read `M14`, which M14 then became something else.

## Track 3: wider domains (direction, not commitment)

Candidates, roughly in order of how well they fit the existing foundation.
None of these is scheduled; the order will be revisited after Track 2.

```text
fluids (SPH, 2D)
thermodynamics / ideal gas
vehicles, as an application of rigid bodies
```

Electromagnetism came off this list and shipped as M18.

Unnumbered on purpose. These were once written as M15 to M18, and M15 then
shipped as gravity for rigid bodies, so the labels pointed at the wrong
work. A candidate gets a number when it is started, not before.

## Track 4: dimension and presentation

This is where the project is ultimately headed, and the reason it is worth
saying so is that it changes how the current code should be read.

```text
3D math                (Vec3, quaternions, inertia tensors)
3D simulation          (the existing domains, in three dimensions)
rendering              (terminal-first holds until then, ADR 0002)
```

Unnumbered, like Tracks 2 and 3, for the same reason: a milestone gets a number
when it starts. These were once M19 to M21, and other work has since shipped
into the numbers below them.

**3D is the destination, not a possibility left open.**
`docs/decisions/0009-three-dimensions-are-the-destination.md` records what that
means for the code today: which 2D decisions are deliberate specializations of
their 3D forms and which would have to be replaced. In short, the scalar angle
and scalar inertia in `malloy_rigid` are the 2D cases of a quaternion and an
inertia tensor, and the integrator, the fixed timestep, the validation and
status model, the scenario dispatch and the impulse formulas carry over almost
unchanged. `malloy_ascii` and the disc-only collision geometry do not.

It is still gated. Nothing 3D is added before its milestone, and the 2D
mechanics it rests on are finished first (`CLAUDE.md` rule 6). Raylib remains
the likely first visualization choice because it gets pixels on screen quickly
without turning this into a graphics project.

## Track 5: separate

```text
M22 - quantum
```

Quantum should start as a separate repo or isolated prototype. Only merge it as
a module if it develops clean boundaries and does not pollute classical
simulation types.
