# 07 - Post-M5 Roadmap

M1-M41 are complete, so this is the live roadmap for what comes next.

Treat each entry as its own milestone: start it only when explicitly asked, and
build one milestone at a time. The all-in-one goal (`CLAUDE.md`,
`docs/decisions/0006-multi-domain-dispatch.md`) does not license building ahead.
Breadth is earned by finishing domains, not by scaffolding for them.

## The bar for "done"

A domain is finished when it has all four (`CLAUDE.md` rule 16):

1. validation of its own settings and state, returning status rather than throwing;
2. an invariant or conserved quantity checked by tests;
3. malformed and boundary input tests;
4. at least one scenario template in `scenarios/`, with a documented expected
   result and at least one `# check` line (`docs/05`).

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
M19 - Vec3 and 3D gravity               [done]
M20 - Quat and torque-free 3D rotation  [done]
M21 - constant applied torque on a 3D body [done]
M22 - a sphere bouncing on a 3D ground plane [done]
M23 - Coulomb friction for 3D contacts     [done]
M24 - sphere against sphere                [done]
M25 - 3D mass properties (inertia tensor)  [done]
M26 - box mass properties and combine      [done]
M27 - mesh mass properties (any shape)     [done]
M28 - a constant applied force (wrench)    [done]
M29 - charged particles in three dimensions [done]
M30 - an oriented box on a ground plane     [done]
M31 - spring networks in three dimensions   [done]
M32 - colliding particles in three dimensions [done]
M33 - box against box (SAT)                 [done]
M34 - box against a sphere                  [done]
M35 - Coulomb friction for particles        [done]
M36 - 2D oriented boxes (SAT)               [done]
M37 - a 2D box on a ground plane            [done]
M38 - a hard-disk gas (kinetic theory)      [done]
M39 - thermal walls (heat flow)             [done]
M40 - per-point iterative contact solver    [done]
M41 - real-time 3D rendering                [done]
```

## What is next: milestone candidates

Nothing is scheduled. The candidates for the next milestones are open GitHub
issues labelled `milestone candidate`, each with a scope sketch; none is
started until it is explicitly asked for, and it gets a number only then.

| Issue | Candidate | Track |
|---|---|---|
| [#74](https://github.com/MalloyTheDev/malloy_sim/issues/74) | immovable 3D bodies (parity with 2D `rigid_static`) | 4 |
| [#48](https://github.com/MalloyTheDev/malloy_sim/issues/48) | box against disc in 2D (until then, a 2D file that mixes the two is refused) | 1 |
| [#75](https://github.com/MalloyTheDev/malloy_sim/issues/75) | a piston (moving wall) to do work on a gas | 3 |
| [#76](https://github.com/MalloyTheDev/malloy_sim/issues/76) | per-wall accommodation coefficients (partially diffuse thermal walls) | 3 |
| [#77](https://github.com/MalloyTheDev/malloy_sim/issues/77) | a local temperature and density profile measurement | 3 |
| [#79](https://github.com/MalloyTheDev/malloy_sim/issues/79) | a uniform-grid broadphase for the particle domain | 3 |
| [#78](https://github.com/MalloyTheDev/malloy_sim/issues/78) | uniform gravity for spring networks | 1 |
| [#64](https://github.com/MalloyTheDev/malloy_sim/issues/64) | one shared `# check` evaluator for the tests and the app | 2 |
| [#81](https://github.com/MalloyTheDev/malloy_sim/issues/81) | a GUI, UI and UX over the rendering | 4 |
| [#82](https://github.com/MalloyTheDev/malloy_sim/issues/82) | package the libraries for reuse by other projects | 2 |
| [#85](https://github.com/MalloyTheDev/malloy_sim/issues/85) | a rigid broadphase and a runtime contact-work budget, so large rigid scenes run | 4 |

Open bugs are issues too; the ones that qualify a claim are noted where the
claim is made. Fixing a bug is not a milestone, but a fix that changes an integrator or a template's numbers
is a deliberate change recorded in `CHANGELOG.md`.

## Track 1: classical mechanics depth (complete)

This was the committed stretch after M5. It turned one working domain into
five, produced the first real template library, and is now complete; its
remaining depth is in the candidates above (#48, #78), and M40 (the per-point
contact solver, ADR 0011) returned to it.

M9 shipped as a support library rather than a domain: `malloy_collide` works on
shapes, not bodies, so it has no world and no scenario template, and rule 16
does not apply to it any more than it does to `malloy_ascii`. M10 is what makes
collision demonstrable.

Track 1 is complete, including the rigid-body contact response that M11
promised and left unnumbered. It became M14.

M14 gave `RigidBody2D` a collision radius and made infinite mass and inertia
mean immovable, which is the representation M11 deliberately deferred on the
grounds that nothing needed statics until contact response did (ADR 0007).
Contacts were disc against disc: oriented boxes and SAT were left to their own
milestones, which became M33 (3D) and M36 (2D).

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
restitution is. Until M40 static holding was not exact: gravity moved a body
before friction acted, so a body inside the friction cone crept down a slope at
`g sin(theta) dt` per unit time. M40 solves contacts before positions move,
which closed issue #44.

What M17 makes visible is the next pressure point. Friction makes stacking look
achievable, and the solver is a single pass with no iteration, so a stack will
not stand. That is where sequential-impulse iteration and then a constraint
solver start pulling, and a constraint solver is where a `Constraint` base class
grows. Rule 12 exists to stop exactly that, so it is named here and deferred
rather than left to drift in. The box contacts of M30 to M37 sharpened it: even
one box on a floor did not fully come to rest. M40 built the per-point solver
without that base class: plain per-point structs local to each world, warm-
started, stopped at a tolerance with a fixed cap (ADR 0011).

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
no such milestone was ever numbered until M14. It was the gap between M9
(collision geometry) and M11 (rigid bodies), and it also needed static,
infinite-mass bodies, which M11 deliberately left unrepresented.

M11 shipped as free motion plus impulses: no persistent forces, no force or
torque accumulators, and no rigid-body contact response. Forces belong with
ballistics, where gravity is the point, and contact response is its own
milestone.

M10 was re-scoped from rigid bodies to colliding particles, so that collision
became demonstrable one milestone sooner and the project gained a second real
domain. Rigid bodies (orientation, angular velocity, torque) moved to M11.

M11's ownership boundary was settled in advance and recorded in
`docs/decisions/0007-rigid-bodies-own-their-state.md`: a dedicated
`RigidBody2D`, not a widened `Body2D`. 2D only, so scalar angle and scalar
inertia; quaternions and inertia tensors stayed with the 3D milestones.

## Track 2: the multi-domain shell

The scenario `type` key and the dispatch switch landed in M10, because that is
when a genuine second world type first existed. It now selects between ten
scenario types across the five domains. What remains is hardening, and it is
deliberately not numbered: a number written down before the work is assigned is
a number that collides with whatever actually ships. This block used to read
`M14`, which M14 then became something else.

Since M14 the shell has gained, without a milestone number: `# check` lines
compared against a real run, documentation counts derived from the code, the
layering test, resource caps on the scenario file, escaped error messages, and
two fuzz targets (`docs/05`, `docs/06`). The open hardening work is tracked as
issues, among them one shared `# check` evaluator so a user's own scenario
checks can run ([#64](https://github.com/MalloyTheDev/malloy_sim/issues/64)),
checks that cannot detect the behaviour they document
([#67](https://github.com/MalloyTheDev/malloy_sim/issues/67)), and template
coverage gaps ([#70](https://github.com/MalloyTheDev/malloy_sim/issues/70)). Continuous
integration has run the builds, every test and a short fuzz run on every push
since issue #66.

## Track 3: wider domains (direction, not commitment)

Candidates, roughly in order of how well they fit the existing foundation.
None of these is scheduled.

```text
thermodynamics, beyond M39: a piston (#75), accommodation coefficients (#76),
    a temperature and density profile (#77), a particle broadphase (#79)
fluids (SPH, 2D)
vehicles, as an application of rigid bodies
```

Electromagnetism came off this list and shipped as M18. Thermodynamics started
as M38, and like ballistics in M12 it was measurement on an existing domain, not
a new one: a box of `malloy_particles` particles read as a gas, with a
temperature, a wall-impulse tally giving the pressure, relaxation statistics,
and a deterministic `lattice` start. The two hard-gas templates match Henderson
(2D) and Carnahan-Starling (3D) once the walls are accounted for.
M39 added thermal walls, which needed the project's first random numbers: a
seeded SplitMix64 owned by the particle world, consumed in a fixed order and
rolled back with a failed step (ADR 0010). A gas now relaxes to its walls'
temperature and conducts heat between a hot and a cold wall, and the
free-molecular heat flow between two thermal walls is tested against its
closed form. What follows needs a moving wall: a piston, and with it work
(#75). Accommodation coefficients (#76) would change the number of random draws
per bounce, so they are an ADR-level change to the draw-order contract.
Fluids (SPH) would be a genuinely new domain, with a neighbour search that the
O(N^2) pair loop cannot afford; the particle broadphase (#79) is the first step
in that direction.

Unnumbered on purpose. These were once written as M15 to M18, and M15 then
shipped as gravity for rigid bodies, so the labels pointed at the wrong
work. A candidate gets a number when it is started, not before.

## Track 4: dimension and presentation

This is where the project is ultimately headed, and the reason it is worth
saying so is that it changes how the current code should be read.

```text
Vec3 and 3D gravity    [done, M19]
quaternions and torque-free 3D rotation [done, M20]
a constant applied torque on a 3D body [done, M21]
sphere against a ground plane          [done, M22]
friction for 3D contacts                [done, M23]
sphere against sphere                   [done, M24]
a general inertia tensor               [done, M25: compound spheres, parallel axis, diagonalized]
a box mass primitive, and combine      [done, M26: three moments, tilt, compound assembly]
mass properties of an arbitrary mesh   [done, M27: signed-tetrahedron volume integrals]
a translational force on a 3D body     [done, M28: constant force setting, F/m, no accumulator; the `force` key since #47]
an oriented box against a plane        [done, M30: penetrating corners; every corner solved since M40]
box against box (SAT, edge-edge)       [done, M33: the clipped manifold; every point solved since M40]
box against a sphere                    [done, M34: the box's nearest point to the sphere centre, the last movable pair in 3D]
3D versions of every domain            [done: nbody M19, rigid M20, charges M29, springs M31, particles M32]
immovable 3D bodies                     [candidate, #74]
a per-point iterative contact solver    [done, M40: resting boxes and short stacks, ADR 0011]
rendering in 3D         [done, M41: malloy_viewer over the run layer and malloy_view3d, ADR 0012]
a GUI, UI and UX        [candidate, #81: after rendering, each its own milestone]
```

M19 took the first step and deliberately took the smallest one. N-body is the
only domain with no contacts and no orientation, so it needed `Vec3` and
nothing else, which let quaternions and inertia tensors stay deferred rather
than be built speculatively.

M20 scoped the same way and got the sharper physics for it. Torque-free
rotation needs no contacts, no collision geometry, no forces and no solver, so
none of those had to be built ahead of their milestone; what it does need is
the one thing 3D has that 2D does not. `RigidBody3D` therefore stores three
principal moments rather than a 3x3 matrix, and still does: M25 to M27 brought
the general tensor as a construction step, diagonalized before it reaches a
body (ADR 0007).

Unnumbered, like Tracks 2 and 3, for the same reason: a milestone gets a number
when it starts. These were once M19 to M21, and other work has since shipped
into the numbers below them.

**3D is the destination, not a possibility left open.**
`docs/decisions/0009-three-dimensions-are-the-destination.md` records what that
means for the code: which 2D decisions are deliberate specializations of their
3D forms and which would have to be replaced. In short, the scalar angle and
scalar inertia in `malloy_rigid` are the 2D cases of a quaternion and an
inertia tensor, and the integrator, the fixed timestep, the validation and
status model, the scenario dispatch and the impulse formulas carried over
almost unchanged, as the amendments to ADR 0009 record. `malloy_ascii`
and the 2D area properties do not carry over; the app draws a 3D world as an
orthographic flattening onto the xy plane until rendering has its milestone.

It is still gated. Nothing 3D is added before its milestone, and the 2D types
stay supported (`CLAUDE.md` rule 6). Raylib remains the likely first
visualization choice because it gets pixels on screen quickly; the GUI, UI and
UX are then built over the rendering, each its own milestone. The simulation
libraries stay free of any graphics dependency, so the rendering consumes them
the way another project would.

## Track 5: separate

```text
quantum (unnumbered)
```

Quantum should start as a separate repo or isolated prototype. Only merge it as
a module if it develops clean boundaries and does not pollute classical
simulation types. This line once read `M22 - quantum`; M22 shipped as a sphere
bouncing on a 3D ground plane, which is why nothing here carries a number.
