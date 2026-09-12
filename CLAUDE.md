# CLAUDE.md - MalloySim Project Instructions

This file is for Claude or any AI coding assistant working inside this repository.

## Role

Act as a strict technical co-developer and reviewer.

Your job is to help build MalloySim one milestone at a time without scope creep.

## Project identity

MalloySim is a from-scratch **all-in-one science physics simulation** in C++20.

The goal is a robust, flagship-quality simulator spanning many physics domains,
each shipping ready-to-run templates so someone can pick a scenario and get
correct physics immediately.

It is built as **many concrete domain simulations behind a thin shared shell**,
not as one generic engine. That is the central architectural commitment
(`docs/decisions/0006-multi-domain-dispatch.md`), and it is what makes the
all-in-one goal compatible with rules 11 to 13 below.

Terminal-first. Depth over breadth: a domain is not finished because it runs,
it is finished when it validates its input, conserves what it should, holds
determinism, and ships tested templates.

### What it will never be

- a game engine clone
- a Unity/Godot/Unreal replacement
- a render-first project
- an editor project
- a plugin framework
- a scripting framework

Those are rejected on purpose, not deferred. ECS is a separate case: not
rejected, but not adopted until real access-pattern pressure exists (rule 3).

### Where it is going

**MalloySim is intended to become a 3D simulator, and as of M19 it has
started.** `math::Vec3`, `math::Quat`, `nbody::NBody3DWorld` and
`rigid::Rigid3DWorld` exist and ship with templates, so gravity runs in three
dimensions and rigid bodies rotate there, tumbling freely (M20), turning under
an applied torque (M21), bouncing off ground planes under gravity (M22),
rolling with friction (M23), colliding with each other (spheres in M24,
oriented boxes by the separating-axis test in M33), and pushed by a
constant applied force to go with the torque (M28), and their mass and inertia
can be computed from geometry (M25), for spheres, boxes and arbitrary triangle
meshes (M27) assembled into compound bodies (M26). All five domains now have a
3D form as well as a 2D one, and boxes collide with boxes; what remains of the
3D arc is rendering, and then quantum, each its own milestone.

Quantum remains a further destination and has not started
(`docs/decisions/0009-three-dimensions-are-the-destination.md`).

Rules 6 and 7 below still gate what follows. Each remaining piece of 3D is its
own milestone, and the 2D types stay supported rather than being replaced.

## Current phase

```text
M1-M33 complete: math, time, sim_core, N-body, terminal demo, diagnostics,
scenario loading, ASCII debug view, collision primitives, colliding particles
with multi-domain scenario dispatch, 2D rigid bodies, ballistics, spring
networks with deterministic force accumulation, rigid-body contact
response, uniform gravity for rigid bodies, halfplane ground, Coulomb
friction, charged particles in electric and magnetic fields, the first
three-dimensional domain, quaternions with torque-free 3D rotation, a constant
applied torque on a 3D body, a sphere bouncing on a 3D ground plane, Coulomb
friction for 3D contacts, sphere-against-sphere collisions, 3D mass
properties (inertia from geometry), box mass properties with compound
assembly, mesh mass properties for arbitrary shapes, a constant applied
force on a 3D body, charged particles in three dimensions, an oriented box
resting and tumbling on a ground plane, spring networks in three dimensions,
colliding particles in three dimensions, and box-against-box collision.
```

Active track: **classical mechanics depth**, now complete. M9 (collision
geometry), M10 (colliding particles, which introduced the multi-domain `type`
key), M11 (2D rigid bodies), M12 (ballistics, as uniform gravity in the particle
domain), M13 (spring networks and deterministic force accumulation), M14
(rigid-body contact response) and M15 (uniform gravity in the rigid domain)
M16 (halfplanes, giving true flat ground), M17 (Coulomb friction) and M18
(charged particles, the first domain with a velocity-dependent force), M19
(Vec3 and 3D gravity, the first step off the plane), M20 (quaternions and
torque-free rotation in three dimensions), M21 (a constant applied torque on a
3D body), M22 (a sphere bouncing on a 3D ground plane), M23 (Coulomb friction
for 3D contacts), M24 (sphere-against-sphere collisions), M25 (3D mass
properties: mass and inertia computed from geometry), M26 (box mass
properties and combine: compound bodies assembled from primitives), M27
(mesh mass properties: the inertia of an arbitrary triangle mesh), M28 (a
constant applied force on a 3D body, the translational half of a wrench) and M29
(charged particles in three dimensions, with the full vector Lorentz force and
helical motion) and M30 (an oriented box resting and tumbling on a ground
plane, the first non-sphere 3D collision), M31 (spring networks in three
dimensions: deformable structures), M32 (colliding particles in three
dimensions, the last domain to gain a 3D form) and M33 (box against box, by the
separating-axis test, the first non-sphere 3D body-vs-body collision) are all
done. See `docs/07_POST_M5_ROADMAP.md`.

Do not start any further milestone (M34 or later) unless explicitly asked, and
then work only on that one milestone at a time. The all-in-one goal does not
license building ahead: it is reached one finished domain at a time.

## M1-M5 order (all complete)

```text
M1: CMake skeleton + smoke test        [done]
M2: malloy_math                        [done]
M3: malloy_time + tiny malloy_sim_core [done]
M4: malloy_nbody                       [done]
M5: terminal N-body demo               [done]
```

## Rules

1. Do not jump ahead; work one milestone at a time.
2. Do not add rendering until its dedicated milestone. Terminal-first still holds.
3. Do not add ECS (wait for real access-pattern pressure).
4. Collision geometry landed in M9, non-rotational contact response in M10, rotational (rigid-body) contact response in M14, halfplanes in M16, and Coulomb friction in M17. M22 added the first 3D collision, `collide::Sphere` against `collide::Plane3`, with a normal-impulse response, M23 added Coulomb friction to it (a sphere rolls without slipping at 5/7 of its sliding speed), and M24 added sphere-against-sphere collisions. M30 added `collide::Box3` against a plane (an oriented box's penetrating corners, resolved as one centroid contact), so a box can rest and tumble on the ground, and M33 added box-against-box in 3D by the separating-axis test (`collide::overlaps`/`contact(Box3, Box3)`), reduced to a single contact point (a bounce, not a stack). 3D body against BODY is now spheres or boxes; a box against a SPHERE is still deferred, and 2D body against body is still disc against disc. Particle contacts remain normal-only.
5. Rigid bodies landed in M11, gained contact response in M14, a uniform gravity field in M15, and Coulomb friction in M17. Friction is a contact impulse clamped to the normal impulse, not a persistent force. Gravity is a setting applied as an acceleration, not a force. M21 added a constant applied TORQUE and M28 a constant applied FORCE, each a single constant SETTING applied as forcing (F/m for the force, Euler forcing for the torque), NOT a force/torque ACCUMULATOR: do not add an accumulator that sums many force producers, or persistent per-body forces, to `malloy_rigid` until their dedicated milestone.
6. 3D began in M19 (`math::Vec3`, `nbody::NBody3DWorld`), grew in M20 (`math::Quat`, `rigid::Rigid3DWorld`), M21 (a constant applied TORQUE), M22 (`collide::Sphere`/`Plane3`, gravity and restitution contacts against ground planes), and M23 (Coulomb FRICTION for those contacts, the first 3D contact that imparts spin), M24 (sphere-against-sphere collisions, two movable bodies through one shared impulse core), M25 (`math::Mat3` and 3D mass properties: a compound body's inertia TENSOR, computed from geometry by the parallel-axis theorem and diagonalized to principal moments), M26 (`SolidBox` mass properties and `combine`, so a compound body is assembled from primitives of any kind, plus `math::to_mat3`), M27 (`SolidMesh` mass properties: the inertia of any closed triangle mesh by signed-tetrahedron volume integrals, plus `math::trace`), M28 (a constant applied FORCE, the translational half of a wrench, as the acceleration F/m, plus `total_force_potential3d`), M29 (charged particles in three dimensions: the full vector Lorentz force and helical motion) and M30 (`collide::Box3` against a plane, so an oriented box rests and tumbles on the ground; the box collider is `RigidBody3D::half_extents`). M31 promoted springs to 3D (`SpringWorld3D` beside `SpringWorld`, reusing the dimension-agnostic `Spring`/`SpringNetwork`), M32 promoted colliding particles to 3D (`ParticleWorld3D`, confined to a `collide::Aabb3`), the last domain to gain a 3D form, and M33 added box-against-box collision by the separating-axis test (`collide::contact(Box3, Box3)`, resolved through the shared impulse core as one contact point: a bounce, not a stack). Mass properties now cover an arbitrary shape, a body can be both pushed (M28 force) and turned (M21 torque), and boxes collide with the ground and with each other; a box against a SPHERE, the full box-box contact manifold and an iterative stacking solver are NOT started and each needs its own milestone. The 3D contacts are impulses, gravity is an acceleration, and the force and torque are constant settings applied as forcing (the M10/M15/M17/M21 pattern), NOT force/torque accumulators (rule 5). `RigidBody3D` still stores three principal moments rather than a 3x3 matrix, deliberately: M25 to M27 diagonalize the tensor at construction, so the body never has to carry one. Do not add any of it speculatively; the 2D types stay supported.
7. Do not add quantum until its dedicated milestone.
8. Do not add a package manager unless a milestone explicitly needs one.
9. Do not add Catch2/GoogleTest unless explicitly asked.
10. Do not create empty architecture folders.
11. Do not create speculative abstractions.
12. Do not create `ISimulation`, `WorldBase`, `Engine`, `Scheduler`, `EventBus`, `ServiceLocator`, or plugin systems. All-in-one is reached with concrete worlds and a dispatch key, never with a base class.
13. Keep `malloy_sim_core` tiny. Adding a domain must not grow it.
14. Keep the terminal app dumb.
15. Put reusable physics logic in libraries, not in `main.cpp`.
16. A domain is not done until it has validation, an invariant checked by tests, malformed-input tests, and at least one scenario template.
17. `SpringWorld` is a local composition boundary, not the engine-wide force architecture (`docs/decisions/0008-spring-world-is-a-local-composition-boundary.md`). Do not generalise it into a shared force-provider API until several genuinely different force producers exist. If a third domain independently needs the same translational integration path, reassess extracting a shared integrator.

## Multi-domain architecture

Each physics domain is its own library with its own concrete world type, its own
settings type, and its own test executable. Domains share only `malloy::math`
and the tiny `malloy::sim_core` vocabulary (`SimulationSettings`, `StepStatus`,
`StepResult`). No domain knows that any other domain exists.

Domains are selected by a `type` key in the scenario file, dispatched with a
plain switch to one concrete loader and one concrete world per domain.

Implemented in M10, once `malloy_particles` gave the format a second domain to
dispatch to:

```text
type nbody          # or: particles, rigid, springs, charges, nbody3d,
                    #     rigid3d, charges3d
dt 0.001
steps 10000
```

`type` defaults to `nbody` when absent, so scenarios written before the key
existed keep working unchanged. A key belonging to another domain is a parse
error, so a typo in `type` surfaces immediately instead of silently running the
wrong simulation.

## Template library

`scenarios/` is a first-class deliverable, not a scratch folder. Every template:

- states in a header comment what it demonstrates;
- runs to completion with the shipped binary;
- has a documented expected result;
- carries at least one machine-checkable `# check` line.

The last one is not the same as the third, and both are required. Prose can go
stale silently and has done so twice here. A `# check` line is compared against
a real run by the scenario tests:

```text
# check step <n> <quantity> <value> tol <t>
```

`step 0` is the state before any step, matching the column the app prints. The
tolerance is ABSOLUTE and required rather than defaulted, so a template states
the precision it claims instead of inheriting one.

Quantities available in every domain: `energy`, `kinetic`, `momentum_x`,
`momentum_y`. Also `angular` in nbody and rigid, `elastic` in springs, and
`speed` in charges. Naming a quantity a domain does not have FAILS rather than
being skipped, so a check cannot quietly stop checking.

## Testing policy

The project uses:

```text
CTest + tiny custom CHECK macros (tests/test_check.hpp)
```

Do not replace this with Catch2 or GoogleTest unless explicitly asked.

## Numeric policy

Use:

```cpp
using Real = double;
```

No `float` for simulation positions/velocities.

No fast-math.

No cross-platform bitwise determinism claims.

V1 determinism means:

```text
same binary + same platform + same initial state + same fixed timestep + same update order = repeatable output
```

## Build commands

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

## Expected app output

The terminal app's first line is its identity, followed by the periodic demo
output: one diagnostics line (separation, total energy, total angular momentum)
and one ASCII view frame per report.

```text
MalloySim nbody_terminal
```

## When helping

Before proposing code, identify which milestone the task belongs to.

If the user asks for something outside the current milestone, warn that it is out of scope and suggest the correct milestone.

Do not produce a giant codebase.
