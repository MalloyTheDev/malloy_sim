# MalloySim C++

MalloySim is a from-scratch **all-in-one science physics simulation** in C++20.

The goal is a robust simulator spanning many physics domains, each shipping
ready-to-run templates so you can pick a scenario and get correct physics
immediately. It is meant to be seen: the destination is a simulator rendered
in 3D, with a GUI, UI and UX on top. It is built simulation-first, because the
system has to exist and be proven correct before an interface is built over
it. The terminal app is only the base it is built and verified on: it will not
be the final product.

It is built as many concrete domain simulations behind a thin shared shell
rather than as one generic engine
(`docs/decisions/0006-multi-domain-dispatch.md`), and its domain libraries are
meant to be reused as the foundation of other science projects.

It is not a game engine, a Unity, Godot or Unreal replacement, an editor, or a
plugin or scripting framework, and it is not trying to become one. Nor is it
render-first: it will be rendered, but over a finished simulation, never the
other way round.

**It is intended to become a 3D simulator, and it has started.** M19 added
`Vec3` and gravity in three dimensions; M20 added quaternions and torque-free
rigid-body rotation, the part of 3D that has no 2D form at all; M21 put a
constant torque on that rotation, so a gyroscope precesses; M22 gave a body a
collision radius and bounced it off a ground plane, the project's first 3D
contact; M23 added friction, so a sliding sphere rolls; M24 made two spheres
collide, conserving momentum; M25 computes a body's mass and inertia from its
geometry, M26 does the same for boxes and assembles compound bodies from
primitives, and M27 does it for an arbitrary closed triangle mesh; M28 pushes a
body with a constant applied force, the translational half of a wrench to go
with M21's torque; M29 takes charged particles into three dimensions, where the
magnetic field becomes a vector and a charge spirals along it; M30 gives a rigid
body a box collider, so an oriented box rests and tumbles on the ground; M31
takes spring networks into three dimensions, so a structure can deform in
space; M32 takes colliding particles into three dimensions, the last of the five
domains to gain a 3D form; M33 collides two oriented boxes by the
separating-axis test; and M34 collides a box with a sphere, the last movable
pair. All five domains now run in space as well as the plane, and any two 3D
primitives collide; what is left of the 3D arc is rendering. 3D is the
destination rather than a possibility left open, and each remaining piece is its
own milestone: see `docs/decisions/0009-three-dimensions-are-the-destination.md`.

Five physics domains are complete and shipping, each in two dimensions and in
three: gravity (N-body), colliding particles, rigid bodies, spring networks, and
charged particles in electric and magnetic fields. Along the way came
ballistics; rigid-body contacts with flat ground and between bodies (every
sphere and box pair in 3D, disc against disc and box against box in 2D); Coulomb
friction for rigid bodies and for particles; a constant applied torque on a 3D
body and, through the library, a constant applied force; mass properties
computed from geometry (spheres, boxes, compound bodies and arbitrary closed
triangle meshes); and the start of thermodynamics: the kinetic theory of a
hard-disk gas (its temperature, wall pressure and relaxation to
Maxwell-Boltzmann) and thermal walls that conduct heat through it, drawing from
a seeded generator (milestones M1-M39).

## Locked baseline

- Language: C++20
- Build system: CMake
- Editor: VS Code + CMake Tools
- Compiler first: MSVC on Windows
- Tests: CTest + tiny custom check macros
- Package manager: none until a milestone needs one
- Graphics: none until the rendering milestone
- Simulation timestep: fixed timestep only
- Numeric type: `double`, and no fast-math
- Determinism: the same binary, platform, initial state, fixed timestep, update
  order and seed repeat a run exactly; there is no cross-platform bitwise claim
- Domain model: one concrete world per domain and dimension, no engine kernel

## Status

All milestones M1-M39 are complete. After the locked M1-M5 roadmap came N-body
diagnostics, scenario loading and the ASCII view (M6-M8); the classical
mechanics track (M9-M37), which took in charged particles (M18), a 3D form for
every domain and collisions between any two 3D primitives (M19-M34), particle
friction (M35) and 2D boxes (M36-M37); and the first two milestones of
thermodynamics (M38-M39). The table under "Milestones complete" lists each one.

The project builds clean under MSVC (`/W4 /WX /permissive-`, with Control Flow
Guard and CET shadow-stack compatibility on every target), and all
14 test executables pass via CTest, along with sixteen checks that run the
terminal app itself and require both its exit code and its output (30 CTest
tests in all). CI runs the Debug and Release builds with every test, and a
short fuzz run, on every push and pull request to `main`. The terminal app runs any of the ten scenario
types from a text file, or two built-in N-body scenarios when given no argument,
and reports each domain's diagnostics alongside an ASCII view of the bodies.

| Module | Type | Provides |
|---|---|---|
| `malloy_math` | INTERFACE | `Real = double`, `Vec2`, `Vec3`, `Quat`, `Mat3` + symmetric eigensolver, vector/scalar helpers |
| `malloy_time` | INTERFACE | `FixedStep` (fixed timestep, tick count, elapsed time) |
| `malloy_sim_core` | STATIC | `SimulationSettings`, `StepStatus`, `StepResult` |
| `malloy_nbody` | STATIC | `Body2D`, `NBodySettings`, `NBodyWorld`, softened gravity, diagnostics; `Body3D` and `NBody3DWorld` |
| `malloy_scenario` | STATIC | parse a scenario text file of any of the ten `type`s into settings and initial state, within caps on steps, entities, ground planes, total work and reports |
| `malloy_ascii` | STATIC | fit a viewport to 2D points, render them as a framed character grid |
| `malloy_collide` | STATIC | 2D `Circle`, `Aabb`, `Obb2`, `Halfplane`; 3D `Sphere`, `Plane3`, `Aabb3`, `Box3`; overlap tests and contacts (normal, depth, point), separating-axis tests for box pairs, and clipped contact manifolds for a box on a plane or on another box |
| `malloy_particles` | STATIC | `Particle2D`/`ParticleWorld` and `Particle3D`/`ParticleWorld3D`: contacts with restitution and Coulomb friction, box walls, gravity, per-wall impulse and heat tallies, thermal walls drawing from a seeded `SplitMix64`; gas statistics (`temperature`, `pressure`, `speed_ratio`, `slow_fraction`) and the `lattice_gas` start in `thermo.hpp` |
| `malloy_rigid` | STATIC | `RigidBody2D`/`RigidWorld`: mass properties, pose integration, impulses, disc-disc and box-box contacts, disc and box contacts with ground halfplanes, uniform gravity, friction; `RigidBody3D`/`Rigid3DWorld`: rotation free or under a constant torque, a constant applied force, sphere and box contacts with ground planes and with each other, and 3D mass properties from geometry (spheres, boxes, meshes, `combine`) |
| `malloy_springs` | STATIC | `Spring`, `SpringNetwork`, deterministic force accumulation, `SpringWorld` and `SpringWorld3D` |
| `malloy_charges` | STATIC | `ChargedParticle2D`/`ChargeWorld` and `ChargedParticle3D`/`Charge3DWorld`: signed Coulomb pairs, uniform E and B fields (B a scalar out of the plane in 2D, a vector in 3D) |
| `malloy_nbody_terminal` | EXECUTABLE | the terminal app, the development base (not the final product): runs any scenario file, or the built-in N-body demo |

Post-M5 work is intentionally gated -- see `docs/07_POST_M5_ROADMAP.md`.

## Domains and templates

MalloySim grows one finished physics domain at a time. Each domain is its own
library with its own concrete world types (one per dimension), its own settings,
and its own test executable. Domains rest only on the shared support libraries
(`malloy_math`, `malloy_time`, and `malloy_collide` for the two that collide) and
the tiny `malloy_sim_core` vocabulary; none of them knows the others exist, and
`malloy_layering_tests` checks that. There is no simulation base class and no
engine kernel, by design.

The libraries are a product in their own right, not internals of the app: they
are meant to be the foundation of other science projects. So a domain can be
linked on its own, its public headers are its API, and nothing a library needs
lives in the app. Packaging them for another project (CMake install and
export) will be its own milestone, when the first consumer needs it
([#82](https://github.com/MalloyTheDev/malloy_sim/issues/82)).

A domain counts as finished only when it has all four of:

1. validation of its own settings and state, returning status rather than throwing;
2. an invariant or conserved quantity checked by tests;
3. malformed and boundary input tests;
4. at least one scenario template in `scenarios/`.

Templates in `scenarios/` are a first-class deliverable: plain text, documented,
and runnable with the shipped binary. Each states what it demonstrates and its
expected result, and carries at least one `# check` line that the scenario tests
compare against a real run (the app reads it as a comment).

The five physics domains each run in two dimensions and in three, and each of
those worlds has its own scenario `type`, so the dispatch has ten.
31 templates ship across those ten types, which the dispatch treats as
10 domains, and every template is parsed, validated and stepped by the test
suite. Both numbers are checked against the directory by the scenario tests, so
neither can go stale.

| `type` | Templates (`scenarios/<name>.scn`) |
|---|---|
| `nbody` (the default) | `two_body`, `three_body_triangle` |
| `nbody3d` | `inclined_orbits` |
| `particles` | `bouncing_particles`, `projectile_arc`, `sliding_friction`, `hard_disk_gas`, `thermal_relaxation`, `heat_conduction` |
| `particles3d` | `particle_collision3d`, `hard_sphere_gas` |
| `rigid` | `spinning_bodies`, `tumbling_impact`, `dropped_bodies`, `flat_ground`, `rolling_and_slipping`, `box_collision2d`, `box_drop2d` |
| `rigid3d` | `intermediate_axis`, `gyroscope`, `bouncing_sphere`, `rolling_sphere`, `sphere_collision`, `box_drop`, `box_collision`, `box_sphere_collision`, `applied_force` |
| `springs` | `spring_chain` |
| `springs3d` | `spring_tetrahedron` |
| `charges` | `cyclotron` |
| `charges3d` | `magnetic_helix` |

A scenario names its domain with a `type` key, dispatched by a plain switch to
one concrete world per `type`. There is no simulation base class:

```text
type particles
dt 0.004
steps 6000
restitution 1.0
bounds -3.0 -3.0 3.0 3.0
particle 1.0  0.30  -2.0 -2.0   1.30 0.90
```

`type` defaults to `nbody` when absent, so scenarios written before it existed
keep working unchanged.

The classical mechanics track is complete, and electromagnetism joined it in
M18. Thermodynamics opened in M38 the way ballistics arrived in M12, as
measurement on an existing domain rather than a new one: a box of colliding
particles read as a gas, with no new `type` key. M39 gave its walls a
temperature, which needed the project's first random numbers, seeded and owned
by the particle world
(`docs/decisions/0010-randomness-is-seeded-and-owned-by-its-domain.md`). Work is
gated one milestone at a time; see `docs/07_POST_M5_ROADMAP.md`.

## Build

From a Developer PowerShell / Developer Command Prompt with MSVC available:

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

The presets use the Visual Studio 2026 generator (x64), so they need CMake 4.2
or newer (`docs/06_CMAKE_AND_VSCODE_WORKFLOW.md`). A `windows-msvc-release`
preset builds the same targets optimised. Two libFuzzer targets, one for the
scenario parser and one for the whole app, build only with clang++ and
`-DMALLOY_BUILD_FUZZERS=ON`; the recipe is in `fuzz/parse_scenario_fuzzer.cpp`.

## Run the demo

With the Visual Studio generator, the debug executable is under:

```powershell
.\out\build\windows-msvc-debug\Debug\malloy_nbody_terminal.exe
```

It runs two hardcoded scenarios and prints, at each report, one line of conserved
system diagnostics (separation, total energy, total angular momentum) followed by
an ASCII view of the bodies. Sample output (diagnostics only, abridged):

```text
MalloySim nbody_terminal

== two-body orbit (normalized units) ==  bodies=2  dt=0.00100000  steps=10000
step      0   sep  1.00000000   E_total  -5.00000000e-07   L_total   1.00000000e-06
...
step  10000   sep  1.00027148   E_total  -4.99999540e-07   L_total   1.00000000e-06

== equilateral three-body (Lagrange) ==  bodies=3  dt=0.00100000  steps=10000
step      0   sep  1.73205081   E_total  -8.66015288e-01   L_total   2.27951999e+00
...
step  10000   sep  1.73143038   E_total  -8.66015093e-01   L_total   2.27951999e+00
```

The three-body triangle holds its shape (separation ~sqrt(3)) while total energy
and angular momentum stay essentially constant -- the conservation you expect
from semi-implicit (symplectic) Euler.

The conserved quantities print in scientific notation so that drift stays
visible at any magnitude. The two-body demo's energy is of order 1e-6, and under
fixed-point formatting every step showed the same value, which would have looked
perfectly conserved even if it were not. Watch `E_total` oscillate rather than
walk in one direction: a bounded oscillation is the symplectic signature, while
a monotone drift would indicate the integrator had been broken.

### The ASCII view

Each report is followed by a frame of the body positions. This is step 1000
of the two-body orbit, with the sun on the left and the planet a sixth of the
way round:

```text
+-------------------------------------------------------------+
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                *                            |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|     *                                                       |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
+-------------------------------------------------------------+
  view x [-0.10000000, 1.10000000]  y [-1.20000000, 1.20000000]
```

The view is fitted to the bodies at the start of a run and afterwards only
grows, so frames share a scale, a near-stationary body keeps its cell, and no
body is ever silently clipped. The printed extents show when it has grown. A 3D
run is drawn as its orthographic projection onto the xy plane, so motion along z
does not show and bodies at different depths can overlap; drawing a 3D scene
properly waits for the rendering milestone.

### Load a scenario from a file

Pass a scenario file as the single argument (the app has no flags):

```powershell
.\out\build\windows-msvc-debug\Debug\malloy_nbody_terminal.exe scenarios\two_body.scn
```

Scenario files are plain text, and `#` starts a comment. This is
`scenarios/two_body.scn` without its comments:

```text
g 1.0
softening 0.000001
dt 0.001
steps 10000
output_every 1000

body 1.0       0.0  0.0   0.0  0.0
body 0.000001  1.0  0.0   0.0  1.0
```

Every template runs the same way. The keys each `type` accepts, and the limits
the parser enforces (a scenario file is the one untrusted input), are documented
in `include/malloy/scenario/scenario.hpp`. A file that does not parse, defines
nothing to simulate, or fails validation is refused with a message on stderr and
exit code 1. With no argument, the app runs the built-in scenarios shown above.

## Known limitations

Open defects that a run of the shipped binary can hit. Each is tracked on
GitHub, and this list shrinks as they close.

- **A box on the ground does not fully settle.** A box that lands tilted now
  rights itself, and in 2D comes to rest, but in 3D it keeps rocking at the
  scale of one step's gravity until the per-point contact solver lands
  ([#41](https://github.com/MalloyTheDev/malloy_sim/issues/41), [#73](https://github.com/MalloyTheDev/malloy_sim/issues/73)). At restitution 1 a nudged box gains
  energy without bound ([#83](https://github.com/MalloyTheDev/malloy_sim/issues/83)).
- **Static friction creeps on slopes.** A rigid body inside the friction cone on
  an incline still slides down it at a steady speed of g sin(theta) dt, because
  contacts are resolved after the position update
  ([#44](https://github.com/MalloyTheDev/malloy_sim/issues/44)).
- **2D box against disc is not implemented.** A `rigid` scenario that declares
  both a box (`rigid_box`) and a disc that collides (`rigid_body` or
  `rigid_static` with a radius above 0) is refused when it is loaded, until the
  pair exists ([#48](https://github.com/MalloyTheDev/malloy_sim/issues/48)).
- **Small header values print as zero.** Run headers and view extents use eight
  fixed decimals, so a value below 5e-9 (for example `dt 1e-9`) prints as 0
  ([#61](https://github.com/MalloyTheDev/malloy_sim/issues/61)).
- **Long paths do not open.** A scenario file whose full path is over 260
  characters (MAX_PATH) cannot be opened, even with long paths enabled in
  Windows
  ([#62](https://github.com/MalloyTheDev/malloy_sim/issues/62)).

Lower-severity defects, including some malformed inputs that load instead of
being refused, are on the
[issue tracker](https://github.com/MalloyTheDev/malloy_sim/issues).

## Milestones complete

| Milestone | Goal | Status |
|---|---|---|
| M1 | CMake skeleton + smoke test | ✅ Done |
| M2 | `malloy_math`: `Real`, `Vec2`, math tests | ✅ Done |
| M3 | `malloy_time` + tiny `malloy_sim_core` | ✅ Done |
| M4 | `malloy_nbody`: `Body2D`, `NBodyWorld`, gravity, tests | ✅ Done |
| M5 | Terminal N-body sun/planet demo | ✅ Done |
| M6 | N-body system diagnostics + three-body demo | ✅ Done |
| M7 | `malloy_scenario`: scenario/config text loading | ✅ Done |
| M8 | `malloy_ascii`: 2D ASCII debug visualization | ✅ Done |
| M9 | `malloy_collide`: 2D collision primitives and contacts | ✅ Done |
| M10 | `malloy_particles`: colliding particles + multi-domain dispatch | ✅ Done |
| M11 | `malloy_rigid`: 2D rigid bodies, mass properties, impulses | ✅ Done |
| M12 | ballistics: uniform gravity in the particle domain | ✅ Done |
| M13 | `malloy_springs`: spring networks + force accumulation | ✅ Done |
| M14 | rigid-body contact response: statics, torque from impacts | ✅ Done |
| M15 | gravity for rigid bodies: uniform field, potential energy | ✅ Done |
| M16 | halfplanes: true flat ground, floors, walls and ramps | ✅ Done |
| M17 | Coulomb friction: rolling, spin-down, static holding | ✅ Done |
| M18 | `malloy_charges`: charged particles, E and B fields | ✅ Done |
| M19 | `Vec3` and 3D gravity: the first step off the plane | ✅ Done |
| M20 | `Quat` and torque-free 3D rotation: the intermediate axis theorem | ✅ Done |
| M21 | a constant torque on a 3D body: the gyroscope precesses | ✅ Done |
| M22 | a sphere bouncing on a 3D ground plane: restitution in 3D | ✅ Done |
| M23 | Coulomb friction for 3D contacts: a sliding sphere rolls at 5/7 v | ✅ Done |
| M24 | sphere-against-sphere collisions: two bodies exchange momentum | ✅ Done |
| M25 | 3D mass properties: mass and inertia computed from geometry | ✅ Done |
| M26 | box mass properties and combine: compound bodies from primitives | ✅ Done |
| M27 | mesh mass properties: inertia of an arbitrary triangle mesh | ✅ Done |
| M28 | a constant applied force on a 3D body: the other half of a wrench | ✅ Done |
| M29 | charged particles in 3D: the vector Lorentz force, helical motion | ✅ Done |
| M30 | an oriented box resting and tumbling on a ground plane | ✅ Done |
| M31 | spring networks in three dimensions: deformable structures | ✅ Done |
| M32 | colliding particles in three dimensions: a box of bouncing spheres | ✅ Done |
| M33 | box against box in 3D: the separating-axis test | ✅ Done |
| M34 | box against a sphere: the last movable pair in 3D | ✅ Done |
| M35 | Coulomb friction for particle contacts, in 2D and 3D | ✅ Done |
| M36 | 2D oriented boxes: box against box by the separating-axis test | ✅ Done |
| M37 | a 2D box resting on a ground plane | ✅ Done |
| M38 | a hard-disk gas: temperature, wall pressure, relaxation | ✅ Done |
| M39 | thermal walls: heat flow through a gas, seeded randomness | ✅ Done |

## What is planned, and what is not

Three different kinds of "not here yet" are separated below, because a reader
has no way to tell them apart otherwise. The first list is a roadmap, the second
waits for a need, and the third is a set of decisions.

### Planned, and gated behind its own milestone

Intended, not yet built, and never added speculatively
(`docs/07_POST_M5_ROADMAP.md`):

- **rendering in 3D, then a GUI, UI and UX**: what the simulation is built to
  be shown through, each its own milestone
  ([#80](https://github.com/MalloyTheDev/malloy_sim/issues/80), [#81](https://github.com/MalloyTheDev/malloy_sim/issues/81)). Rendering is the remaining piece of
  the 3D arc; M19 to M34 built the rest:
  `Vec3`, `Quat` and `Mat3`, a 3D form of every domain, mass properties from
  geometry, a constant torque and force, and collisions between any two 3D
  primitives (`docs/decisions/0009-three-dimensions-are-the-destination.md`).
  Until then the terminal app is the working base, and the M8 debug view
  stays ASCII text
- **quantum**, further out still, and meant to start as a separate repository
  or an isolated prototype
- an iterative, per-point contact solver, which a resting box or a stack needs
  (box contacts are resolved today as one contact at the penetration-weighted
  centroid of the manifold: a bounce, not a stack), and box against disc in 2D
- persistent forces and force/torque accumulators (gravity, the 3D torque and
  the 3D force are each a single constant setting, not a registered force
  producer), including persistent per-particle forces or springs on particles
  (springs are their own domain)
- the rest of thermodynamics: work done by a moving piston, and partially
  diffuse walls (M38 shipped the kinetic theory of a hard-disk gas and M39
  thermal walls with heat flow)
- vehicles and fluids, which are candidates rather than commitments

### Not until a milestone needs it

Neither on the roadmap nor rejected: each waits for a milestone that genuinely
needs it (`docs/09_MISTAKES_TO_AVOID.md`):

- a package manager (`CLAUDE.md` rule 8)
- threading
- an asset manager, serialization, or a logging framework
- a CLI parser (the app takes one positional path)
- ECS, which is not rejected but is not adopted until real access-pattern
  pressure exists rather than an expectation of it (rule 3)

### Not planned at all

Rejected rather than deferred, because adopting them would change what the
project is:

- a game engine clone, or a Unity, Godot or Unreal replacement
- a render-first project: MalloySim is meant to be rendered, but over a
  finished simulation, so each rendering milestone comes after the physics it
  shows
- an editor
- a plugin framework or a scripting framework
- a simulation base class, engine kernel, scheduler, event bus or service
  locator (rule 12: all-in-one is reached with concrete worlds and a dispatch
  key)

## Before committing

Keep the build and tests green:

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

For the full plan and history, read:

- `CLAUDE.md`
- `CHANGELOG.md`
- `docs/00_START_HERE.md`
- `docs/02_MILESTONE_ROADMAP_M1_M5.md`
- `docs/07_POST_M5_ROADMAP.md`
- `docs/decisions/` (the architecture decision records)
