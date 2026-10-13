# 03 - Module Boundaries

> All modules below are implemented (M2-M41). These boundaries are in force
> in the shipped code; keep them when extending the project.

> New physics domains follow the same shape: one library, one concrete world
> type per dimension, its own settings type, one test executable. A domain
> depends on `malloy_math`, `malloy_sim_core` and `malloy_time`, plus a support
> library such as `malloy_collide` where it needs geometry. No domain depends on
> another and none derives from a shared base
> (`docs/decisions/0006-multi-domain-dispatch.md`). `malloy_layering_tests`
> enforces this from the build files and the sources (`docs/05`).

The module graph at M41 (an arrow means "links"):

```text
malloy_nbody_terminal (app)
  -> malloy_scenario, malloy_ascii, malloy_nbody, malloy_charges
malloy_viewer (app, built only with MALLOY_BUILD_VIEWER)
  -> malloy_scenario, malloy_view3d, raylib
malloy_scenario
  -> every domain, malloy_collide, malloy_sim_core, malloy_math
malloy_view3d
  -> every domain, malloy_collide, malloy_math
malloy_particles, malloy_rigid
  -> malloy_collide, malloy_sim_core, malloy_time, malloy_math
malloy_nbody, malloy_springs, malloy_charges
  -> malloy_sim_core, malloy_time, malloy_math
malloy_collide, malloy_ascii, malloy_time, malloy_sim_core
  -> malloy_math
```

## `malloy_math`

Responsible for `Real`, `Vec2`, `Vec3` (M19), `Quat` (M20), `Mat3` with a
symmetric eigensolver (`eigen_symmetric`), the conversions `to_quat` and
`to_mat3`, `outer` and `trace` (M25 to M27), scalar/vector helpers, approximate
equality, and finite and squarable checks (`is_squarable`, which rejects a
vector whose square would overflow).

Not responsible for mass, force, velocity semantics, acceleration semantics, time, gravity, collision, rendering, or random numbers (the one generator lives in `malloy_particles` until a second domain needs randomness; ADR 0010).

## `malloy_time`

Responsible for fixed timestep representation, tick count, elapsed simulation time, and positive dt validation.

Used by every world: all ten (`NBodyWorld`, `NBody3DWorld`, `ParticleWorld`, `ParticleWorld3D`, `RigidWorld`, `Rigid3DWorld`, `SpringWorld`, `SpringWorld3D`, `ChargeWorld`, `Charge3DWorld`) hold a `FixedStep` for tick counting and elapsed time rather than reimplementing a counter per domain. `SimulationSettings` remains the public place `dt` is configured; `FixedStep` is the runtime counter built from it.

Not responsible for wall-clock time, frame pacing, sleeping, render interpolation, or physics formulas.

## `malloy_sim_core`

Responsible for tiny shared simulation vocabulary: `SimulationSettings`, `StepResult`, and `StepStatus`. It depends on `malloy_math` and nothing else, and `malloy_layering_tests` fails if that changes.

Not responsible for polymorphic simulation interfaces, generic world base classes, entity systems, schedulers, event buses, service locators, plugin systems, scene graphs, or engine kernels.

## `malloy_nbody`

Responsible for `Body2D`, `NBodySettings`, `NBodyWorld`, pairwise gravity, softening, semi-implicit Euler integration, deterministic update order, N-body-specific validation, and the system diagnostics (energy, momentum, angular momentum). Since M19 it also owns `Body3D` and `NBody3DWorld`, the same physics in three dimensions with the same `NBodySettings`: the dimension is not a domain (ADR 0009).

Not responsible for printing, rendering, app loop policy, config files, input, or GUI.

## `malloy_scenario`

Responsible for parsing a scenario/config text file into a `Scenario` for whichever domain its `type` key names (ten types), and for reporting syntax errors with a line number. It knows every domain's key set, which is the price of a plain dispatch switch and is cheaper than the abstraction it replaces. Since M38 one key pair, `lattice`/`lattice3`, turns a single line into many particles: it is capped against the entity limit before anything is allocated, and expanded only after the whole file is read, so the result does not depend on whether `bounds` comes first. Since M39 it also reads `wall_temperature` and `seed` for both particle types.

It also bounds what a file may ask for, because the scenario file is the one untrusted input: steps, entities, ground planes, total work (steps and reports times a step's work, each pair test weighted by its measured worst-case cost under the M40 solver, issues #46 and #84), potential contact points (at most 2000000, which bounds the solver's memory), reports, and line length are all capped, each refusal naming the line where it can (the limits are listed in `scenario.hpp`). A `type rigid` file that declares both a box and a disc that collides is refused, because box against disc is not implemented ([#48](https://github.com/MalloyTheDev/malloy_sim/issues/48)). Its error messages are printable ASCII, with any token it repeats escaped (issue #28); a UTF-8 byte-order mark is skipped and UTF-16 is refused (issue #32). `entity_count` tells the app whether a parsed file defines anything to simulate. It never throws.

Not responsible for semantic validation (that stays in each domain's own `World::validate()`), physics, printing, rendering, or CLI argument handling. Lines beginning `# check` are comments to the parser; only the scenario tests evaluate them (`docs/05`).

Since M41 it also holds the RUN LAYER (`include/malloy/scenario/run.hpp`): a `make_<type>_world(const Scenario&)` per scenario type, a plain readout struct and `readout()` per world, and `report_line()`, the exact line the terminal app prints. Both front ends call it, each with its own switch over the type (ADR 0006, ADR 0012); there is no base class, variant or registry.

## `malloy_view3d`

The 3D render model (M41), the sibling of `malloy_ascii`. Responsible for turning a 3D world's state into final world-space primitives (spheres, oriented boxes with their corners, plane patches, segments, points, arrows) and triangle meshes, fitting a camera to a frame (grow-only, with context such as walls and ground never enlarging the fit), projecting points for tests, and pacing fixed steps against a wall clock. Not responsible for any graphics library, window, input or drawing, and it never changes a world (a lock-step test pins that).

## `malloy_viewer`

The real-time 3D viewer (M41, ADR 0012), an executable and the only target that uses a graphics library (raylib, in `apps/viewer/window.cpp` alone). Draws the five 3D types, refuses 2D ones, stops at the scenario's `steps`, and restarts by rebuilding the world from the parsed scenario. Not responsible for physics, geometry (that is `malloy_view3d`) or world building (the run layer).

## `malloy_ascii`

Responsible for turning 2D points into a character grid: viewport fitting (`fit_viewport`) and framed grid rendering (`render`).

Not responsible for physics, body types, simulation state, terminal control sequences, color, animation, or any graphics API. It never sees a third dimension: the app flattens 3D positions onto the xy plane before calling it.

## `malloy_collide`

Responsible for 2D collision primitives (`Circle`, `Aabb`, since M16 the `Halfplane`, and since M36 the oriented `Obb2`, whose separating-axis contact with another `Obb2` is the 2D sibling of the `Box3` pair, with `contacts(Obb2, Obb2)` giving the clipped manifold of up to two points and, since M37, `contacts(Obb2, Halfplane)` the corners of a box inside a ground plane), overlap tests, contact data (normal, penetration depth, contact point) with a documented deterministic answer wherever the contact normal is geometrically undefined, and the area properties of a shape: area, centroid, and polar second moment of area about the centroid. Since M22 it also carries the 3D primitives `Sphere` and `Plane3`, since M30 the oriented `Box3`, whose contact with a plane is the manifold of its penetrating corners (a list, not a single point), since M32 the axis-aligned `Aabb3` (the walls a 3D particle box confines to), since M33 the separating-axis contact between two oriented boxes (`overlaps`/`contact(Box3, Box3)`), plus `contacts(Box3, Box3)`, the clipped contact manifold (the incident face clipped to the reference face, up to eight points; a single point for edge against edge), and since M34 the box-against-sphere contact (`contact(Box3, Sphere)` and its reverse), by the box's nearest point to the sphere centre. There is no box-against-circle contact in 2D ([#48](https://github.com/MalloyTheDev/malloy_sim/issues/48)).

Not responsible for density, mass, inertia, bodies, velocity, contact response, integration, broadphase acceleration, or scenario loading. The area properties stop at geometry: turning them into mass properties belongs to whichever domain owns bodies (ADR 0007). Like `malloy_ascii` it works on shapes, not on simulation types, so it never sees a `Body2D`.

## `malloy_particles`

Responsible for `Particle2D`, `ParticleSettings`, `ParticleWorld`, non-rotational contact response (positional correction, an impulse along the contact normal and, since M35, a Coulomb friction impulse along the tangent, clamped to the normal one, on both particle and wall contacts; a particle has no orientation, so it imparts no spin), wall containment, uniform gravity (M12, ballistics), and its own validation and diagnostics including potential and total energy and angular momentum. Since M32 it also has the 3D siblings `Particle3D`, `ParticleSettings3D` and `ParticleWorld3D`: spheres colliding (the M22/M24 geometry) and confined to a 3D `collide::Aabb3` box. It is the last of the five domains to gain a 3D form.

Since M38 it also reads a box of particles as a gas (`thermo.hpp`): each world keeps a per-wall tally of the normal impulse its walls have delivered (`wall_impulse()`, bookkeeping that changes no trajectory and is rolled back with a rolled-back step), and pure functions turn the state and that tally into a temperature, a wall pressure and two relaxation statistics, plus `lattice_gas`/`lattice_gas3d`, a deterministic grid start with no random numbers.

Since M39 its walls can be THERMAL (`wall_temperature` in the settings): a thermal wall re-emits a particle from its temperature's distribution (`emission_normal_speed`, `emission_tangent`), drawing from the world's own seeded `SplitMix64` (`splitmix64.hpp`, ADR 0010), and a per-wall heat tally (`wall_heat()`) records the kinetic energy each thermal wall gave the gas; a specular wall's entry is 0 by definition. The generator is rolled back with a rolled-back step, and a box without a thermal wall never draws.

Not responsible for collision geometry (that is `malloy_collide`), orientation, angular velocity, torque, any random numbers but its own seeded generator (no `<random>` engines or distributions, no global or time-seeded generator; ADR 0010), broadphase acceleration (every pair is visited; a uniform grid is a milestone candidate, [#79](https://github.com/MalloyTheDev/malloy_sim/issues/79)), or scenario loading. It carries its own body type rather than widening `nbody::Body2D`, because each domain owns its concrete state.

## `malloy_rigid`

Responsible for `RigidBody2D` (pose plus mass distribution), mass-property construction from a shape and a density, the parallel-axis theorem, world/local conversions, pose integration, impulse application at a point, and rigid-body diagnostics.

Also responsible for rigid contact response since M14: disc against disc contacts (and, since M36, box against box, when a `RigidBody2D` carries a `half_extents` box collider and collides through `collide::Obb2`, and since M37 box against a `Halfplane` ground; since M40 every manifold point is solved, ADR 0011), impulses that generate torque because they act away from the centre of mass, and immovable bodies represented as infinite mass and inertia. Infinity is per quantity: infinite mass alone is a body that can spin but not translate, infinite inertia alone one that can translate but not spin, and `is_static()` means both. Box against disc is not built: a 2D scenario that mixes a box with a disc that collides is refused at load ([#48](https://github.com/MalloyTheDev/malloy_sim/issues/48)).

Since M16 a world also owns immovable ground planes, carried in `RigidSettings`. They are `collide::Halfplane` values, so the geometry stays in `malloy_collide` and `malloy_rigid` only resolves against it.

Since M17 contacts also carry Coulomb friction: a tangential impulse clamped to the friction coefficient times the normal impulse. It is a contact impulse rather than a persistent force, so it needs no force accumulator and rule 5 is untouched. `malloy_particles` gained the same clamp in M35, on its particle and wall contacts, with no lever arm since a particle has no orientation. Since M40 friction is accumulated per contact point and clamped to that point's accumulated normal impulse, and contacts are solved before positions are integrated, so static friction holds (issue #44).

Since M15 it also owns a uniform gravity field, carried in `RigidSettings` and applied as an acceleration before the position update. Static bodies are skipped, and gravitational potential energy is reported alongside kinetic.

Since M40 each world owns a per-point iterative contact solver (ADR 0011): contacts detected at the start of the step with a speculative margin, an impact phase for restitution, a warm-started rest phase stopped at a tolerance with a fixed cap, and a nonlinear Gauss-Seidel position projection. Its contact points are plain structs local to the world's source file; there is no Constraint type and no solver interface, and the warm-start impulses are world state, rolled back with the bodies.

Not responsible for shape geometry (that is `malloy_collide`), persistent forces, force or torque accumulators, or scenario loading. In `RigidBody2D` the inertia is a scalar and the orientation is a scalar angle; those are the 2D CASES of an inertia tensor and a quaternion, not alternatives to them (`docs/decisions/0009-three-dimensions-are-the-destination.md`).

### The 3D rigid world

Since M20 the module also owns `RigidBody3D` and `Rigid3DWorld`: rotation in three dimensions, with a quaternion orientation and three principal moments of inertia. It lives here rather than in a library of its own because rigid-body dynamics is one domain and the dimension is not a domain, the same reasoning M19 used for 3D gravity. `malloy_rigid` links `malloy_collide` for the sphere, plane and box primitives, the same as the 2D path. Milestone by milestone:

- **M20**: torque-free rotation (Euler's equations, so a body with three distinct moments tumbles).
- **M21**: a constant world-frame TORQUE as a setting (the M15 pattern, applied as one kick per step rather than through a force/torque accumulator; issue #88 split the kick from the free rotation that follows the contact solve), which makes a gyroscope precess.
- **M22**: gravity (an acceleration) and restitution contacts against immovable ground planes, resolved with a normal impulse.
- **M23**: Coulomb friction, a tangential impulse clamped to the normal one, which DOES have a lever arm and so spins a sliding sphere up until it rolls without slipping (at 5/7 of its sliding speed).
- **M24**: sphere-against-sphere collisions, two movable bodies exchanging momentum through the SAME impulse core the ground uses, with the plane expressed as a zero-inverse-mass participant so the formula is written once (ADR 0008).
- **M25**: 3D mass properties. `mass_properties_3d` computes a compound body's mass, centre of mass and principal moments from solid spheres and density, via the parallel-axis theorem and the diagonalization of the resulting inertia tensor. It is a construction layer, the 3D sibling of the 2D `mass_properties`, not a domain: no world, no template.
- **M26**: `SolidBox` as a second mass primitive (three distinct moments and a real orientation, so a single tilted box already has a non-diagonal lab-frame tensor) and `combine`, which assembles a compound body from mass-property sets of any kind by reconstructing each tensor (`math::to_mat3`), shifting it by parallel axis and diagonalizing the sum. `box_inertia` gives a box's moments in its own axis order (issue #21).
- **M27**: `SolidMesh`, the general primitive: the mass properties of a solid bounded by a closed, consistently wound triangle mesh (checked exactly, issue #34), computed as signed-tetrahedron volume integrals, so an axis-aligned box mesh reproduces the M26 box exactly. Its tetrahedra are spanned from the mean of the vertices, so a mesh far from the origin keeps full precision (issue #45).
- **M28**: a constant applied FORCE in `Rigid3DSettings`, the translational half of a wrench whose rotational half is M21's torque. It enters as the acceleration F/m before the position update and acts through the centre of mass, so it makes no torque; like the torque it is a single constant setting, not a force accumulator (rule 5). The `force` scenario key sets it (issue #47), and `applied_force.scn` demonstrates it.
- **M30**: a box collider (`half_extents`, zero meaning the sphere `radius` is used instead), resolved against a ground plane through the manifold of its penetrating corners (since M40 every corner is solved, so a tilted box comes to rest).
- **M33**: box against box, through the same impulse core (ADR 0008), through the clipped manifold (since M40 every point of it is solved, so boxes stack; before that, one contact at its centroid).
- **M34**: box against sphere through the same core (the box's nearest point to the sphere centre), so every 3D movable pair now collides.

`RigidBody3D` still stores three principal moments and an orientation, never a 3x3 tensor: M25 to M27 diagonalize the tensor at construction (ADR 0007). Every `RigidBody3D` is movable; the only fixed 3D geometry is a ground plane ([#74](https://github.com/MalloyTheDev/malloy_sim/issues/74)). Since M40 every manifold point is solved iteratively (ADR 0011), so a box comes to rest and short stacks stand; stacks taller than about five boxes may need more than the solver's sweep cap.

## `malloy_springs`

Responsible for `Spring`, `SpringNetwork`, `SpringBody2D`, the pure `accumulate_spring_forces` kernel, and `SpringWorld`, which composes that kernel with the minimal translational integration needed to make the domain runnable, including a per-spring timestep stability check in `validate()`. Since M31 it also has the 3D siblings `SpringBody3D` and `SpringWorld3D` (and a `Vec3` overload of the kernel); `Spring` and `SpringNetwork` are dimension-agnostic (ids and scalars) and are shared unchanged, since the dimension is not a domain (ADR 0009).

Not responsible for rigid-body attachment points or torque, collision, gravity (a milestone candidate, [#78](https://github.com/MalloyTheDev/malloy_sim/issues/78)), or a generic engine-wide force-provider API. `accumulate_spring_forces` neither integrates nor mutates a body, so a later force producer can reuse the pipeline (ADR 0008).

## `malloy_charges`

Responsible for `ChargedParticle2D`, `ChargeSettings` and `ChargeWorld` (M18): charged particles under the pairwise Coulomb interaction (signed charges, with the same squared softening as gravity), a uniform electric field, and a uniform magnetic field, out of plane and so a scalar in 2D. Since M29 it also has `ChargedParticle3D`, `Charge3DSettings` and `Charge3DWorld`, where the magnetic field is a vector and a charge spirals along it. It owns its validation and diagnostics: momentum, kinetic energy, and the Coulomb and field potential energy (a magnetic field does no work, so it has none).

The step is half the electric kick, an exact rotation of each velocity by the magnetic field, the other half of the kick, then the move (`docs/04`, issue #42): in crossed fields the drift is exactly perpendicular to E and the energy error is a bounded oscillation.

Not responsible for collisions, induced or time-varying fields, radiation, or scenario loading. Like every domain it knows no other domain exists (ADR 0006).

## `malloy_nbody_terminal`

The development base: the app the simulation is built and verified through.
Since M41 the viewer (`malloy_viewer`) is a second front end; neither is the
final product, so nothing a later GUI would need may live only here.

Responsible for hardcoded demo setup (the built-in scenarios it runs with no argument), reading a single positional scenario path from argv, refusing a file that fails to parse or defines nothing to simulate, dispatching with its own plain switch over the scenario type to the run layer, which builds each world and formats the report line, calling library APIs, fixed number of steps, terminal output (one diagnostics line per report, with the seed in the run header when a wall is thermal), view framing policy (which viewport each frame shows, and the orthographic flattening of 3D positions onto the xy plane), and returning nonzero on a parse, validation or step failure.

Not responsible for reusable physics logic, flags or a CLI parser, or evaluating `# check` lines ([#64](https://github.com/MalloyTheDev/malloy_sim/issues/64)).

## Fuzz targets

`fuzz/` holds two libFuzzer harnesses, `malloy_parse_scenario_fuzzer` and `malloy_app_fuzzer`. They are not modules: they are built only with `MALLOY_BUILD_FUZZERS=ON` and clang, link the scenario loader (and, for the app fuzzer, the app's own `main.cpp` and `malloy_ascii`), and add nothing to any library (`docs/06`).
