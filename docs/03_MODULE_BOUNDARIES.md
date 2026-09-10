# 03 - Module Boundaries

> All modules below are implemented (M2-M15). These boundaries are in force
> in the shipped code; keep them when extending the project.

> New physics domains follow the same shape: one library, one concrete world
> type, one settings type, one test executable, depending only on
> `malloy_math` and `malloy_sim_core`. No domain depends on another and none
> derives from a shared base (`docs/decisions/0006-multi-domain-dispatch.md`).

## `malloy_math`

Responsible for `Real`, `Vec2`, scalar/vector helpers, approximate equality, and finite checks.

Not responsible for mass, force, velocity semantics, acceleration semantics, time, gravity, collision, or rendering.

## `malloy_time`

Responsible for fixed timestep representation, tick count, elapsed simulation time, and positive dt validation.

Used by every world: `NBodyWorld` and `ParticleWorld` each hold a `FixedStep` for tick counting and elapsed time rather than reimplementing a counter per domain. `SimulationSettings` remains the public place `dt` is configured; `FixedStep` is the runtime counter built from it.

Not responsible for wall-clock time, frame pacing, sleeping, render interpolation, or physics formulas.

## `malloy_sim_core`

Responsible for tiny shared simulation vocabulary: `SimulationSettings`, `StepResult`, and `StepStatus`.

Not responsible for polymorphic simulation interfaces, generic world base classes, entity systems, schedulers, event buses, service locators, plugin systems, scene graphs, or engine kernels.

## `malloy_nbody`

Responsible for `Body2D`, `NBodySettings`, `NBodyWorld`, pairwise gravity, softening, semi-implicit Euler integration, deterministic update order, and N-body-specific validation.

Not responsible for printing, rendering, app loop policy, config files, input, or GUI.

## `malloy_scenario`

Responsible for parsing a scenario/config text file into a `Scenario` for whichever domain its `type` key names, and for reporting syntax errors with a line number. It knows every domain's key set, which is the price of a plain dispatch switch and is cheaper than the abstraction it replaces.

Not responsible for semantic validation (that stays in each domain's own `World::validate()`), physics, printing, rendering, or CLI argument handling.

## `malloy_ascii`

Responsible for turning 2D points into a character grid: viewport fitting (`fit_viewport`) and framed grid rendering (`render`).

Not responsible for physics, body types, simulation state, terminal control sequences, color, animation, or any graphics API.

## `malloy_collide`

Responsible for 2D collision primitives (`Circle`, `Aabb`), overlap tests, contact data (normal, penetration depth, contact point) with a documented deterministic answer wherever the contact normal is geometrically undefined, and the area properties of a shape: area, centroid, and polar second moment of area about the centroid.

Not responsible for density, mass, inertia, bodies, velocity, contact response, integration, broadphase acceleration, or scenario loading. The area properties stop at geometry: turning them into mass properties belongs to whichever domain owns bodies (ADR 0007). Like `malloy_ascii` it works on shapes, not on simulation types, so it never sees a `Body2D`.

## `malloy_particles`

Responsible for `Particle2D`, `ParticleSettings`, `ParticleWorld`, non-rotational contact response (positional correction plus an impulse along the contact normal), wall containment, uniform gravity, and its own validation and diagnostics including potential and total energy.

Not responsible for collision geometry (that is `malloy_collide`), orientation, angular velocity, torque, gravity, or scenario loading. It carries its own body type rather than widening `nbody::Body2D`, because each domain owns its concrete state.

## `malloy_rigid`

Responsible for `RigidBody2D` (pose plus mass distribution), mass-property construction from a shape and a density, the parallel-axis theorem, world/local conversions, pose integration, impulse application at a point, and rigid-body diagnostics.

Also responsible for rigid contact response since M14: disc against disc contacts, impulses that generate torque because they act away from the centre of mass, and immovable bodies represented as infinite mass and inertia. Infinity is per quantity: infinite mass alone is a body that can spin but not translate, infinite inertia alone one that can translate but not spin, and `is_static()` means both.

Since M15 it also owns a uniform gravity field, carried in `RigidSettings` and applied as an acceleration before the position update. Static bodies are skipped, and gravitational potential energy is reported alongside kinetic.

Not responsible for shape geometry (that is `malloy_collide`), persistent forces, force or torque accumulators, oriented-box contacts, orientation in 3D, or scenario loading. Inertia is a scalar and orientation is a scalar angle: quaternions and inertia tensors are 3D concerns deferred to M19.

## `malloy_springs`

Responsible for `Spring`, `SpringNetwork`, `SpringBody2D`, the pure `accumulate_spring_forces` kernel, and `SpringWorld`, which composes that kernel with the minimal translational integration needed to make the domain runnable.

Not responsible for rigid-body attachment points or torque, collision, gravity, or a generic engine-wide force-provider API. `accumulate_spring_forces` neither integrates nor mutates a body, so a later force producer can reuse the pipeline (ADR 0008).

## `malloy_nbody_terminal`

Responsible for hardcoded demo setup, reading a single positional scenario path from argv and choosing between a file and the built-in scenarios, calling library APIs, fixed number of steps, formatted terminal output, view framing policy (which viewport each frame shows), and returning nonzero on validation or step failure.

Not responsible for reusable physics logic.
