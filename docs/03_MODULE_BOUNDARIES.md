# 03 - Module Boundaries

> All modules below are implemented (M2-M31). These boundaries are in force
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

Responsible for 2D collision primitives (`Circle`, `Aabb`), overlap tests, contact data (normal, penetration depth, contact point) with a documented deterministic answer wherever the contact normal is geometrically undefined, and the area properties of a shape: area, centroid, and polar second moment of area about the centroid. Since M22 it also carries the 3D primitives `Sphere` and `Plane3`, and since M30 the oriented `Box3`, whose contact with a plane is the manifold of its penetrating corners (a list, not a single point).

Not responsible for density, mass, inertia, bodies, velocity, contact response, integration, broadphase acceleration, or scenario loading. The area properties stop at geometry: turning them into mass properties belongs to whichever domain owns bodies (ADR 0007). Like `malloy_ascii` it works on shapes, not on simulation types, so it never sees a `Body2D`.

## `malloy_particles`

Responsible for `Particle2D`, `ParticleSettings`, `ParticleWorld`, non-rotational contact response (positional correction plus an impulse along the contact normal), wall containment, uniform gravity, and its own validation and diagnostics including potential and total energy.

Not responsible for collision geometry (that is `malloy_collide`), orientation, angular velocity, torque, gravity, or scenario loading. It carries its own body type rather than widening `nbody::Body2D`, because each domain owns its concrete state.

## `malloy_rigid`

Responsible for `RigidBody2D` (pose plus mass distribution), mass-property construction from a shape and a density, the parallel-axis theorem, world/local conversions, pose integration, impulse application at a point, and rigid-body diagnostics.

Also responsible for rigid contact response since M14: disc against disc contacts, impulses that generate torque because they act away from the centre of mass, and immovable bodies represented as infinite mass and inertia. Infinity is per quantity: infinite mass alone is a body that can spin but not translate, infinite inertia alone one that can translate but not spin, and `is_static()` means both.

Since M16 a world also owns immovable ground planes, carried in `RigidSettings`. They are `collide::Halfplane` values, so the geometry stays in `malloy_collide` and `malloy_rigid` only resolves against it.

Since M17 contacts also carry Coulomb friction: a tangential impulse clamped to the friction coefficient times the normal impulse. It is a contact impulse rather than a persistent force, so it needs no force accumulator and rule 5 is untouched. `malloy_particles` has no friction.

Since M15 it also owns a uniform gravity field, carried in `RigidSettings` and applied as an acceleration before the position update. Static bodies are skipped, and gravitational potential energy is reported alongside kinetic.

Not responsible for shape geometry (that is `malloy_collide`, which since M22 also has `Sphere` and `Plane3`), persistent forces, force or torque accumulators, oriented-box contacts, or scenario loading. In `RigidBody2D` the inertia is a scalar and the orientation is a scalar angle; those are the 2D CASES of an inertia tensor and a quaternion, not alternatives to them (`docs/decisions/0009-three-dimensions-are-the-destination.md`).

Since M20 the module also owns `RigidBody3D` and `Rigid3DWorld`: rotation in three dimensions, with a quaternion orientation and three principal moments of inertia. It lives here rather than in a library of its own because rigid-body dynamics is one domain and the dimension is not a domain, the same reasoning M19 used for 3D gravity. M20 was torque-free; M21 added a constant world-frame torque as a setting (the M15 pattern, applied as Euler forcing rather than through a force/torque accumulator), which is what makes a gyroscope precess; M22 added gravity (an acceleration) and restitution contacts against immovable ground planes, resolved with a normal impulse; M23 added Coulomb friction, a tangential impulse clamped to the normal one, which DOES have a lever arm and so spins a sliding sphere up until it rolls without slipping (at 5/7 of its sliding speed). M24 added sphere-against-sphere collisions: two movable bodies exchanging momentum through the SAME impulse core the ground uses, with the plane expressed as a zero-inverse-mass participant so the formula is written once (ADR 0008). M25 added 3D mass properties: `mass_properties_3d` computes a compound body's mass, centre of mass and principal moments from solid spheres and density, via the parallel-axis theorem and the diagonalization of the resulting inertia tensor (`math::Mat3` and a symmetric eigensolver, new in `malloy_math`). It is a construction layer, the 3D sibling of the 2D `mass_properties`, not a domain: no world, no template, the same standing M9's collision primitives had. The inertia is still a diagonal rather than a general 3x3 tensor. M26 added `SolidBox` as a second mass primitive (three distinct moments and a real orientation, so a single tilted box already has a non-diagonal lab-frame tensor) and `combine`, which assembles a compound body from mass-property sets of any kind by reconstructing each tensor (`math::to_mat3`, the inverse of M25's `to_quat`), shifting it by parallel-axis and diagonalizing the sum. M27 added `SolidMesh`, the general primitive: the mass properties of a solid bounded by a closed triangle mesh, computed as signed-tetrahedron volume integrals (the divergence theorem, `math::trace` turning the covariance into an inertia tensor), so an axis-aligned box mesh reproduces the M26 box exactly and any closed shape now has mass properties. M28 added a constant applied FORCE to `Rigid3DSettings`, the translational half of a wrench whose rotational half is M21's torque: it enters as the acceleration F/m before the position update (so it scales with mass, unlike gravity) and acts through the centre of mass (so it makes no torque), and like the torque it is a single constant setting applied as forcing, not a force accumulator (rule 5). M30 gave `RigidBody3D` a box collider (`half_extents`, zero meaning the sphere `radius` is used instead) and resolves it against a ground plane as the manifold of its penetrating corners reduced to one centroid contact, so an oriented box rests and tumbles without a spurious torque on a symmetric landing; box against box is deferred. `malloy_rigid` gained a link on `malloy_collide` for the sphere, plane and box primitives, the same as the 2D path.

## `malloy_springs`

Responsible for `Spring`, `SpringNetwork`, `SpringBody2D`, the pure `accumulate_spring_forces` kernel, and `SpringWorld`, which composes that kernel with the minimal translational integration needed to make the domain runnable. Since M31 it also has the 3D siblings `SpringBody3D` and `SpringWorld3D` (and a `Vec3` overload of the kernel); `Spring` and `SpringNetwork` are dimension-agnostic (ids and scalars) and are shared unchanged, since the dimension is not a domain (ADR 0009).

Not responsible for rigid-body attachment points or torque, collision, gravity, or a generic engine-wide force-provider API. `accumulate_spring_forces` neither integrates nor mutates a body, so a later force producer can reuse the pipeline (ADR 0008).

## `malloy_nbody_terminal`

Responsible for hardcoded demo setup, reading a single positional scenario path from argv and choosing between a file and the built-in scenarios, calling library APIs, fixed number of steps, formatted terminal output, view framing policy (which viewport each frame shows), and returning nonzero on validation or step failure.

Not responsible for reusable physics logic.
