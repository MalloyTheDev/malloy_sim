# ADR 0011 - Contacts Are Solved Per Point, Iteratively, Inside Each Rigid World

## Status

Accepted, recorded at M40.

## Context

Until M40 the rigid worlds resolved contacts in ONE pass, after the position
update. Each box contact manifold (box on the ground, box against box) was
reduced to a single contact at a centroid, and a positional correction moved
the bodies out by the deepest penetration. Three defects followed from that
design rather than from any one bug:

- **#41:** a box that landed anything but flat rocked and gained energy. A
  penetration-weighted centroid (the review pass after M39) helped, but one
  contact cannot stop several corners at once, so a 3D box kept rocking at
  the scale of one step's gravity.
- **#83:** at restitution 1, a nudged box on the floor gained energy without
  bound, 442 times its start in 3D.
- **#44:** gravity moved a body before friction acted, so a body inside the
  friction cone crept down a slope at g sin(theta) dt per unit time.

Stacks could not stand: a staggered stack of three boxes fell within about 13
seconds, and an aligned stack reported motion that was not happening, because
the correction pushed it back up every step. Rendering (#80) would put all of
this on screen.

## Decision

Each rigid world (`RigidWorld`, `Rigid3DWorld`) solves its contacts PER POINT,
iteratively, with its own solver in its own source file.

**The step order changes:**

- **A. Detect** contacts at the START-of-step poses, every pair in a fixed
  order (body pairs ascending, then body against plane). The queries are
  speculative: a pair still apart by up to 1e-4 of the smaller body's size is
  kept at a negative penetration, so a resting contact cannot drop out on
  rounding.
- **B. Impact phase,** only when the restitution is above 0. Points
  approaching faster than a threshold of 4 dt (|g| + |F| / m_min) at the
  middle of the step's forcing get Newton's restitution. Points that were not
  approaching take no part, so the phase never raises the kinetic energy.
- **Forcing** is applied exactly as before, so free flight is bit-identical.
- **C. Rest phase,** every point: sequential impulses, accumulated and clamped
  (the normal impulse never below 0, the friction impulse within the disc of
  radius mu times that point's accumulated normal impulse).
  - In 3D the friction disc is projected in the metric of the point's tangent
    effective mass, which is Coulomb's maximum-dissipation law.
  - The phase is warm-started from the last step's impulses.
- **Integrate positions** as before, from the solved velocities.
- **D. Project positions** (nonlinear Gauss-Seidel, 3 passes, a fifth of the
  penetration beyond the slop per pass). The projection keeps the body's world
  angular momentum through a turn, so it changes no momentum and adds no kinetic
  energy.

**Iteration stops at a tolerance, with a fixed cap.** The rest phase stops at
1e-8 of the velocity scale, the impact phase at 1e-12, and each phase at 100
sweeps at most. The count is set by the state, so a run repeats exactly.

**The warm-start impulses are world state.** They are snapshotted with the
bodies and restored by a rolled-back step. A world rebuilt from `bodies()`
alone starts without them, so it differs for its first steps, deterministically.

**The code shape keeps the project's rules:**

- The contact points are plain structs, private to each world's source file,
  in a list rebuilt every step. There is no Constraint type, no solver
  interface and no shared 2D/3D template (rules 11 and 12), and no force
  accumulator (rule 5).
- Each world keeps a single impulse core, as ADR 0008 asks.
- `malloy_sim_core` does not change (rule 13, ADR 0004).
- Every parameter is a named constant with its derivation beside it, stated in
  the world header's step contract. No scenario key and no settings field was
  added.

**The collide layer gains margin-aware queries.** Each contact query has an
overload that keeps points within a gap, with the plain query unchanged. The 3D
box-box separating-axis test prefers a face axis over an edge axis within 5 per
cent.

- That preference was tried and removed once before (issue #22), because it put
  contact points outside the boxes.
- It is adopted now only with a randomized proof in the tests: 100,000
  configurations at scales from 1e-6 to 1e6, every manifold point inside both
  boxes to rounding.
- Without it, a twisted five-box stack drifted 4.5 mm in 10 s; with it, 0.09 mm.

## Rejected alternatives

- **Keeping the old step order.** A per-point solve with the old order still
  crept down slopes at exactly g sin(theta) dt. #44 is the order, not the
  solver.
- **A fixed iteration count.** Ten or twenty sweeps let a five-box stack drift
  millimetres and broke the symmetry of a frictional landing. A tolerance with
  a cap is cheaper (no sweeps at rest) and no less repeatable.
- **No warm start.** It took 77 sweeps a step instead of 2, and a box on one
  100 times lighter sank through it.
- **Baumgarte stabilisation.** It fed the correction into the velocity: a box
  started 2 cm inside the floor was launched with ten times the energy the
  overlap held.
- **Split impulse.** It broke mirror symmetry and gained energy on elastic
  bounces. With a rotational part it made friction stacks creep sideways.
- **Scaling the friction impulse radially onto the cone.** That is Coulomb's law
  only for an isotropic effective mass; it turned a box sliding diagonally.
- **New scenario keys for the solver parameters.** Each key would have cost
  validation, parse and fuzz surface, and the derived constants already hold
  across the time steps tested (dt 0.001 to 1/60).

## Consequences

- **The defects are fixed.**
  - A box that lands tilted comes to rest.
  - A nudged box at restitution 1 no longer gains energy.
  - Static friction holds: creep falls from centimetres to about 1e-11 m in 10
    s.
  - Short stacks (three to five boxes, aligned, offset and yawed) stand, and so
    does a box on one 100 times lighter.
  - A resting body no longer micro-bounces.
- **Templates with a contact moved.** Where a closed form exists the checks now
  pin it exactly: a body on a ramp has momentum m g sin cos t, and a rolling
  sphere has 34/3. Contact-free templates are bit-identical.
- **The cost grows with the points and sweeps.**
  - **Work.** The parse-time work budget charges every rigid pair and plane
    test at its measured worst case, with every point solved at the 100-sweep
    cap (a 3D box pair 45860 units, about 10 ns each).
  - **Memory.** The parser caps the potential contact points at 2000000, which
    bounds the solver's memory at about 2 GiB.
  - **Large scenes.** The worst-case charge limits large rigid scenes: 1000
    spheres get about 74 steps. A broadphase and a runtime contact-work budget
    are proposed in #85.
- **Limits.**
  - Stacks taller than about five boxes may need more than the 100-sweep cap.
  - A box much heavier than the one under it (1e4 times) still sinks, slowly
    and without gaining energy.
  - Box against disc in 2D is still not implemented (#48).
