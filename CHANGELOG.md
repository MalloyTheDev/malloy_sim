# Changelog

All notable changes to MalloySim are recorded here. The format follows
`docs/10_CHANGELOG_TEMPLATE.md`; dates are ISO-8601.

## [Unreleased]

### Added

- `docs/decisions/0007-rigid-bodies-own-their-state.md`: settles the M11
  ownership boundary before M11 starts. A dedicated `RigidBody2D` rather than a
  widened `Body2D`, following the M10 precedent that each domain owns the state
  its model requires. Also records the three-way split between shape geometry,
  mass-property construction, and integration; the angle-canonicalization
  position; and that 2D means scalar angle and scalar inertia, with quaternions
  and inertia tensors deferred to M19.
- `docs/decisions/0006-multi-domain-dispatch.md`: records how the all-in-one
  goal is reached without a simulation base class. Each domain is a concrete
  library selected by a scenario `type` key and a plain dispatch switch. The
  `type` key is introduced only once a second domain exists.
- `CLAUDE.md` sections for multi-domain architecture and the template library,
  plus rule 16 defining when a domain counts as finished (validation, a tested
  invariant, malformed-input tests, and a scenario template).
- `README.md` section describing domains and templates.

### Changed

- Project identity amended from "simulation-first C++20 engine project" to
  "all-in-one science physics simulation". Terminal-first still holds, the
  anti-framework rules are unchanged, and work stays gated one milestone at a
  time.
- `docs/07_POST_M5_ROADMAP.md` restructured into tracks. The committed next
  stretch is classical mechanics depth (collision, rigid bodies, ballistics,
  springs and oscillators); later tracks are direction rather than commitment.
- `docs/00_START_HERE.md`, `docs/01_V3_ARCHITECTURE_DECISION.md`,
  `docs/03_MODULE_BOUNDARIES.md` and `README.md` updated to match.
- `README.md` milestone table extended through M8 and retitled, since M6-M8
  were complete but only M1-M5 were listed.
- `NBodyWorld` and `ParticleWorld` hold a `time::FixedStep` instead of a private
  `tick_count_` each. `malloy_time` had been built and tested since M3 while
  linked by nothing but its own test, and the two worlds had begun duplicating
  what it already does. `SimulationSettings` remains where `dt` is configured;
  `FixedStep` is the runtime counter built from it.
- Both worlds gained `elapsed_time()`, computed as ticks * dt rather than
  accumulated, so it cannot drift.

### Fixed

- `NBodyWorld::step()` now validates the state after the update as well as
  before it, and rolls back rather than reporting success when a step produced
  non-finite values (#3). Overflow from an extreme mass or dt reaches this point
  even when every input was finite and passed validation. The pre-step copy is
  held in a reused member so the common case does not allocate per step:
  measured on a 12-body, 200k-step run, a fresh vector each step cost about 24%
  and the reused buffer costs about 7%.
- The terminal app drives its step loop with an `int64` counter (#4). With
  `steps` at `INT_MAX` an `int` counter reached `INT_MAX`, passed the loop test,
  and overflowed on increment, which is undefined behavior and in practice never
  terminated.
- The scenario parser rejects a negative `steps` or `output_every` with a
  line-numbered error instead of accepting them (#4). Zero stays legal for both.
- The scenario parser rejects trailing tokens instead of silently discarding
  them (#6). The motivating case is a body line written with 3D fields, which
  used to be truncated to the first five and run as a different simulation than
  the file described. Comments and trailing whitespace, including CRLF, are
  still accepted.
- Scenario files are opened in binary mode, so a `0x1A` byte no longer silently
  truncates a file mid-parse under the MSVC CRT, and a stream that goes bad is
  reported as a parse failure rather than a success (#7).
- `render` skips a point whose normalized coordinates are not finite, and clamps
  the computed cell before indexing (#5). Every comparison against NaN is false,
  so a NaN escaped the clip test and reached `std::lround`, whose NaN result is
  unspecified: 0 on MSVC, `LONG_MIN` on GCC, which then cast to a huge index.
  `fit_viewport` no longer returns a non-finite rectangle, which was how a NaN
  arose from finite input.
- `CMakePresets.json` declares the CMake 4.2 it actually requires (#9). The
  `Visual Studio 18 2026` generator was added in 4.2, while the presets claimed
  3.21, so the documented build command failed on any CMake below 4.2. The
  project's own `CMakeLists.txt` still works at 3.21, and the docs now
  distinguish the two.
- `scenarios/three_body_triangle.scn` carries its velocities at full precision
  rather than rounded, so it reproduces the built-in three-body demo exactly
  (#12). Both templates now document their expected final report.
- The terminal demo prints the conserved quantities in scientific notation
  (#10). The two-body demo's energy is of order 1e-6, and under fixed(8) every
  step printed the identical `-0.00000050`, so the output would have looked
  perfectly conserved even if it were not. `E_total` now visibly oscillates
  within a bounded window, which is the symplectic signature; a monotone walk
  would indicate a broken integrator. Separation stays fixed-point. README
  sample output and both template expected-result comments were regenerated.
- Documentation corrected where it contradicted the code (#11). `docs/01` no
  longer gates graphics, the package manager, or the CMake layout on M5, which
  shipped. `malloy_time` is documented as currently linked by nothing but its
  own test, and removed from the dependency chain it was never part of;
  adopting `FixedStep` in `NBodyWorld` is left as a design change for its own
  milestone rather than papered over. ADRs 0002 and 0003 are restated against
  the milestone that owns them instead of expiring at M5. The `docs/00` read
  order includes `docs/09` and `docs/10`. The `docs/03` app boundary records
  that the app reads a scenario path from argv. The `type nbody` snippets in
  `CLAUDE.md` and ADR 0006 are labelled as planned and not parsed today.
- Removed every em dash (U+2014) from tracked content (#1). Nineteen were
  heading separators, where a hyphen reads correctly; one was prose in
  `sim_core.hpp`, where a colon does.
- Coincident bodies with `softening == 0` no longer produce NaN (#2). The
  inverse-cube law is undefined at zero separation: `1/sqrt(0)` is `+inf`, and
  `inf` times a delta that is exactly zero is NaN. `compute_accelerations` now
  skips such a pair, which is also the correct answer when `g == 0`. Both
  preconditions passed validation before, and `softening` defaults to 0 in the
  scenario format, so any scenario omitting that line was exposed.
- `total_potential_energy` carries the matching guard, so force and energy stay
  consistent as `docs/04` requires. Verified numerically: the central-difference
  gradient of the potential still matches the acceleration.
- `specific_orbital_energy` no longer evaluates 0/0 when `g == 0` and the bodies
  coincide. Its header now states that it is the unsoftened Kepler energy and is
  therefore not the quantity the softened integrator conserves.
- `-ffp-contract=off` is now set for GCC and Clang. The MSVC branch has always
  set `/fp:precise`, but the other branch set no floating-point flag at all, and
  both compilers default to `-ffp-contract=fast`, which fuses `a*b+c` into an
  FMA and silently changes results. `docs/04` forbids fast-math, so the policy
  was being enforced on one compiler only.
- `docs/05` and ADR 0006 claimed per-module test executables catch dependency
  leaks. They catch LINK-time leaks only: every module puts the whole `include/`
  tree on its consumers' include path, so using another module's header without
  linking it still compiles. Both documents now say so.
- The M10 changelog recorded the "resolve each pair twice" mutation as an
  equivalent mutant on the strength of reasoning. Measuring it disproved that,
  and the M10 entry has been corrected.

### Removed

### Tests

- Closed four gaps in `malloy_nbody_tests` where a plausible physics regression
  would have shipped undetected (#8):
  - an exact single-step assertion that distinguishes semi-implicit from
    explicit Euler, which the previous tolerance-based orbit tests could not;
  - exact acceleration magnitudes, pinning both the `G * m / r^3` law and that
    softening enters the denominator squared. The existing force-symmetry check
    survived a dropped `G`, a wrong inverse power, and a sign flip;
  - an angular-momentum case where `position.y * velocity.x` is nonzero. Every
    prior case left that term at zero, so half the formula was never exercised;
  - direct tests of `Body2D::is_valid` and `NBodySettings::is_valid` for
    infinite and NaN mass, G, and softening. Because `inf >= 0` is true, the
    `is_finite` clauses were dead.
- Regression tests for the coincident-body guard, including `g == 0`.
- Regression tests for every fix above: post-step validation and rollback,
  negative and zero run lengths, trailing tokens on each directive type, a
  stream that goes bad, an infinite viewport, and `fit_viewport` overflow.
- Every shipped scenario template is now parsed, validated, and stepped by
  `malloy_scenario_tests` (#12), and `parse_scenario_file` has coverage for the
  documented missing-file behavior. The templates were previously the
  least-tested files in the repository. Confirmed live by breaking a template
  and observing the test fail.
- Every new assertion was confirmed by mutation testing: eight broken
  implementations (explicit Euler, dropped G, inverse-square, unsquared
  softening, deleted angular-momentum term, both removed `is_finite` clauses,
  and a removed coincident guard) each make the suite fail.
- Closed the `malloy_math` and `malloy_time` coverage gaps found during the M8
  review and never filed: `distance` and `distance_squared` were only tested
  against the origin, where `a - b` and `a + b` are identical; `dot` was tested
  only as `dot(a, a)` and on a perpendicular pair; `approx_equal(Vec2)` never
  had x as the sole difference, so a y-only implementation passed everything;
  the `<=` boundary was untested in both overloads; `FixedStep::elapsed_time`'s
  no-drift guarantee was untested because `dt = 0.25` over four ticks is
  bit-identical either way; and negative infinity was untested in
  `FixedStep::create`.
- `elapsed_time()` on both worlds, including that an unusable `dt` leaves it at
  zero rather than reporting a bogus time.
- `docs/05` records the rule that produced both of the coverage gaps above:
  prefer asymmetric test configurations, because a symmetric setup cancels the
  very errors it is meant to catch.
- Verified by mutation testing: six broken implementations (a+b in distance,
  y-only vector equality, a sign-flipped dot term, a strict `<` boundary, and
  accumulated elapsed time in both `FixedStep` and the worlds) each fail.

### Architecture Notes

- The all-in-one goal raises rather than retires the engine-kernel risk named
  in `docs/01`. ADR 0006 is the standing answer: adding a domain must leave
  `malloy_sim_core` unchanged.
- The charter amendment changed no code. `malloy_sim_core` still contains
  exactly three types and `NBodyWorld` is still concrete.
- The #2 fix chooses to skip a degenerate pair rather than fail the step. That
  keeps the simulation running and matches the documented role of softening,
  but it does mean a modelling mistake stays silent. Reporting non-finite state
  loudly is tracked separately (#3).

## [M11] - 2026-09-09  (2D rigid bodies)

### Added

- `malloy_rigid` (STATIC, `malloy::rigid`): `RigidBody2D` (pose plus mass
  distribution), `MassProperties` and `mass_properties` from a shape and a
  density, the explicit parallel-axis `shift_inertia`, world/local conversions,
  `velocity_at`, a concrete `RigidWorld` with pose integration and
  `apply_impulse_at`, and total linear momentum, angular momentum and kinetic
  energy diagnostics.
- Area properties in `malloy_collide`: `area`, `centroid`, and
  `second_moment_of_area` about the centroid, for `Circle` and `Aabb`. Pure
  geometry, with no density and no mass anywhere.
- `malloy_rigid_tests`.
- `type rigid` in the scenario format, with a `rigid_body` key carrying mass,
  inertia, the local centre-of-mass offset, pose, and both velocities.
- `scenarios/spinning_bodies.scn`.

### Architecture Notes

- Implements ADR 0007. `RigidBody2D` is its own type rather than a widened
  `nbody::Body2D`: the two obey different equations of motion, and gravity has
  no use for an angle. `Body2D` and `Particle2D` are untouched.
- The three responsibilities ADR 0007 required to stay separable are separate:
  shape area properties live in `malloy_collide` with no density; mass-property
  construction and the parallel-axis step live in `malloy_rigid`; integration
  lives in `RigidWorld`. A centroid mistake therefore fails a geometry test
  with an exact closed-form expected value, never a rotational one.
- `inertia` is documented as being about the centre of mass, and `r` is measured
  from the centre of mass everywhere it appears. One documented reference point
  is what stops a parallel-axis term being applied twice or not at all.
- `step()` translates the CENTRE OF MASS and derives the body origin from it, so
  a body whose origin is offset orbits its centre rather than spinning about the
  wrong point.
- The angle is not wrapped by integration. Canonicalisation stays a separate
  concern with one owner (ADR 0007), and a test pins that the angle keeps
  accumulating past 2*pi.
- Impulse-only: no persistent forces, no force or torque accumulators, and no
  rigid-body contact response. Accumulators would introduce a determinism
  surface (accumulation order across callers is observable, as pair order was in
  M10) and nothing needs them yet. Forces belong with ballistics.
- Static or infinite-mass bodies are not represented. Nothing needs them until
  contact response does, and inventing a representation early would be
  speculative (rule 11).
- Third domain, third concrete runner, one more `switch` arm. `malloy_sim_core`
  is still exactly three types.

### Tests

- Every rotational test uses one deliberately asymmetric body: origin away from
  the centre of mass, nonzero angle, inertia not 1, mass not 1, both velocity
  components nonzero, and nonzero spin. Per `docs/05`, a body missing any of
  these hides a whole class of mistake.
- Mass properties against closed-form values with a density that is not 1, plus
  a cross-check against the textbook `I = mR^2/2` and `I = m(w^2+h^2)/12`
  relations, so a formula and its cross-check cannot both be wrong the same way.
- Parallel-axis with a nonzero distance, and a bad distance returning the input
  rather than a plausible wrong answer.
- Pose round-tripping, a quarter turn pinning counter-clockwise as positive, and
  a centre of mass that is the origin plus the ROTATED offset.
- `velocity_at` evaluated at the centre of mass, one unit away from it, and at
  the body origin, which is the case that distinguishes an arm measured from the
  centre of mass from one measured from the origin.
- An offset spinning body whose centre of mass must not move while its origin
  swings a quarter turn around it.
- Invariants: with no impulses, linear momentum, angular momentum and kinetic
  energy hold over a 5000-step run.
- Impulses: through the centre of mass changes linear velocity only; off centre
  changes angular velocity by exactly `(r x J)/I` with both components of `r`
  and `J` nonzero; mirroring the arm mirrors the spin, so both torque signs
  occur; and total momentum changes by exactly `J` and `p x J`.
- Verified by mutation testing: twelve broken implementations each make the
  suite fail. They include the torque arm measured from the body origin, a
  flipped cross-product sign, a clockwise `rotate`, a wrong-signed `omega x r`,
  `step` moving the origin instead of the centre of mass, angular momentum
  missing either its spin or its orbital term, parallel-axis using `m*d`, mass
  properties dropping density from the inertia, and both second-moment constants
  being wrong. All twelve were caught on the first attempt, which is the
  asymmetric-body rule from `docs/05` doing its job.

## [M10] - 2026-09-06  (colliding particles + multi-domain dispatch)

### Added

- `malloy_particles` (STATIC, `malloy::particles`): `Particle2D`,
  `ParticleSettings`, and a concrete `ParticleWorld` with non-rotational contact
  response (positional correction split by inverse mass, plus an impulse along
  the contact normal), wall containment, and `total_momentum` /
  `total_kinetic_energy` diagnostics.
- `malloy_particles_tests`.
- A `type` key in the scenario format selecting the domain, with `restitution`,
  `bounds`, and `particle` keys for the particles domain.
- `scenarios/bouncing_particles.scn`, the first non-gravity template.

### Changed

- `Scenario` carries a `ScenarioType` tag and both domains' field sets. `type`
  defaults to `NBody` when absent, so the format change is additive and both
  existing templates still parse and run byte-identically.
- `apps/nbody_terminal/main.cpp` dispatches on that tag to one concrete runner
  per domain, and reports the quantities each domain actually conserves.
- `print_view` takes a point list rather than a body list, so both domains share
  the ASCII view without either knowing the other's type.

### Architecture Notes

- This is the first milestone with two domains, and it is what ADR 0006 was
  waiting for: the `type` key was introduced only once a second world type
  actually existed. The whole dispatch mechanism is one switch and one concrete
  runner per domain. No base class, no virtual step, no registry, and
  `malloy_sim_core` is unchanged at three types.
- `Particle2D` is a separate type rather than a widened `nbody::Body2D`.
  Gravity has no use for a radius, and `Body2D` is documented as pure
  simulation state, so each domain owns its own concrete state instead. This
  also settles the M9 question of where a collision radius lives.
- M10 was re-scoped from rigid bodies to colliding particles, so collision
  became demonstrable one milestone sooner. Rigid bodies (orientation, angular
  velocity, torque) moved to M11.
- Walls are immovable, so they carry momentum away. Momentum is a conserved
  quantity only for a run where nothing touches a wall, and the tests and the
  template both say so rather than asserting a conservation that does not hold.

### Tests

- Invariants: particle/particle collisions conserve momentum at every
  restitution, restitution 1 also conserves kinetic energy exactly, and anything
  below it strictly removes energy.
- Exact hand-computed results: equal-mass head-on elastic collision exchanges
  velocities; restitution 0 leaves no separation at all; free flight moves by
  exactly `velocity * dt`.
- All four walls, each a separate branch, with position and both velocity
  components asserted.
- Containment over a 4000-step run in a crowded box, and no energy created from
  nothing over a 3000-step elastic run.
- Determinism to exact equality, and validation of settings, particles, and a
  particle too large for its box.
- Multi-domain parsing: absent `type` still means nbody, a key from the wrong
  domain is a line-numbered error, an unknown type is rejected by name, `type`
  after a domain key is rejected, and every new key's field count is checked.
- Verified by mutation testing. Nine mutations were applied; two escaped and
  exposed real coverage gaps, both since closed:
  - re-impulsing an already separating pair conserves both momentum and energy,
    so no invariant could see it; it needed a direct velocity assertion.
  - splitting the positional correction by the wrong mass is invisible with
    equal masses and never touches velocity; it needed an unequal-mass pair.
  A third mutation (resolving each pair twice) is NOT caught. It was originally
  recorded here as an equivalent mutant on the strength of reasoning; measuring
  it afterwards disproved that. Diffing a full 6000-step run shows the results
  diverge from about step 4500, in the eighth significant digit and growing,
  because the doubled pass applies an extra rounding-level positional
  correction. It is a change in update order, not a physics error, and no
  portable test can catch it: the only assertion that would is a hardcoded
  golden value, which asserts the cross-toolchain bitwise determinism docs/04
  explicitly declines to claim.

## [M9] - 2026-09-06  (collision primitives)

### Added

- `malloy_collide` (STATIC, `malloy::collide`): `Circle` and `Aabb` shapes with
  their own validation, cheap `overlaps` tests, and `contact` queries returning
  `std::optional<Contact>` (normal, penetration depth, contact point) for all
  three shape pairs.
- `malloy_collide_tests`.

### Architecture Notes

- `malloy_collide` depends only on `malloy::math` and knows nothing about
  `Body2D` or any simulation type, exactly as `malloy_ascii` knows points rather
  than bodies. Pairing a body with a shape is the caller's job, so
  `malloy_nbody` stays free of collision concerns and `Body2D` keeps its
  documented scope of position, velocity, and mass.
- Shipped as a support library, not a domain: it has no world and no scenario
  template, so rule 16 does not apply to it. M10 (rigid bodies) is what will
  make collision demonstrable.
- Contact conventions are fixed and documented: the normal points from `a` to
  `b` and is unit length, penetration is never negative, and moving `b` by
  `+normal * penetration` separates the pair exactly. Shapes that merely touch
  produce a contact with zero penetration rather than none, because the boundary
  is exactly representable and therefore testable.
- Every configuration where the normal is geometrically undefined has a fixed,
  documented answer rather than a NaN: coincident circle centers and an exactly
  centered circle inside a box both give +x, and an equal box overlap on both
  axes resolves along x. Fixed choices keep degenerate cases repeatable
  (docs/04 determinism).
- No broadphase. Pair testing is the caller's loop, matching the existing O(n^2)
  gravity loop. A spatial hash is the right answer when N actually hurts, and
  building it now would be speculative infrastructure.
- No new `StepStatus` values were needed, so `malloy_sim_core` is unchanged and
  still contains exactly three types (ADR 0004, ADR 0006).

### Tests

- Separated, exactly touching, and overlapping cases for all three shape pairs,
  with hand-computed penetration depths and unit normals including diagonal
  ones.
- Degenerate and hostile input: coincident centers, zero-radius circles,
  zero-area boxes, inverted boxes, NaN and infinite coordinates and radii, a
  circle centered exactly inside a box, and coordinates large enough that both
  the squared distance and the squared radius sum overflow.
- A round trip asserting the convention actually works: moving `b` by
  `+normal * penetration` leaves the pair exactly touching.
- A 61x61 sweep asserting `overlaps` and `contact().has_value()` agree, since
  they are separate code paths that could silently diverge.
- Verified by mutation testing: nine broken implementations (flipped normal,
  dropped touching case, removed coincident fallback, deepest-axis separation,
  forgotten radius in the inside case, changed tie-break, both dropped
  validation clauses, and removed overflow guards) each make the suite fail.
  The overflow mutation initially survived, which showed the first version of
  that test never reached the guard; the test was corrected rather than the
  claim being softened.

## [M8] - 2026-09-06  (simple 2D debug visualization)

### Added

- `malloy_ascii` (STATIC, `malloy::ascii`): `Viewport`, `fit_viewport` (bounding
  box of the finite points plus a margin, widened when an extent is zero), and
  `render` (points into a framed character grid, y up, with out-of-view and
  non-finite points clipped).
- `malloy_ascii_tests`.

### Changed

- `apps/nbody_terminal/main.cpp` prints an ASCII view of the body positions
  after each diagnostics line. One viewport is held for the whole run, seeded
  from the initial state and only ever grown, so successive frames share a
  scale and a near-stationary body keeps its cell.

### Tests

- Grid placement (center, y-up orientation, custom characters), clipping of
  out-of-view points, empty input, and `fit_viewport` margin and
  degenerate-extent behavior.
- Boundary and malformed input: non-finite coordinates, a zero-extent viewport,
  and non-positive grid sizes, each of which must degrade to an empty or
  clamped grid rather than index out of range.

### Architecture Notes

- `malloy_ascii` depends only on `malloy::math`: it knows about points, not
  bodies, so the physics libraries stay free of presentation code.
- Framing policy (which viewport a frame shows) lives in the app; the library
  only fits and draws.
- Still no graphics API, no terminal control sequences, and no CLI parser.

## [M7] - 2026-06-19  (scenario/config loading)

### Added

- `malloy_scenario` (STATIC, `malloy::scenario`): a minimal line-based text
  parser turning a scenario file into `{SimulationSettings, NBodySettings,
  bodies, steps, output_every}`. Parsing returns a result (line number +
  message) and never throws; semantic validation stays in `NBodyWorld`.
- Example scenarios in `scenarios/` (`two_body.scn`, `three_body_triangle.scn`).
- `malloy_scenario_tests`.

### Changed

- `apps/nbody_terminal/main.cpp` loads and runs a scenario file given as its
  single argument; with no argument it runs the built-in scenarios. The report
  handles any body count (N<2 omits the separation column).

### Architecture Notes

- No CLI parser and no serialization library: a single positional path and a
  hand-rolled plain-text format. Config parsing lives in its own library, not
  in the dumb app or the physics library.

## [M6] - 2026-06-19  (N-body diagnostics)

### Added

- Whole-system N-body diagnostics in `malloy_nbody`
  (`include/malloy/nbody/diagnostics.hpp`): `total_kinetic_energy`,
  `total_potential_energy` (softened), `total_energy`, `total_momentum`,
  `center_of_mass`, `total_angular_momentum`.
- A second hardcoded demo scenario: a rotating equilateral three-body (Lagrange)
  configuration, exercising N>2.

### Changed

- `apps/nbody_terminal/main.cpp` refactored to run multiple scenarios via a
  shared helper and report system diagnostics (separation, total energy, total
  angular momentum).

### Tests

- Known-value diagnostics test and an N=3 conservation test (linear momentum to
  machine precision; energy and angular momentum bounded over a run).

### Architecture Notes

- The app stays dumb: still no CLI/config; all diagnostics live in the library.

## [M5] - 2026-06-19  (terminal N-body demo)

### Added

- `apps/nbody_terminal/main.cpp` now runs the hardcoded normalized sun + planet
  orbit (G=1, dt=0.001, softening=1e-6, 10000 steps), printing radius and
  specific orbital energy every 1000 steps and returning nonzero on
  validation/step failure.
- `malloy::nbody::specific_orbital_energy` (`include/malloy/nbody/diagnostics.hpp`).

### Changed

- `malloy_nbody_terminal` now links `malloy::nbody`.

### Tests

- Specific-orbital-energy known-value and approximate-conservation tests.

### Architecture Notes

- The app stays dumb: all physics (including the energy formula) lives in the
  library, not `main.cpp`.

## [M4] - 2026-06-19  (malloy_nbody)

### Added

- `malloy_nbody` (STATIC, `malloy::nbody`): `Body2D`, `NBodySettings`,
  `NBodyWorld` with softened pairwise gravity and semi-implicit Euler in the
  documented update order.

### Tests

- `malloy_nbody_tests`: all 11 required M4 cases (settings/state validation,
  single-body zero acceleration, softening no-NaN, force symmetry, determinism,
  near-circular two-body orbit).

### Architecture Notes

- `NBodyWorld` is concrete; validation returns `StepStatus` and never throws.

## [M3] - 2026-06-19  (malloy_time + malloy_sim_core)

### Added

- `malloy_time` (INTERFACE, `malloy::time`): `FixedStep`.
- `malloy_sim_core` (STATIC, `malloy::sim_core`): `SimulationSettings`,
  `StepStatus`, `StepResult`.

### Tests

- `malloy_time_tests`, `malloy_sim_core_tests`.

### Architecture Notes

- `malloy_sim_core` kept tiny: no engine kernel, no polymorphism (ADR 0004).

## [M2] - 2026-06-19  (malloy_math)

### Added

- `malloy_math` (INTERFACE, `malloy::math`): `Real = double`, `Vec2`, and
  vector/scalar helpers (`dot`, lengths, distances, `normalize`, `approx_equal`,
  `is_finite`).
- `MALLOY_CHECK_VEC2_NEAR` test helper.

### Tests

- `malloy_math_tests`.

## [M1] - 2026-06-19  (CMake skeleton and smoke test)

### Added

- Initial CMake skeleton and `CMakePresets.json` (Visual Studio 2026 MSVC
  debug/release).
- `malloy_nbody_terminal` app stub and `malloy_smoke_tests` CTest target.
- MIT `LICENSE`; project planning docs.

### Tests

- Smoke test passes through CTest.

### Architecture Notes

- Locked C++20 + CMake + VS Code + MSVC baseline. No physics code yet.
