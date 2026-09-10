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

- The terminal app drives every world through one `drive()` template instead of
  four copies of the same loop, and extracts positions through one
  `positions_of()` template instead of three. `rigid_points` stays separate
  because it emits two points per body and is genuinely a different function.
  414 lines down to 380, and output is byte-identical across the built-in demo
  and all six templates.
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

- Three stale milestone claims corrected: `README.md` and
  `docs/07_POST_M5_ROADMAP.md` both still said M1-M8 were the completed set, and
  `docs/08_AI_HANDOFF_PROMPT.md` said M1-M11 and omitted ballistics. A handoff
  prompt that understates what exists is worse than no prompt.
- `docs/07_POST_M5_ROADMAP.md` said M1-M12 were the completed set, two
  milestones stale. This is the second time that line has gone stale: it was
  corrected in the pre-M13 audit and drifted again through M13 and M14. It is a
  hand-maintained count with nothing guarding it.
- The `restitution` scenario key wrote to both `particle_settings.restitution`
  and `rigid_restitution` regardless of the declared type, so a rigid scenario
  silently carried a particle setting it never reads. It now writes only the
  field its own domain uses. Harmless before, but it read as a mistake to
  anyone inspecting a parsed scenario.
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

- The `simulation_settings()`, `nbody_settings()` and `particle_settings()`
  accessors were public API with no test at all. Now covered on all three
  worlds, with every field a distinct value so a getter returning the wrong
  member shows, and confirmed by a mutation that makes one return a default.
- The scenario tests enumerate `scenarios/` with `std::filesystem` instead of a
  hardcoded list. The list happened to be in sync, but the structure was the
  risk: a template added without editing the test file would have been silently
  untested, quietly breaking rule 16. Verified by dropping a broken template
  into the directory and watching it fail with no test edit.
- Rule 16's other half is now enforced too. Every template must contain a
  documented expected result, checked by reading the file rather than trusted to
  review. Verified by adding a template with no `# Expected` block and watching
  it fail.
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

## [M15] - 2026-09-10  (gravity for rigid bodies)

### Added

- A uniform gravity field in `malloy_rigid`, carried in the new `RigidSettings`
  and applied to every non-static body before the position update, which keeps
  the integration semi-implicit (symplectic) Euler.
- `RigidSettings`, holding restitution and gravity with an `is_valid()` that
  rejects a restitution outside [0, 1] and a non-finite gravity.
- `total_potential_energy(bodies, gravity)` and `total_energy(bodies, gravity)`
  in `malloy::rigid`. The diagnostics line now reports kinetic plus potential,
  so it stays near flat during free fall instead of climbing.
- `gravity <gx> <gy>` is accepted for `type rigid`, alongside the `type
  particles` use it was written for in M12.
- `scenarios/dropped_bodies.scn`.
- `include/malloy/rigid/rigid_settings.hpp`, so a translation unit that needs
  only the settings does not have to include the whole world class. This follows
  the `malloy/particles/particle_settings.hpp` precedent, and `scenario.hpp`
  now includes the smaller header.
- A README template-count guard in the scenario tests. The stated count has gone
  stale twice, so it is now read out of README.md and compared against the
  directory, and a wrong number fails the suite.

### Changed

- BREAKING, C++ API: `RigidWorld` takes a `RigidSettings` where it took a bare
  `math::Real restitution`. M14 added restitution positionally and M15 would
  have added gravity beside it, which is exactly where a transposed pair of
  reals stops being a compile error and starts being a wrong simulation. The
  default is unchanged, so `RigidWorld{sim, bodies}` still means restitution 1
  and no gravity. The `restitution()` accessor is replaced by `settings()`.
- The rigid diagnostics column is labelled `E` rather than `KE`, because it now
  reports kinetic plus potential energy. With no gravity the two are identical,
  so the figures documented in `scenarios/spinning_bodies.scn` and
  `scenarios/tumbling_impact.scn` are unchanged; the wording in the latter was
  updated to match the new label.
- `MALLOY_CHECK_NEAR` and `MALLOY_CHECK_VEC2_NEAR` now FAIL on a NaN. See below.

### Fixed

- `MALLOY_CHECK_NEAR` and `MALLOY_CHECK_VEC2_NEAR` silently passed any NaN.
  Both compared `std::abs(difference) > epsilon`, and every comparison involving
  a NaN is false, so a NaN difference took the passing branch. That affected
  every one of the roughly 370 near-equality assertions in the project at once:
  a test asserting that some quantity was 0 would pass just as happily if the
  code under test returned NaN.

  Both macros now use the negated form `!(difference <= epsilon)`, which is true
  for a NaN and therefore fails. Rebuilding and running all eleven suites with
  the stricter macro produced no new failures, so nothing in the project was
  passing only because of the hole; it was a latent hole rather than a live one.

  Two consequences worth knowing. A NaN tolerance is now rejected as well, so a
  bad epsilon cannot quietly accept every value. And two infinities are no
  longer "near", because infinity minus infinity is NaN: use `MALLOY_CHECK_EQ`
  for a value expected to be infinite.

  Found while mutation testing M15, when a mutation that made
  `total_potential_energy` return NaN for a static body was not caught by an
  assertion written specifically to catch it.

### Architecture Notes

- Gravity is an acceleration, not a force. It is added to a velocity directly
  and never divided by a mass, which is why bodies of different mass fall
  identically. A static body is skipped explicitly: infinity times zero is NaN,
  and without the skip an immovable wall would accelerate downwards forever
  while still refusing to be pushed by anything.
- Rotation under gravity is EMERGENT, not implemented. The field acts through
  the centre of mass and so generates no torque however long a body falls, but a
  contact away from the centre of mass does have a moment arm. Dropping a body
  whose centre of mass is offset from its collision disc makes it rock, and that
  is M14's contact arm driven by M15's field with no new machinery between them.
- `malloy_rigid` still has no force or torque accumulators. Gravity is a setting
  applied in one pass, not a registered force producer.
- Statics are excluded from potential energy, because infinite mass times a
  position is not a number.

### What conserves, and what does not

Total energy is NOT conserved under a constant field, and the amount lost is
exact rather than approximate. Semi-implicit Euler sheds
`(1/2) * (sum of m) * |g|^2 * dt^2` per free-flight step, the same constant
derived for the particle domain in M12 because it is the same integrator on the
same kind of field. `scenarios/dropped_bodies.scn` measures 1.563840e-01 over
its first 500 steps of free fall against a predicted 1.563837e-01.

Linear and angular momentum are not conserved either, and for a more ordinary
reason: gravity and an immovable floor are both external to the bodies.

### Tests

- Zero gravity reproduces the pre-M15 trajectory exactly, position and angle, so
  every scenario written before the setting existed is unaffected.
- The integrator is checked to be semi-implicit rather than explicit: with
  g = (0, -10) and dt = 0.5 one step gives velocity -5 and a centre-of-mass
  displacement of -2.5, not the 0 that explicit Euler would give.
- Gravity accelerates a 0.5 kg body and a 500 kg body identically, and moves an
  immovable one not at all.
- The free-fall energy loss matches the derived constant to 1e-10.
- Potential energy is measured at the centre of mass, not the body origin, and
  an offset centre of mass changes it accordingly.
- A body dropped on a floor lands, stops, and does not sink into it.
- The emergent case: a body with an offset centre of mass resting on a floor
  rocks, which it would not do if the centre of mass were at the disc centre.
- The claim `scenarios/dropped_bodies.scn` is built on is asserted directly:
  two bodies differing tenfold in mass hold the same height, angle and vertical
  velocity to 1e-12 across a fall and a bounce. A companion test breaks their
  shared inertia-to-mass ratio and confirms the agreement then stops, which is
  what shows the bounce depends on `I/m` rather than on mass.
- The check macros now have their own tests, in the smoke executable because it
  links no libraries: ordinary near-equality, NaN on either side, NaN in either
  Vec2 component, a NaN tolerance, infinities, that a failure is reported and
  that a pass is silent.
- Verified by mutation testing. Five mutations to the gravity work were applied
  and four were caught at once: gravity applied after the position update
  instead of before, gravity scaled by mass, the static skip removed, and the
  sign of the potential energy flipped. The fifth, dropping the static skip from
  `total_potential_energy`, escaped, and chasing it found the NaN hole in the
  check macros recorded above. With the macros fixed it is caught. The fixed
  macros were themselves mutated back to the old comparison, which the new
  smoke tests catch, and the README count guard was mutated to a wrong count,
  which the scenario tests catch.

## [M14] - 2026-09-10  (rigid-body contact response)

### Added

- Rigid contact response in `malloy_rigid`. Bodies with a nonzero radius collide
  as discs, and the impulse acts at the contact point rather than at the centre
  of mass, so it generates torque and the bodies tumble.
- `RigidBody2D::radius`, a collision disc centred on the BODY ORIGIN.
- Static bodies: infinite mass and inertia mean immovable, with
  `is_static()`, `inverse_mass()` and `inverse_inertia()`.
- `RigidWorld` restitution, validated to lie in [0, 1].
- `rigid_static <radius> <px> <py> <angle>` in the scenario format, and a
  `restitution` key now accepted for `type rigid` as well as `type particles`.
- `scenarios/tumbling_impact.scn`.

### Changed

- BREAKING, scenario format: `rigid_body` takes eleven fields rather than ten.
  The radius is third, after mass and inertia. `scenarios/spinning_bodies.scn`
  was updated with radius 0, meaning no collision, and its output is unchanged.
  Any scenario file outside this repository using `rigid_body` needs the extra
  field.
- `RigidBody2D::is_valid()` now ACCEPTS infinite mass and inertia. M11 rejected
  them on the grounds that nothing needed statics until contact response did
  (ADR 0007), and a test asserted that rejection. Contact response needs them,
  so the contract changed deliberately and the test changed with it. NaN is
  still rejected, because NaN > 0 is false.
- The diagnostics skip immovable bodies. Infinite mass times zero velocity is
  NaN, and an immovable body carries no momentum or kinetic energy to report.

### Architecture Notes

- The disc is centred on the body origin, not on the centre of mass, and that
  offset is the entire mechanism. A contact normal on a disc centred at the
  centre of mass passes straight through it, giving a zero moment arm, so such a
  body can never be spun by an impact. Offsetting them is what makes contacts
  rotational.
- Immovability is infinity rather than a flag, so `inverse_mass()` and
  `inverse_inertia()` come out as exactly zero and every impulse formula handles
  a static body with no branch: it contributes nothing to the effective mass and
  receives nothing from the impulse.
- Contacts are disc against disc. Oriented boxes need SAT, which
  `malloy_collide` does not have, and that is its own milestone.
- Still no persistent forces and no force or torque accumulators in
  `malloy_rigid`.

### What conserves, and what does not

The contact IMPULSE conserves both linear and angular momentum exactly. It is
equal and opposite and acts at one shared point, so the two angular
contributions are `-(p x J)` and `+(p x J)`.

The positional CORRECTION does not conserve angular momentum. It moves positions
without changing velocities, which shifts the orbital term `m (r x v)` by
`c x (v_b - v_a)` where `c` is the correction. This was found by a failing test
and then isolated by disabling the correction, which made angular momentum hold
to 1e-8. The perturbation is proportional to penetration depth: the long-run
test drifts by about 5e-4, while a single contact with penetration nine orders
smaller drifts by about 5e-9. Linear momentum is unaffected either way, because
the correction never touches a velocity.

This is the usual cost of resolving penetration by moving bodies. It is recorded
in the header and in the tests rather than hidden behind a loose tolerance.

### Tests

- Statics: infinite mass and inertia are valid, report `is_static()`, give
  exactly zero inverses, and cannot be moved by any impulse.
- A head-on contact between discs whose centres of mass are at their disc
  centres exchanges velocities and produces no spin at all, which is the
  control case.
- The milestone's point: an off-centre contact against an immovable body
  produces spin the body did not arrive with, and mirroring the offset mirrors
  the spin, so both torque signs occur.
- A zero radius means the body passes through others untouched.
- Momentum invariants at three restitutions, plus a zero-penetration case that
  isolates the impulse from the correction artefact.
- Restitution 1 conserves kinetic energy across a contact including its
  rotational term; below 1 strictly removes some.
- Verified by mutation testing. Nine mutations were applied and seven were
  caught immediately: rotational terms dropped from the effective mass, arms
  measured from the body origin, a flipped angular impulse sign, and inverse
  mass used where inverse inertia belongs. Two escaped and exposed real gaps,
  both since closed:
  - the relative velocity at a contact ignoring `omega x r`, which needed a
    spinning body with zero linear velocity whose surface is nonetheless closing;
  - a separating pair being impulsed again, which conserves momentum and so
    needed a direct velocity assertion.
  A tenth mutation, removing the infinity special case from `inverse_mass()`, is
  an equivalent mutant: `1.0 / inf` is exactly `0.0` in IEEE 754, verified rather
  than assumed, so the branch is readability and not correctness.

## [M13] - 2026-09-09  (spring networks and deterministic force accumulation)

- The app refactor was triggered by a measurement rather than by taste. The
  four driving loops differed only in a comment, and that comment was the
  rationale for the int64 loop counter: with `steps == INT_MAX` an int counter
  reaches INT_MAX, passes the loop test, and overflows on increment, which is
  the confirmed infinite loop from #4. The fix was replicated four times and
  explained once, so three copies carried a correctness-critical detail with
  nothing saying why it mattered. That is the real cost of the duplication, not
  the line count.
- `drive()` and `positions_of()` are function templates over duck-typed worlds
  and bodies, not base classes. The four worlds share nothing but a `step()`
  returning `StepResult`, and giving them a common base is exactly what rule 12
  forbids. At M10 a comment justified keeping two copies of the position
  extractor; at three copies the project's own threshold says evidence rather
  than coincidence (ADR 0008), and that comment has been replaced rather than
  left contradicting the code.
- The app still contains no physics: no timestep arithmetic, no velocity
  updates, no trigonometry, no division by mass. Rule 14 holds in substance.
### Added

- `malloy_springs` (STATIC, `malloy::springs`): `Spring`, `SpringNetwork`,
  `SpringBody2D`, the pure `accumulate_spring_forces` kernel, and `SpringWorld`,
  plus momentum, kinetic and elastic-energy diagnostics.
- `malloy_springs_tests`.
- `type springs` in the scenario format, with `spring_body` and `spring` keys.
- `scenarios/spring_chain.scn`.
- `docs/decisions/0008-spring-world-is-a-local-composition-boundary.md`.
- `CLAUDE.md` rule 17, so the ADR 0008 boundary is enforced by the rule list and
  not only by the ADR.

### Architecture Notes

- Springs got their own domain rather than going into `ParticleWorld`. Gravity
  was environmental configuration, one setting applying to everything, so M12
  folded it in. A spring network is interaction topology, which is new state.
  Admitting it into `ParticleWorld` would have admitted cloth, rods, distance
  constraints, cables and breakable joints by the same argument, until the
  particle domain became the general physics engine.
- This is the first many-to-one force pipeline here: a body may be an endpoint
  of several springs and must receive every contribution.
- `accumulate_spring_forces` is a pure kernel. It reads body state, adds into a
  force buffer, and neither integrates nor mutates a body, so a later force
  producer can reuse the pipeline. `SpringWorld` composes it with integration.
- The step contract is fixed and documented: clear the accumulator, evaluate
  springs in stable network order, accumulate equal and opposite endpoint
  forces, integrate bodies in stable index order. Both orders are observable
  rather than incidental, because floating-point addition is not associative
  (docs/04).
- `SpringWorld` is a LOCAL composition boundary, not the engine-wide force
  architecture (ADR 0008). Generalising it into a shared force-provider API is
  deferred until several genuinely different producers exist, because designing
  that interface from one example would be choosing between `step(dt, forces)`,
  `apply_forces`, `accumulator()`, `ForceGenerator` and similar on no evidence.
- The five lines of translational integration duplicated from `ParticleWorld`
  are intentional, with an explicit trigger: if a third domain independently
  needs the same path, reassess extracting a shared integrator. Two copies can
  be coincidence, three are evidence.
- `SpringBody2D` is its own type rather than a reused `Body2D` or `Particle2D`,
  so `malloy_springs` depends on no other domain.
- Rigid-body spring attachment points and the torque they generate stay
  deferred. They would add attachment transforms, torque accumulation and
  angular integration on top of spring topology and damping, which is too many
  concepts for one milestone.

### The spring-damper formulation

For a spring from body a to body b:

```text
d   = x_b - x_a
L   = |d|
n   = d / L                 unit vector, a toward b
e   = L - rest_length       positive when stretched
v_r = (v_b - v_a) . n       positive when separating

F   = (stiffness * e + damping * v_r) * n
F_a = +F
F_b = -F
```

Signs check on both cases. Stretched puts F along +n so the pair contracts;
compressed puts it along -n so they push apart; separating gives a force that
opposes the separation. `F_a + F_b` is exactly zero by construction, which is
why a spring network cannot change total momentum at any stiffness or damping.

Damping uses the AXIAL relative velocity, not the full relative speed. A damper
resists motion along its own axis, and the full magnitude would lose the sign
that makes it oppose rather than drive the motion.

A spring whose endpoints are coincident has no axis and contributes nothing; any
direction would be equally arbitrary.

### Tests

- Hooke's law with exact hand-computed values, in extension and in compression,
  and zero force at rest length with no relative motion.
- Damping opposing relative axial motion in both directions, and a case where
  the relative velocity has a perpendicular component that must not contribute:
  axial 3 of a relative (3,4), so a full-magnitude implementation gives 10
  instead of 6 and a perpendicular leak shows as a nonzero y.
- Many-to-one accumulation, with two springs of DIFFERENT stiffness on one body
  so neither contribution can be mistaken for the total, and the whole-network
  sum being zero whatever the topology.
- That the kernel adds rather than assigns, by pre-filling the buffer.
- Removing one spring changing only that spring's contribution.
- A self-spring, out-of-range endpoints, coincident endpoints, a mismatched
  force buffer, and non-finite values, each handled explicitly.
- End-to-end behaviour a force-only library could not test: a stretched spring
  actually contracts, an undamped network holds kinetic plus elastic within the
  integrator's drift, and a damped one loses more than half.
- Verified by mutation testing: eleven broken implementations each make the
  suite fail, including the four called for by name (accumulation replaced with
  assignment, only the first incident spring processed, the same sign applied to
  both endpoints, and damping using the full relative speed) plus a flipped
  axis, a flipped damping sign, an ignored rest length, an uncleared
  accumulator, and each of the three validation guards removed.

## [M12] - 2026-09-09  (ballistics)

### Added

- Uniform `gravity` in `ParticleSettings`, applied before the position update so
  the integration stays semi-implicit (symplectic) Euler, matching `NBodyWorld`.
- `total_potential_energy` and `total_energy` in `malloy_particles`.
- A `gravity` key for `type particles` scenarios.
- `scenarios/projectile_arc.scn`.

### Changed

- The terminal app reports total energy rather than kinetic energy for particle
  scenarios. With zero gravity the two are identical, so
  `bouncing_particles.scn` still reports 6.82500000e+00 and only the column
  label changed.

### Architecture Notes

- Ballistics shipped as gravity inside `malloy_particles`, not as a separate
  library. A 2D projectile is a colliding particle under gravity, and a
  dedicated domain would have duplicated almost all of `malloy_particles` to add
  one setting. ADR 0006's "duplication beats abstraction" was about a few small
  fields, not a whole world.
- Gravity is world configuration, not body state, so `Particle2D` is unchanged
  and ADR 0006 is untouched. It defaults to zero, which is asserted to reproduce
  the pre-gravity behaviour exactly rather than merely closely.
- No force or torque accumulators were introduced. A uniform field is applied
  directly, so the accumulation-order determinism surface flagged during M11
  still does not exist.

### The energy invariant is NOT conservation

Semi-implicit Euler under a constant field does not conserve energy. Deriving
the step:

```text
v' = v + g dt
x' = x + v' dt
E  = (1/2) m |v|^2 - m (g . x)
E' - E = -(1/2) m |g|^2 dt^2
```

so energy falls by exactly `(1/2) m |g|^2 dt^2` every free-flight step. This is
a constant secular drift, not the bounded oscillation a symplectic integrator
gives on a bounded orbit, because projectile motion is unbounded.

The tests therefore assert the exact per-step change rather than approximate
conservation, and the same for momentum: `dp = (sum m) g dt` exactly per step.
An exact statement of what the integrator does is a stronger test than a loose
statement of what it nearly does, and it does not require claiming something
false. `scenarios/projectile_arc.scn` documents the figure for its own values.

### Tests

- Zero gravity reproduces the pre-gravity behaviour to exact equality over 400
  steps, which is the backward-compatibility guard for every earlier scenario.
- A single hand-computable step pinning that gravity is applied before the
  position update: `g = (0,-10)`, `dt = 0.5` gives velocity -5 and position
  -2.5, distinguishing semi-implicit Euler from both explicit Euler (0) and the
  exact half-a-t-squared solution (-1.25).
- Gravity accelerates equally regardless of mass, with both components nonzero.
- The exact momentum change and the exact energy drift over long runs, each
  asserted against a hand-computed figure, with a check that the drift is far
  larger than the tolerance so the assertion is not vacuous.
- Potential energy sign pinned on both axes, not just the one gravity usually
  uses.
- A dropped ball bounces and settles on the floor over 20000 steps.
- Non-finite gravity is rejected by validation.
- Verified by mutation testing: six broken implementations each make the suite
  fail, including gravity applied after the position update (explicit Euler
  order), integrated with `dt` squared, scaled by mass as though it were a force
  rather than an acceleration, potential energy with a flipped sign or using
  velocity instead of position, and validation dropping the finite check.

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
