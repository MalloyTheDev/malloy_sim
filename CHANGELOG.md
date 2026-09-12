# Changelog

All notable changes to MalloySim are recorded here. The format follows
`docs/10_CHANGELOG_TEMPLATE.md`; dates are ISO-8601.

## [Unreleased]

### Fixed

- M20 review pass. Four defects found reviewing the milestone after it shipped,
  all in the documentation and diagnostics rather than the integrator, and all
  corrected here.

  `rigid::total_angular_momentum3d` summed only the SPIN term `R (I omega)` and
  omitted the orbital `m (r x v)`, so a body flying past the world origin
  reported no angular momentum about it. The 2D `total_angular_momentum` has
  always summed both and warns against dropping either; the 3D one now matches.
  The per-body accessor is renamed `spin_angular_momentum` to say what it is, a
  mutation dropping the orbital term is caught, and the orbital term's exact
  constancy for a free body is asserted.

  The intermediate-axis growth rate in `scenarios/intermediate_axis.scn` was
  written `W sqrt((I2 - I1)(I2 - I3) / (I1 I3))`, whose second bracket is
  NEGATIVE for the intermediate axis. Corrected to `(I2 - I1)(I3 - I2)`, with
  both stable axes' rates spelled out too, since that leading sign is the whole
  theorem. The test comment carried the same slip and is fixed; the asserted
  numbers were already right.

  The renormalization comment in `rigid3d.cpp` still said normalizing perturbs
  world-frame angular momentum. It does not: a nonzero quaternion and any
  scaling of it represent the same rotation. Normalizing is required only
  because `rotate` uses the optimized unit-only form, which returns
  `p + k^2 (R p - p)` for a quaternion of norm k; that property is now its own
  test.

  The zero-drift condition was stated as "a principal axis, or a sphere". The
  exact condition is that `I omega` is parallel to `omega`, i.e. omega lies in
  an EIGENSPACE of the inertia, which for an axisymmetric body includes every
  direction in the degenerate plane and not only its two principal axes. That
  case is now tested directly, and the prose and the drift-law comment say
  eigenspace.

  Also: the angular-velocity step is named FORWARD Euler wherever it is
  described, since the exact O(dt^2) drift is that method's error and the word
  "semi-implicit" belongs to the velocity-before-orientation coupling; the
  finite-horizon nature of the stable-axis result is stated in the template
  (forward Euler amplifies even an oscillatory mode by sqrt(1 + (mu dt)^2) per
  step, so "stable" means at this dt over this run); and the quaternion tests
  are split out of `tests/math/vec2_tests.cpp` into `tests/math/quat_tests.cpp`,
  a fourteenth test executable.

- M19 shipped without scenario parse tests for `type nbody3d` or the `body3`
  key, while every other domain has them. Rule 16 requires malformed-input
  tests and the parser is where malformed input arrives, so the milestone was
  not finished.

  Now covered: all seven `body3` fields with distinct nonzero values, so a
  transposed or dropped one shows and the z components in particular are not
  hidden by zeros; six fields refused, since a 2D body line is not a 3D one
  with a component missing; `body` refused inside `nbody3d` so a scenario
  cannot half-convert; `body3` refused in every other domain; `g` and
  `softening` shown to reach the shared `NBodySettings` and not to leak into
  the charge domain; another domain's keys refused inside `nbody3d`; and a
  softening whose square overflows refused by the same bound the 2D world has.

  Four mutations, all caught: `body3` accepted anywhere, its y and z read in
  the wrong order, `g` no longer serving `nbody3d`, and `softening` routed to
  the wrong settings.

- docs/09's scope-creep watchlist still listed collision, rigid bodies and 3D
  among things not to add. All three have had their milestone, and unlike
  docs/02 that file is live guidance rather than a historical specification, so
  the list is now separated from what has shipped. The note added with it is
  the durable part: the rule is the "until a milestone calls for it" clause,
  not the list of names, because a list goes stale as milestones land and the
  discipline does not.


### Added

- Issue #20: `malloy_particles` gains `total_angular_momentum`, which it had no
  way to report. It runs the same position-only correction `malloy_rigid` does,
  so it carries the same `c x (v_b - v_a)` artefact, and the quantity that
  artefact perturbs was not measurable in that domain at all.

- The wall clamp's energy contribution is derived and asserted rather than
  merely observed.

  Clamping a particle back to a wall changes its position without changing its
  velocity, and against gravity that ADDS potential energy: `m |g| delta` once
  per event, with `delta` bounded by how far it overshot in one step, so
  `0 <= delta < |v| dt`.

  The project had characterized exactly one integrator error term, the
  `(1/2) m |g|^2 dt^2` free-flight loss. This one is comparable and opposite in
  sign: at dt = 4e-3 an elastic ball on a floor gains +0.94 from clamping
  against a free-flight loss of -1.6.

  The new test puts a whole step in closed form, with restitution 1 so the
  impulse changes no kinetic energy:

  ```text
  gravity    v1 = -6 - 1   = -7
  position   y1 = 1 - 0.7  =  0.3
  clamp      surface 0.5, so delta = 0.2

  E0 = 56,  E1 = 59,  change = +3
     = free flight -1  +  clamp +4
  ```

  Both terms are asserted, and so is the fact that the clamp dominates the
  characterized one by a factor of four here. It is not a rounding-level
  correction.

  The rigid domain's version of the same term is `m |g| pen` per event, by the
  same derivation, and the header for `ParticleWorld::step` now states the
  contract rather than leaving it implied.

  Three mutations, all caught.


### Fixed

- Issue #16: `is_valid` accepted state that every diagnostic then reported as
  `inf`.

  Everything which squares a vector reaches infinity above sqrt(DBL_MAX), about
  1.34e154, while the vector itself is still finite and every component is a
  normal number. Validation tested only `is_finite`, which accepts up to
  1.8e308, so a window covering the entire upper half of the exponent range
  reported a world as sound while every energy it printed was infinite.

  `math::is_squarable` is the exact predicate, and all five body types now
  require it on position and velocity, plus `local_center_of_mass` for a rigid
  body since `center_of_mass` squares that too. Five mutations, one per domain,
  all caught.

- Issue #19: `tests/nbody/nbody_world_tests.cpp` asserted a radial deviation at
  0.02 against a real 4.21e-04, and that 4.21e-04 is not physics.

  The stored velocity in semi-implicit Euler is half a step behind the
  position, so a tangential initial velocity is short of the discrete circular
  orbit by a radial component and seeds an epicycle of eccentricity
  `dt * Omega / 2`. For that world it is 5.0e-04 against a physical
  eccentricity of about 1e-6, so the wobble is roughly 500 times the physics.

  The run only reaches part of it, and the shortfall is part of the derivation
  rather than slack: the deviation grows as `R e sin(Omega t)`, and 1000 steps
  is `Omega t = 1.0` rad, about a sixth of an orbit. Predicted maximum
  `5.000002e-04 * sin(1.0000005) = 4.207358e-04`, measured 4.210280e-04, which
  is 0.07 per cent.

  Now bounded from BOTH sides. Correcting the initial conditions with a
  half-step kick would drop this into the 1e-6 range, and that changes what an
  initial velocity MEANS across the whole N-body domain and every template
  figure with it, so it should fail here and be updated deliberately rather
  than pass unnoticed. That remains a milestone rather than a patch.

- Issue #18: the angle accumulates by repeated addition, and how far it drifts
  is now written down and pinned.

  `malloy_time` computes elapsed time as `tick_count * dt` specifically so it
  cannot drift that way. The angle is the same pattern and cannot get the same
  treatment, because angular velocity is changed by contacts and there is no
  constant increment to multiply.

  Summing N terms accumulates at most `u * theta * N / 2` with `u = 2^-53`. A
  test now compares 100000 steps of constant spin against the exactly computed
  `N * omega * dt` and holds inside that bound. ADR 0007 is amended with the
  measurement and with why the two alternatives, a per-body tick count and
  compensated summation, are not taken: both add state to `RigidBody2D` for a
  quantity nothing reads at that precision.

### Added

- A coordinate magnitude budget in docs/04, which had no such convention. Two
  limits, and the smaller is the one that matters: `position += velocity * dt`
  stops advancing once the increment falls below half an ulp, so at
  `|x| ~ 1e13` with metre-scale motion at millisecond steps a body travels
  1.95x too far and at 1e14 it freezes, while every diagnostic reads perfectly
  conserved because they are all computed from velocity. Budget: keep `|x|`
  below about 1e12, scaled by `v * dt`.

- The half-step velocity stagger is documented in docs/04, and
  `scenarios/two_body.scn` says which part of its wobble is numerical.


### Added

- Issue #13: a layering test, `malloy_layering_tests`, which checks the module
  graph against the build files and the sources.

  The graph was enforced at LINK time and nowhere else. Every module declares
  `target_include_directories(<mod> PUBLIC include)`, so the whole include tree
  is on every consumer's search path, and a module can include another's header
  without linking it. The build succeeds, and the rule it breaks is the one the
  architecture rests on: domains share only `malloy::math` and the sim_core
  vocabulary, and no domain knows that any other exists (ADR 0006).

  Three things are checked. Every `#include <malloy/X/...>` must be covered by
  a declared link, which is the check the link step cannot make. No physics
  domain may link or include another, stated directly rather than inferred, so
  adding the link would not make it legal. And `malloy_sim_core` may depend on
  `malloy::math` and nothing else, which is CLAUDE.md rule 13 written down as a
  test.

  The declared graph is parsed out of CMakeLists.txt rather than kept in a
  table here, so the two cannot drift. Both spellings of
  `target_link_libraries` are handled, since INTERFACE targets put it on one
  line and STATIC ones spread it over several.

  Verified by mutation: a domain including another domain's header, a domain
  declaring a link on another domain, and `sim_core` growing a dependency are
  all caught, each with a message naming the file and the rule.

### Fixed

- `malloy_scenario` used `collide::Halfplane` and `collide::Aabb` without
  declaring `malloy::collide`, found by the test above on its first run. It
  linked only because collide arrives transitively through particles and rigid,
  so the day either of those stopped depending on it, scenario would have
  broken with an error pointing somewhere else entirely.


### Added

- Issue #14: scenario templates now carry machine-checkable assertions, and the
  scenario tests evaluate them against a real run.

  A template's documented figures are a deliverable, and until now nothing
  compared them against anything. The words `# Expected` were enforced; the
  numbers under them were not. That is how two figures shipped wrong in this
  repository and were corrected by hand, and how the `rolling_and_slipping`
  table had to be regenerated by hand after the lever-arm fix moved it.

  The form is one line, ignored by the parser like any other comment:

  ```text
  # check step <n> <quantity> <value> tol <t>
  ```

  `step 0` is the state before any step, matching the column the app prints.
  The tolerance is ABSOLUTE and required rather than defaulted, so a template
  states the precision it is claiming instead of inheriting one.

  Quantities available in every domain are `energy`, `kinetic`, `momentum_x`
  and `momentum_y`, plus `angular` in nbody and rigid, `elastic` in springs and
  `speed` in charges. Naming a quantity a domain does not have FAILS rather
  than being skipped, so a check cannot quietly stop checking.

  Every one of the eleven templates now carries at least one, and a template
  with none fails the suite. Evaluation is a function template over the world
  type rather than a base class, the same approach the app uses for `drive()`,
  since the five worlds share no ancestor (ADR 0006).

  Verified six ways, all caught: a documented value that drifts, a tolerance
  tighter than the value's own precision, a template that loses all its check
  lines, a malformed line, a check naming a step past the end of the run, and a
  quantity the domain does not have.

  Cost: the scenario suite goes from 0.22s to 0.28s, and the whole suite runs
  in 0.71s.

### Fixed

- The scenario format reference still listed four domain types after M18 added
  a fifth.


### Fixed

- Issue #17: a softening value large enough that its SQUARE overflows was
  accepted, and silently switched the pairwise interaction off.

  `softening` enters the denominator squared, so a finite value above
  sqrt(DBL_MAX), about 1.34e154, squares to infinity. Every separation then
  becomes infinite, every force is exactly zero and every potential is exactly
  zero. The failure is not merely undetected, it looks BETTER than a correct
  run, because a simulation in which nothing happens conserves everything
  perfectly.

  Both domains that have a softening now require its square to be finite, which
  is the exact condition and needs no magic constant.

  Worth recording how this was found. The issue was filed against `malloy_nbody`
  after an audit. Checking whether M18 had inherited it showed that it had: the
  `softening` added to `malloy_charges` one commit earlier had the identical
  gap, and two like charges that should have flown apart sat still for 2000
  steps while the energy column read `0.00000000e+00` against a correct 0.5. A
  filed issue caught a fresh instance of itself.

- Issue #15, both halves, in `malloy_springs`.

  A spring whose endpoints are more than about 1.34e154 apart has a separation
  whose square overflows, and the force kernel could not tell that from
  coincident endpoints: the guard documented for the coincident case swallowed
  it, so the spring silently ceased to exist while the bodies coasted apart
  forever and every step still returned Ok. Measured: kinetic energy of exactly
  `0.00000000e+00` for 2000 steps with a stiffness of 1.

  And nothing bounded the timestep against a spring's stiffness. Symplectic
  Euler on a spring pair diverges geometrically past a bounded dt, and until
  the separation reached the range above, nothing reported that either. Just
  under the limit the "conserved" energy of an undamped network climbed by a
  factor of 25000 while momentum still read exactly zero; just over it, the run
  printed `inf` for 6000 steps with every step Ok. The first defect is what
  made the second silent: once the separation overflowed, the spring was
  dropped, the bodies coasted, and an exponential blow-up became a permanently
  finite and permanently wrong state.

  `SpringWorld::validate()` now rejects both. The separation is checked every
  step rather than only at load, because one that starts representable can grow
  into that range during a run, which is exactly what a diverging network does.

  The stability bound is `dt < 2(sqrt(gamma^2 + omega^2) - gamma) / omega^2`
  with `omega^2 = k/mu` and `gamma = c/(2 mu)`, reducing to `dt < 2/omega` when
  undamped. Damping TIGHTENS it rather than helping, which the undamped formula
  would get wrong, so a case that is inside the undamped bound and outside the
  real one is asserted.

  The check is per spring, which is NECESSARY but not sufficient for a network:
  a network's stiffness matrix is a sum of positive-semidefinite per-spring
  terms, so its largest eigenvalue is at least that of any single spring, and a
  spring unstable alone is unstable in company. The converse does not hold, and
  the header says so rather than implying a guarantee it cannot make.

  The shipped `scenarios/spring_chain.scn` runs with a 321x margin, so no
  template is affected.

  Verified by mutation: removing the overflow check, removing the stability
  gate, and dropping damping from the bound are all caught.



### Fixed

- A statement-ordering bug in `apply_contact` made every contact lever arm
  `R + depth/2` instead of `R - depth/2`. The positional correction lifted the
  centre of mass by the full penetration, and only then was the arm measured
  against a `hit->point` still holding its pre-correction value, so the two
  quantities came from different frames. Of the three consistent answers, that
  is the one that errs outward, and nobody chose it.

  All contact arithmetic now takes one arm, fixed before the correction, in the
  frame the contact data actually describes. A new `velocity_at_arm` helper
  takes the arm directly, because `velocity_at` re-derives it from the body's
  CURRENT centre of mass and so silently measured two different arms depending
  on where it was called. The friction tangent was the third site and was still
  reading the post-correction frame after the first two were fixed.

  The normal impulse is provably unaffected: the arm changed purely along the
  normal, and `cross(k*n, n)` is zero, so `cross(arm, normal)` is identical
  either way. Confirmed by measurement, and it is why exactly one of the ten
  templates moved. `scenarios/rolling_and_slipping.scn` is the only one with
  friction, and its figures are updated.

  The correction lands on the derived value: with `a = R - depth/2` the closed
  form predicts a final momentum of 11.333309898 and a final energy of
  56.93001145, and the binary now prints `1.13333099e+01` and `5.69300114e+01`,
  which is exact at the nine significant figures it reports.

- Three explanations that were wrong while their numbers were right, all in
  work shipped in M16 and M17, and all found by independently re-deriving the
  physics rather than by re-reading the code:

  - The M17 lever arm above. The documented `R - depth/2` predicted a residual
    of the correct MAGNITUDE and the wrong SIGN, and the assertion guarding it
    was a magnitude bound, so nothing could catch the discrepancy. That
    assertion now derives a signed expectation from the run's own gravity and
    holds to 1e-9, and reverting the ordering fix is caught by it.
  - `scenarios/rolling_and_slipping.scn` claimed friction exerts no torque
    about the world origin and blamed the observed angular-momentum drift on
    gravity and the normal impulses "nearly" cancelling. Backwards on both
    counts. The contact point sits half a penetration off the origin line, so
    friction is the only net torque, and it accounts for the drift as
    `(depth/2) * (p - p0)`. Gravity, the push-out and the normal impulse cancel
    EXACTLY rather than nearly. The file's own table was the disproof all
    along: L goes bit-constant the moment friction switches off at step 3000,
    which a residue from gravity could not do.
  - `scenarios/flat_ground.scn` attributed its final energy change to the wall.
    It is the ramp: across that event the body accelerates at
    `(+4.7088, -3.5316)`, the frictionless-ramp signature the same file derives
    two points earlier, and its horizontal momentum RISES, which a wall normal
    of `(-1, 0)` cannot do.

- Two numbers offered as evidence in the M14 angular-momentum test were wrong.
  The short-run drift is `-1.20527e-09`, not "around 5e-9", and the two
  penetrations differ by 5.9 orders, not nine. The comment also contradicted
  itself: 5e-4 reduced by nine orders is 5e-13, and back-solving its own
  figures implied a penetration of 2.0 between discs of radius 0.5. Both
  measured values are now quoted, along with the reason the ratio is not the
  penetration ratio alone: the perturbation is `c x (v_b - v_a)`, so it scales
  with the relative velocity too.

- Two overstated precision claims. `flat_ground.scn` claimed agreement "to
  9e-16" from values printed to nine significant figures, which can only
  support about 1e-8; the true double-precision figure is 8.9e-15.
  `dropped_bodies.scn` presented `1.563837e-01 against a measured 1.563840e-01`
  as a discrepancy, when the print resolution on those values is about 5e-7 and
  the prediction sits inside it. There was nothing to explain.

- The scalloped-floor numbers described two different surfaces. The 0.0546
  ripple belongs to the floor's own surface and goes with a normal tilt of
  24.62 degrees; a radius-0.4 body rides a locus that ripples by 0.0318, and
  that is the one whose normal swings by 14.48 degrees. No body ever rides a
  0.0546 ripple. The judgement that the normal swing is the number that matters
  was right; the two figures quoted beside it were not a matched pair.

### Changed

- Two N-body assertions were far looser than the errors they measure, and one
  could not fail at all.

  Angular momentum is EXACT for this integrator, to all orders in dt: for
  pairwise central forces the (i,j) and (j,i) torque contributions cancel
  identically, and softening changes only the magnitude. The only possible
  error is rounding, measured at 1.8e-15 on a quantity of order 2.3, about
  8 ulp. It was asserted at 0.02, roughly 1e13 times too loose, which would
  have passed for an integrator that was not symplectic at all. Now 1e-12.

  Energy oscillates about the shadow Hamiltonian with bounded amplitude rather
  than drifting; measured at 7.42e-08 and asserted at 0.02. Now 1e-6.



### Added

- `docs/decisions/0009-three-dimensions-are-the-destination.md`. Every mention
  of 3D in the documentation was a prohibition, and it sat in the same list as
  plugins, an editor and a scripting layer, which are rejected rather than
  deferred. A reader had no way to tell the two apart, and the natural reading
  of README.md was that MalloySim will always be 2D. It will not: 3D is the
  destination.

  The ADR records what that means for the code as it stands, so a 2D decision
  can be checked against it instead of guessed at. The fixed timestep, the
  sim_core status and validation model, semi-implicit Euler, the multi-domain
  dispatch of ADR 0006, the determinism policy and the testing approach carry
  over. Vec2, the scalar angle, the scalar inertia, the scalar cross product and
  the impulse formulas SPECIALIZE, meaning the 2D form is the 3D one with a
  dimension removed, which is the justification ADR 0007 already relied on when
  it chose a scalar angle. `malloy_ascii` and the area-based collision
  properties do not carry over at all.

  It also records the part that is genuinely new, and is why 3D is a milestone
  rather than a widening: in 2D the inertia is a scalar and the angular velocity
  lies along a fixed axis, so a body's inertia never changes in world space and
  there is no gyroscopic term. In 3D the tensor rotates with the body and
  Euler's equations carry omega x (I omega), so a spinning body precesses and
  can tumble about its intermediate axis under no torque at all. That is physics
  the current code does not contain in any form.

  Three things are deliberately left undecided: whether 3D types are separate or
  the existing ones become templates, whether 3D domains sit beside the 2D ones
  or replace them, and whether 2D stays supported.

- An ADR-range guard in the scenario tests, since docs/00_START_HERE.md states
  which ADRs exist and that number went stale the moment 0009 was written. The
  highest number in `docs/decisions/` is now the truth, and a stale claim fails
  the suite. This is the fourth documented count derived rather than maintained,
  after the template count, the test-executable count and the milestone range.

### Changed

- README.md and CLAUDE.md now separate what is PLANNED and gated from what is
  REJECTED. 3D and quantum were listed among plugins, editors and scripting
  layers, which made a roadmap item look like a refusal. Nothing about the
  gates changed: 3D still arrives as its own milestone and nothing 3D is added
  before then (rules 6 and 7).
- docs/03, docs/07 and ADR 0007 no longer defer 3D work "to M19". That number
  was assigned before M14, M15 and M16 shipped into the sequence below it, and
  Track 4 is now unnumbered like Tracks 2 and 3.
- docs/08_AI_HANDOFF_PROMPT.md no longer instructs a reader not to add collision
  or rigid bodies, both of which shipped in M9 and M11.


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

- Immovable bodies with a velocity are now documented and tested as KINEMATIC
  bodies: moving platforms that push everything and are pushed by nothing.
  Gravity does not accelerate one, no impulse can slow it, and it carries no
  momentum in the diagnostics, so it acts on the world without being part of
  it. Nothing specified this in either direction, so adding a static skip to
  the pose loop passed the suite just as well. Behaviour is unchanged; the
  contract is now written down in the step description and pinned by a test.

### Fixed

- `RigidBody2D::is_static()` was an OR over two independently inverted
  quantities, so a body with finite mass and infinite inertia reported itself
  static while `inverse_mass()` remained non-zero. Contacts moved it and the
  pose loop integrated it, while gravity and all four diagnostics skipped it
  entirely. The momentum handed to such a body vanished from the report, which
  silently falsified the exact linear-momentum conservation that
  `rigid_world.hpp` claims for every contact.

  Infinity is now per quantity throughout, which is what the M14 note already
  described: infinite mass alone is a body that can spin but not translate,
  infinite inertia alone one that can translate but not spin, and `is_static()`
  is an AND meaning immovable in both senses. `has_infinite_mass()` and
  `has_infinite_inertia()` are the questions the consumers now ask. Gravity
  tests the mass, so a body that merely cannot spin still falls. Angular
  momentum and kinetic energy guard their translational and rotational terms
  separately, so the real half of a half-infinite body is reported and only the
  infinite half is dropped.

  BREAKING for anyone relying on the old meaning: `is_static()` now returns
  false for a half-infinite body, and the M14 assertion that pinned the OR
  behaviour changed with it, deliberately and with the reason recorded beside
  it.

  Found by an audit for assertions that cannot fail, and confirmed the same way
  as the rest: BOTH directions of the inconsistency escaped the suite, so the
  tests specified neither semantics. Six mutations now cover it, including each
  diagnostic reverting to the blanket check, and all six are caught. One of them
  is caught by the NaN-rejecting macro fixed earlier in this entry, because
  reverting the momentum guard lets `inf * 0` reach the assertion.

  Also newly tested: two bodies of infinite mass and finite inertia are no
  longer short-circuited by the static early-out, so they reach the impulse
  code with an effective mass of exactly zero. The existing guard returns
  rather than dividing, but nothing had exercised that path.

- Scenario templates are now run for their FULL documented number of steps by
  the test suite, not one step. A single step cannot tell a template that runs
  from one that diverges on step 900, and `step()` reports `InvalidState` as
  soon as any body stops being finite, so running to the end is also the NaN
  check. Costs about 0.2 seconds across all eight in a Debug build.

  Verified by inverting the sign of the rigid positional correction, so contacts
  push bodies together instead of apart: the template run fails, where a
  single-step check passed. Comparing the documented FIGURES against a run is a
  separate and larger job, filed as issue #14.

### Fixed

- Documentation claims contradicted by the code, all confirmed against the
  source rather than by review:
  - README.md listed "rigid bodies" and "rigid-body contact response" under
    **Out of scope (gated)**, a section headed "never added speculatively",
    while the same file's milestone table marks both Done. They shipped in M11
    and M14. Replaced with what is genuinely still out of scope: oriented-box
    contacts and SAT, and friction, since every contact in the project is
    normal-only and bodies therefore slide forever.
  - docs/00_START_HERE.md said the app "runs all three domains" and named
    three, omitting springs. It contradicted its own line 44. There are four.
  - docs/07_POST_M5_ROADMAP.md still headed "M1-M14 are complete" while its own
    Complete block listed M15 as done.
  - docs/07 numbered its unscheduled candidates M14 through M18, colliding with
    the real M14 and M15. Those blocks are now unnumbered, because a number
    assigned before the work starts is a number that collides with whatever
    actually ships.
  - docs/03_MODULE_BOUNDARIES.md and the matching comment in scenario.hpp both
    said semantic validation "stays in NBodyWorld". Since M10 to M13 there are
    four domains, each with its own validate().
  - ADR 0007 said M11 "has not started" with no date on the claim, unlike ADR
    0002 which dates the identical kind of statement.
  - CMakeLists.txt still described the project as "a simulation-first C++20
    engine project" although the CHANGELOG records that identity being amended.
  - The parse error for a bare `type` key named two of the four accepted
    domains.

### Added

- A milestone-status guard in the scenario tests. The "M1-M<n> complete" line
  has gone stale three times: corrected once in the pre-M13 audit, again
  through M13 and M14, and again through M15. Correcting it a fourth time would
  only reset the clock, so it is now derived. CHANGELOG.md is the source of
  truth, since a milestone is complete exactly when it has an entry, and every
  such claim across CLAUDE.md, README.md and docs 00, 07 and 08 is checked
  against it. M1-M5 is exempt, being the locked original roadmap and a fixed
  historical range. A document that stops stating the range at all also fails,
  so the check cannot go vacuous.

  Verified three ways: a doc left at the old number fails, a doc that drops the
  claim fails, and adding a new milestone entry without updating the docs fails.

### Fixed

- Six coverage gaps found by auditing the test suite for assertions that cannot
  fail. Each was confirmed by applying the mutation and watching the whole suite
  stay green, then confirmed closed by watching the new assertion catch it.

  - The spring force kernel was only ever tested with a horizontal spring. Every
    assertion that depended on the value or direction of a spring force used
    endpoints like (0,0) and (3,0), so replacing `dot(delta, delta)` with
    `delta.x * delta.x` passed the entire suite. The two genuinely 2D tests
    could not see it: momentum conservation is structural, since the code writes
    `forces[a] += force; forces[b] -= force;`, and the determinism test compares
    two worlds carrying the same defect. Under that mistake every VERTICAL
    spring would hit the coincident-endpoint guard and produce exactly zero
    force forever. Now covered by an off-axis, off-origin spring with an exact
    3-4-5 separation, and by a purely vertical one.
  - Wall restitution was unpinned on the left and top walls. The four-wall table
    ran at restitution 1, where the multiplier is a no-op, so deleting
    `* restitution` from either branch changed nothing. Only the right wall was
    ever checked at a value other than 1. Now all four run at 0.4, each with a
    different tangential speed so a swapped component shows as well.
  - The softening term in the N-body POTENTIAL energy was never pinned, though
    the force side was. Every call used either softening 0 or a softening so
    much smaller than the separation that dropping it moved the result by about
    3e-7 under tolerances of 0.02. Now checked with softening comparable to the
    separation, where the answer is exactly 5 rather than 3.
  - Two of the four circle-in-box faces, left and top, were never reached with a
    distinguishing configuration, so their normals were unconstrained.
  - Box against box picks a sign from the relative centres, and every existing
    case put the second box to the +x or -y of the first, leaving two of the
    four sign branches unasserted.
  - `MALLOY_CHECK_TRUE(grid.size() > 0)` in the ASCII tests could not fail:
    `render` always emits a border. Replaced with the exact frame geometry.

- `fit_viewport`'s header contract omitted its overflow fallback. Points near
  the limits of double make the extent overflow to infinity even though every
  input is finite, and the function abandons the fit and returns the default
  box. The consequence was undocumented and untested: that box does NOT contain
  the points it was built from, so they are not drawn at all. Both are now
  stated in the header and pinned by the test, which previously asserted only
  that the viewport was finite.

## [M32] - 2026-09-12  (colliding particles in three dimensions)

### Added

- `particles::Particle3D`, `particles::ParticleSettings3D` and
  `particles::ParticleWorld3D`, the 3D sibling of the M10 particle domain:
  spheres that collide and are confined to a box, with the diagnostics
  `total_momentum3d`, `total_kinetic_energy3d`, `total_potential_energy3d`,
  `total_angular_momentum3d` and `total_energy3d`.
- `collide::Aabb3`, an axis-aligned 3D box, the walls the particles are confined
  to (the 3D sibling of `Aabb`).
- A `particles3d` scenario type (`particle3`, `bounds3`, and the shared
  `restitution` and `gravity3` keys) and the `particle_collision3d` template.
  Twenty-one templates now, across ten domains.

### The last domain to reach three dimensions

Every one of the five domains now has a 3D form; M32 lifts the last of them.
Particles collide as spheres, which is the geometry M22 and M24 already built,
so the contact response is the 2D one with a Vec3 normal: a positional
correction plus an equal-and-opposite impulse along the line of centres, which
imparts no spin (a particle has no orientation). Wall containment is a per-axis
clamp-and-reflect against an `Aabb3`, with a third axis added. Gravity and the
diagnostics lift the same way.

### Invariants

An elastic head-on collision conserves both total momentum (the impulse is equal
and opposite) and total kinetic energy (a perfect bounce loses none), and equal
masses swap velocities; an inelastic one keeps the momentum but sheds energy. A
particle confined to the box never escapes on any axis and, with elastic walls,
keeps its speed exactly. And a run confined to the z = 0 plane matches the 2D
`ParticleWorld` step for step, including gravity, walls and the divide-by-mass,
which ties the new world to the proven one. The shipped `particle_collision3d`
template sends two spheres head-on along the space diagonal and holds momentum
and energy through the bounce.

### Mutation testing

Twelve mutations, all caught: the collision impulse's equal-and-opposite split,
its restitution factor and closing-velocity guard, the wall clamp, the gravity
kick, the particle validation, the momentum, kinetic and potential diagnostics,
`Aabb3`'s min-less-than-max check, and the scenario loader's field order and
gravity routing.

## [M31] - 2026-09-11  (spring networks in three dimensions)

### Added

- `springs::SpringBody3D` and `springs::SpringWorld3D`, the 3D siblings of the
  M13 spring domain, with a Vec3 overload of `accumulate_spring_forces` and the
  diagnostics `total_momentum3d`, `total_kinetic_energy3d` and
  `total_elastic_energy3d`.
- A `springs3d` scenario type (`spring_body3` for a body; the existing `spring`
  key is reused) and the `spring_tetrahedron` template. Twenty templates now,
  across nine domains.

### The dimension is not a domain

`Spring` and `SpringNetwork` carry only body ids and the scalar rest length,
stiffness and damping, none of which has a dimension, so they are shared with
the 2D world UNCHANGED: M31 duplicates only the body and the force. Hooke's law
and the damper are the same in space as in the plane (the force acts along the
spring's axis, a Vec3 now), so the kernel, the symplectic integrator and the
per-spring stability bound are a direct Vec2-to-Vec3 lift. What it buys is
structures that deform in three dimensions: a chain, a cloth, a lattice.

### Invariants

Total momentum is conserved exactly, at any stiffness and any initial motion,
because every spring's endpoint forces are equal and opposite. Undamped, kinetic
plus elastic energy is conserved up to the symplectic drift; damping removes it
on purpose, which is testable. And a network confined to the z = 0 plane matches
the 2D `SpringWorld` step for step (including the damping term and the
divide-by-mass), which ties the new integrator to the proven one. The shipped
`spring_tetrahedron` template releases four masses on six stretched springs from
rest: it breathes in and out, its momentum stays exactly zero, and its energy
holds near the value the initial stretch stores.

### Mutation testing

Twelve mutations, all caught: the force's equal-and-opposite split, its
stiffness, extension, axis and damping terms, the divide-by-mass in the
integrator, the body validation (mass and squarable state), the momentum,
kinetic and elastic diagnostics, and the scenario loader routing a spring to the
wrong network.

## [M30] - 2026-09-11  (an oriented box on a ground plane)

### Added

- `collide::Box3`, an oriented box (centre, half-widths, orientation), with
  `overlaps(Box3, Plane3)` and `contacts(Box3, Plane3)`, the first contact query
  that returns a MANIFOLD (the list of the box's penetrating corners) rather
  than a single point.
- `rigid::RigidBody3D::half_extents`, a box collider: when set, the body
  collides as an oriented box against ground planes instead of as a sphere. Zero
  (the default) keeps the sphere behaviour, so every pre-M30 body is unchanged.
- A `rigid_box3d` scenario key, whose inertia is computed from the box and mass
  (M26) rather than typed in, and the `box_drop` template. Nineteen templates
  now, across eight domains.

### The first non-sphere 3D collision

Until now a 3D body only ever collided as a sphere, even though M25 to M27 could
compute the inertia of a box or an arbitrary mesh. M30 closes that loop for a box
against a plane. It needs no separating-axis search: a plane has one normal, so
the box is inside the solid exactly where its corners are, and the eight corners
tested against that normal are the whole story. Up to four touch at once (a flat
face), so the contact is a manifold.

That manifold is reduced to a SINGLE contact at the centroid of the penetrating
corners, resolved through the same impulse core the sphere uses. For a flat face
the centroid lies directly below the centre of mass, so the normal impulse
passes through it and makes no torque: a flat drop stays flat and its angular
momentum is conserved exactly, where a per-corner sequential solver would have
spun it up out of nothing. For a corner or edge landing the centroid is that
corner or edge, the correct lever arm, so the box tips and tumbles. A full
per-corner manifold with the iterative solver a stack needs is a later
milestone, the same boundary the 2D solver drew.

A box's inertia and its collision shape describe the one box, because
`rigid_box3d` computes the inertia from the geometry (M26). Box against box,
which needs a separating-axis test and an edge-edge case, is deferred; a box
collides with ground planes only.

### Invariants

A flat symmetric drop induces no spin and no horizontal drift, exactly, through
the whole bounce, and settles resting with the centre of mass at half a side
height. A corner-first landing does spin the box up and never sinks through the
floor. A box placed at rest stays put. The geometry is checked against closed
forms: an axis-aligned cube reports its four bottom corners at the right depth
(and their x*y products cancel, so all four sign combinations are present), and
a box tilted 45 degrees reports exactly the two corners that dip below, at
0.5*sqrt(2) - h.

### Mutation testing

Twelve mutations. Ten caught (one after tightening the corner test so a
collapsed corner set no longer passes): the corner enumeration, the penetration
depth and its sign, the contact normal in both the geometry and the resolution,
the deepest-penetration push-out, the box-collider validation, and the scenario
loader's density, half-extents and computed inertia. Two survive: reading a zero
half-extent as a box is a true equivalent (the resulting `Box3` is invalid and
yields no contact either way), and not averaging the manifold centroid is a
coverage gap, since it changes the lever arm only on an off-centre landing whose
exact spin the symmetric tests here do not pin.

## [M29] - 2026-09-11  (charged particles in three dimensions)

### Added

- `charges::ChargedParticle3D`, `charges::Charge3DSettings` and
  `charges::Charge3DWorld`, the 3D sibling of the M18 charges domain, with a full
  VECTOR magnetic field and helical motion, plus the diagnostics
  `total_momentum3d`, `total_kinetic_energy3d`, `total_potential_energy3d` and
  `total_energy3d`.
- A `charges3d` scenario type (keys `coulomb`, `efield3`, `bfield3`, `softening`,
  `charge3`) and the `magnetic_helix` template. Eighteen templates now ship
  across eight domains.

### The physics two dimensions could not express

M18's magnetic field was a scalar, out of the plane, and it turned a particle in
a flat circle. In three dimensions the field is a vector with a direction, and
q(v x B) turns only the velocity ACROSS the field, leaving the component ALONG
it untouched. So a particle circles across the field and drifts along it at
once: a helix. That combination has no 2D form, the same way M20's torque-free
tumbling did not.

The magnetic turn is still the exact rotation M18 insisted on, now about the
field DIRECTION rather than a fixed axis:
`rotate(from_axis_angle(B, -(q |B| / m) dt), v)`. A magnetic force does no work,
so this cannot change speed, and because it is a rotation about B it carries the
along-field component to the last bit. Applying it as a kick would inflate the
speed by sqrt(1 + (q |B| dt / m)^2) every step, exactly as in 2D. The electric
and Coulomb terms are the 2D ones with a third component.

### Invariants

The speed is preserved exactly under a magnetic field, even a tilted one and
even with a velocity component along it (a 3-4-5 velocity keeps |v| = 5 to every
digit); the energy is then purely kinetic and constant. The along-field velocity
is carried bit for bit (a rotation about z leaves z untouched), the drift half
of the helix. A full cyclotron period returns the velocity, and the turn's sense
and rate match the field. The electric field is a force per charge divided by
mass, so a heavier charge accelerates less and an opposite one the other way.
And a run confined to the z = 0 plane, with B along z, matches the 2D
`ChargeWorld` step for step, which ties the new integrator to the proven one.

### Mutation testing

Fourteen mutations, all caught: the magnetic angle's sign, its division by mass
and its factor of |B|, and the rotation axis; the pairwise force's
equal-and-opposite split, the electric field's division by mass, and the Coulomb
sign; the particle and settings validation; and the kinetic, momentum and
potential diagnostics.

## [M28] - 2026-09-11  (a constant applied force on a 3D body)

### Added

- A constant applied FORCE on `Rigid3DSettings`, in the world frame: the
  translational half of a wrench, where M21's torque is the rotational half.
- `rigid::total_force_potential3d`, the potential energy of that force,
  -(F . r) summed over bodies, so energy still balances when a force acts.

### The other half of the wrench

M21 put a torque on a 3D body; M28 puts a force. It enters the step as the
acceleration F/m it produces, in the same place gravity does (before the
position update, which keeps the world semi-implicit), but unlike gravity it
scales with 1/mass. That is the whole difference between a force and an
acceleration, and it is why the force reads the body's mass, now computed from
geometry (M25 to M27). It acts through the centre of mass, so on its own it
makes no torque: a centre-of-mass force has no lever arm.

Together the force and torque settings express any constant wrench at the centre
of mass, which is how an off-centre push is represented (a force plus the couple
r x F). So M28 needs no application point and no lever-arm machinery: the couple
is the torque that is already there.

It is a single constant SETTING applied as forcing, not the force accumulator
rule 5 defers. M28 is that rule's dedicated milestone for a translational force,
and it adds one force, not a pipeline that sums many; the rule now reads as a
ban on accumulators rather than on forces.

### Invariants

Newton's second law, discretely exact in momentum: a body from rest under a
constant force reaches velocity (F/m) t and momentum F t, because each step adds
the impulse F dt exactly. The force scales with mass, so two spheres of the same
size and different density gain the same momentum while the lighter goes faster,
in inverse proportion to mass, which reads the mass computed in M25 to M27. A
centre-of-mass force makes no torque: a tumbling body pushed by a force has
bit-for-bit the same angular velocity and orientation as the same body left
free, while its velocity differs. The complete wrench is orthogonal: with force
and torque both set, the angular state matches a torque-only run and the
velocity a force-only run, both exactly. And kinetic plus the force potential
sheds exactly (1/2) |F|^2 / m dt^2 per free-flight step, the force's analogue of
the gravity drift.

### Mutation testing

Six mutations, all caught: the force having no effect, the force not divided by
mass, its sign flipped, the force validation dropped, and the force potential's
sign flipped or wrongly mass-scaled.

## [M27] - 2026-09-11  (mesh mass properties)

### Added

- `rigid::SolidMesh` and `mass_properties_3d(const SolidMesh&)`, the general
  mass primitive: a solid uniform body bounded by a closed, outward-wound
  triangle mesh. Its mass, centre of mass and inertia are the volume integrals
  over the enclosed solid.
- `math::trace`, the diagonal sum, which turns a covariance tensor into an
  inertia tensor (I = trace(C) I - C).

### The general shape, not a third special case

M25 did spheres and M26 did boxes; both are shapes a mesh can also describe, so
M27 is the general one. A closed triangle mesh's volume integrals are computed
without meshing the interior, by the divergence theorem: each triangle spans a
signed tetrahedron with the origin, and summing the per-tetrahedron integrals
gives the integral over the enclosed solid wherever the origin sits relative to
the body, because the parts outside the solid cancel between overlapping signed
tetrahedra. The volume is the signed sum of det/6 (det = a . (b x c)); the
centre of mass is the signed sum of the tetrahedra centroids; and the second
moments come from the canonical-tetrahedron covariance mapped to each
tetrahedron, (det/120)(a a^T + b b^T + c c^T + (a+b+c)(a+b+c)^T). The covariance
about the origin becomes an inertia tensor (trace(C) I - C), is shifted to the
centre of mass by the parallel-axis theorem, and diagonalized by the same
`finalize` M25 and M26 use, so a mesh drops into a `RigidBody3D` exactly as they
do.

Winding carries the sign: a mesh wound outward has positive volume, one wound
inward has negative volume and is rejected, along with a degenerate or
non-closed mesh (non-positive volume) and the usual malformed input.

### Validated against closed forms, then run through the dynamics

The anchor is exact: an axis-aligned box mesh reproduces the M26 `SolidBox` down
to the full tensor, which ties the general integral to a known result. A cube
placed far from the origin pins the centre of mass and the large parallel-axis
cancellation; a tetrahedron with one vertex at the origin pins the volume and
centroid while exercising the degenerate (zero-signed-volume) faces; a tilted
box mesh matches the M26 box turned the same way; a mesh result composes with
`combine`; and a box mesh, having three distinct principal moments, is run
through M20's dynamics as an intermediate-axis object that flips. Like M25 and
M26 this is a construction layer, not a domain (rule 16 does not apply): no
world, no scenario key, no template.

### Mutation testing

Eighteen mutations. Sixteen caught: the signed-volume triple product, the
volume, centroid and covariance constants, the covariance's outer-of-sum term,
the covariance-to-inertia conversion, the parallel-axis shift, the density
factor on both mass and inertia, the density and volume-positivity gates, the
face-count minimum, the triangle-index bounds check (whose removal faults the
out-of-range mesh under the checked runtime instead of returning zero), and
`trace` reading the wrong entries. Two are documented equivalents: dropping the
vertex-count minimum and dropping the per-vertex finiteness check, because a
body of fewer than four vertices or one with a non-finite vertex has no
positive, finite volume, so the volume gate rejects it either way.

## [M26] - 2026-09-11  (box mass properties and combine)

### Added

- `rigid::SolidBox` and `mass_properties_3d(const SolidBox&)`, the second mass
  primitive: a solid uniform box with three distinct principal moments,
  (1/3) m (h_j^2 + h_k^2) about each of its own axes, carried into the lab frame
  by its orientation.
- `rigid::combine`, which merges two mass-property sets into one compound body,
  so a body can be assembled from primitives of any kind (spheres, boxes, and
  whatever comes next) one at a time.
- `math::to_mat3`, a unit quaternion as its rotation matrix, the inverse of
  M25's `to_quat`. It is what lets `combine` rebuild a part's inertia tensor
  from the orientation the part stores.

### A second primitive, and assembly from primitives

M25 computed mass properties for a body made of solid spheres. A sphere is the
easy case: it is isotropic, so its inertia is one number and it has no
orientation. A box is the first primitive that is neither, so it is what forces
the general machinery to be right. Its three moments differ, and a tilted box
already has an inertia tensor that is not diagonal in the lab frame, which is
exactly the case M25's diagonalization was built for but could not yet exercise
from a single part.

`combine` is the assembly step ADR 0007 always implied. Masses add, the centre
of mass is their mass-weighted mean, and each part's inertia tensor is
reconstructed from its stored principal moments and orientation
(R diag(moments) R^T, which is where `to_mat3` comes in), shifted to the shared
centre by the parallel-axis theorem, summed, and diagonalized again. It is
commutative and associative up to the ordering of the eigenvalues, so the order
parts are folded in does not matter, and a zero-mass body is its identity, so a
fold over parts can start from nothing.

### Validated against closed forms, then pinned against the list path

The box moments are checked against (1/3) m (h_j^2 + h_k^2) with the half-extents
chosen so the principal frame is the identity and each moment lands on a known
axis; a tilted box then checks that the full lab-frame tensor, reconstructed
from the result, equals the one built independently from the closed-form moments
and the tilt, so a builder that ignored the orientation is caught.

`combine` is pinned against M25's list path, which never calls it: two spheres
combined one at a time must equal the same two built together in one list, which
ties combine's reconstruction and its parallel-axis shift to an independent
route. Commutativity, associativity and the zero-mass identity are checked
directly. A box already has three distinct principal moments, so it is run
through M20's dynamics as an intermediate-axis object: spun about its middle
axis it flips, spun about an extreme one it holds.

Like M25, this is a construction layer rather than a domain (rule 16 does not
apply): no world, no scenario key, no template. The M20 integration is its
demonstration.

### Mutation testing

Twenty mutations. Nineteen caught: the box mass coefficient, the (1/3) factor,
the moment-to-axis pairing, the orientation dropped from the box tensor and its
transpose removed; in combine, the mass sum, the mass-weighted centre, the
per-part reconstruction, the parallel-axis shift and its d d^T term, and both
zero-mass identity branches; the box validation of each half-extent, of
finiteness, density and unit orientation; and `to_mat3`'s columns and its sense
of rotation. The twentieth is a documented equivalent: flipping the sign of the
offset in the parallel-axis shift changes nothing, because m (|d|^2 I - d d^T)
is even in d.

## [M25] - 2026-09-11  (3D mass properties: inertia from geometry)

### Added

- `math::Mat3`, a concrete 3x3 matrix stored by columns, with the operations
  mass properties needs: matrix-vector and matrix-matrix products, transpose,
  the outer product, a symmetric eigensolver (`eigen_symmetric`, cyclic Jacobi),
  and a rotation-matrix-to-quaternion conversion (`to_quat`, Shepperd's method).
- `rigid::mass_properties_3d`, which computes a compound body's mass, centre of
  mass and principal moments of inertia from a set of solid uniform spheres and
  their densities, plus `rigid::rigid_body_from` to build a `RigidBody3D` from
  the result.

### Physics from geometry, not typed in

M20 through M24 took a body's mass and inertia as given. M25 computes them from
what the body is made of: mass is density times volume, the centre of mass is
the mass-weighted mean of the parts, and the inertia tensor about that centre is
each sphere's own (2/5) m r^2 shifted to the common centre by the parallel-axis
theorem, m (|d|^2 I - d d^T). This is the 3D sibling of the 2D `mass_properties`
and ADR 0007's step 2 (geometry, then mass properties, then integration), and it
is the piece a future generated shape will lean on to carry real physics.

A compound body's inertia tensor is generally not diagonal in the frame its
parts were given in, so M25 diagonalizes it: the eigenvalues are the principal
moments and the eigenvectors are the principal axes. That is exactly what
`RigidBody3D` stores, three moments and an orientation, so the result drops
straight in and the body never has to carry a full tensor. The general tensor
ADR 0007 and 0009 deferred has arrived as a construction step, not as a change
to how a body stores its inertia.

### The math is a support library, tested against closed forms

`malloy_math` gained a `Mat3` and a symmetric eigensolver, which are new kinds
of code for this project (linear algebra, not just vectors), so they are tested
on their own before anything rests on them: the eigensolver reconstructs its
input (`V diag Vᵀ`), returns a right-handed frame with ascending eigenvalues,
handles a repeated eigenvalue, and round-trips through `to_quat`.

The mass properties are then validated against derived closed forms:

- a single solid sphere: mass = density (4/3) pi r^3, isotropic moment
  (2/5) m r^2;
- a dumbbell (two spheres on an axis): the axial moment is the two spheres' own
  with no shift, the transverse pair each gains m d^2;
- a diagonal dumbbell: the smallest-moment principal axis points along the
  dumbbell, which pins the orientation the diagonalization returns;
- different densities: the centre of mass sits nearer the denser part.

The strongest test runs the whole chain through M20's dynamics: a cross of
spheres with three distinct principal moments, built by `mass_properties_3d`,
tumbles about its intermediate principal axis and holds about the extreme two,
exactly the intermediate-axis theorem. The derived moments, their ordering, the
orientation and `rigid_body_from` all have to be right for the flip to land on
the middle axis.

Like M9's collision primitives, this is a construction layer rather than a
domain, so rule 16 does not apply: no world, no scenario key, no template. The
M20 integration test is its demonstration.

### Mutation testing

Fifteen mutations, all caught: the outer product transposed, a dropped matrix
column, the eigensolver's sort reversed, its right-handedness fix removed, its
eigenvector accumulation dropped, `to_quat`'s branch inverted (caught by a
180-degree rotation, where the other branch divides by zero); and in the mass
properties, the sphere's (2/5) coefficient, the volume constant, the
mass-weighting of the centre, the parallel-axis sign, its |d|^2 term, the own
inertia not summed, the radius and density validation, and `rigid_body_from`
dropping the orientation. Several needed sharper tests than the first pass had:
a 180-degree rotation for `to_quat`, a valid-plus-invalid part pair for the
per-part validation (a single bad part is caught by the positive-mass guard
instead), and a non-identity orientation for `rigid_body_from` (the cross's
frame is the identity).

## [M24] - 2026-09-11  (sphere-against-sphere collisions)

### Added

- `collide::Sphere` against `collide::Sphere`: `overlaps` and `contact`, the 3D
  sibling of circle against circle, with the same +x fallback for coincident
  centres.
- Two-body contact response in `Rigid3DWorld`: a body-body pass in `step()`
  over every pair in a fixed order, alongside the existing ground pass.
- `scenarios/sphere_collision.scn`.

### One impulse core, two contact types

The contact response is now a single two-body core, `resolve_pair`, that takes
two participants and exchanges an equal and opposite impulse between them. The
ground is expressed as a participant with zero inverse mass, so a sphere against
a plane and a sphere against a sphere run the same formula; it is never written
twice. This is the 3D echo of the 2D `resolve_ground` stand-in body (ADR 0008,
"the second copy that becomes a third"), and refactoring M22/M23's
`resolve_ground3d` onto it left every existing sphere-plane result bit for bit
unchanged, which the M22 and M23 tests confirm.

### The headline is momentum

An immovable plane is a momentum sink: it absorbs whatever it must, so
sphere-plane contacts say nothing about momentum. Two movable bodies do, and the
sharp invariant is that TOTAL linear momentum is conserved, because the impulse
is equal and opposite. It is asserted to hold to rounding across a general
collision (unequal masses and radii, off-axis velocities, spin, and friction),
over four hundred steps.

The velocity outcomes are exact where they have closed forms:

- equal masses, elastic, head-on: the velocities are exchanged;
- equal masses, one at rest: the mover stops and the target leaves with the
  whole velocity (a Newton's cradle of two);
- unequal masses, elastic: the exact 1D elastic result
  `va' = ((ma-mb) va + 2 mb vb) / (ma+mb)`, and likewise for b;
- the relative normal velocity reverses by -e, off-axis.

Energy is conserved across a single elastic impulse (the head-on cases, exact),
and across a frictionless elastic multi-body collision to the scheme's contact
drift (~1e-7, the same class as M14's positional-correction artefact). With
friction or e < 1 it strictly falls. The positional correction splits by inverse
mass, so the heavier body of an overlapping pair moves less (a 1 : 3 mass pair
moves 3 : 1).

### Mutation testing

Fourteen mutations. Thirteen caught: the sphere-sphere normal (unnormalized and
reversed), penetration sign, overlap test, and coincident-centre fallback; both
contact arms (b's sign and a's radius, pinned by a glancing unequal collision's
induced spins); the radius guard on the second body; the impulse not reaching b,
both bodies taking the same sign, friction not reaching b; the positional
correction not splitting by mass; and the body-body pass dropped.

One is EQUIVALENT and uncatchable: removing the "two immovable things never
resolve" guard. A RigidBody3D always has positive finite mass, and there is no
plane-against-plane path, so the inverse-mass sum is always positive and the
guard never fires. It stays as defensive intent.

## [M23] - 2026-09-11  (Coulomb friction for 3D contacts)

### Added

- `friction` on `Rigid3DSettings` (default 0, frictionless), and a tangential
  friction impulse in the sphere-plane contact, clamped to the Coulomb cone.
- `type rigid3d` gains a `friction` key (shared spelling with the 2D `rigid`
  key, writing the 3D settings).
- `scenarios/rolling_sphere.scn`.

### The first 3D contact that imparts spin

M22 was rotationally trivial on purpose: a normal impulse on a centred sphere
passes through the centre of mass, r x n = 0, so it could not turn the body.
Friction is tangential, r x t is not zero, and it is the first contact here that
does work through the rotational effective-mass term ADR 0009 named,
`n . (I^-1 (r x n)) x r`, now put to use for the tangent. A sliding sphere is
dragged at its contact, slows, and spins up until the contact point is no longer
sliding: it rolls without slipping.

The response mirrors the 2D friction of M17: recompute the relative velocity
after the normal impulse, take the tangent along the residual slide, size the
impulse to arrest it, and clamp it to `friction` times the normal impulse
(Coulomb's cone). No new dynamics infrastructure: the world-frame inverse
inertia is the same body-frame bridge Euler's equations use for the torque
(M21), R I^-1 R^T, now applied to a contact impulse.

### The headline is a derived constant

A solid sphere sliding at v0 rolls without slipping at

```text
v_roll = v0 / (1 + I/(m R^2)) = (5/7) v0,
```

independent of the friction coefficient and of gravity, which set only how
quickly the slide phase ends. It is the 3D echo of the 2D rolling ratio (a disc
rolls at 2/3 of its sliding speed). The tests assert it two ways:

- Exactly, in one step, with a large coefficient that arrests the slide fully:
  v_roll = 5/7 v0 and the contact point comes to rest, to rounding, on a
  DIAGONAL slide so the spin axis is not a coordinate axis.
- Asymptotically under gravity, and INDEPENDENT of the coefficient: two
  different coefficients settle at the same 5/7 v0 to a part in a million, while
  each is far from the starting speed.

Also asserted: the Coulomb clamp (a small coefficient drops the speed by exactly
mu (1+e) |vz|, not the full 2/7 v0 a complete arrest would give); friction does
nothing in mid-air, since no normal impulse means no tangential one; and with
friction off the world is M22 bit for bit.

### The world-frame inverse inertia is rotated, and a test proves it

For a sphere the inertia is isotropic and I^-1 commutes with the orientation, so
the rotation into and out of the body frame is a no-op and could be dropped
unnoticed. It is not a no-op for a body whose principal moments differ. A test
tilts such a body 90 degrees so its third principal axis lies on the roll axis
and checks that it rolls at v0 / (1 + I3/(m R^2)), the I3 value, not the
body-frame I2 that ignoring the rotation would give. The two are far apart, so
the rotation is pinned.

### Mutation testing

Twenty-two mutations. Twenty caught on the first pass; a further one after adding
the tilted non-isotropic test above (the escape was the dropped body-frame
rotation, invisible to a sphere). The catches span the friction early-out, the
rotational term in the tangent effective mass, the angular-velocity update
(dropped and sign-flipped), the friction velocity update sign, the Coulomb clamp
(removed and mis-signed), the contact arm sign, the contact velocity ignoring
spin, the inverse inertia's conjugate rotation and its division by the moments,
and the friction validation.

One mutation is EQUIVALENT and cannot be caught: removing the "not sliding,
invent no tangent" guard. When the slide is exactly zero the tangent is 0/0, and
the next guard (tangential effective mass not positive) returns for the NaN
anyway, so the two produce identical behaviour. The guard stays as the clean
place to handle it rather than relying on that backstop.

## [M22] - 2026-09-11  (a sphere bouncing on a 3D ground plane)

### Added

- `collide::Sphere` and `collide::Plane3`, the 3D siblings of `Circle` and
  `Halfplane`, in new headers inside `malloy_collide` (2D untouched). Plus
  `overlaps` and a `contact(Sphere, Plane3)` returning a `Contact3` (normal,
  penetration, point). Like circle against halfplane, this pair has no
  degenerate case: the normal is the plane's own.
- A collision `radius` on `RigidBody3D`. Zero (the default) means the body does
  not collide, so every pre-M22 body is unchanged.
- `restitution`, a `gravity` acceleration, and a list of ground `Plane3`s on
  `Rigid3DSettings`, all defaulting to the pre-M22 behaviour (elastic, no
  gravity, no ground).
- `Rigid3DWorld::step()` gains a gravity kick and a contact-resolution pass; a
  `total_potential_energy3d` diagnostic.
- `type rigid3d` gains a `radius` field on `rigid_body3d` (18 fields now), and
  the keys `gravity3`, `plane3` and `restitution`.
- `scenarios/bouncing_sphere.scn`.

### The first 3D contact, and why it is the M10 echo, not M14

M20 gave a rigid body free rotation and M21 gave it a torque, but nothing it
could hit. M22 gives it a radius, drops it under a gravity acceleration, and
bounces it off an immovable plane. It is the smallest slice that reaches a new
capability, chosen the way every 3D milestone has been: it needs no inertia
tensor, no translational force, no solver, and no sphere-sphere geometry, so
none of those were built.

It reaches new physics without new dynamics infrastructure because a contact
here is a specialization, not the genuinely-new tier ADR 0009 reserved for the
tumble. A uniform sphere's inertia is isotropic, and a CENTRED sphere's contact
point lies on the line from its centre straight to the plane, so the arm from
the centre of mass to the contact is parallel to the normal and r x n = 0. A
normal impulse therefore imparts no spin and never touches the inertia: the
effective mass is just 1/m against the immovable plane. That makes M22 the 3D
echo of M10's colliding particles, which were also translational, rather than
M14's rigid contacts, which needed the disc centred off the centre of mass to
get a moment arm. Spin from a contact needs a tangential (friction) component,
r x t != 0, and that is a later milestone.

The ground plane is stored as geometry in the settings, not as a body, so
`RigidBody3D` never has to represent an infinite mass and its validity still
rejects one. Gravity is an acceleration applied before the position update (the
M15 pattern, rule 5), not a force, and the contact is an impulse, so no
force or torque accumulator was added (rules 5, 6, 17).

### Derived, then measured

The tests assert derived constants, not loosened tolerances:

- Velocity restitution, EXACT: across the impulse the normal component of the
  velocity reverses to -e times itself, to rounding, on a tilted plane whose
  normal has all three components. The tangential component is untouched (no
  friction) and the spin is untouched (no lever arm), both asserted.
- Energy across a straight-down bounce: the kinetic energy scales by exactly
  e^2, since the speed scales by e.
- Rebound height e^2 h for a drop from rest, which follows from the velocity
  law and is shown in the template.
- Free-flight energy shed = (1/2)(sum m)|g|^2 dt^2 per step, the SAME derived
  constant as M12 and M15, now in 3D and with a non-axis-aligned g. An
  equality, asserted every step, matched to about 1e-13.
- The contact is rotationally INERT, which is the milestone's scope claim
  written as a test: a spinning body under a torque, dropped onto a floor, has
  bit-identical orientation and angular velocity to the same body with no floor.
  The bounce changes where it is, never how it spins.

### Mutation testing

Twenty-one mutations, all caught: the sphere-plane normal not negated, the
depth sign flipped, the overlap test inverted, the plane's unit-normal check
and the sphere's radius check dropped, the signed distance flipped; the contact
response's closing-speed sign, its restitution factor, its impulse direction,
its positional correction (dropped and reversed), its radius guard and its
separating guard; the gravity kick dropped and flipped; the whole contact pass
dropped; the settings validation weakened on restitution, gravity and ground
planes; the body validity ignoring the radius; and the potential energy sign.

One first appeared to escape: removing the "radius 0 does not collide" guard
changed no test, because no test dropped a zero-radius body onto a plane. That
is a real gap, now closed with a test that a zero-radius body falls straight
through a plane, bit-identical to the same body with no ground.

### Fixed

- A latent brace malformation in `tests/scenario/scenario_tests.cpp`, left by
  the M21 insertion: a stray `{` wrapped the whole torque section and was
  balanced only by the following block borrowing its close. It compiled, but
  the next inserted block would have unbalanced it. Removed while adding the
  M22 parse tests.

## [M21] - 2026-09-11  (a constant applied torque on a 3D body)

### Added

- `rigid::Rigid3DSettings`, a settings type for the 3D rigid world, holding a
  single constant `torque`. Valid when the torque is SQUARABLE, the same bound
  positions, velocities and M18's softening carry, because it enters |u|^2 in
  the drift law.
- `Rigid3DWorld` gains that settings argument (defaulted, so a world built
  without it is the torque-free M20 world bit for bit) and applies the torque
  each step through Euler's equations.
- `type rigid3d` gains a `torque <tx> <ty> <tz>` key: one world-frame torque for
  the whole world, the way the 2D rigid domain has one gravity.
- `scenarios/gyroscope.scn`, a fast top under a torque across its spin.

### What torque adds, and why it is the smallest next step

M20 was torque-free. A body spun and tumbled, but nothing acted on the
rotation, so a body on a principal axis stayed there forever. M21 makes
something act on it, and does so the way M15 made gravity act on a 2D body: a
constant field in the settings, applied as forcing, NOT a force/torque
accumulator (rule 5, rule 17). No contacts, no collision geometry, no
translational force and no solver had to be built to reach it.

The physics it reaches is the one a torque does that has no 2D analogue: it
moves the spin AXIS, not just the spin rate. Push a still body and it moves the
way you push; push a fast gyroscope and the axis travels at right angles to the
push. In two dimensions a torque can only speed a spin up or slow it down,
because there is only one rotation axis.

### The torque is world-frame, and that is what makes the invariant clean

The setting is a torque in the WORLD frame, an external couple fixed in the lab.
Euler's equations are diagonal only in the body frame, so it is rotated in by
the body's current orientation, T_body = R^-1 T_world, before the forward-Euler
angular step. A body-fixed torque (a thruster bolted to the body, constant in
the body frame) is a different thing and is deferred; it is not reachable from
this one without knowing the body.

World-frame is the choice that gives the headline invariant. The physical law
dL/dt = torque holds in the WORLD frame whatever the body does, so under a
constant torque the world-frame angular momentum grows along a straight line:

```text
L(t) = L(0) + T_world t
```

exactly in the continuum, and tracked to first order by the scheme. The
gyroscope template reads that off a single column: L starts at (0, 0, 15), the
torque is (0.4, 0, 0), and L_x grows as 0.4 t (1.2, 2.4, 3.6, 4.8 at the four
reports) while L_z holds at 15, so |L| climbs from 15 as the axis tilts.

### One case is exact, not merely first order

A body at rest, unrotated, with the torque along a principal axis, spins up
EXACTLY. The body only ever rotates about that axis, a rotation about an axis
fixes it, so the body-frame torque stays equal to the world torque bit for bit;
the gyroscopic term is zero because omega stays on the axis; and what remains is
I wx' = T, which forward Euler integrates without error because the right side
is constant. So

```text
wx(n) = n dt T / Ix,   wy = wz = 0,   L_world = (T t, 0, 0)
```

to rounding, with the transverse components not merely small but zero and the
axis not nearly fixed but fixed. That is the test's anchor: an exact invariant,
not a tolerance.

### The forced discrete laws, derived and asserted every step

With L = I omega and u = L x omega + T_body, the forward-Euler step
omega' = omega + dt I^-1 u gives, exactly,

```text
|L'|^2 = |L|^2 + 2 dt (L . T_body) + dt^2 |u|^2
T'     = T     +     dt (omega . T_body) + dt^2 (I^-1 u) . u / 2
```

The cross term is orthogonal to both L and omega, so it drops out of the
first-order part and only the torque survives there. With T_body = 0 these are
exactly the M20 laws, which is why the torque-free world is recovered bit for
bit. Both are checked every step in the tests against a torque and a tilted,
tumbling body, to about 1e-13. Note L . T_body = L_world . T_world and
omega . T_body = omega_world . T_world: the work a torque does is frame
independent.

### Mutation testing

Every torque path caught: the torque dropped from Euler's equations; applied in
the world frame without rotating to the body; rotated by the orientation
instead of its conjugate (backwards); its components swapped; negated; the
settings validation weakened to accept a non-squarable or a merely-finite
torque; the world's validate no longer checking the settings; the parser
reading two components or dropping the third; and two mutations of the
template's own `# check` values.

A methodology note: one header-only mutation (weakening the inline
`Rigid3DSettings::is_valid`) first appeared to escape, but that was a stale test
binary, not a gap. Changing an inline function in a header did not trigger a
rebuild of a test translation unit that only includes it, so the harness ran the
old binary. Rebuilding the test unit, the mutation is caught. The mutation
harness now touches the test units before building so a header-only change
propagates.

## [M20] - 2026-09-10  (quaternions and torque-free rotation in three dimensions)

### Added

- `math::Quat`, a quaternion, used only ever as a UNIT quaternion representing
  an orientation. Hamilton product, conjugate, norm, `normalize`, `rotate`,
  `from_axis_angle`, and an `is_unit` that is deliberately stricter than
  `is_finite`.
- `rigid::RigidBody3D`: position, velocity, a quaternion orientation, a
  BODY-frame angular velocity, mass, and three PRINCIPAL moments of inertia.
- `rigid::Rigid3DWorld`, torque-free rotation, inside `malloy_rigid` rather
  than a library of its own, on the same reasoning M19 used for 3D gravity:
  rigid-body dynamics is one domain and the dimension is not a domain.
- `rigid::angular_momentum` and `rigid::rotational_energy` for a single body,
  plus `total_linear_momentum3d`, `total_angular_momentum3d` and
  `total_kinetic_energy3d` for a world.
- `type rigid3d` in the scenario format, with a seventeen-field
  `rigid_body3d` line. The orientation is given as an axis and an angle rather
  than as four quaternion components, because a hand-written quaternion is
  almost never a unit one and `from_axis_angle` cannot produce anything else.
- `scenarios/intermediate_axis.scn`, and a seventh branch of the app's dispatch
  switch.

### What is actually new here

This is the first 3D milestone whose physics has no 2D form at all. M19 widened
gravity: the same equations with one more component. Rotation does not widen.

In two dimensions the inertia is a scalar and the angular velocity lies along a
fixed axis, so a body's inertia never changes in world space and the term
`omega x (I omega)` is identically zero. A free body spins at a constant rate
forever, and there is nothing else it could do. In three dimensions Euler's
equations carry that term:

```text
Ix wx' = (Iy - Iz) wy wz
Iy wy' = (Iz - Ix) wz wx
Iz wz' = (Ix - Iy) wx wy
```

Every right-hand side is a DIFFERENCE of principal moments, which is why a
sphere still cannot do anything but spin and a body with three distinct moments
can tumble under no torque at all.

The milestone is built around the intermediate axis theorem. With moments
(1, 2, 3) and an identical one-part-in-two-thousand nudge, a body spun about
the largest or the smallest moment holds that nudge to within a part in a
hundred over twenty thousand steps, and a body spun about the one in between
turns completely over. The same code, the same perturbation, the only
difference is the axis.

### Scope, and what was deliberately not built

M20 scoped the way M19 did, by choosing the slice that needs nothing else
first. Torque-free rotation has no contacts, no collision geometry, no forces
and no solver, so none of those had to be built ahead of their milestone.

`RigidBody3D` stores three principal moments and NOT a 3x3 matrix. That is a
choice of axes rather than a restriction: a symmetric tensor is always
diagonalizable, so every rigid body has such a frame. A general tensor is
needed once bodies are built from composed shapes and the parallel-axis step
moves inertia off the principal axes, which arrives with 3D mass properties and
3D contacts. Until then it would be three extra zeros (rule 11). ADR 0007's
scope boundary is amended with what came due and what was deferred again.

The angular velocity is stored in the BODY frame, which is what makes Euler's
equations diagonal. In the world frame the inertia rotates with the body and
the equations would need a full tensor at every step.

### What conserves, and what does not

Neither angular momentum nor energy is exactly conserved by the discrete
scheme, and the size of the difference is an equality rather than a bound.
Writing `u = L x omega`, so the update is `omega' = omega + dt I^-1 u`, both
`L . u` and `omega . u` vanish identically. The first-order terms cancel and
exactly the second-order ones survive:

```text
|L'|^2 = |L|^2 + dt^2 |u|^2
T'     = T     + dt^2 (I^-1 u) . u / 2
```

Both corrections are positive, so this is a slow gain and not a loss. Both are
checked every step in the tests, and they explain the cases that ARE exact:
`u` is zero exactly when `I omega` is parallel to `omega`, that is when omega
lies in an eigenspace of the inertia. For three distinct moments that means a
principal axis; for a sphere any direction; for an axisymmetric body any
direction in the degenerate plane. In those configurations the angular
velocity, the world-frame spin angular momentum and the rotational energy are
bit-identical after five thousand steps.

`total_angular_momentum3d` is the SPIN of each body plus its orbital
`m (r x v)`, matching the 2D diagnostic; the orbital term is exactly constant
since no force acts. Its magnitude is a vector sum and so does not itself
inherit the per-body `|L|^2` law, unlike the kinetic energy, which is a scalar
sum and does. The drift itself comes from the FORWARD-Euler angular velocity
step; "semi-implicit" here names only the velocity-before-orientation coupling,
as in every other world.

That also makes the template's diagnostics readable rather than mysterious. The
E and |L| columns sit still, jump, and sit still again, and the jump is exactly
where the flip is: the drift term is zero for a principal-axis spin and large
only while the body is turning over.

Linear momentum IS exact, since nothing here applies a force.

### Renormalization costs nothing

Integrating a unit quaternion does not keep it unit, and the growth is exact
rather than approximate, because the quaternion norm is multiplicative and the
update is a right multiplication:

```text
|q'| = |q| |(1, (dt/2) omega)| = |q| sqrt(1 + |omega|^2 dt^2 / 4)
```

Unlike the positional correction in the 2D contact code, which this was first
documented as resembling, renormalizing is not a projection that costs
anything. Scaling a quaternion does not change the rotation it represents, so
it moves no body and perturbs no conserved quantity. It restores the magnitude
`rotate` assumes, which would otherwise scale every rotated vector by |q|^2,
and it stops the factor above from compounding.

### Derived, then measured

Following the practice the rest of the project uses, the tests assert derived
constants rather than loosened tolerances:

- the discrete growth rate about the intermediate axis. The transverse map has
  eigenvalues `1 +/- (W/sqrt(3)) dt`, so from a nudge of `e` the growing
  component is exactly `(e/2)[(1 + r dt)^n + (1 - r dt)^n]`. Measured against
  that to fifteen significant figures.
- the bound the two stable axes respect. Explicit Euler on a rotation grows the
  amplitude by exactly `sqrt(1 + (mu dt)^2)` per step, the same polygon factor
  as M18's cyclotron, with `mu = W/sqrt(3)` about the smallest moment and
  `mu = W` about the largest. Two different numbers, so this checks the
  coefficients rather than repeating one case twice.
- the peak angular speed during the tumble. Starting on the separatrix fixes
  `|L|^2 = (I2 W)^2` and `2T = I2 W^2`; where the body passes through wy = 0
  those give `wz^2 = W^2/3` and `wx^2 = W^2`, so `|omega| = 2W/sqrt(3)`.
  Measured 2.3107 against 2.3094, the gap being the accumulated per-step drift
  above.
- the symmetric top's precession, which is the cyclotron polygon again: `wz` is
  bit-constant because `(Ix - Iy)` is exactly zero, and the transverse pair
  turns by exactly `atan(Omega dt)` per step, never by `Omega dt`.
- the orientation after a principal-axis spin, which is that polygon once more,
  a step at a time in the (w, axis) plane at half the angular rate.

Two properties are separated by measuring them at two step sizes. Under
torque-free rotation the WORLD-frame angular momentum is fixed while the
BODY-frame vector `I omega` tumbles with the body. The former is integration
error and falls with dt, the latter is the physics and does not: at dt = 1e-3
the world-frame deviation is 2.925e-3 and at dt = 1e-4 it is 2.923e-4, a ratio
of 10.006 that pins the scheme's first order, while the body-frame swing is
5.10 at both. That ratio is asserted, because "small" alone would also be
satisfied by an implementation that never rotated anything.

### Mutation testing

Twenty-two mutations, all caught: each Euler coefficient's sign, a difference
turned into a sum, a transposed component, a division by the wrong moment, the
renormalization removed, the quaternion multiplied on the wrong side, the
half-step factor dropped, angular momentum left in the body frame, the missing
half in the energy, the translation dropped, three validity checks weakened,
and six mutations of the quaternion algebra itself.

One escaped on the first pass: using the stale angular velocity for the
orientation, which makes the scheme fully explicit instead of semi-implicit.
The ordering was documented and untested, and it is otherwise invisible,
because swapping the two changes the result only at O(dt), the same order as
the scheme's own error, so no convergence test can see it. It is visible as a
DIRECTION: starting from the identity the quaternion increment points along
whichever omega was used, and the two differ by `rate * dt`. That is now
asserted, and the mutation is caught.

Eight further mutations cover the new scenario key and the template's own
`# check` lines.

### Fixed

- Three documentation claims that were stale before this milestone started, all
  of them now derived from the repository instead of restated.

  README.md still said 3D "has not started" and that everything shipped was 2D.
  That had been wrong since M19.

  README.md claimed the templates covered "five domains" when M19 had made it
  six. Only the template COUNT was guarded; the domain count was prose, so it
  drifted silently. The scenario tests now count distinct scenario types across
  the shipped templates and check that number too.

  docs/08_AI_HANDOFF_PROMPT.md claimed "ADRs 0001-0008" after ADR 0009 was
  written. The guard for exactly this claim existed but read only docs/00, so
  fixing the guarded copy let the unguarded one drift through two milestones.
  The guard now reads both.

- The M20 headers claimed, while it was being built, that angular momentum and
  energy were exactly conserved and that renormalization perturbed world-frame
  angular momentum. The tests disproved all three before the milestone shipped.
  The derivations above replace them.

## [M19] - 2026-09-10  (Vec3 and gravity in three dimensions)

### Added

- `math::Vec3`, a concrete 3D vector beside `Vec2`, with a `cross` that returns
  a VECTOR.
- `nbody::Body3D` and `nbody::NBody3DWorld`, Newtonian gravity in three
  dimensions, inside `malloy_nbody` rather than a library of its own.
- `nbody::total_angular_momentum(const std::vector<Body3D>&)`, which returns a
  Vec3 where the 2D world returns a scalar.
- `type nbody3d` in the scenario format, with `body3 <mass> <px> <py> <pz> <vx>
  <vy> <vz>`. The `g` and `softening` keys now serve both gravity worlds,
  because they mean the same thing in either dimension.
- `scenarios/inclined_orbits.scn`.
- An orthographic projection in the terminal app, so a 3D scenario still draws.

### The three questions ADR 0009 left open, answered

That ADR deliberately deferred three decisions until a concrete case existed.
Building M19 settled all three, and the amendment there records them.

**Separate types, not templates.** `Vec3` is its own type. The two are not one
algebra with a different component count: `cross` is a SCALAR in 2D and a
VECTOR in 3D, so a template would need specializing for the one operation that
matters most to rigid-body dynamics, and `perp` and rotation by a scalar angle
have no 3D form at all. A template would have unified the component-wise
arithmetic, which is the easy half.

**3D worlds sit inside their domain's library.** Gravity is one domain; the
dimension is not a domain. Splitting them would put the same equations in two
libraries and require keeping them in step by hand. The dispatch switch grew by
one branch, exactly what ADR 0006 says a new world should cost.

**2D stays supported**, which follows from the above rather than being a
separate decision. The 2D types are untouched and all twelve templates run.

### Architecture Notes

- This is the smallest honest step into 3D, and that was the point. N-body is
  the only domain with no contacts and no orientation, so it needed `Vec3` and
  nothing else. Quaternions and inertia tensors stay deferred rather than being
  built speculatively for a user that does not exist yet (rule 11).
- `malloy_ascii` does not carry over, as ADR 0009 said it would not. The app
  projects orthographically along z and nothing cleverer: motion along z is
  invisible, a tilted circular orbit reads as an ellipse, and two bodies at
  different depths can overlap. Making that a camera would turn this into a
  graphics project, which rule 2 prevents until rendering has its milestone.
  The template says so rather than leaving a reader to wonder.
- `malloy_sim_core` is unchanged, and the layering test confirms the new code
  adds no dependency: `malloy_nbody` still links math, sim_core and time.

### What conserves, and what does not

Momentum and angular momentum are both EXACT for pairwise central forces, to
all orders in dt, because the (i,j) and (j,i) torque contributions are
`x_i x x_j` and `x_j x x_i` and cancel identically. Softening changes only the
magnitude and does not break it. The only error possible is rounding.

Angular momentum being a vector is what lets three dimensions state something
two cannot: because its DIRECTION is fixed, every orbit is confined forever to
the plane it started in. Two orbits in different planes stay in different
planes. `scenarios/inclined_orbits.scn` is built around that, and the test
asserts it directly by checking an orbiting body's position stays perpendicular
to L at every step.

Energy is not conserved and is not meant to be. It oscillates about a nearby
quantity with bounded amplitude and no secular drift, which is what a symplectic
method gives, and the template's E column is left visibly wandering in the ninth
digit with a note saying so.

### Tests

- `Vec3` arithmetic with every component distinct, so a dropped or transposed
  one shows, and a 2-3-6 triple so lengths are exact rather than nearly so.
- The cross product: right-handed on all three axis pairs, anticommutative,
  zero on a parallel pair, perpendicular to both inputs, satisfying
  `|a x b|^2 + (a . b)^2 = |a|^2 |b|^2`, and one hand-computed value off every
  axis. The handedness assertions matter because a sign error there is
  invisible to every magnitude (docs/05).
- Every 3D configuration is deliberately NOT confined to a coordinate plane,
  and the angular momentum in the conservation test has all three components
  nonzero, checked before the conservation is. A planar configuration would
  conserve two components trivially by both being zero.
- The orbital-plane invariant, asserted at every step of a tilted orbit.
- The strongest check available: a 2D configuration embedded at z = 0 runs in
  both worlds and they agree BIT FOR BIT, with the third component staying
  exactly zero. That puts the new code against an implementation tested since
  M4 rather than against a fresh derivation that could share a mistake with it.
- The softening contract, pinned at the same number the 2D world uses because
  it is the same denominator.
- Verified by mutation testing. Nine mutations were applied and all nine were
  caught: two different wrong cross products, a dot product missing its z term,
  `is_squarable` weakened to `is_finite`, a pair force applied to one member
  only, an unsquared softening, the integrator turned into explicit Euler, the
  mass dropped from angular momentum, and the body validation weakened.

## [M18] - 2026-09-10  (charged particles)

### Added

- `malloy_charges`, a fifth domain: point charges in electric and magnetic
  fields.
- `ChargedParticle2D`, whose charge is SIGNED and may be zero. No other body
  model in this project has a quantity that decides whether an interaction
  attracts or repels.
- `ChargeSettings`: Coulomb constant, uniform electric field, uniform magnetic
  field, and softening in the same sense as the N-body domain.
- `ChargeWorld`, with `total_momentum`, `total_kinetic_energy`,
  `total_potential_energy` and `total_energy`.
- `type charges` in the scenario format, with `coulomb`, `efield`, `bfield`,
  `charge <mass> <q> <px> <py> <vx> <vy>`, and `softening`, which now serves
  two domains for the same reason: both have a 1/r^2 pair term that goes to
  infinity for a coincident pair.
- `scenarios/cyclotron.scn`.

### Architecture Notes

- This is the first domain in the project whose force depends on VELOCITY.
  Everything before it was position-dependent (gravity, Coulomb, Hooke) or an
  instantaneous impulse.
- The magnetic term is applied as a ROTATION, not as another kick, and that is
  a physical statement rather than an optimisation. A magnetic force is always
  perpendicular to the velocity, so it does no work and cannot change a
  particle's speed. Applying it the way every other force here is applied,
  `v += (q/m)(v x B) dt`, adds a vector at right angles to v, so the two are
  the legs of a right triangle and the speed becomes `|v| sqrt(1 + (q b dt/m)^2)`
  every step, without bound. A cyclotron orbit integrated that way spirals
  outward: at the template's parameters it grows 22 per cent in ten turns.

  In two dimensions the magnetic sub-problem has an exact solution, because a
  uniform out-of-plane field rotates the velocity at a constant rate and does
  nothing else. So it is integrated exactly and the speed is preserved to the
  last bit, which is stronger than the Boris push manages: Boris preserves the
  speed exactly but carries an O(dt^2) frequency error, and an exact rotation
  has neither.

  The cost is operator splitting. The electric and magnetic parts are applied
  in sequence rather than together, which is first order in dt in their
  interaction, exactly as the single kick in every other domain is.
- The magnetic field is a SCALAR, and that is a consequence of two dimensions
  rather than a simplification: only the out-of-plane component of B produces
  an in-plane force. In 3D it becomes a vector again and the force stops being
  expressible this way (ADR 0009).
- `malloy_sim_core` is unchanged, and the dispatch switch grew by exactly one
  branch, which is the whole cost ADR 0006 said a new domain should have.

### The rule-17 trigger fired, and the answer is still no

ADR 0008 said to reassess extracting a shared translational integrator once a
third domain independently needed the same path. There are now five. The
reassessment is recorded in that ADR, and M18 turned out to be evidence AGAINST
extraction rather than for it: `charges` puts a rotation between the kick and
the move, so it does not share the path at all. What all five have in common is
one line, `x += v * dt`, and extracting one line would create a coupling point
that four domains would use differently and the fifth would not fit.

### What conserves, and what does not

Under a magnetic field alone, the speed of every particle is constant and so is
the total energy. Not to some order in dt: exactly, up to the rounding of
`cos^2 + sin^2`, which is about 2 ulp per step. Over 10000 steps that bounds the
drift in `|v|^2` at 4.4e-12 relative, and the measured figure is -5.6e-13.

The pairwise Coulomb interaction conserves momentum exactly, because both
members of a pair are written in the same visit. A uniform field does not, and
should not: it is external.

Two discretisation results, both derived and then measured rather than bounded
by a tolerance:

- The discrete cyclotron orbit is a regular polygon, not a circle. Its
  circumradius exceeds the true radius `m|v|/(|q| b)` by `(theta/2)/sin(theta/2)`,
  about `1 + theta^2/24`. Measured by the separation of two vertices half a turn
  apart, which is exactly the diameter wherever the centre happens to be.
- The discrete E x B drift is the continuum drift `|E|/b` ROTATED by half a
  step's angle and inflated by the same polygon factor. Both corrections carry
  `theta = -q b dt/m`, so the continuum independence from charge and mass
  survives only to order dt^2. The two particles in the test gyrate at rates
  differing by a factor of three and their drifts differ in the seventh digit.

Worth noting what the continuum guiding centre is NOT: `position + (m/qb)(v.y, -v.x)`
is not the polygon's centre, it is off by about `|v| dt / 2`. Measuring the
radius from it gives 2.99994 rather than 3.0000049, which looks like a failing
radius and is really a wrong centre.

### Tests

- Validation, including that a charge of zero and a negative charge are both
  legal while a mass of zero is not.
- The speed is unchanged over 10000 steps, to a bound derived from the rounding
  of `cos^2 + sin^2` rather than picked.
- The orbit closes: after exactly one turn, chosen as a whole number of steps,
  both position and velocity return to their starting values.
- The period does not depend on speed: two particles with the same charge-to-
  mass ratio and a ninefold speed difference return together.
- The radius does depend on mass and charge, checked at three combinations
  against the derived polygon circumradius.
- Handedness: a positive charge in a positive field turns clockwise and a
  negative one turns the other way, mirror images to the last bit. Invisible to
  both speed and radius, so it needs its own assertion (docs/05).
- E x B drift, asserted per particle against the derived discrete value, with
  the continuum claim stated separately and correctly qualified.
- Coulomb is signed: like charges repel, unlike attract, a neutral particle
  does neither, and the potential energy of a like pair is POSITIVE, which
  gravity's never is.
- A uniform electric field is not gravity: three particles with different
  charge-to-mass ratios accelerate by different amounts and in different
  directions.
- Softening enters the denominator squared, pinned at the same number the
  N-body domain uses because it is the same contract.
- Verified by mutation testing. Ten mutations were applied and nine were caught
  at once: the turn applied as a kick, its sense flipped, the mass dropped from
  the rotation rate, Coulomb's sign removed, the pair force applied to one
  member only, the electric field treated as an acceleration, both potential
  signs, and validation dropped. The tenth, un-squaring the softening, escaped
  and exposed a real gap: the softening test asserted only that the result was
  finite, never what it was. Now pinned, and now caught.

## [M17] - 2026-09-10  (Coulomb friction)

### Added

- Coulomb friction for rigid contacts. After the normal impulse, a tangential
  impulse is applied along the slide, clamped to `friction` times the normal
  impulse.
- `RigidSettings::friction`, defaulting to 0. Validated as non-negative and
  finite, and deliberately NOT capped at 1: a coefficient above 1 is physically
  real, and clamping it would silently change a caller's model.
- `friction <mu>` in the scenario format, for `type rigid`.
- `scenarios/rolling_and_slipping.scn`.

### Architecture Notes

- Friction is a contact impulse, not a persistent force. It is computed and
  consumed inside a contact exactly as restitution is, so `malloy_rigid` still
  has no force or torque accumulators and rule 5 is untouched. It also adds no
  state to `RigidBody2D`.
- The clamp is a MINIMUM of two quantities: what it would take to stop the slide
  outright, and the cone limit. Taking the cone limit unconditionally is the
  classic way friction ends up adding energy, because it overshoots a slide that
  was about to stop and reverses it.
- Friction is gated on the normal impulse by construction, since the limit is a
  multiple of it. No normal impulse means no friction, automatically, so a body
  in mid-air cannot be dragged sideways however large the coefficient.
- There is no fallback direction for a zero tangent, and there must not be.
  Inventing one would push a body that is not sliding, which is how a resting
  arrangement drifts.
- This milestone depends on M16. The invariant below assumes a contact arm of
  constant length R, which a floor built from overlapping discs does not
  provide: its normal swings by up to 14.48 degrees, so there is no constant arm
  and the closed-form result does not hold.

### What conserves, and what does not

A tangential impulse at the contact point changes velocity and spin together in
a fixed ratio, so it leaves `I*omega - m*R*v` unchanged however large the
impulse is and whether or not the cone clamped it. That single fact gives every
result in the template in closed form.

For a uniform disc launched at `v0` with no spin, the constant is `-m*R*v0`, and
rolling without slipping means `omega = -v/R`, which makes it `-(3/2)*m*R*v`.
So:

    v_roll = (2/3) * v0

and exactly one third of the launch kinetic energy is gone. Neither figure
depends on the friction coefficient or on gravity: those decide only how long
the slipping phase lasts, never what it ends at. The test asserts it across
three different combinations of the two for that reason, because a result that
cannot be tuned is worth more than one that can.

Once rolling, friction switches itself off: there is no tangential relative
velocity left, so no impulse is applied, and rolling dissipates nothing. In the
template, E, L and momentum are identical to every printed digit from step 3000
to step 6000.

The result is not exact in the code, and the residual is derived rather than
tuned. The contact POINT is reported midway through the overlap, so the arm is
`R - depth/2` rather than `R` while the invariant assumes `R`. In steady contact
the penetration is `g*dt^2/(1 + restitution)`, which predicts a rolling residual
of `-v * (depth/2) / R`. For the headline test that is
`-2.0 * (9.81 * 0.0005^2 / 2) / 0.5 = -4.905e-6`.

NOTE, added later: as shipped in M17 the code did NOT do this. A statement
ordering bug made the effective arm `R + depth/2`, so the residual was
`+4.905e-6` and this paragraph described a behaviour the code did not have. The
magnitude was right, which is why nothing caught it: the assertion was a
magnitude bound. Both the code and the claim are corrected in the Unreleased
section above.

Kinetic energy is strictly non-increasing across every contact, at every
coefficient, which is asserted at every step of a long run rather than at the
end.

### Tests

- Validation: negative and non-finite coefficients refused, above 1 accepted.
- Zero friction reproduces the pre-M17 trajectory exactly.
- The 2/3 rolling ratio at three different (friction, gravity) pairs, with the
  rolling condition checked physically as "the contact point is not moving"
  rather than through a sign convention.
- Exactly one third of the launch energy gone.
- Backspin reversal: a disc launched forward with backspin rolls out, stops and
  comes back. Nothing in the project could do this before, since it needs spin
  converted into translation.
- Spin-down: a body dropped with spin and no translation loses the spin and
  gains motion from it.
- Static holding on a slope, using a body of infinite inertia so it slides
  rather than rolls: it holds when the coefficient reaches `tan(theta)` and
  otherwise accelerates at exactly `g(sin - mu cos)`. That body is only
  expressible because infinity is per quantity, so infinite inertia means
  "cannot spin" rather than "immovable".
- Friction cannot act in mid-air, however large the coefficient.
- Verified by mutation testing. Nine mutations were applied and seven were
  caught immediately: the impulse sign flipped, the cone removed, the cone
  applied unconditionally rather than as a minimum, the gating on the normal
  impulse removed, the rotational terms dropped from the tangential effective
  mass, the friction torque dropped, and validation dropped.

  One escaped and exposed a real gap. The step contract states that friction is
  computed from the velocity REMAINING after the normal impulse, and that this
  is observable, but nothing tested it: on a flat floor the normal impulse does
  not change the tangential velocity, so both orderings agree. The case that
  distinguishes them is a body whose centre of mass is offset from its disc,
  dropped straight down with no spin. Its tangential velocity before the normal
  impulse is exactly zero, and the impulse then spins it, which is what gives
  friction something to resist. Now tested, and now caught.

  The ninth is an equivalent mutant, and it was measured rather than argued.
  Removing the explicit zero-slide guard changes nothing, because a zero slide
  makes the tangent NaN, which makes the tangential effective mass NaN, which
  fails the negated `!(effective > 0)` check and returns. Confirmed by running
  all ten templates with the guard removed and diffing: byte-identical. The
  guard is kept because relying on NaN propagation reaching a negated
  comparison is fragile, not because the code needs it today.

### Not in this milestone

`malloy_particles` has no friction; particle contacts remain normal-only, and
`scenarios/projectile_arc.scn` still documents that correctly. The solver is a
single pass with no iteration, so a stack of bodies will not stand. That is
where sequential-impulse iteration and then a constraint solver start pulling,
and a constraint solver is where a `Constraint` base class grows, which rule 12
exists to prevent. Named here and deferred rather than left to drift in.

## [M16] - 2026-09-10  (halfplanes: true flat ground)

### Added

- `collide::Halfplane`: everything on one side of an infinite straight line,
  stored as a unit normal and an offset. The normal points OUT of the solid, so
  the signed distance `dot(normal, p) - offset` is positive in free space.
- `overlaps(Circle, Halfplane)` and `contact(Circle, Halfplane)`.
- Ground planes in `RigidSettings`, resolved after body-against-body contacts
  in ascending body then plane order.
- `ground <nx> <ny> <offset>` in the scenario format, for `type rigid`, and it
  may appear more than once.
- `scenarios/flat_ground.scn`: a ramp, a floor and a wall.

### Changed

- The rigid impulse formula moved into `apply_contact`, called by both
  `resolve_contact` (disc against disc) and the new `resolve_ground`. Writing
  it a second time for the plane path is exactly the second copy that becomes a
  third (ADR 0008). No behaviour changed: the existing suite passed unaltered
  across the extraction.

### Architecture Notes

- The reason this exists is the contact NORMAL, not the surface height. A floor
  built from overlapping discs has a normal that points at whichever disc centre
  is nearest, so it SWINGS as a body moves along it: up to 14.48 degrees for
  the radius-0.4 bodies in `scenarios/dropped_bodies.scn`, which documents its
  height ripple as about 0.05 and thereby understates the problem considerably.
  (That 0.05 is the FLOOR's own surface ripple. The locus a radius-0.4 body's
  centre actually rides ripples by 0.0318. The two were quoted together as
  though they described one surface, which is corrected in Unreleased above.) Getting that
  under one degree with discs needs roughly 160 of them. A plane has one normal
  everywhere.
- Circle against halfplane is the only pair in `malloy_collide` with NO
  degenerate case. Every other query has to infer a direction from two centres
  and falls back to +x when they coincide (`contact.hpp`). The plane supplies
  the direction itself, so even a circle whose centre lies exactly on the line
  gets the correct normal rather than a documented arbitrary one.
- The plane stands in as a body of infinite mass and inertia, which is what
  keeps one copy of the impulse formula. Its inverse mass and inverse inertia
  are exactly zero, so it absorbs the impulse without moving, takes no share of
  the positional correction, and its arm cannot matter because that term is
  multiplied by zero. No static early-out is needed: a body that cannot move
  leaves the effective mass at zero and the existing guard returns.
- Normalization happens once at the input boundary, in the scenario loader, so
  the geometry type keeps a strict unit-normal invariant instead of normalizing
  on every query. A direction that cannot be normalized is refused rather than
  silently turned into one, since `(0, 0)` and `(1e-300, 0)` would otherwise
  both look acceptable while meaning different things.
- There is deliberately no `area`, `centroid` or `second_moment_of_area` for a
  halfplane. All three are infinite, and returning 0 the way the invalid-shape
  path does would be indistinguishable from a real answer.
- Body against body is still disc against disc. Oriented boxes and SAT remain
  their own milestone; a halfplane is a signed-distance test, not a
  separating-axis loop.

### What conserves, and what does not

Nothing new. A ground plane is immovable and therefore external to the system,
exactly like a static body, so neither momentum nor energy is conserved in its
presence and neither is meant to be.

The property worth stating is the one that is now exact: a body sliding along a
plane keeps its tangential velocity to the last bit, and picks up no spin at
all, because the normal never turns. The test asserts that with a tolerance of
zero rather than a small number.

One rounding caveat, recorded rather than hidden. On a SLANTED plane the spin a
centred body picks up is about 7e-18 rather than exactly zero. The arm is
parallel to the impulse, so the true cross product vanishes, but `cross`
computes `arm.x * impulse.y - arm.y * impulse.x`, and on a slanted normal those
two products multiply the same pair of components in opposite orders and round
differently. On an axis-aligned plane it is exactly zero, because one factor in
each product is a literal zero.

### Tests

- Validation: a unit normal is required, and a zero, non-unit, infinite or NaN
  one is refused. Ground planes are validated with the rest of the settings, so
  a malformed one is refused before it can produce a NaN normal.
- Contact geometry off-origin and off-axis, including a 3-4-5 plane where both
  components matter, a ceiling that pushes the other way, exact touching
  reported with zero penetration, and the documented separation contract:
  moving the circle by `-normal * penetration` leaves it exactly touching.
- The case that would be degenerate for circle against circle: a centre lying
  exactly on the line still gets the plane's own normal, with no fallback.
- The milestone's point: a body sliding along a flat floor for 3000 steps keeps
  its horizontal velocity EXACTLY and its spin at exactly zero.
- A frictionless body on a 3-4-5 incline, derived in closed form for one step:
  gravity gives `(0, -0.1)`, the contact removes the normal component, and the
  result is `(0.048, -0.036)`, whose magnitude is `g dt sin(theta)` and whose
  direction is down the slope.
- An off-centre contact against a plane DOES generate spin, and mirroring the
  offset mirrors it, so the zero above is a property of the geometry rather
  than of planes being inert.
- A corner made of two planes, so a loop that resolved only the first would
  fail; a zero-radius body falling straight through; and an empty ground list
  reproducing the pre-M16 trajectory exactly.
- Verified by mutation testing. Seven mutations were applied and all seven were
  caught: a non-unit normal accepted, the offset dropped from the signed
  distance, the contact normal not negated, the penetration sign inverted, the
  half-depth dropped from the contact point, only the first ground plane
  resolved, and ground planes left unvalidated.

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
- Documentation count guards in the scenario tests. Counts stated in prose go
  stale, and this one had twice, so the numbers are now read back out of the
  documents and compared against reality: the template count in README.md
  against the `scenarios/` directory, and the test-executable count in both
  README.md and docs/00_START_HERE.md against the number CMake actually
  registered, which it passes in as `MALLOY_TEST_EXECUTABLE_COUNT`. A wrong
  number, or a missing one, fails the suite.

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
