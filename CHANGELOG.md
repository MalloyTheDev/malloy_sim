# Changelog

All notable changes to MalloySim are recorded here. The format follows
`docs/10_CHANGELOG_TEMPLATE.md`; dates are ISO-8601.

## [Unreleased]

### Added

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

### Fixed

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
- Every new assertion was confirmed by mutation testing: eight broken
  implementations (explicit Euler, dropped G, inverse-square, unsquared
  softening, deleted angular-momentum term, both removed `is_finite` clauses,
  and a removed coincident guard) each make the suite fail.

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
