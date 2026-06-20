# Changelog

All notable changes to MalloySim are recorded here. The format follows
`docs/10_CHANGELOG_TEMPLATE.md`; dates are ISO-8601.

## [Unreleased]

### Added

### Changed

### Fixed

### Removed

### Tests

### Architecture Notes

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
