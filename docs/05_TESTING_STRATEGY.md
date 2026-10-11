# 05 - Testing Strategy

## Final decision

M1-M5 use:

```text
CTest + tiny custom CHECK macros
```

Do not use Catch2 or GoogleTest in M1-M5.

## Why

- zero external dependency friction
- no FetchContent/network failures
- no package manager
- CTest is already available
- simple tests are enough for M1-M5
- per-module test executables catch LINK-time dependency leaks

Note the limit: every module puts the whole `include/` tree on its consumers'
include path, so a header-only leak (using another module's header without
linking it) still compiles. The guard is real but it is link-time only.

## Test executable layout

All fourteen exist and pass (one per module, plus `malloy_quat_tests` for the
quaternion half of `malloy_math` and `malloy_layering_tests`, which checks the
build graph against the module boundaries):

```text
malloy_smoke_tests
malloy_math_tests
malloy_quat_tests
malloy_time_tests
malloy_sim_core_tests
malloy_nbody_tests
malloy_scenario_tests
malloy_ascii_tests
malloy_collide_tests
malloy_particles_tests
malloy_rigid_tests
malloy_springs_tests
malloy_charges_tests
malloy_layering_tests
```

Each test executable should return 0 on success, return nonzero on failure, print useful failure messages, and link only the module it tests.

## Test helpers

Use `MALLOY_CHECK_TRUE`, `MALLOY_CHECK_FALSE`, `MALLOY_CHECK_EQ`,
`MALLOY_CHECK_NEAR`, and `MALLOY_CHECK_VEC2_NEAR` (all in `tests/test_check.hpp`).


## Prefer asymmetric configurations

A test whose setup is symmetric can pass against a wrong implementation, because
the symmetry cancels the error. This has already happened twice in this project,
both times found by mutation testing rather than by review:

- `distance(a, b)` was only ever tested with `b` at the origin, where `a - b`
  and `a + b` are identical. An implementation using the wrong operator passed.
- `approx_equal(Vec2)` had no case where x was the sole difference, so an
  implementation checking only y passed everything.

The same shape of mistake hides sign errors, swapped operands, reversed normals,
and reversed impulses. When writing a physics test, deliberately break symmetry:

- unequal masses, unequal radii, unequal extents;
- positions and velocities that are off-axis and off-origin;
- values distinct in every field, so a swapped pair changes the result;
- for anything rotational: an off-centre application point, a body whose origin
  is not its centre of mass, and a nonzero initial angular velocity.

A useful check on a new test is to ask what wrong implementation would still
pass it, and then either add a case that would not, or confirm by mutation.

## Statistical checks (M38)

A gas is chaotic, but a run is deterministic (docs/04), so a check on a
statistical quantity returns the same number every time on the same binary and
platform. That makes it tempting to set the tolerance from the observed value.
Do not: a tolerance fitted to one run proves nothing about the physics.

Instead, every tolerance on a statistical quantity is DERIVED and the derivation
is written beside the check:

- three standard deviations of the quantity's own scatter, from its sampling
  distribution where one is known (binomial for a fraction, Beta for an energy
  share), or measured as the run-to-run spread when the timestep is perturbed
  slightly (which reshuffles every collision) for a time average like pressure;
- plus an explicit budget for model error when the comparison is with an
  approximate theory (a first-order wall correction, a finite step);
- and it must EXCLUDE the obvious wrong answer. Say what it excludes: the ideal
  gas pressure, the unrelaxed start, the starting energy ratio.

Exact quantities keep tight tolerances even in a statistical template: a
conserved temperature, or a statistic of the deterministic start.

Since M39 some quantities are random, drawn from a seeded generator. A
template states its seed, and its tolerances come from the spread measured
across many seeds (40 to 60 for the M39 templates), never from the one seed it
ships with. The generator itself is pinned by known answers rather than
statistics, at compile time and at run time, and the draw order by exact
comparison with the public samplers.

A bound can replace a standard deviation when the error has one. The point
particle ideal-gas test computes its own tolerance from the run (whole bounces
counted, plus the clamp lag), and also asserts that tolerance is tight enough to
mean something.

## When to migrate to Catch2

Only consider Catch2 after M5 when test count, diagnostics, parameterization, or CI reporting become painful.
