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

All eleven exist and pass (one per module):

```text
malloy_smoke_tests
malloy_math_tests
malloy_time_tests
malloy_sim_core_tests
malloy_nbody_tests
malloy_scenario_tests
malloy_ascii_tests
malloy_collide_tests
malloy_particles_tests
malloy_rigid_tests
malloy_springs_tests
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

## When to migrate to Catch2

Only consider Catch2 after M5 when test count, diagnostics, parameterization, or CI reporting become painful.
