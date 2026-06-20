# 05 — Testing Strategy

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
- per-module test executables catch dependency leaks

## Test executable layout

All five exist and pass (one per module):

```text
malloy_smoke_tests
malloy_math_tests
malloy_time_tests
malloy_sim_core_tests
malloy_nbody_tests
```

Each test executable should return 0 on success, return nonzero on failure, print useful failure messages, and link only the module it tests.

## Test helpers

Use `MALLOY_CHECK_TRUE`, `MALLOY_CHECK_FALSE`, `MALLOY_CHECK_EQ`,
`MALLOY_CHECK_NEAR`, and `MALLOY_CHECK_VEC2_NEAR` (all in `tests/test_check.hpp`).

## When to migrate to Catch2

Only consider Catch2 after M5 when test count, diagnostics, parameterization, or CI reporting become painful.
