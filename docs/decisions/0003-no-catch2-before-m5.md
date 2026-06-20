# ADR 0003 — No Catch2 Before M5

## Status

Accepted.

## Decision

M1-M5 use CTest plus tiny custom check macros.

Catch2 is delayed until after M5.

## Rationale

The first milestones should avoid external dependency friction.

CTest is already available through CMake. The early tests are small enough that `tests/test_check.hpp` is sufficient.
