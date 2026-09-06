# ADR 0003 - No Catch2 Before M5

## Status

Accepted.

## Decision

The project uses CTest plus tiny custom check macros.

Catch2 and GoogleTest are not added unless explicitly asked for (`CLAUDE.md`
rule 9). This was originally scoped to M1-M5; the constraint outlived that
window and now stands on its own.

## Rationale

The first milestones should avoid external dependency friction.

CTest is already available through CMake. The early tests are small enough that `tests/test_check.hpp` is sufficient.
