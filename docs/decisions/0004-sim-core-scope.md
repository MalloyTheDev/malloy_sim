# ADR 0004 — Tiny `malloy_sim_core`

## Status

Accepted.

## Decision

`malloy_sim_core` exists, but stays tiny and non-polymorphic.

It contains only shared simulation vocabulary: `SimulationSettings`, `StepResult`, and `StepStatus`.

It does not contain generic engine architecture.
