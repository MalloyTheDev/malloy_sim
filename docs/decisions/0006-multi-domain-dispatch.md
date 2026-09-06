# ADR 0006 - Multi-Domain Dispatch Without Abstraction

## Status

Accepted.

## Context

MalloySim's goal is an all-in-one science physics simulation: many domains,
each with ready-to-run templates. The obvious way to build that is a generic
engine with a simulation interface every domain implements.

That is exactly the failure mode this project was set up to avoid.
`docs/01_V3_ARCHITECTURE_DECISION.md` names a speculative engine framework in
`malloy_sim_core` as the biggest risk, `docs/09_MISTAKES_TO_AVOID.md` lists
turning `malloy_sim_core` into an engine kernel as mistake 1, and `CLAUDE.md`
rules 11 and 12 forbid `ISimulation`, `WorldBase`, `Engine`, and friends.

The two goals look contradictory. They are not.

## Decision

Each physics domain is its own library with its own concrete world type, its own
settings type, and its own test executable. Domains share only `malloy::math`
and the tiny `malloy::sim_core` vocabulary. No domain references another, and
no domain derives from anything.

A scenario file selects its domain with a `type` key. The loader dispatches on
that key with a plain switch, calling one concrete parser and constructing one
concrete world.

This is planned, not implemented: the shipped parser rejects `type` as an
unknown key, and there is no dispatch switch yet. The intended shape is:

```text
type nbody          # planned, not parsed today
dt 0.001
steps 10000
body 1.0  0.0 0.0  0.0 0.0
```

There is no base class, no virtual `step()`, no registry, and no plugin
mechanism. Adding a domain means adding a library, a parser branch, and a test
executable. It does not mean touching `malloy_sim_core`.

The `type` key is introduced only when a second domain actually exists. Until
then the scenario format stays as M7 shipped it.

## Rationale

All-in-one science simulators in the wild are not single generic engines. They
are collections of concrete simulations behind a thin shell. The generic engine
is what a project builds when it wants to look finished before any physics
works.

Concrete-per-domain keeps every property this project already has:

- `malloy_sim_core` stays three types.
- `NBodyWorld` stays concrete and owns its own `step()`.
- Each domain validates its own input and returns status rather than throwing.
- Each domain is independently testable, and a test executable links only the
  module it tests, so dependency leaks stay visible.
- The dependency direction still points one way, with nothing pointing upward.

The cost is real and accepted: a switch statement grows by one branch per
domain, and shared behavior is duplicated rather than inherited. That
duplication is cheaper than the abstraction it replaces, and it can be
factored later against evidence from domains that actually exist.

## Consequences

- Breadth is bounded by finished domains, not by framework capability. This is
  intended: depth over breadth.
- A domain is not done until it has validation, an invariant checked by tests,
  malformed-input tests, and at least one scenario template (`CLAUDE.md` rule 16).
- If a genuine shared abstraction ever becomes obvious from three or more real
  domains, it may be extracted then, with this ADR superseded rather than
  quietly ignored.
