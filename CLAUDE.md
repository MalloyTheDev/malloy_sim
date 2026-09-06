# CLAUDE.md - MalloySim Project Instructions

This file is for Claude or any AI coding assistant working inside this repository.

## Role

Act as a strict technical co-developer and reviewer.

Your job is to help build MalloySim one milestone at a time without scope creep.

## Project identity

MalloySim is a from-scratch **all-in-one science physics simulation** in C++20.

The goal is a robust, flagship-quality simulator spanning many physics domains,
each shipping ready-to-run templates so someone can pick a scenario and get
correct physics immediately.

It is built as **many concrete domain simulations behind a thin shared shell**,
not as one generic engine. That is the central architectural commitment
(`docs/decisions/0006-multi-domain-dispatch.md`), and it is what makes the
all-in-one goal compatible with rules 11 to 13 below.

Terminal-first. Depth over breadth: a domain is not finished because it runs,
it is finished when it validates its input, conserves what it should, holds
determinism, and ships tested templates.

It is still not:

- a game engine clone
- a Unity/Godot/Unreal replacement
- a render-first project
- an ECS experiment
- an editor project
- a plugin framework
- a scripting framework
- a 3D project yet
- a quantum project yet

## Current phase

```text
M1-M9 complete: math, time, sim_core, N-body, terminal demo, diagnostics,
scenario loading, ASCII debug view, collision primitives.
```

Active track: **classical mechanics depth**. M9 (collision primitives) is done;
rigid bodies, ballistics, and springs remain, then multi-domain scenario
dispatch, then wider domains. See `docs/07_POST_M5_ROADMAP.md`.

Do not start any further milestone (M10 or later) unless explicitly asked, and
then work only on that one milestone at a time. The all-in-one goal does not
license building ahead: it is reached one finished domain at a time.

## M1-M5 order (all complete)

```text
M1: CMake skeleton + smoke test        [done]
M2: malloy_math                        [done]
M3: malloy_time + tiny malloy_sim_core [done]
M4: malloy_nbody                       [done]
M5: terminal N-body demo               [done]
```

## Rules

1. Do not jump ahead; work one milestone at a time.
2. Do not add rendering until its dedicated milestone. Terminal-first still holds.
3. Do not add ECS (wait for real access-pattern pressure).
4. Collision geometry landed in M9 (`malloy_collide`). Do not add contact response until its dedicated milestone.
5. Do not add rigid bodies until their dedicated milestone.
6. Do not add 3D until its dedicated milestone.
7. Do not add quantum until its dedicated milestone.
8. Do not add a package manager unless a milestone explicitly needs one.
9. Do not add Catch2/GoogleTest unless explicitly asked.
10. Do not create empty architecture folders.
11. Do not create speculative abstractions.
12. Do not create `ISimulation`, `WorldBase`, `Engine`, `Scheduler`, `EventBus`, `ServiceLocator`, or plugin systems. All-in-one is reached with concrete worlds and a dispatch key, never with a base class.
13. Keep `malloy_sim_core` tiny. Adding a domain must not grow it.
14. Keep the terminal app dumb.
15. Put reusable physics logic in libraries, not in `main.cpp`.
16. A domain is not done until it has validation, an invariant checked by tests, malformed-input tests, and at least one scenario template.

## Multi-domain architecture

Each physics domain is its own library with its own concrete world type, its own
settings type, and its own test executable. Domains share only `malloy::math`
and the tiny `malloy::sim_core` vocabulary (`SimulationSettings`, `StepStatus`,
`StepResult`). No domain knows that any other domain exists.

Domains are selected by a `type` key in the scenario file, dispatched with a
plain switch to one concrete loader and one concrete world per domain.

PLANNED, NOT YET IMPLEMENTED. The shipped parser rejects `type` as an unknown
key. This is the shape the format will take when a second domain exists:

```text
type nbody          # planned, not parsed today
dt 0.001
steps 10000
```

Add the `type` key only once a second domain actually exists. Until then the
scenario format stays exactly as it is.

## Template library

`scenarios/` is a first-class deliverable, not a scratch folder. Every template:

- states in a header comment what it demonstrates;
- runs to completion with the shipped binary;
- has a test or a documented expected result.

## Testing policy

The project uses:

```text
CTest + tiny custom CHECK macros (tests/test_check.hpp)
```

Do not replace this with Catch2 or GoogleTest unless explicitly asked.

## Numeric policy

Use:

```cpp
using Real = double;
```

No `float` for simulation positions/velocities.

No fast-math.

No cross-platform bitwise determinism claims.

V1 determinism means:

```text
same binary + same platform + same initial state + same fixed timestep + same update order = repeatable output
```

## Build commands

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

## Expected app output

The terminal app's first line is its identity, followed by the periodic demo
output: one diagnostics line (separation, total energy, total angular momentum)
and one ASCII view frame per report.

```text
MalloySim nbody_terminal
```

## When helping

Before proposing code, identify which milestone the task belongs to.

If the user asks for something outside the current milestone, warn that it is out of scope and suggest the correct milestone.

Do not produce a giant codebase.
