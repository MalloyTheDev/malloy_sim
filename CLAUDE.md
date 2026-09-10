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

### What it will never be

- a game engine clone
- a Unity/Godot/Unreal replacement
- a render-first project
- an editor project
- a plugin framework
- a scripting framework

Those are rejected on purpose, not deferred. ECS is a separate case: not
rejected, but not adopted until real access-pattern pressure exists (rule 3).

### Where it is going

**MalloySim is intended to become a 3D simulator.** It is 2D today and every
milestone so far is 2D, but three dimensions are the destination rather than a
door left open, and quantum is a further one
(`docs/decisions/0009-three-dimensions-are-the-destination.md`).

Neither has started, and rules 6 and 7 below still gate them: 3D arrives as its
own milestone, once the 2D mechanics it rests on are finished. The point of
saying it here is that several decisions already made only make sense as
preparation for it, and someone reading the rules should know which of today's
2D choices are stepping stones and which would have to change.

## Current phase

```text
M1-M18 complete: math, time, sim_core, N-body, terminal demo, diagnostics,
scenario loading, ASCII debug view, collision primitives, colliding particles
with multi-domain scenario dispatch, 2D rigid bodies, ballistics, spring
networks with deterministic force accumulation, rigid-body contact
response, uniform gravity for rigid bodies, halfplane ground, Coulomb
friction, and charged particles in electric and magnetic fields.
```

Active track: **classical mechanics depth**, now complete. M9 (collision
geometry), M10 (colliding particles, which introduced the multi-domain `type`
key), M11 (2D rigid bodies), M12 (ballistics, as uniform gravity in the particle
domain), M13 (spring networks and deterministic force accumulation), M14
(rigid-body contact response) and M15 (uniform gravity in the rigid domain)
M16 (halfplanes, giving true flat ground), M17 (Coulomb friction) and M18
(charged particles, the first domain with a velocity-dependent force) are all
done. See
`docs/07_POST_M5_ROADMAP.md`.

Do not start any further milestone (M19 or later) unless explicitly asked, and
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
4. Collision geometry landed in M9, non-rotational contact response in M10, rotational (rigid-body) contact response in M14, halfplanes in M16, and Coulomb friction in M17. Body against body is still disc against disc; oriented boxes and SAT are not implemented. Friction is in `malloy_rigid` only: particle contacts remain normal-only.
5. Rigid bodies landed in M11, gained contact response in M14, a uniform gravity field in M15, and Coulomb friction in M17. Friction is a contact impulse clamped to the normal impulse, not a persistent force. Gravity is a setting applied as an acceleration, not a force: do not add persistent forces or force/torque accumulators to `malloy_rigid` until their dedicated milestone.
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
17. `SpringWorld` is a local composition boundary, not the engine-wide force architecture (`docs/decisions/0008-spring-world-is-a-local-composition-boundary.md`). Do not generalise it into a shared force-provider API until several genuinely different force producers exist. If a third domain independently needs the same translational integration path, reassess extracting a shared integrator.

## Multi-domain architecture

Each physics domain is its own library with its own concrete world type, its own
settings type, and its own test executable. Domains share only `malloy::math`
and the tiny `malloy::sim_core` vocabulary (`SimulationSettings`, `StepStatus`,
`StepResult`). No domain knows that any other domain exists.

Domains are selected by a `type` key in the scenario file, dispatched with a
plain switch to one concrete loader and one concrete world per domain.

Implemented in M10, once `malloy_particles` gave the format a second domain to
dispatch to:

```text
type nbody          # or: particles, rigid, springs, charges
dt 0.001
steps 10000
```

`type` defaults to `nbody` when absent, so scenarios written before the key
existed keep working unchanged. A key belonging to another domain is a parse
error, so a typo in `type` surfaces immediately instead of silently running the
wrong simulation.

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
