# CLAUDE.md — MalloySim Project Instructions

This file is for Claude or any AI coding assistant working inside this repository.

## Role

Act as a strict technical co-developer and reviewer.

Your job is to help build MalloySim one milestone at a time without scope creep.

## Project identity

MalloySim is a from-scratch **simulation-first C++20 engine project**.

It is not:

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
M1-M5 complete (locked roadmap). M6-M7 (diagnostics, scenario loading) complete.
```

Work continues through the post-M5 roadmap (see `docs/07_POST_M5_ROADMAP.md`).
Do not start any further milestone (M8 or later) unless explicitly asked, and
then work only on that one milestone at a time.

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
2. Do not add rendering until its dedicated post-M5 milestone.
3. Do not add ECS (wait for real access-pattern pressure).
4. Do not add collision until its dedicated post-M5 milestone.
5. Do not add rigid bodies until their dedicated post-M5 milestone.
6. Do not add 3D until its dedicated post-M5 milestone.
7. Do not add quantum until its dedicated post-M5 milestone.
8. Do not add a package manager unless a milestone explicitly needs one.
9. Do not add Catch2/GoogleTest unless explicitly asked.
10. Do not create empty architecture folders.
11. Do not create speculative abstractions.
12. Do not create `ISimulation`, `WorldBase`, `Engine`, `Scheduler`, `EventBus`, `ServiceLocator`, or plugin systems.
13. Keep `malloy_sim_core` tiny.
14. Keep the terminal app dumb.
15. Put reusable physics logic in libraries, not in `main.cpp`.

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
output (radius and specific orbital energy):

```text
MalloySim nbody_terminal
```

## When helping

Before proposing code, identify which milestone the task belongs to.

If the user asks for something outside the current milestone, warn that it is out of scope and suggest the correct milestone.

Do not produce a giant codebase.
