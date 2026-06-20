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

The repository starts at:

```text
M1: CMake skeleton and smoke test
```

Do not implement M2 or later unless explicitly asked.

## Locked M1-M5 order

```text
M1: CMake skeleton + smoke test
M2: malloy_math
M3: malloy_time + tiny malloy_sim_core
M4: malloy_nbody
M5: terminal N-body demo
```

## Rules

1. Do not jump ahead.
2. Do not add rendering.
3. Do not add ECS.
4. Do not add collision before M5.
5. Do not add rigid bodies before M5.
6. Do not add 3D before M5.
7. Do not add quantum before M5.
8. Do not add a package manager before M5.
9. Do not add Catch2 before M5.
10. Do not create empty architecture folders.
11. Do not create speculative abstractions.
12. Do not create `ISimulation`, `WorldBase`, `Engine`, `Scheduler`, `EventBus`, `ServiceLocator`, or plugin systems.
13. Keep `malloy_sim_core` tiny.
14. Keep the terminal app dumb.
15. Put reusable physics logic in libraries, not in `main.cpp`.

## Testing policy

M1-M5 use:

```text
CTest + tiny custom CHECK macros
```

Do not replace this with Catch2 or GoogleTest unless explicitly asked after M5.

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

## Expected M1 app output

```text
MalloySim nbody_terminal
```

## When helping

Before proposing code, identify which milestone the task belongs to.

If the user asks for something outside the current milestone, warn that it is out of scope and suggest the correct milestone.

Do not produce a giant codebase.
