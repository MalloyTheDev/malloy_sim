# MalloySim C++

MalloySim is a from-scratch **simulation-first C++20 engine project**.

It starts as a terminal-only 2D classical simulation project. It is not a Unity, Godot, Unreal, ECS, editor, rendering, plugin, or scripting project.

The first serious goal is a deterministic terminal-based 2D N-body gravity demo.

## Locked baseline

- Language: C++20
- Build system: CMake
- Editor: VS Code + CMake Tools
- Compiler first: MSVC on Windows
- Tests: CTest + tiny custom check macros for M1-M5
- Package manager: none for M1-M5
- Graphics: none for M1-M5
- Simulation timestep: fixed timestep only
- Numeric type: `double`

## Current milestone

This repository begins at **M1: CMake skeleton and smoke test**.

M1 proves only:

- CMake configures.
- MSVC builds C++20.
- One app target builds.
- One smoke test target builds.
- CTest discovers and runs the test.
- The terminal app prints one line.

M1 does not prove physics yet.

## Build

From a Developer PowerShell / Developer Command Prompt with MSVC available:

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

## Run the M1 app

With the Visual Studio generator, the debug executable should be under:

```powershell
.\out\build\windows-msvc-debug\Debug\malloy_nbody_terminal.exe
```

Expected output:

```text
MalloySim nbody_terminal
```

## M1-M5 roadmap

| Milestone | Goal |
|---|---|
| M1 | CMake skeleton + smoke test |
| M2 | `malloy_math`: `Real`, `Vec2`, math tests |
| M3 | `malloy_time` + tiny `malloy_sim_core` |
| M4 | `malloy_nbody`: `Body2D`, `NBodyWorld`, gravity, tests |
| M5 | Terminal N-body sun/planet demo |

## What not to build yet

Do not add these before M5:

- rendering
- Raylib / SDL / GLFW / SFML
- collision
- rigid bodies
- vehicles
- 3D
- quantum
- ECS
- plugins
- editor
- scripting
- asset manager
- threading
- package manager
- config/scenario files
- CLI parser

## First commit target

Commit message:

```text
M1: CMake skeleton and smoke test
```

Required before commit:

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

For the full plan, read:

- `CLAUDE.md`
- `docs/00_START_HERE.md`
- `docs/01_V3_ARCHITECTURE_DECISION.md`
- `docs/02_MILESTONE_ROADMAP_M1_M5.md`
