# MalloySim C++

MalloySim is a from-scratch **simulation-first C++20 engine project**.

It starts as a terminal-only 2D classical simulation project. It is not a Unity, Godot, Unreal, ECS, editor, rendering, plugin, or scripting project.

Its first serious goal, a deterministic terminal-based 2D N-body gravity demo, is complete (milestones M1-M5).

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

## Status

Milestones **M1-M5 are complete** -- the locked initial roadmap is done. The
project builds clean under MSVC (`/W4 /permissive-`), all five test executables
pass via CTest, and the terminal app runs a stable two-body orbit.

| Module | Type | Provides |
|---|---|---|
| `malloy_math` | INTERFACE | `Real = double`, `Vec2`, vector/scalar helpers |
| `malloy_time` | INTERFACE | `FixedStep` (fixed timestep, tick count, elapsed time) |
| `malloy_sim_core` | STATIC | `SimulationSettings`, `StepStatus`, `StepResult` |
| `malloy_nbody` | STATIC | `Body2D`, `NBodySettings`, `NBodyWorld`, softened gravity |
| `malloy_nbody_terminal` | EXECUTABLE | the terminal N-body demo |

Post-M5 work is intentionally gated -- see `docs/07_POST_M5_ROADMAP.md`.

## Build

From a Developer PowerShell / Developer Command Prompt with MSVC available:

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

## Run the demo

With the Visual Studio generator, the debug executable is under:

```powershell
.\out\build\windows-msvc-debug\Debug\malloy_nbody_terminal.exe
```

It runs the normalized sun + planet orbit and prints radius and specific
orbital energy every 1000 steps. Sample output (abridged):

```text
MalloySim nbody_terminal
two-body orbit (normalized units): G=1 dt=0.001 steps=10000
step      0   radius   1.00000000   specific_energy    -0.50000100
step   1000   radius   0.99957897   specific_energy    -0.50000089
...
step  10000   radius   1.00027148   specific_energy    -0.50000054
```

The radius stays close to 1.0 and the specific energy is conserved to ~1e-7
over 10,000 steps -- the expected behavior of semi-implicit (symplectic) Euler.

## M1-M5 roadmap (complete)

| Milestone | Goal | Status |
|---|---|---|
| M1 | CMake skeleton + smoke test | ✅ Done |
| M2 | `malloy_math`: `Real`, `Vec2`, math tests | ✅ Done |
| M3 | `malloy_time` + tiny `malloy_sim_core` | ✅ Done |
| M4 | `malloy_nbody`: `Body2D`, `NBodyWorld`, gravity, tests | ✅ Done |
| M5 | Terminal N-body sun/planet demo | ✅ Done |

## Out of scope (gated)

These remain out of scope until taken up as their own dedicated post-M5
milestone (`docs/07_POST_M5_ROADMAP.md`) -- never added speculatively:

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

## Before committing

Keep the build and tests green:

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

For the full plan and history, read:

- `CLAUDE.md`
- `CHANGELOG.md`
- `docs/00_START_HERE.md`
- `docs/02_MILESTONE_ROADMAP_M1_M5.md`
- `docs/07_POST_M5_ROADMAP.md`
