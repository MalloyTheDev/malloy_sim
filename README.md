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

The locked **M1-M5 roadmap is complete**, plus the first post-M5 milestone
(**M6: N-body diagnostics**). The project builds clean under MSVC
(`/W4 /permissive-`), all five test executables pass via CTest, and the terminal
app runs two N-body scenarios (a two-body orbit and a three-body triangle)
reporting conserved system diagnostics.

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

It runs two hardcoded scenarios and prints conserved system diagnostics
(separation, total energy, total angular momentum) periodically. Sample output
(abridged):

```text
MalloySim nbody_terminal

== two-body orbit (normalized units) ==  bodies=2  dt=0.00100000  steps=10000
step      0   sep  1.00000000   E_total   -0.00000050   L_total    0.00000100
...
step  10000   sep  1.00027148   E_total   -0.00000050   L_total    0.00000100

== equilateral three-body (Lagrange) ==  bodies=3  dt=0.00100000  steps=10000
step      0   sep  1.73205081   E_total   -0.86601529   L_total    2.27951999
...
step  10000   sep  1.73143038   E_total   -0.86601509   L_total    2.27951999
```

The three-body triangle holds its shape (separation ~sqrt(3)) while total energy
and angular momentum stay essentially constant -- the conservation you expect
from semi-implicit (symplectic) Euler.

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
