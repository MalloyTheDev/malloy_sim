# MalloySim C++

MalloySim is a from-scratch **all-in-one science physics simulation** in C++20.

The goal is a robust simulator spanning many physics domains, each shipping
ready-to-run templates so you can pick a scenario and get correct physics
immediately. It is terminal-first, and it is built as many concrete domain
simulations behind a thin shared shell rather than as one generic engine
(`docs/decisions/0006-multi-domain-dispatch.md`).

It is not a Unity, Godot, Unreal, ECS, editor, rendering, plugin, or scripting
project.

Its first domain, a deterministic terminal 2D N-body gravity simulation, is
complete and shipping (milestones M1-M8).

## Locked baseline

- Language: C++20
- Build system: CMake
- Editor: VS Code + CMake Tools
- Compiler first: MSVC on Windows
- Tests: CTest + tiny custom check macros
- Package manager: none until a milestone needs one
- Graphics: none until the rendering milestone
- Simulation timestep: fixed timestep only
- Numeric type: `double`
- Domain model: one concrete world per domain, no engine kernel

## Status

The locked **M1-M5 roadmap is complete**, plus the first post-M5 milestones
(**M6: N-body diagnostics**, **M7: scenario loading**, **M8: ASCII debug view**).
The project builds clean under MSVC (`/W4 /permissive-`), all seven test
executables pass via CTest, and the terminal app runs N-body scenarios --
built-in, or loaded from a text file -- reporting conserved system diagnostics
alongside an ASCII view of the bodies.

| Module | Type | Provides |
|---|---|---|
| `malloy_math` | INTERFACE | `Real = double`, `Vec2`, vector/scalar helpers |
| `malloy_time` | INTERFACE | `FixedStep` (fixed timestep, tick count, elapsed time) |
| `malloy_sim_core` | STATIC | `SimulationSettings`, `StepStatus`, `StepResult` |
| `malloy_nbody` | STATIC | `Body2D`, `NBodySettings`, `NBodyWorld`, softened gravity, diagnostics |
| `malloy_scenario` | STATIC | parse a scenario/config text file into bodies + settings |
| `malloy_ascii` | STATIC | fit a viewport to 2D points, render them as a framed character grid |
| `malloy_nbody_terminal` | EXECUTABLE | the terminal N-body demo |

Post-M5 work is intentionally gated -- see `docs/07_POST_M5_ROADMAP.md`.

## Domains and templates

MalloySim grows one finished physics domain at a time. Each domain is its own
library with its own concrete world type, its own settings, and its own test
executable. Domains share only `malloy_math` and the tiny `malloy_sim_core`
vocabulary; none of them knows the others exist. There is no simulation base
class and no engine kernel, by design.

A domain counts as finished only when it has all four of:

1. validation of its own settings and state, returning status rather than throwing;
2. an invariant or conserved quantity checked by tests;
3. malformed and boundary input tests;
4. at least one scenario template in `scenarios/`.

Templates in `scenarios/` are a first-class deliverable: plain text, documented,
and runnable with the shipped binary. Gravity ships two today.

The active track is classical mechanics depth (collision, rigid bodies,
ballistics, springs and oscillators). Work is gated one milestone at a time; see
`docs/07_POST_M5_ROADMAP.md`.

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

It runs two hardcoded scenarios and prints, at each report, one line of conserved
system diagnostics (separation, total energy, total angular momentum) followed by
an ASCII view of the bodies. Sample output (diagnostics only, abridged):

```text
MalloySim nbody_terminal

== two-body orbit (normalized units) ==  bodies=2  dt=0.00100000  steps=10000
step      0   sep  1.00000000   E_total  -5.00000000e-07   L_total   1.00000000e-06
...
step  10000   sep  1.00027148   E_total  -4.99999540e-07   L_total   1.00000000e-06

== equilateral three-body (Lagrange) ==  bodies=3  dt=0.00100000  steps=10000
step      0   sep  1.73205081   E_total  -8.66015288e-01   L_total   2.27951999e+00
...
step  10000   sep  1.73143038   E_total  -8.66015093e-01   L_total   2.27951999e+00
```

The three-body triangle holds its shape (separation ~sqrt(3)) while total energy
and angular momentum stay essentially constant -- the conservation you expect
from semi-implicit (symplectic) Euler.

The conserved quantities print in scientific notation so that drift stays
visible at any magnitude. The two-body demo's energy is of order 1e-6, and under
fixed-point formatting every step showed the same value, which would have looked
perfectly conserved even if it were not. Watch `E_total` oscillate rather than
walk in one direction: a bounded oscillation is the symplectic signature, while
a monotone drift would indicate the integrator had been broken.

### The ASCII view

Each report is followed by a frame of the body positions. This is step 1000
of the two-body orbit, with the sun on the left and the planet a sixth of the
way round:

```text
+-------------------------------------------------------------+
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                *                            |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|     *                                                       |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
|                                                             |
+-------------------------------------------------------------+
  view x [-0.10000000, 1.10000000]  y [-1.20000000, 1.20000000]
```

The view is fitted to the bodies at the start of a run and afterwards only
grows, so frames share a scale, a near-stationary body keeps its cell, and no
body is ever silently clipped. The printed extents show when it has grown.

### Load a scenario from a file

Pass a scenario file as the single argument:

```powershell
.\out\build\windows-msvc-debug\Debug\malloy_nbody_terminal.exe scenarios\two_body.scn
```

Scenario files are plain text (`#` starts a comment):

```text
g 1.0
softening 0.000001
dt 0.001
steps 10000
output_every 1000
body 1.0       0.0 0.0   0.0 0.0
body 0.000001  1.0 0.0   0.0 1.0
```

Example scenarios live in `scenarios/`. With no argument, the app runs the
built-in scenarios shown above.

## Milestones complete

| Milestone | Goal | Status |
|---|---|---|
| M1 | CMake skeleton + smoke test | ✅ Done |
| M2 | `malloy_math`: `Real`, `Vec2`, math tests | ✅ Done |
| M3 | `malloy_time` + tiny `malloy_sim_core` | ✅ Done |
| M4 | `malloy_nbody`: `Body2D`, `NBodyWorld`, gravity, tests | ✅ Done |
| M5 | Terminal N-body sun/planet demo | ✅ Done |
| M6 | N-body system diagnostics + three-body demo | ✅ Done |
| M7 | `malloy_scenario`: scenario/config text loading | ✅ Done |
| M8 | `malloy_ascii`: 2D ASCII debug visualization | ✅ Done |

## Out of scope (gated)

These remain out of scope until taken up as their own dedicated post-M5
milestone (`docs/07_POST_M5_ROADMAP.md`) -- never added speculatively:

- graphical rendering (the M8 debug view is ASCII text only)
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
