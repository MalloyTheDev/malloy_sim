# 01 - V3 Architecture Decision

## Final verdict

This architecture is fully realized: every target below was built across M1-M5
and the dependency direction holds in the shipped code.

The V3 architecture is:

```text
C++20 + CMake + VS Code + MSVC
terminal-first
simulation-first
2D first
no graphics before the rendering milestone
no ECS before real pressure
no package manager until a milestone needs one
```

## Core identity

MalloySim is an all-in-one science physics simulation in C++, built
simulation-first and terminal-first.

It started as a terminal 2D N-body gravity simulation. It grows one finished
physics domain at a time: collision, 2D rigid bodies, ballistics, springs and
oscillators, then fluids, 3D, and quantum. Each domain is a concrete library
with its own world type, settings type, and tests, selected by a scenario
`type` key rather than through a generic engine interface
(`docs/decisions/0006-multi-domain-dispatch.md`).

The measure of the project is depth, not coverage: a domain counts when it
validates its input, holds its invariants, and ships tested templates.

## Final locked decisions

| Area | Decision |
|---|---|
| Language | C++20 |
| Build system | CMake |
| Minimum CMake | 3.21 for the project, 4.2 for the shipped presets |
| Editor | VS Code |
| VS Code workflow | CMake Tools |
| Compiler first | MSVC on Windows |
| Package manager | None until a milestone needs one |
| Test framework | CTest + tiny custom CHECK macros |
| Test executable style | One test executable per module |
| CMake layout | Single root `CMakeLists.txt` |
| Numeric type | `malloy::math::Real = double` |
| Timestep | Fixed timestep only |
| Graphics | Excluded until the rendering milestone |
| ECS | Excluded until real access-pattern pressure |
| First app | `malloy_nbody_terminal` |
| First simulation | Normalized two-body orbit |

## Final target list by M5 (all built)

```text
malloy_project_options      INTERFACE
malloy_project_warnings     INTERFACE

malloy_math                 INTERFACE
malloy_time                 INTERFACE
malloy_sim_core             STATIC
malloy_nbody                STATIC

malloy_nbody_terminal       EXECUTABLE

malloy_smoke_tests          EXECUTABLE
malloy_math_tests           EXECUTABLE
malloy_time_tests           EXECUTABLE
malloy_sim_core_tests       EXECUTABLE
malloy_nbody_tests          EXECUTABLE
```

## Final dependency direction

```text
app
 ↓
malloy_nbody
 ↓
malloy_sim_core + malloy_math
 ↓
foundation
```

No dependency may point upward.

`malloy_time` is built and tested but is currently linked by nothing except its
own test executable: `NBodyWorld` keeps its own tick count and
`SimulationSettings` holds its own `dt`. It is shown outside the chain above on
purpose. Adopting `FixedStep` in `NBodyWorld` would make the chain read
`malloy_sim_core + malloy_time + malloy_math`, but that is a design change for
its own milestone, not a documentation fix.

## Biggest risk

The biggest risk is `malloy_sim_core` becoming a speculative engine framework.

It must stay tiny. It should contain only shared vocabulary and status/result types.

It must not contain `ISimulation`, `virtual step()`, `WorldBase`, `Engine`, `Entity`, `Component`, `Scene`, `Scheduler`, `EventBus`, `ServiceLocator`, `Plugin`, `Registry`, or `SystemManager`.

`NBodyWorld` remains concrete and owns its own `step()` behavior.

As of M8 this holds: `malloy_sim_core` contains only `SimulationSettings`,
`StepStatus`, and `StepResult`, and `NBodyWorld` is a concrete class. Keep it
that way.

The all-in-one goal raises this risk rather than retiring it, because "support
every domain" is the exact argument that produces an engine kernel. ADR 0006 is
the standing answer: many concrete worlds and a dispatch key, never a base
class. Adding a domain must leave `malloy_sim_core` unchanged.
