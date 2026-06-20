# 01 — V3 Architecture Decision

## Final verdict

This architecture is fully realized: every target below was built across M1-M5
and the dependency direction holds in the shipped code.

The V3 architecture is:

```text
C++20 + CMake + VS Code + MSVC
terminal-first
simulation-first
2D first
no graphics before M5
no ECS before real pressure
no package manager before M5
```

## Core identity

MalloySim is a simulation-first C++ project.

It starts with a terminal 2D N-body gravity simulation and later may grow into collision, 2D rigid-body physics, ballistics, vehicles, a visual sandbox, 3D, and quantum as a separate domain.

## Final locked decisions

| Area | Decision |
|---|---|
| Language | C++20 |
| Build system | CMake |
| Minimum CMake | 3.21 |
| Editor | VS Code |
| VS Code workflow | CMake Tools |
| Compiler first | MSVC on Windows |
| Package manager | None for M1-M5 |
| Test framework | CTest + tiny custom CHECK macros |
| Test executable style | One test executable per module |
| CMake layout | Single root `CMakeLists.txt` through M5 |
| Numeric type | `malloy::math::Real = double` |
| Timestep | Fixed timestep only |
| Graphics | Excluded until after M5 |
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
malloy_sim_core + malloy_time + malloy_math
 ↓
foundation
```

No dependency may point upward.

## Biggest risk

The biggest risk is `malloy_sim_core` becoming a speculative engine framework.

It must stay tiny. It should contain only shared vocabulary and status/result types.

It must not contain `ISimulation`, `virtual step()`, `WorldBase`, `Engine`, `Entity`, `Component`, `Scene`, `Scheduler`, `EventBus`, `ServiceLocator`, `Plugin`, `Registry`, or `SystemManager`.

`NBodyWorld` remains concrete and owns its own `step()` behavior.

As of M5 this held: `malloy_sim_core` contains only `SimulationSettings`,
`StepStatus`, and `StepResult`, and `NBodyWorld` is a concrete class. Keep it
that way.
