# 02 — Milestone Roadmap M1-M5

## M1 — CMake Skeleton and Smoke Test

Goal: prove the toolchain and workflow.

M1 must prove:

- CMake configures.
- MSVC builds C++20.
- App target builds.
- Test target builds.
- CTest discovers/runs test.
- App prints one line.

M1 does not include `Vec2`, `Real`, physics, N-body, Catch2, rendering, CI, or empty module folders.

## M2 — `malloy_math`

Add target:

```text
malloy_math INTERFACE
malloy_math_tests EXECUTABLE
malloy::math ALIAS
```

Files:

```text
include/malloy/math/real.hpp
include/malloy/math/vec2.hpp
include/malloy/math/math.hpp
tests/math/vec2_tests.cpp
```

Contents:

- `Real = double`
- `Vec2`
- constructors
- arithmetic operators
- compound assignment operators
- unary minus
- `dot`
- `length_squared`
- `length`
- `distance_squared`
- `distance`
- `normalize`
- `approx_equal(Real, Real, epsilon)`
- `approx_equal(Vec2, Vec2, epsilon)`
- `is_finite(Real)`
- `is_finite(Vec2)`

Do not add `Vec3`, `Mat4`, `Quat`, templates, SIMD, GLM, Eigen, or a units library.

## M3 — `malloy_time` + tiny `malloy_sim_core`

Add targets:

```text
malloy_time INTERFACE
malloy_sim_core STATIC
malloy_time_tests EXECUTABLE
malloy_sim_core_tests EXECUTABLE
```

`malloy_time` contains `FixedStep`, `dt`, `tick_count`, `elapsed_time`, `advance()`, and positive dt validation.

`malloy_sim_core` contains `SimulationSettings`, `StepStatus`, and `StepResult`.

No real-time loop, wall-clock timer, sleep, frame pacing, render interpolation, scheduler, `ISimulation`, `virtual step()`, `WorldBase`, or engine kernel.

## M4 — `malloy_nbody`

Add targets:

```text
malloy_nbody STATIC
malloy_nbody_tests EXECUTABLE
malloy::nbody ALIAS
```

Types:

- `Body2D`
- `NBodySettings`
- `NBodyWorld`

Integration order:

1. clear accelerations
2. compute pairwise accelerations
3. update velocities using `acceleration * dt`
4. update positions using `velocity * dt`
5. increment tick
6. return `StepResult`

Required M4 tests:

- `dt <= 0` invalid
- `G < 0` invalid
- `softening < 0` invalid
- `mass <= 0` invalid
- non-finite position invalid
- non-finite velocity invalid
- one body alone gets zero acceleration
- zero/tiny distance does not produce NaN with softening
- pairwise force symmetry: `mass_a * accel_a ≈ -(mass_b * accel_b)`
- same initial state repeated gives equivalent final state
- near-circular two-body radius stays within tolerance for a short run

## M5 — Terminal N-body Demo

Update `apps/nbody_terminal/main.cpp`.

Behavior:

- hardcoded sun + planet
- fixed timestep
- fixed step count
- formatted terminal output
- periodic output
- print radius
- print specific orbital energy
- return nonzero on validation/step failure

No CLI parser, config files, logging framework, graphics, input loop, real-time sleep, rendering library, collision, rigid bodies, 3D, quantum, ECS, or serialization.

## M5 constants

```text
G = 1.0
sun_mass = 1.0
planet_mass = 0.000001
radius = 1.0
sun_position = (0.0, 0.0)
sun_velocity = (0.0, 0.0)
planet_position = (1.0, 0.0)
planet_velocity = (0.0, 1.0)
softening = 0.000001
dt = 0.001
steps = 10000
output_every = 1000
```
