# 04 — Numeric and Physics Conventions

> These conventions are in force in the shipped M1-M5 code (the softening
> formula and validation rules below are implemented in `malloy_nbody`).

Use:

```cpp
using Real = double;
```

Location:

```text
include/malloy/math/real.hpp
```

Namespace:

```cpp
malloy::math::Real
```

`Vec2` stores `Real`. Do not template `Vec2` in v1.

## Floating point rules

- Use `double` for positions and velocities.
- Do not use `float` for simulation state.
- Do not enable fast-math.
- Do not use `==` for floating-point simulation comparisons.
- Use explicit tolerances.

## Determinism scope

Do not claim cross-platform bitwise determinism.

V1 determinism means same binary, same platform, same initial state, same fixed timestep, same update order, and same result within expected tolerance.

## Units

Use normalized demo units first. Do not use real astronomical SI units in the first demo.

## N-body softening

Use softening from the first N-body implementation.

Formula:

```text
delta = other.position - current.position
r2 = dot(delta, delta) + softening * softening
inv_r = 1 / sqrt(r2)
inv_r3 = inv_r * inv_r * inv_r
acceleration += G * other.mass * delta * inv_r3
```

## Validation rules

Validate `dt > 0`, `G >= 0`, `softening >= 0`, `mass > 0`, finite position, and finite velocity.

Use status/result returns for normal validation failures. Do not throw for normal validation failures in M1-M5.
