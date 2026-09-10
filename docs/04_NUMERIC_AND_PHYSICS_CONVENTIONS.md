# 04 - Numeric and Physics Conventions

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

## The stored velocity is half a step behind

Worth knowing before reading any orbit result.

Eliminating the velocity from the two update lines of semi-implicit Euler gives

```text
x[n+1] - 2 x[n] + x[n-1] = dt^2 a(x[n])
```

which is exactly the Stormer-Verlet position recurrence. The trajectory is
therefore second-order accurate, but the STORED velocity,
`v[n+1] = (x[n+1] - x[n]) / dt`, is a backward difference: it approximates the
true velocity at `t[n+1] - dt/2`, not at `t[n+1]`.

Two consequences, both first order in dt.

Supplying a physically self-consistent `(x, v)` pair at t=0 is a half-step
inconsistency. A tangential velocity on a circular orbit is short of the
discrete circle by a radial component of `(dt/2)|a|`, which seeds an epicycle of
eccentricity `dt * Omega / 2`. In `scenarios/two_body.scn` that is 5e-4 against
a physical eccentricity of about 1e-6, so the radial wobble in that template is
roughly 500 times the physics and is entirely a startup artefact. It halves when
dt halves.

And every reported energy carries an offset of `(dt/2) dV/dt`, which is why
energy oscillates about a mean rather than sitting on it. The oscillation is
bounded and shows no secular drift, as a symplectic method should.

Neither is a defect to fix in passing. Correcting the initial conditions with a
half-step kick would change what an initial velocity MEANS across the whole
N-body domain and every template figure with it, so it belongs to a milestone
rather than a patch (issue #19).

## Coordinate magnitude budget

Two distinct limits, neither of them `DBL_MAX`, and the smaller one is the one
worth remembering.

### The one you can actually reach: about 1e12

`position += velocity * dt` stops advancing once `|v * dt|` falls below half an
ulp of the coordinate, and rounds to a FULL ulp when it sits just above half
one. The threshold is `|x| > 2 v dt / eps`, which is `9.0e15 * (v dt)`.

Measured, one body drifting at 1 m/s with `dt = 1e-3` for 5000 steps, so it
should travel exactly 5.0:

| coordinate | travelled | error |
|---|---|---|
| 1e12 | 4.8828125 | 2.3 per cent short |
| 1e13 | 9.765625 | 1.95x too far |
| 1e14 | 0.0 | frozen |

The middle row is the one to fear: the body moves at almost exactly twice its
stated speed. And nothing in the output can show it, because every diagnostic
is computed from `velocity`, which the rounding never touches. Energy and
momentum read as perfectly conserved throughout.

**Budget: keep `|x|` below about 1e12 for metre-scale motion at millisecond
steps.** Scale it with `v * dt` if either differs. This is a documented
convention rather than a validated one, because the useful limit depends on how
much precision a given simulation needs, not on what is representable.

### The one that is validated: about 1.34e154

Everything that squares a vector reaches infinity above `sqrt(DBL_MAX)`, while
the vector itself is still finite. `is_valid` used to test only `is_finite`,
which accepts up to 1.8e308, so a window covering the entire upper half of the
exponent range reported a world as sound while every energy it printed was
`inf`.

Every body type now requires `math::is_squarable` on its position and velocity,
so that window is closed. The same reasoning bounds `softening` in the N-body
and charge domains, and the endpoint separation of a spring.

## Validation rules

Validate `dt > 0`, `G >= 0`, `softening >= 0`, `mass > 0`, finite position, and finite velocity.

Use status/result returns for normal validation failures. Do not throw for normal validation failures in M1-M5.
