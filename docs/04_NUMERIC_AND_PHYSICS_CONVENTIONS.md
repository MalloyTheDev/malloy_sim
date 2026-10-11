# 04 - Numeric and Physics Conventions

> These conventions are in force in every domain of the shipped code, through
> M39. Where a rule belongs to one domain (the softening formula to
> `malloy_nbody` and `malloy_charges`, the gas definitions to
> `malloy_particles`), its section says so.

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

`Vec2`, `Vec3`, `Quat` and `Mat3` all store `Real`. They are separate concrete
types, not templates on dimension: `cross` is a scalar in two dimensions and a
vector in three, so a template would need specializing exactly where it matters
(ADR 0009, amendment at M19).

## Floating point rules

- Use `double` for positions and velocities.
- Do not use `float` for simulation state.
- Do not enable fast-math. MSVC builds use `/fp:precise`, and GCC and Clang
  builds use `-ffp-contract=off`, so no compiler fuses `a*b+c` into an FMA and
  changes a result (`CMakeLists.txt`).
- Do not use `==` for floating-point simulation comparisons.
- Use explicit tolerances.

## Determinism scope

Do not claim cross-platform bitwise determinism.

V1 determinism means same binary, same platform, same initial state, same fixed timestep, same update order, same seed (M39, ADR 0010), and same result within expected tolerance.

Update order is part of each world's documented step contract, because
floating-point addition is not associative: particle and body pairs in
ascending `(i, j)` order, springs in declaration order (ADR 0008), particle
walls particle by particle with x before y (before z), contacts before ground
planes, and random draws in a fixed order and count (ADR 0010). Changing any of
them changes results and is a deliberate act.

## Units

Every domain uses normalized units (ADR 0005), and Boltzmann's constant is 1
(M38). There is no unit system and no SI conversion; do not introduce real
astronomical or SI magnitudes into a template without a milestone that needs
them.

## Integration and contacts

Every world advances with semi-implicit (symplectic) Euler: velocities first,
then positions from the updated velocities. The specializations:

- **Charges** (2D and 3D): the electric kick, then an EXACT rotation of the
  velocity by the magnetic field, then the move. A magnetic field alone
  preserves every speed to rounding error, with no error of any order in dt.
  Known defect: this kick-then-rotate (Lie) splitting makes motion in crossed E
  and B fields drift against the
  electric force, gaining energy without bound at about `(qE)^2 dt / (2m)` per
  unit time; no shipped template uses an electric field, so none shows it
  ([#42](https://github.com/MalloyTheDev/malloy_sim/issues/42)).
- **3D rigid bodies**: the angular velocity advances by FORWARD Euler on Euler's
  equations in the body frame, the orientation from the updated angular
  velocity, and the quaternion is renormalized every step.
- **2D rigid bodies**: the angle accumulates by addition and is not wrapped;
  its drift is bounded and pinned (ADR 0007).

Gravity is an acceleration, not a force, so it does not scale with mass. The 3D
rigid domain's applied force (M28) and torque (M21) are single constant
settings applied as forcing, not force or torque accumulators (`CLAUDE.md`
rule 5). Only `malloy_springs` accumulates forces, through its own network
(ADR 0008).

Contacts are impulses, resolved after the position update, in a fixed order and
in ONE pass with no iteration: a positional correction along the normal, a
normal impulse that reverses the closing speed scaled by the restitution, and
then a Coulomb friction impulse computed from the velocity that remains,
clamped to the friction coefficient times the normal impulse. A box contact
manifold (box on ground, box against box) is reduced to one contact at its
centroid with the deepest penetration.

Three fixed-step consequences follow from that design. The first is derived and
documented; the other two are open defects.

- **Micro-bounce.** There is no restitution threshold, so a body resting on the
  ground keeps a vertical speed of `e g dt / (1 + e)`: each step's gravity is
  bounced back at restitution `e` (2.67e-3 in `box_drop.scn`).
- **Static friction creeps.** Gravity moves a body before friction acts, so a
  rigid body on a slope inside the friction cone still creeps down it, at a steady
  `g sin(theta) dt` per unit time
  ([#44](https://github.com/MalloyTheDev/malloy_sim/issues/44)).
- **Boxes do not settle.** The centroid reduction and a translation-only
  positional correction together make a box that lands anything but flat rock
  indefinitely and gain energy, even at restitution 0, in 2D and 3D
  ([#41](https://github.com/MalloyTheDev/malloy_sim/issues/41)). The symmetric
  box-drop templates do not show it. The real fix is a per-point iterative
  contact solver, a milestone candidate
  ([#73](https://github.com/MalloyTheDev/malloy_sim/issues/73)).

## Kinetic theory (M38)

A box of colliding particles is read as a gas (`include/malloy/particles/thermo.hpp`).
Each definition below is the one the code computes and the templates check.

- **Boltzmann's constant is 1**, in keeping with the normalized units, so a
  temperature is an energy.
- **Temperature** is the mean kinetic energy per degree of freedom,
  `T = 2 KE / (d N)` with `d` = 2 or 3. It is measured in the lab frame, so a
  gas drifting as a whole counts its drift as heat; a gas in a box has none once
  the walls have absorbed it.
- **Pressure** is the total NORMAL impulse the walls have delivered since the run
  began, divided by the elapsed time and by the total length of wall (the
  perimeter `2 (W + H)`) or, in 3D, the total wall area. It is a time average
  over the whole run, what a gauge covering every wall would read. Before any
  time has passed it does not exist (an empty value, not 0). Friction's
  tangential impulse is along a wall, not into it, and is not counted.
- **The ideal gas** is the point-particle limit: particles of radius 0 never
  collide, and obey `P A = N T` in 2D and `P V = N T` in 3D. Two effects bound
  the measured error: the wall clamp lands a particle ON the wall and
  wastes part of a step, lowering each bounce rate by at most `v dt / L`; and a
  finite run counts whole bounces, so each particle's tally is within one bounce
  of its mean.
- **Finite radii raise the pressure** by the excluded-volume factor `Z` of the
  bulk gas (Henderson's equation for hard disks, Carnahan-Starling for hard
  spheres). In a box this small two wall effects matter as much: a centre cannot
  come within `r` of a wall, so the gas fills the box shrunk by `r` on every
  side and is denser against the walls than in the bulk (first-order adsorption
  `(2/3) rho^2 sigma^3` per unit length in 2D, `pi rho^2 sigma^4 / 4` per unit
  area in 3D, with `sigma = 2r`); and particles touch only the part of each wall
  `r` in from its ends, while the gauge spreads that force over the whole wall.
  `hard_disk_gas.scn` works the numbers through.
- **`speed_ratio`** is `<v^4> / <v^2>^2`: exactly 1 when every speed is equal,
  2 for a 2D Maxwell-Boltzmann gas, 5/3 in 3D. **`slow_fraction`** is the share
  of particles with `v^2` below the mean: `1 - 1/e = 0.632` in 2D, 0.608 in 3D.
  Both are unweighted by mass, so those targets are for a gas of equal masses.
  `slow_fraction` is decided by rounding when every speed is equal, so it is read
  after relaxation, and `speed_ratio` describes the start.
- **A deterministic start.** A gas starts from `lattice`: a grid with every
  speed equal and directions in exactly opposite pairs, pair `k` at `k` golden
  angles (2D) or on two independent equidistributed sequences (3D). Collisions
  do the randomizing.

## Thermal walls and randomness (M39)

A wall at temperature `T_w > 0` is THERMAL; 0 is the ordinary specular wall.

- **Emission.** A thermal wall fully accommodates: it forgets how a particle
  arrived and re-emits it from the flux of a Maxwell-Boltzmann gas at `T_w`
  leaving it. The normal speed is flux-weighted, a Rayleigh law
  `sqrt(-2 (T_w/m) ln u)`, and each tangential component is Gaussian with
  variance `T_w/m` (Box-Muller). A half-Gaussian normal law would leave a
  colliding gas at `d/(d+1)` of `T_w` (collisionless points in a 2D box at
  about half of it). The particle stays on the wall, and restitution and
  friction do not apply there.
- **Tallies.** The impulse tally records the true normal momentum change,
  `m (|v_n,in| + v_n,out)`. The heat tally records the kinetic energy after
  minus before, positive into the gas. It records only thermal re-emission, so
  a specular wall's entry is 0 by definition, not because the wall exchanges no
  energy: at restitution below 1 a specular wall absorbs energy that no tally
  shows. With restitution 1, no friction and no gravity, the gas's kinetic
  energy changes by exactly the summed heat; pair and wall friction dissipate
  even at restitution 1, so the identity needs all three conditions.
- **The generator.** SplitMix64 (`include/malloy/particles/splitmix64.hpp`),
  one stream per world, seeded from the `seed` setting. Deviates are strictly
  inside (0, 1), so no logarithm is ever infinite. A failed step restores the
  stream. Draws come only at thermal bounces, three per bounce, by particle
  index then wall axis. Order and count are observable behaviour (ADR 0010).
- **Determinism.** The integer stream is the same everywhere. The doubles derived
  from it pass through `log`, `cos` and `sin`, so bitwise agreement across
  platforms is still not claimed.
- **Two fixed-step biases.**
  - A wall clamp wastes part of a step per bounce, which raises a point gas's
    temperature by about `dt / (2 d E[tau])`, with `E[tau]` the mean flight
    time between walls. That is 0.4 per cent in a unit box at `dt = 0.01`.
  - The same lag adds `dt` to each round trip in free-molecular heat flow. The
    closed form between a hot and a cold wall, per particle, is
    `q = ((d+1)/2) (T_h - T_c) / (L sqrt(pi m / 2) (T_h^-1/2 + T_c^-1/2) + dt)`.
- **Choosing dt.** An emitted normal speed can reach `8.57 sqrt(T_w/m)` (the
  smallest deviate is 2^-53), and the whole emitted speed, normal and
  tangential together, about `12.1 sqrt(T_w/m)`. Keep that times `dt` well
  under a particle diameter.

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

Two neighbouring windows are still open. A rigid body's mass and inertia are
not bounded the same way, so a huge one passes validation while every energy
prints `inf` ([#51](https://github.com/MalloyTheDev/malloy_sim/issues/51)); and
a particle world's wall tallies can overflow while its particles stay valid,
so the step is not rolled back
([#50](https://github.com/MalloyTheDev/malloy_sim/issues/50)).

## Validation rules

Every world has a `validate()` that returns a `sim_core::StepStatus`: `Ok`,
`InvalidSettings` (a setting is unusable: `dt <= 0`, a restitution outside
[0, 1], a negative friction, `G` or softening, a wall temperature that is
negative or not finite, a timestep past a spring's stability limit) or
`InvalidState` (a body is unusable: a non-positive mass, a position or velocity
that is not squarable). Each world's header lists exactly what its
`validate()` checks.

`step()` validates before and after the update. On a failure it restores the
pre-step state (in the particle worlds that includes the wall tallies and the
random stream) and returns the failing status. Nothing throws for a normal
validation failure, and the scenario parser returns a status rather than
throwing too. The parser checks only syntax and resource caps; semantic rules
live in exactly one place, the world's own `validate()`.

Infinite mass and inertia are allowed in one place: a 2D rigid body, where
they mean immovable, per quantity (`docs/03`, `malloy_rigid`). The 3D rigid
domain has no
immovable bodies, only ground planes
([#74](https://github.com/MalloyTheDev/malloy_sim/issues/74)).
