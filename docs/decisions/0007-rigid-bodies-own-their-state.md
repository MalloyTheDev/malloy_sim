# ADR 0007 - Rigid Bodies Own Their State

## Status

Accepted. Recorded ahead of M11, which had not started as of this writing.

## Context

M11 adds 2D rigid bodies. The obvious shortcut is to widen `nbody::Body2D` with
an angle, angular velocity, and a moment of inertia, since it already carries
position, velocity, and mass.

## Decision

M11 introduces a dedicated `RigidBody2D` in its own domain library.
`nbody::Body2D` and `particles::Particle2D` are left exactly as they are.

`RigidBody2D` owns:

```text
position / linear velocity
angle / angular velocity
mass / inverse mass
inertia / inverse inertia
local center of mass
```

It does not inherit from either existing body type. Composition, or duplicating
a few trivially small fields, is preferred over an inheritance hierarchy across
two different physical models.

## Rationale

A translational body and a rigid body do not differ by a few extra fields. They
obey different equations of motion and expose different invariants. Widening
`Body2D` would give every consumer rotational semantics whether it needs them or
not, and gravity has no use for an angle.

This follows the precedent M10 set when `Particle2D` was introduced rather than
adding a radius to `Body2D`: each domain owns the state its model requires
(ADR 0006).

The distinction that justifies the separate type: a rigid body owns a pose and a
mass distribution, not a point with an angle attached.

## Responsibility split

Three responsibilities stay separable, so that an error in one cannot masquerade
as an error in another:

1. shape-local centroid and second moment of area (`malloy_collide`: pure
   geometry, no density and no mass, testable against closed-form values);
2. body mass-property construction, including the parallel-axis theorem
   `I = I_com + m * d^2` when the body origin is not the centroid;
3. rigid-body integration.

The parallel-axis step is explicit rather than something that emerges from shape
code, so a centroid mistake fails a geometry test rather than a rotational one.

## Angle representation

Integration is `theta_next = theta + omega * dt`. Canonicalization is a separate
concern: wrapping is useful for a bounded public representation and for
trigonometric evaluation, but forcing it into the state transition introduces
discontinuities and complicates comparison. If canonical angles are wanted they
are to be a documented invariant with one owner, not an `fmod` called wherever
convenient.

## Scope boundary

2D only. Orientation is a scalar angle and inertia is a scalar that is
rotation-invariant, so there are no quaternions and no inertia tensors at M11.
Quaternion normalization and inertia tensor transforms are 3D concerns and stay
deferred to M19 with the rest of 3D (`docs/04`, `CLAUDE.md` rule 6).

## Testing consequence

The first rotational tests must break symmetry on every axis at once: body
origin away from the center of mass, application point away from the center of
mass, a force with both components nonzero, a nonzero initial angle, an inertia
that is not 1, and both signs of torque.

The reason is specific. Torque in 2D is `tau = rx*Fy - ry*Fx`, with `r` measured
from the center of mass. If the origin is the center of mass, an implementation
that wrongly measures from the body origin passes. If `ry` or `Fx` is zero, half
the determinant vanishes and another class of mistake passes
(`docs/05`, "Prefer asymmetric configurations").

## Numerical note

`-ffp-contract=off` settles the compiler-contract question but not conditioning.
`rx*Fy - ry*Fx` is a difference of products, so when the two products are close
the relative error in a small torque can be large even though both products are
individually accurate.

No compensated arithmetic is to be added preemptively. Establish the scale
assumptions first and measure whether the loss matters in the supported numeric
range. Near-cancellation cases belong in numerical-characterization tests that
check sign stability, bounded magnitude, and finite results, not in golden
values pinning one bit pattern, which would assert the cross-toolchain bitwise
determinism `docs/04` declines to claim.
