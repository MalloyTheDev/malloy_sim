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
deferred to the 3D milestone (`docs/04`, `CLAUDE.md` rule 6,
`docs/decisions/0009-three-dimensions-are-the-destination.md`). The scalar forms
here are those quantities in two dimensions, not a different model, which is why
this is a deferral rather than a decision to be undone later.

### Amendment at M20: the deferral came due, and held

M20 built `RigidBody3D` beside `RigidBody2D`, so the claim above can be checked
rather than asserted. It held in both halves.

The scalar angle became a quaternion and the scalar inertia became three
principal moments, and neither 2D form had to be undone: `RigidBody2D` is
unchanged, and every 2D template still runs. Quaternion normalization did
arrive with the 3D body, exactly where this ADR said it would, and it turned out
to cost nothing physically: scaling a quaternion does not change the rotation it
represents, so renormalizing restores a magnitude rather than correcting a pose.

The inertia is the more interesting half. A general 3x3 tensor is still NOT
here, and M20 did not build one. Three principal moments are not a
simplification of a tensor: a symmetric tensor is always diagonalizable, so
storing the diagonal is a choice of axes rather than a restriction. A general
tensor is needed once bodies are built from composed shapes and the
parallel-axis step moves inertia off the principal axes, which arrives with 3D
mass properties and 3D contacts. Until then it would be three extra zeros.

So the scope boundary above is now the boundary between `RigidBody2D` and
`RigidBody3D` rather than between this project and a later one, and one clause
of it, the inertia tensor, is deferred again on its own terms.

### Amendment at M25: the tensor deferral came due, on those terms

The M20 amendment said a general tensor "becomes necessary only once bodies are
built from composed shapes and the parallel-axis step moves inertia off the
principal axes." M25 built exactly that: `mass_properties_3d` composes solid
spheres, shifts each by the parallel-axis theorem, and forms the compound
body's inertia tensor, which is not diagonal in general. `math::Mat3` and a
symmetric eigensolver arrived to hold and diagonalize it.

But the prediction that `RigidBody3D` need not store a tensor held. The tensor
is diagonalized at construction into principal moments and an orientation, which
is what the body already carries, so the body's dynamics are unchanged and no
`RigidBody3D` gained a 3x3. The tensor is a construction intermediate, not a
stored state: exactly the split ADR 0009 predicted.

### Amendment at M26: the assembly step, step 2 finished

Step 2 of the responsibility split above is "body mass-property construction,
including the parallel-axis theorem". M25 did it for solid spheres. M26 finished
it. `SolidBox` is the second mass primitive, and the first that is not
isotropic: a box has three distinct moments and a real orientation, so a single
tilted box already produces an inertia tensor that is not diagonal in the lab
frame. That is what forces the general machinery, rather than leaving it
exercised only by multi-part compounds as it was at M25.

`combine` is the assembly the split named but neither M25 nor this ADR had yet
written: it merges two mass-property sets into one, reconstructing each part's
tensor from its stored moments and orientation (`math::to_mat3`, the inverse of
M25's `to_quat`), shifting each to the shared centre by the parallel-axis
theorem, summing, and diagonalizing again. Because the reconstruction and the
shift are explicit, a mistake in either fails against the M25 list path, which
computes the same compound without ever calling `combine`, exactly the
"an error in one cannot masquerade as an error in another" separation this ADR
asked for. The prediction still holds: nothing about M26 made a body store a
tensor. The result of `combine` is diagonalized, and what drops into a
`RigidBody3D` is still three moments and an orientation.

### Amendment at M27: the general shape, step 1 finished too

Step 1 of the split, "shape-local centroid and second moment", was written for
2D polygons. M27 is its 3D form: `SolidMesh` computes the mass properties of a
solid bounded by a closed triangle mesh, by the divergence theorem, so any shape
a mesh can describe now has an inertia tensor. The sphere (M25) and the box
(M26) become special cases of it, which is the honest reading of "the general
one": a box mesh reproduces the `SolidBox` result down to the full tensor, and
that equality is a test.

The split still holds at the seam it was drawn for. The mesh integrals are pure
geometry and density (step 1 and step 2), producing a tensor about the centre of
mass; that tensor is diagonalized by the same `finalize` the sphere and box use,
into the three moments and orientation the body stores (step 3's input). No mesh
is kept on a `RigidBody3D`, no tensor is stored, and the diagonalization is still
a construction step rather than state. The prediction from M20, restated at M25
and M26, has now survived the general case: three principal moments and an
orientation are enough, whatever the shape.

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

## Amendment: the cost of accumulating the angle by addition

The decision above integrates the angle as `theta_next = theta + omega*dt` and
leaves canonicalization out. That still stands, and the code, the header
contract and the tests all agree with it. What was never written down is the
price.

`malloy_time` computes elapsed time as `tick_count * dt` specifically so it
"cannot drift the way repeated floating-point addition would". The angle is the
same pattern and does not get the same treatment, and it cannot: angular
velocity is changed by contacts, so there is no constant increment to multiply.

The size of it: summing N terms accumulates at most `u * theta * N / 2` with
`u = 2^-53`. Measured against exactly summed values at omega = 1.4, dt = 0.004
and N = 1e7, the drift is -1.13e-05 rad, about 2.33 arcsec, against a bound of
3.1e-05. The bound is tight to within 3x, and a test now pins it.

Stagnation is not reachable. `theta += c` stops advancing at
`theta > c * 2^53`, needing more than 9e15 steps against a parser cap of 2.1e9.

The alternatives are a per-body accumulated tick count or compensated
summation, both of which add state to `RigidBody2D`. Neither is taken. The
drift is bounded, derivable and pinned, which is the honest resolution for a
quantity nothing in the project reads at that precision (issue #18).
