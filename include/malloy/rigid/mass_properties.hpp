#pragma once

#include <malloy/collide/shapes.hpp>
#include <malloy/math/vec2.hpp>

namespace malloy::rigid
{
// The mass distribution of a shape, ready to drop into a RigidBody2D.
//
// This is step 2 of the three-way split in ADR 0007. Step 1 is the pure
// geometry in malloy_collide (area, centroid, second moment of area, with no
// density anywhere), and step 3 is integration. Keeping them separate is what
// stops a centroid mistake from masquerading as a rotational-integration bug.
struct MassProperties
{
    math::Real mass{0.0};
    // About the centre of mass, matching RigidBody2D::inertia.
    math::Real inertia{0.0};
    // Centre of mass, in the same frame the shape was given in.
    math::Vec2 center_of_mass{};
};

// Mass properties of a uniform lamina: mass = density * area, and inertia =
// density * (polar second moment of area about the centroid).
//
// Returns all zeros for an invalid shape or a non-positive density, rather than
// a meaningless number.
MassProperties mass_properties(const collide::Circle& shape, math::Real density);
MassProperties mass_properties(const collide::Aabb& shape, math::Real density);

// Parallel-axis theorem: I = I_com + m * d^2.
//
// Explicit rather than something that emerges from shape code, so that moving a
// reference point is always a visible step (ADR 0007). `distance` is measured
// from the centre of mass to the new reference point. A negative or non-finite
// distance returns the input unchanged, since squaring would hide the error.
math::Real shift_inertia(math::Real inertia_about_com, math::Real mass,
                         math::Real distance);
} // namespace malloy::rigid
