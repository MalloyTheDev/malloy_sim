#pragma once

#include <malloy/math/vec2.hpp>

namespace malloy::collide
{
// A circle in 2D. A zero radius is allowed and means a point.
//
// This module deliberately knows nothing about Body2D or any simulation type:
// it works on shapes, the way malloy_ascii works on points
// (docs/03_MODULE_BOUNDARIES.md). Pairing a body with a shape is the caller's
// job, which keeps collision geometry reusable and keeps malloy_nbody free of
// collision concerns.
struct Circle
{
    math::Vec2 center{};
    math::Real radius{0.0};

    // Valid when the radius is non-negative and every value is finite.
    bool is_valid() const;
};

// An axis-aligned box, stored as its lower and upper corners.
struct Aabb
{
    math::Vec2 min{};
    math::Vec2 max{};

    // Valid when min <= max on both axes and every value is finite. A zero-area
    // box (min == max) is allowed and means a point.
    bool is_valid() const;
};

// An oriented box in 2D (M36): a centre, half-widths along its OWN two axes, and
// an orientation angle in radians. The 2D sibling of `Box3`. Its local axes are
// (cos, sin) and (-sin, cos), so a zero angle is an `Aabb` of the same extents.
// It is what lets a 2D rigid body collide as a box rather than a disc.
struct Obb2
{
    math::Vec2 center{};
    math::Vec2 half_extents{};
    math::Real orientation{0.0};

    // Valid when both half-extents are strictly positive and finite, and the
    // centre and orientation are finite. A zero or negative extent is not a box,
    // matching `Box3`.
    bool is_valid() const;
};

// A halfplane: everything on one side of an infinite straight line. Used for
// ground and walls, where a disc is the wrong shape and an Aabb is a lie
// (a floor is not 40 units thick, and a body that tunnels past its far face
// should not pop out the bottom).
//
// The line is the set of points p with dot(normal, p) == offset, so the signed
// distance from p to it is dot(normal, p) - offset. `normal` points OUT of the
// solid side, into free space, so that distance is positive for a point in
// free space and negative for one inside the solid. A floor at y = -2 is
// normal (0, 1) with offset -2; a left wall at x = 3 that keeps bodies to its
// left is normal (-1, 0) with offset -3.
//
// Unlike every other pair in this module, a circle against a halfplane has NO
// degenerate case: the normal is the plane's own and never has to be inferred
// from the relative position of two centres, so none of the documented
// fallbacks in contact.hpp apply to it. That is the point of the primitive.
struct Halfplane
{
    math::Vec2 normal{0.0, 1.0};
    math::Real offset{0.0};

    // Valid when the normal is finite, the offset is finite, and the normal is
    // a UNIT vector.
    //
    // Unit length is required rather than normalized on use. Normalizing
    // silently would make (0, 0) and (1e-300, 0) both look acceptable while
    // meaning different things, and the signed distance above is only a
    // distance at all when the normal is unit. Callers that accept arbitrary
    // input, such as the scenario loader, normalize once at the boundary and
    // reject what cannot be normalized.
    bool is_valid() const;
};

// --- Area properties ---
//
// There is deliberately no area, centroid or second moment for a Halfplane.
// All three are infinite, and returning 0 the way the invalid-shape path does
// would be indistinguishable from a real answer.
//
// Pure geometry: no density and no mass appear here, so these are testable
// against closed-form values with nothing physical involved. Converting them
// into mass properties belongs to whichever domain owns bodies
// (docs/decisions/0007-rigid-bodies-own-their-state.md).
//
// All three return 0 for an invalid shape rather than a meaningless number.

math::Real area(const Circle& c);
math::Real area(const Aabb& box);

// The centroid, which for both of these shapes is their geometric center.
math::Vec2 centroid(const Circle& c);
math::Vec2 centroid(const Aabb& box);

// Polar second moment of area about the centroid: the integral of r^2 dA, with
// r measured from the centroid. For a disc of radius R this is pi*R^4/2, and
// for a w by h box it is w*h*(w^2 + h^2)/12.
//
// It is "about the centroid" on purpose. Moving it to another reference point
// is the parallel-axis step, and that belongs with the body, not the shape.
math::Real second_moment_of_area(const Circle& c);
math::Real second_moment_of_area(const Aabb& box);
} // namespace malloy::collide
