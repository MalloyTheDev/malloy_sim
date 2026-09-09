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

// --- Area properties ---
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
