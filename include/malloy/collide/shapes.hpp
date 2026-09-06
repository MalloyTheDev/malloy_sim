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
} // namespace malloy::collide
