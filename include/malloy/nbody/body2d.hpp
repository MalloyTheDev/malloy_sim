#pragma once

#include <malloy/math/vec2.hpp>

namespace malloy::nbody
{
// A point mass in 2D. This is pure simulation state: position, velocity and
// mass. Acceleration is a per-step derived quantity owned by the integrator
// (NBodyWorld), not stored here.
struct Body2D
{
    math::Vec2 position{};
    math::Vec2 velocity{};
    math::Real mass{1.0};

    // Valid when mass is strictly positive and every value is finite.
    bool is_valid() const;
};
} // namespace malloy::nbody
