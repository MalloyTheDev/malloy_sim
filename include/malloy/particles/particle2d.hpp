#pragma once

#include <malloy/math/vec2.hpp>

namespace malloy::particles
{
// A circular point mass that can collide.
//
// This is deliberately not nbody::Body2D. That type is documented as pure
// simulation state (position, velocity, mass) and gravity has no use for a
// radius, so rather than widening it, this domain carries its own body type.
// Each domain owning its own concrete state is the point of ADR 0006.
struct Particle2D
{
    math::Vec2 position{};
    math::Vec2 velocity{};
    math::Real mass{1.0};
    math::Real radius{0.5};

    // Valid when mass is strictly positive, radius is non-negative, and every
    // value is finite. A zero radius is allowed and means a point particle,
    // which can still collide with the walls.
    bool is_valid() const;
};
} // namespace malloy::particles
