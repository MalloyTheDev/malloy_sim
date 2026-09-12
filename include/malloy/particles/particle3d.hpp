#pragma once

#include <malloy/math/vec3.hpp>

namespace malloy::particles
{
// A spherical point mass that can collide, the 3D sibling of Particle2D.
// Position and velocity gain a third component; mass and radius are unchanged.
//
// Its own concrete type rather than a widened nbody::Body3D, on the same
// reasoning M10 used in 2D: gravity has no use for a radius (ADR 0006). It
// lives beside the 2D particle because the dimension is not a domain (ADR 0009).
struct Particle3D
{
    math::Vec3 position{};
    math::Vec3 velocity{};
    math::Real mass{1.0};
    math::Real radius{0.5};

    // Valid when mass is strictly positive, radius is non-negative, and
    // position and velocity are SQUARABLE (docs/04). A zero radius is a point
    // particle, which can still collide with the walls.
    bool is_valid() const;
};
} // namespace malloy::particles
