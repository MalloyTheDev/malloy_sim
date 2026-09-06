#pragma once

#include <malloy/collide/shapes.hpp>
#include <malloy/math/real.hpp>

namespace malloy::particles
{
// Parameters for a colliding-particle simulation. The fixed timestep is not
// here: that is generic and lives in sim_core::SimulationSettings
// (docs/03_MODULE_BOUNDARIES.md).
struct ParticleSettings
{
    // Bounciness of every collision, both particle/particle and particle/wall.
    // 1 is perfectly elastic and conserves kinetic energy; 0 is perfectly
    // inelastic, so the pair stops separating along the contact normal.
    math::Real restitution{1.0};

    // The box the particles are confined to. Walls are immovable, which is why
    // they carry momentum away: only the no-wall case conserves it.
    collide::Aabb bounds{math::Vec2{-1.0, -1.0}, math::Vec2{1.0, 1.0}};

    // Valid when restitution is in [0, 1] and finite, and bounds is a valid box.
    bool is_valid() const;
};
} // namespace malloy::particles
