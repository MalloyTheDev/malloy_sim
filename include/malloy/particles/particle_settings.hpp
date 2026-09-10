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

    // Uniform acceleration applied to every particle, in world units. Defaults
    // to zero, so every scenario written before this existed behaves exactly as
    // it did. A setting rather than per-body state: it is a property of the
    // world, not of a particle, so no body type changes (ADR 0006).
    math::Vec2 gravity{};

    // Valid when restitution is in [0, 1] and finite, bounds is a valid box,
    // and gravity is finite.
    bool is_valid() const;
};
} // namespace malloy::particles
