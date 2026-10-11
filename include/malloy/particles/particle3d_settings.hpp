#pragma once

#include <cstdint>

#include <malloy/collide/shapes3d.hpp>
#include <malloy/math/real.hpp>
#include <malloy/math/vec3.hpp>

namespace malloy::particles
{
// The temperature of each of the six walls (M39), the 3D sibling of
// WallTemperature2D: a wall at T > 0 is thermal, 0 (the default) is the
// ordinary specular wall. Valid when every entry is non-negative and finite.
struct WallTemperature3D
{
    math::Real min_x{0.0};
    math::Real max_x{0.0};
    math::Real min_y{0.0};
    math::Real max_y{0.0};
    math::Real min_z{0.0};
    math::Real max_z{0.0};

    bool is_valid() const;
};

// Parameters for a 3D colliding-particle simulation, the 3D sibling of
// ParticleSettings. The fixed timestep is not here (it lives in
// sim_core::SimulationSettings).
struct ParticleSettings3D
{
    // Bounciness of every collision, particle/particle and particle/wall. 1 is
    // perfectly elastic; 0 is perfectly inelastic.
    math::Real restitution{1.0};

    // The box the particles are confined to, now a 3D `Aabb3`. Walls are
    // immovable, which is why they carry momentum away: only the no-wall case
    // conserves it.
    collide::Aabb3 bounds{math::Vec3{-1.0, -1.0, -1.0}, math::Vec3{1.0, 1.0, 1.0}};

    // Uniform acceleration applied to every particle. Defaults to zero. A
    // setting rather than per-body state (ADR 0006).
    math::Vec3 gravity{};

    // Coulomb friction for every contact, particle/particle and particle/wall:
    // a tangential impulse clamped to `friction` times the normal impulse, which
    // damps sliding rather than reversing it. A particle carries no orientation,
    // so this imparts no spin. Defaults to zero, and it is the LAST field so a
    // `{restitution, bounds}` or `{restitution, bounds, gravity}` brace still
    // compiles and a scenario written before it existed behaves bit-for-bit as
    // it did. Like restitution it is consumed on the spot, not a persistent
    // force (the 3D sibling of ParticleSettings).
    math::Real friction{0.0};

    // M39: which walls are thermal, and at what temperature; restitution and
    // friction do not apply at a thermal wall. As in 2D, added after `friction`.
    WallTemperature3D wall_temperature{};

    // M39: the seed of the world's random number generator, which only a
    // thermal wall draws from. Every value is valid. The LAST field.
    std::uint64_t seed{0};

    // Valid when restitution is in [0, 1] and finite, friction is non-negative
    // and finite, bounds is a valid box, gravity is finite, and every wall
    // temperature is non-negative and finite.
    bool is_valid() const;
};
} // namespace malloy::particles
