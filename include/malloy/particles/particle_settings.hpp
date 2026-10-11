#pragma once

#include <cstdint>

#include <malloy/collide/shapes.hpp>
#include <malloy/math/real.hpp>

namespace malloy::particles
{
// The temperature of each wall of the box (M39), named as the impulse tally
// names them. A wall at T > 0 is THERMAL: it forgets how a particle arrived and
// re-emits it as a wall at temperature T would, so heat flows between the wall
// and the gas (docs/04). 0, the default, is the ordinary specular wall, so a
// box with no thermal wall behaves exactly as it did before M39 and draws no
// random number. Valid when every entry is non-negative and finite.
struct WallTemperature2D
{
    math::Real min_x{0.0};
    math::Real max_x{0.0};
    math::Real min_y{0.0};
    math::Real max_y{0.0};

    bool is_valid() const;
};

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

    // Coulomb friction for every contact, both particle/particle and
    // particle/wall: a tangential impulse clamped to `friction` times the normal
    // impulse, which damps sliding rather than reversing it. A particle carries
    // no orientation, so this imparts no spin (unlike the rigid domain, where the
    // same clamp has a lever arm). Defaults to zero, and it is the LAST field so
    // that a `{restitution, bounds}` or `{restitution, bounds, gravity}` brace
    // still compiles and a scenario written before this existed behaves
    // bit-for-bit as it did: with no friction a contact only ever touched the
    // normal direction. It is a contact impulse consumed on the spot, not a
    // persistent force, so no body type changes and no accumulator is added.
    math::Real friction{0.0};

    // M39: which walls are thermal, and at what temperature. At a thermal wall
    // restitution and friction do not apply: the wall fully accommodates, so
    // what leaves it does not depend on what arrived. They still act on pair
    // contacts and on the specular walls. Added after `friction` so that every
    // earlier brace still compiles and behaves as before.
    WallTemperature2D wall_temperature{};

    // M39: the seed of the world's random number generator, which only a
    // thermal wall draws from, so without one it changes nothing. Every value is
    // valid. The LAST field, for the same reason as wall_temperature.
    std::uint64_t seed{0};

    // Valid when restitution is in [0, 1] and finite, friction is non-negative
    // and finite, bounds is a valid box, gravity is finite, and every wall
    // temperature is non-negative and finite.
    bool is_valid() const;
};
} // namespace malloy::particles
