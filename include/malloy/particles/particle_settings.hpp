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

    // Valid when restitution is in [0, 1] and finite, friction is non-negative
    // and finite, bounds is a valid box, and gravity is finite.
    bool is_valid() const;
};
} // namespace malloy::particles
