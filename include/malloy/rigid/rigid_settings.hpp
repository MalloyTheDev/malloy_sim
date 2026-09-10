#pragma once

#include <vector>

#include <malloy/collide/shapes.hpp>
#include <malloy/math/real.hpp>
#include <malloy/math/vec2.hpp>

namespace malloy::rigid
{
// Parameters for a rigid-body simulation. The fixed timestep is not here: that
// is generic and lives in sim_core::SimulationSettings.
//
// A struct rather than more constructor parameters. M14 added restitution
// positionally and M15 would have added gravity beside it, which is exactly
// where a transposed pair of reals stops being a compile error.
struct RigidSettings
{
    // Bounciness of every contact. 1 is perfectly elastic and conserves kinetic
    // energy; 0 is perfectly inelastic, so the pair stops separating along the
    // contact normal.
    math::Real restitution{1.0};

    // Uniform acceleration applied to every non-static body. Defaults to zero,
    // so a world written before this existed behaves exactly as it did.
    //
    // Gravity acts through the centre of mass, so it generates no torque on its
    // own. A body still ends up rotating under it, because a CONTACT away from
    // the centre of mass does have a moment arm: that is M14's machinery being
    // switched on rather than anything new here.
    math::Vec2 gravity{};

    // Coulomb friction coefficient for every contact. 0 is frictionless, which
    // is how the project behaved before M17 and remains the default.
    //
    // The tangential impulse is clamped to `friction` times the normal impulse,
    // which is what makes it self-limiting: no normal impulse means no friction
    // at all, automatically, so a body in mid-air cannot be accelerated
    // sideways by it.
    //
    // Not capped at 1. A coefficient above 1 is physically real (rubber on
    // rubber), and clamping it would silently change a caller's model.
    math::Real friction{0.0};

    // Immovable ground planes: floors, walls, ramps. Each resolves exactly as a
    // body of infinite mass and inertia would, with one difference that is the
    // whole reason the primitive exists. A floor built from discs has a contact
    // normal that TURNS as a body rolls across it, by as much as 14 degrees for
    // the floor in dropped_bodies.scn, because the normal points at whichever
    // disc centre is nearest. A plane's normal is its own and never turns.
    //
    // Empty by default, so a world written before this existed behaves exactly
    // as it did.
    std::vector<collide::Halfplane> ground;

    // Valid when restitution is in [0, 1] and finite, friction is non-negative
    // and finite, gravity is finite, and every ground plane is itself valid.
    bool is_valid() const;
};
} // namespace malloy::rigid
