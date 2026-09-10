#pragma once

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

    // Valid when restitution is in [0, 1] and finite, and gravity is finite.
    bool is_valid() const;
};
} // namespace malloy::rigid
