#pragma once

#include <malloy/math/vec3.hpp>

namespace malloy::charges
{
// A point charge in three dimensions: the 3D sibling of ChargedParticle2D,
// exactly as Body3D is to Body2D. The scalar mass and signed charge are
// unchanged; only position and velocity gain a third component.
//
// It lives beside the 2D type rather than replacing it, on the same reasoning
// M19 and M20 used: the dimension is not a domain (ADR 0009), so the 2D form
// stays supported. Charge is SIGNED, the one thing no other body model here
// carries, so this is still its own concrete type (ADR 0006), not a widened
// Body3D.
struct ChargedParticle3D
{
    math::Vec3 position{};
    math::Vec3 velocity{};
    math::Real mass{1.0};

    // Signed, and zero is meaningful: a neutral particle is carried by the
    // fields of others while contributing nothing to them. Like charges repel
    // and unlike charges attract, which is the sign of `charge_a * charge_b`.
    math::Real charge{0.0};

    // Valid when mass is strictly positive and finite, the charge is finite
    // (not necessarily positive or nonzero, the whole difference from gravity),
    // and position and velocity are SQUARABLE (docs/04): a magnetic turn takes
    // |v|, and the free-flight energy drift squares the electric acceleration.
    bool is_valid() const;
};
} // namespace malloy::charges
