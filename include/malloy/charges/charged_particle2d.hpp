#pragma once

#include <malloy/math/vec2.hpp>

namespace malloy::charges
{
// A point charge: a point mass that also carries an electric charge.
//
// Deliberately not `nbody::Body2D` and not `particles::Particle2D`. Charge is
// SIGNED, and no other body model in this project has a signed quantity that
// determines whether an interaction attracts or repels. Widening one of those
// types to carry it would make gravity and collisions answer questions they
// have no use for. Each domain owning its own concrete state is the point of
// ADR 0006.
struct ChargedParticle2D
{
    math::Vec2 position{};
    math::Vec2 velocity{};
    math::Real mass{1.0};

    // Signed, and zero is meaningful: a neutral particle is carried by the
    // fields of others while contributing nothing to them. Like charges repel
    // and unlike charges attract, which is the sign of `charge_a * charge_b`.
    math::Real charge{0.0};

    // Valid when mass is strictly positive and every value is finite.
    //
    // The charge is NOT required to be positive or nonzero. That is the whole
    // difference from gravity, where mass is strictly positive and the force
    // is therefore unconditionally attractive.
    bool is_valid() const;
};
} // namespace malloy::charges
