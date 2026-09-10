#pragma once

#include <malloy/math/real.hpp>
#include <malloy/math/vec2.hpp>

namespace malloy::charges
{
// Parameters for a charged-particle simulation. The fixed timestep is not
// here: that is generic and lives in sim_core::SimulationSettings.
struct ChargeSettings
{
    // Coulomb constant. Normalized rather than SI, for the reason ADR 0005
    // gives for the gravitational constant: a demo written in coulombs and
    // metres spends its whole dynamic range on unit conversion.
    math::Real k{1.0};

    // Uniform electric field. A charge feels q*E, so it accelerates as
    // (q/m)*E.
    //
    // That is NOT how gravity behaves, and the difference is testable: gravity
    // is an acceleration and every body falls at the same rate, while a heavier
    // charge in the same electric field accelerates less, and a charge of the
    // opposite sign accelerates the other way.
    math::Vec2 electric{};

    // Uniform magnetic field, perpendicular to the plane.
    //
    // A scalar, not a vector, and that is a consequence of two dimensions
    // rather than a simplification: only the out-of-plane component of B
    // produces an in-plane force. With B = (0, 0, b) and v = (vx, vy, 0),
    // q(v x B) = q*b*(vy, -vx), which stays in the plane.
    //
    // In 3D this becomes a vector again, and the force stops being expressible
    // this way (docs/decisions/0009-three-dimensions-are-the-destination.md).
    math::Real magnetic{0.0};

    // Softening, in the same sense and with the same units as the N-body
    // domain (docs/04): it enters the denominator SQUARED, so a coincident
    // pair produces a large finite force rather than infinity.
    //
    // Charges need it more than masses do, because two like charges released
    // near each other accelerate apart without bound, and nothing in the model
    // stops them meeting in the first place.
    math::Real softening{0.0};

    // Valid when k, the electric field, the magnetic field and the softening
    // are all finite, and the softening is non-negative. k may be zero, which
    // switches the pairwise interaction off and leaves only the fields.
    bool is_valid() const;
};
} // namespace malloy::charges
