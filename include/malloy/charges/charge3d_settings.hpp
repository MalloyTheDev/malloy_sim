#pragma once

#include <malloy/math/real.hpp>
#include <malloy/math/vec3.hpp>

namespace malloy::charges
{
// Parameters for a 3D charged-particle simulation, the 3D sibling of
// ChargeSettings. The fixed timestep is not here: that is generic and lives in
// sim_core::SimulationSettings.
struct Charge3DSettings
{
    // Coulomb constant, normalized rather than SI (ADR 0005). May be zero,
    // which switches the pairwise interaction off and leaves only the fields.
    math::Real k{1.0};

    // Uniform electric field. A charge feels q*E, so it accelerates as
    // (q/m)*E: a heavier charge accelerates less, and an opposite sign the
    // other way, which is the whole difference from gravity.
    math::Vec3 electric{};

    // Uniform magnetic field, now a full VECTOR rather than the 2D scalar.
    //
    // In two dimensions only the out-of-plane component of B produced an
    // in-plane force, so a scalar sufficed. In three dimensions the field has a
    // direction, and q(v x B) turns the velocity about that direction while
    // leaving the component ALONG it untouched. That is what makes the motion
    // helical rather than circular, and it is the physics that has no 2D form
    // (docs/decisions/0009-three-dimensions-are-the-destination.md).
    math::Vec3 magnetic{};

    // Softening, entering the denominator SQUARED as in the N-body domain
    // (docs/04), so a coincident pair produces a large finite force rather than
    // infinity.
    math::Real softening{0.0};

    // Valid when k is finite, the electric and magnetic fields are SQUARABLE
    // (the magnetic length squares its components, and the electric field's
    // square enters the free-flight energy drift), the softening is
    // non-negative and finite, and the SQUARE of the softening is finite.
    //
    // The squared tests are not redundant with finiteness. A finite value above
    // sqrt(DBL_MAX), about 1.34e154, squares to infinity, which for the field
    // would make a magnetic turn NaN and for the softening would make every
    // separation infinite: a simulation in which nothing happens conserves
    // everything perfectly, the most convincing wrong answer there is.
    bool is_valid() const;
};
} // namespace malloy::charges
