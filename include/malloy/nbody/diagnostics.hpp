#pragma once

#include <malloy/math/real.hpp>
#include <malloy/nbody/body2d.hpp>

namespace malloy::nbody
{
// Specific orbital energy of the two-body relative motion (energy per unit
// reduced mass): v_rel^2 / 2 - G * (m_a + m_b) / r, where r and v_rel are the
// separation and relative speed of the two bodies.
//
// The exact dynamics conserve this, so it is a useful eyeball check on how
// well the integrator behaves. Requires r > 0 (the two bodies not coincident).
math::Real specific_orbital_energy(const Body2D& a, const Body2D& b, math::Real g);
} // namespace malloy::nbody
