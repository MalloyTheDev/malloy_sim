#pragma once

#include <vector>

#include <malloy/math/vec2.hpp>
#include <malloy/nbody/body2d.hpp>

namespace malloy::nbody
{
// --- Two-body diagnostic ---

// Specific orbital energy of the two-body relative motion (energy per unit
// reduced mass): v_rel^2 / 2 - G * (m_a + m_b) / r, where r and v_rel are the
// separation and relative speed of the two bodies. Requires r > 0.
math::Real specific_orbital_energy(const Body2D& a, const Body2D& b, math::Real g);

// --- Whole-system diagnostics (any number of bodies) ---
//
// For an isolated N-body system these are conserved by the exact dynamics, so
// they are the natural way to judge how well the integrator behaves over a run.

// Total kinetic energy: sum of (1/2) m |v|^2.
math::Real total_kinetic_energy(const std::vector<Body2D>& bodies);

// Total potential energy with softening: sum over unique pairs of
// -G m_i m_j / sqrt(r^2 + softening^2). Uses the same softening as the
// integrator so that total_energy is a quantity the dynamics conserve.
math::Real total_potential_energy(const std::vector<Body2D>& bodies, math::Real g,
                                  math::Real softening);

// Total mechanical energy: kinetic + potential.
math::Real total_energy(const std::vector<Body2D>& bodies, math::Real g,
                        math::Real softening);

// Total linear momentum: sum of m v.
math::Vec2 total_momentum(const std::vector<Body2D>& bodies);

// Center of mass: (sum m x) / (sum m). Returns the origin if total mass is 0.
math::Vec2 center_of_mass(const std::vector<Body2D>& bodies);

// Total angular momentum about the origin (a scalar in 2D):
// sum of m (x_x v_y - x_y v_x).
math::Real total_angular_momentum(const std::vector<Body2D>& bodies);
} // namespace malloy::nbody
