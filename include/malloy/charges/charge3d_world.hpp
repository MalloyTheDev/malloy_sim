#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/charges/charge3d_settings.hpp>
#include <malloy/charges/charged_particle3d.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::charges
{
// A concrete 3D electrostatics and magnetostatics simulation, the 3D sibling of
// ChargeWorld. Concrete and derived from nothing (ADR 0006); it shares
// malloy::math and the sim_core vocabulary and knows about no other domain.
class Charge3DWorld
{
public:
    Charge3DWorld(sim_core::SimulationSettings simulation_settings,
                  Charge3DSettings charge_settings,
                  std::vector<ChargedParticle3D> particles);

    const Charge3DSettings& charge_settings() const { return charge_settings_; }

    // Ok, or InvalidSettings (dt, k, either field, softening) / InvalidState
    // (a particle).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step:
    //
    //   1. add HALF the ELECTRIC acceleration times dt to every velocity, which
    //      is the pairwise Coulomb term plus the uniform field, evaluated once
    //      from the positions at the start of the step;
    //   2. ROTATE every velocity about the magnetic-field direction by the
    //      angle -(q |B| / m) dt;
    //   3. add the other half of the same electric kick;
    //   4. move every particle by the updated velocity * dt.
    //
    // A particle that nothing turns (B == 0, or no charge) gets the whole kick
    // at once, which keeps such a run bit-identical to plain semi-implicit
    // Euler. Velocities before positions, so this is the same semi-implicit
    // Euler as every other domain (docs/04), with the turn sandwiched in the
    // kick.
    //
    // Step 2 is the 3D form of the exact magnetic turn. A magnetic force is
    // always perpendicular to the velocity, so it does no work and cannot
    // change speed; applying it as `v += (q/m)(v x B) dt` would inflate the
    // speed by sqrt(1 + (q |B| dt / m)^2) every step. In three dimensions a
    // uniform field rotates the velocity about the field DIRECTION at a
    // constant rate, leaving the component along the field untouched, so the
    // exact update is a rotation about that axis. Implemented with a
    // quaternion, `rotate(from_axis_angle(B, -(q |B| / m) dt), v)`, so the
    // speed is preserved to rounding error (a few ulp per step), and the
    // component along B is carried exactly: circular motion across the field, uniform drift along it, which
    // together are a helix. This is the physics two dimensions cannot express.
    //
    // As in 2D this is operator splitting, and the symmetric order is what
    // keeps crossed fields right (issue #42): across B the fixed point of the
    // velocity map is exactly perpendicular to E, the continuum E x B / |B|^2
    // scaled by (theta/2) cot(theta/2) with theta = -(q |B| / m) dt, where
    // kick-then-rotate turned it by theta/2 and so crept against qE and gained
    // energy without bound. Along B the turn does nothing, and the two halves
    // add up to one kick, as they would with no field at all.
    //
    // On validation failure, before or after the update, the state is left
    // unchanged and the failing status is returned; this never throws.
    sim_core::StepResult step();

    // The acceleration on each particle from the pairwise Coulomb interaction
    // plus the uniform electric field. Exposed for testing, and it excludes the
    // magnetic term because that one is a rotation, not an acceleration applied
    // for a time. A pair whose coupling k q_i q_j is exactly zero contributes
    // nothing at any separation, as in 2D (issue #49).
    std::vector<math::Vec3> compute_electric_accelerations() const;

    const std::vector<ChargedParticle3D>& particles() const { return particles_; }
    std::uint64_t tick_count() const { return step_ ? step_->tick_count() : 0; }

    math::Real elapsed_time() const
    {
        return step_ ? step_->elapsed_time() : math::Real{0};
    }

    const sim_core::SimulationSettings& simulation_settings() const
    {
        return simulation_settings_;
    }

private:
    sim_core::SimulationSettings simulation_settings_;
    Charge3DSettings charge_settings_;
    std::vector<ChargedParticle3D> particles_;
    std::vector<ChargedParticle3D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Diagnostics ---

math::Vec3 total_momentum3d(const std::vector<ChargedParticle3D>& particles);

math::Real total_kinetic_energy3d(const std::vector<ChargedParticle3D>& particles);

// Coulomb potential energy of every pair (+k q_i q_j / r, positive for like
// charges, which want to fly apart) plus the potential of the uniform electric
// field (-q (E . r), the shape of gravity's -m (g . r)). The magnetic field
// contributes nothing: a magnetic force does no work, so there is no potential
// to associate with it. Uses the same softening as the integrator, so
// total_energy3d is a quantity the dynamics conserve rather than approach.
math::Real total_potential_energy3d(const std::vector<ChargedParticle3D>& particles,
                                    const Charge3DSettings& settings);

// Kinetic plus potential. Conserved to rounding when the magnetic field is the
// only thing acting, because the rotation is exact and changes no speed at all.
// With a uniform electric field PERPENDICULAR to B it is bounded for the whole
// run, an oscillation of amplitude about |q| |E|^2 dt / (2 |B|) for a charge
// released from rest, with no secular term (issue #42; see the 2D
// total_energy for the closed form). A component of E along B, or an electric
// field with nothing turning the particle, is the projectile case of every
// other domain: semi-implicit Euler under a constant force, which loses
// m |a|^2 dt^2 / 2 a step, with a the part of the acceleration that is not
// turned. Pairwise forces make it drift at the same order as every other
// domain here.
math::Real total_energy3d(const std::vector<ChargedParticle3D>& particles,
                          const Charge3DSettings& settings);
} // namespace malloy::charges
