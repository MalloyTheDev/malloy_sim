#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/charges/charge_settings.hpp>
#include <malloy/charges/charged_particle2d.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::charges
{
// A concrete electrostatics and magnetostatics simulation.
//
// Concrete and derived from nothing, like every other world here (ADR 0006).
// It shares malloy::math and the sim_core vocabulary with the other domains
// and knows about none of them.
class ChargeWorld
{
public:
    ChargeWorld(sim_core::SimulationSettings simulation_settings,
                ChargeSettings charge_settings,
                std::vector<ChargedParticle2D> particles);

    const ChargeSettings& charge_settings() const { return charge_settings_; }

    // Ok, or InvalidSettings (dt, k, either field, softening) / InvalidState
    // (a particle).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step:
    //
    //   1. add HALF the ELECTRIC acceleration times dt to every velocity, which
    //      is the pairwise Coulomb term plus the uniform field, evaluated once
    //      from the positions at the start of the step;
    //   2. ROTATE every velocity by the magnetic angle -(q*b/m)*dt;
    //   3. add the other half of the same electric kick;
    //   4. move every particle by the updated velocity * dt.
    //
    // A particle that nothing turns (b == 0, or no charge) gets the whole kick
    // at once instead of two halves: the same thing in exact arithmetic, and
    // it keeps such a run bit-identical to plain semi-implicit Euler.
    // Velocities before positions, so this is the same semi-implicit Euler as
    // every other domain (docs/04), with the turn sandwiched in the kick.
    //
    // Step 2 is a rotation rather than another kick, and that is a physical
    // statement rather than an optimisation. A magnetic force is always
    // perpendicular to the velocity, so it does no work and cannot change a
    // particle's speed. Applying it as `v += (q/m)(v x B) dt` breaks that: the
    // added vector is perpendicular to v, so the two are the legs of a right
    // triangle and the new speed is |v| * sqrt(1 + (q*b*dt/m)^2), which grows
    // every step without bound. A cyclotron orbit integrated that way spirals
    // outward forever.
    //
    // In two dimensions the magnetic sub-problem has an exact solution, because
    // a uniform out-of-plane field rotates the velocity at a constant rate and
    // nothing else. So it is integrated exactly, and the speed is preserved to
    // rounding error (a few ulp per step) rather than to some order in dt.
    //
    // The electric and magnetic parts are still applied in sequence rather
    // than together (operator splitting), and the ORDER of the sequence is
    // what decides the motion in crossed fields. Kick-then-rotate, used before
    // issue #42, has a velocity fixed point equal to the E x B drift turned by
    // half a step's angle, so the guiding centre also crept against qE at
    // |qE/m| dt / 2 and the energy grew by (qE)^2 dt / (2m) per unit time,
    // without bound. The symmetric order here (the order of the Boris pusher,
    // with the rotation exact) has its fixed point exactly perpendicular to E:
    // the continuum drift scaled by (theta/2) cot(theta/2), theta = -q b dt/m,
    // the same for either sign of charge, and with no creep along the field.
    //
    // On validation failure, before or after the update, the state is left
    // unchanged and the failing status is returned; this never throws.
    sim_core::StepResult step();

    // The acceleration on each particle from the pairwise Coulomb interaction
    // plus the uniform electric field. Exposed for testing, and it excludes
    // the magnetic term because that one is not an acceleration applied for a
    // time, it is a rotation.
    //
    // A pair whose coupling k q_i q_j is exactly zero (k == 0, or a neutral
    // particle) contributes nothing at ANY separation: it is skipped before
    // 1/r^3 is formed, which overflows to +inf below about 1.8e-103 apart
    // (issue #49). So does a coincident pair with no softening.
    std::vector<math::Vec2> compute_electric_accelerations() const;

    const std::vector<ChargedParticle2D>& particles() const { return particles_; }
    std::uint64_t tick_count() const { return step_ ? step_->tick_count() : 0; }

    // Accumulated simulation time, computed as ticks * dt by FixedStep.
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
    ChargeSettings charge_settings_;
    std::vector<ChargedParticle2D> particles_;
    std::vector<ChargedParticle2D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Diagnostics ---

math::Vec2 total_momentum(const std::vector<ChargedParticle2D>& particles);

math::Real total_kinetic_energy(const std::vector<ChargedParticle2D>& particles);

// Coulomb potential energy of every pair, plus the potential of the uniform
// electric field.
//
// The pair term is `+k q_i q_j / r`, and its SIGN carries the physics: it is
// positive for like charges, which is a configuration that wants to fly apart,
// and negative for unlike ones. Gravity has no equivalent, since `-G m_i m_j/r`
// is unconditionally negative.
//
// Uses the same softening as the integrator, so total_energy is a quantity the
// dynamics conserve rather than one they merely approach.
//
// The magnetic field contributes nothing, and that is not an omission: a
// magnetic force does no work, so there is no potential to associate with it.
math::Real total_potential_energy(const std::vector<ChargedParticle2D>& particles,
                                  const ChargeSettings& settings);

// Kinetic plus potential.
//
// Conserved to rounding when the magnetic field is the only thing acting,
// because the rotation is exact and changes no speed at all.
//
// In a uniform electric field across the magnetic one it is not conserved
// but it is BOUNDED, for the whole run: for a charge released from rest the
// error is exactly -m |a|^2 (dt/2)^2 cot(theta/2) sin(n theta) after n steps,
// with a = qE/m and theta = -q b dt/m, an oscillation at the cyclotron
// frequency of amplitude about |q| |E|^2 dt / (2 b) and no secular term
// (issue #42; the order before it gained energy without bound).
//
// With a uniform electric field and nothing turning the particle it drifts,
// as a projectile under gravity does in every other domain: semi-implicit
// Euler under a constant force loses m |a|^2 dt^2 / 2 a step. With pairwise
// forces it is not conserved either, for the same reason and by the same
// mechanism as every other domain here.
math::Real total_energy(const std::vector<ChargedParticle2D>& particles,
                        const ChargeSettings& settings);
} // namespace malloy::charges
