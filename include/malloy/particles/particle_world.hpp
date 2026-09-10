#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec2.hpp>
#include <malloy/particles/particle2d.hpp>
#include <malloy/particles/particle_settings.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::particles
{
// A concrete 2D colliding-particle simulation.
//
// Like NBodyWorld it is concrete, owns its own fixed-step update, and derives
// from nothing (ADR 0004, ADR 0006). The two worlds share only malloy_math and
// the sim_core vocabulary; neither knows the other exists.
class ParticleWorld
{
public:
    ParticleWorld(sim_core::SimulationSettings simulation_settings,
                  ParticleSettings particle_settings,
                  std::vector<Particle2D> particles);

    // Check the settings and every particle. Returns Ok, or the first failure
    // category found: InvalidSettings (dt/restitution/bounds) or InvalidState
    // (a particle, or a particle too large to fit inside the bounds).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step, in this order:
    //
    //   1. add gravity * dt to every velocity;
    //   2. move every particle by the UPDATED velocity * dt, which makes this
    //      semi-implicit (symplectic) Euler, the same order NBodyWorld uses;
    //   3. resolve particle/particle contacts in ascending pair order;
    //   4. resolve wall contacts.
    //
    // The pair order is fixed so a run is repeatable (docs/04). On validation
    // failure, before or after the update, the state is left unchanged and the
    // failing status is returned; this never throws.
    sim_core::StepResult step();

    const std::vector<Particle2D>& particles() const { return particles_; }
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
    const ParticleSettings& particle_settings() const { return particle_settings_; }

private:
    sim_core::SimulationSettings simulation_settings_;
    ParticleSettings particle_settings_;
    std::vector<Particle2D> particles_;
    // Scratch copy of the pre-step state, used only to roll back a step that
    // produced non-finite values. A member so the common case reuses its
    // allocation rather than allocating once per step.
    std::vector<Particle2D> previous_;
    // Tick counting and elapsed time come from malloy_time rather than being
    // reimplemented per domain. Empty when dt is unusable.
    std::optional<time::FixedStep> step_;
};

// --- Whole-system diagnostics ---
//
// Momentum is conserved exactly by particle/particle collisions at any
// restitution, because the impulse is equal and opposite. Walls are immovable
// and therefore do change it, so momentum is only a conserved quantity for a
// run where nothing touches a wall.
math::Vec2 total_momentum(const std::vector<Particle2D>& particles);

// Kinetic energy is conserved only when restitution is exactly 1 AND gravity is
// zero. Below that, every collision removes some, which is what makes
// restitution testable.
math::Real total_kinetic_energy(const std::vector<Particle2D>& particles);

// Gravitational potential energy in the uniform field: the sum of -m (g . r).
// Zero when gravity is zero, so total_energy reduces to the kinetic term.
math::Real total_potential_energy(const std::vector<Particle2D>& particles,
                                  const math::Vec2& gravity);

// Kinetic plus potential.
//
// This is NOT conserved under a constant field. Semi-implicit Euler loses
// exactly (1/2) * (sum of m) * |g|^2 * dt^2 per free-flight step, a constant
// secular drift rather than a bounded oscillation, because the motion is
// unbounded. The exact per-step change is what the tests assert; claiming
// conservation here would be claiming something false.
math::Real total_energy(const std::vector<Particle2D>& particles,
                        const math::Vec2& gravity);
} // namespace malloy::particles
