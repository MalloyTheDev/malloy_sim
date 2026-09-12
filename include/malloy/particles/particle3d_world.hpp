#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec3.hpp>
#include <malloy/particles/particle3d.hpp>
#include <malloy/particles/particle3d_settings.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::particles
{
// A concrete 3D colliding-particle simulation, the 3D sibling of ParticleWorld.
// Concrete and derived from nothing (ADR 0006). Particles collide as spheres
// (the M22/M24 sphere geometry) and are confined to a 3D box.
class ParticleWorld3D
{
public:
    ParticleWorld3D(sim_core::SimulationSettings simulation_settings,
                    ParticleSettings3D particle_settings,
                    std::vector<Particle3D> particles);

    // Ok, or InvalidSettings (dt/restitution/bounds) / InvalidState (a
    // particle, or one too large to fit inside the bounds).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step, in the same order as the 2D world:
    //
    //   1. add gravity * dt to every velocity;
    //   2. move every particle by the UPDATED velocity * dt (semi-implicit
    //      Euler);
    //   3. resolve particle/particle contacts in ascending pair order;
    //   4. clamp every particle back inside the box.
    //
    // The contact response is a positional correction plus an equal-and-opposite
    // impulse along the line of centres, so it imparts no spin (a particle has
    // no orientation) and conserves momentum. The wall clamp is a reflection,
    // not an exchange, so walls carry momentum away. The pair order is fixed so
    // a run repeats (docs/04). On validation failure the state is left unchanged
    // and the failing status returned; never throws.
    sim_core::StepResult step();

    const std::vector<Particle3D>& particles() const { return particles_; }
    std::uint64_t tick_count() const { return step_ ? step_->tick_count() : 0; }
    math::Real elapsed_time() const
    {
        return step_ ? step_->elapsed_time() : math::Real{0};
    }
    const sim_core::SimulationSettings& simulation_settings() const
    {
        return simulation_settings_;
    }
    const ParticleSettings3D& particle_settings() const { return particle_settings_; }

private:
    sim_core::SimulationSettings simulation_settings_;
    ParticleSettings3D particle_settings_;
    std::vector<Particle3D> particles_;
    std::vector<Particle3D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Whole-system diagnostics ---
//
// Momentum is conserved exactly by particle/particle collisions at any
// restitution (equal and opposite impulses); walls change it, so momentum is
// conserved only for a run where nothing touches a wall.
math::Vec3 total_momentum3d(const std::vector<Particle3D>& particles);

// Kinetic energy is conserved only when restitution is exactly 1 AND gravity is
// zero; below that every collision removes some, which is what makes
// restitution testable.
math::Real total_kinetic_energy3d(const std::vector<Particle3D>& particles);

// Total angular momentum about the world origin, sum of m (r x v), a VECTOR in
// three dimensions (the 2D version was a scalar). Not conserved with walls, the
// same as in 2D: an immovable wall is an external agent.
math::Vec3 total_angular_momentum3d(const std::vector<Particle3D>& particles);

// Gravitational potential energy in the uniform field, sum of -m (g . r). Zero
// when gravity is zero.
math::Real total_potential_energy3d(const std::vector<Particle3D>& particles,
                                    const math::Vec3& gravity);

// Kinetic plus potential. NOT conserved under a constant field: semi-implicit
// Euler sheds (1/2)(sum m)|g|^2 dt^2 per free-flight step, the same secular
// drift the 2D world documents.
math::Real total_energy3d(const std::vector<Particle3D>& particles,
                          const math::Vec3& gravity);
} // namespace malloy::particles
