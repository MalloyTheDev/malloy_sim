#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec3.hpp>
#include <malloy/particles/particle3d.hpp>
#include <malloy/particles/particle3d_settings.hpp>
#include <malloy/particles/splitmix64.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::particles
{
// The total normal impulse each of the six walls has delivered since the run
// began (M38), the 3D sibling of WallImpulse2D: non-negative, the normal
// momentum change only (m (1 + e) |v_n| at a specular wall, m (|v_n,in| +
// v_n,out) at a thermal one), and divided by the elapsed time and the wall area
// it is the pressure.
struct WallImpulse3D
{
    math::Real min_x{0.0};
    math::Real max_x{0.0};
    math::Real min_y{0.0};
    math::Real max_y{0.0};
    math::Real min_z{0.0};
    math::Real max_z{0.0};
};

// The total heat each of the six walls has delivered to the gas (M39), the 3D
// sibling of WallHeat2D: kinetic energy after minus before, at every bounce off
// a thermal wall, positive into the gas; a specular wall's entry stays 0 by
// definition. The energy identity needs restitution 1, no friction and no
// gravity, as in 2D.
struct WallHeat3D
{
    math::Real min_x{0.0};
    math::Real max_x{0.0};
    math::Real min_y{0.0};
    math::Real max_y{0.0};
    math::Real min_z{0.0};
    math::Real max_z{0.0};
};

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
    //   4. clamp every particle back inside the box, particle by particle in
    //      ascending index, each particle's x walls, then y, then z.
    //
    // A THERMAL wall (M39) re-emits instead of reflecting, drawing exactly three
    // uniforms per bounce in a fixed order: the normal speed, then the two
    // tangential components (thermo.hpp), which go to the wall's two tangent
    // axes in ascending order (an x wall sets y then z, a y wall x then z, a z
    // wall x then y). No other step draws.
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

    // The impulse each wall has delivered so far; read-only bookkeeping that
    // changes no trajectory and is rolled back with a rolled-back step.
    const WallImpulse3D& wall_impulse() const { return wall_impulse_; }

    // The heat each wall has delivered to the gas so far (M39), rolled back
    // with a rolled-back step.
    const WallHeat3D& wall_heat() const { return wall_heat_; }

    // The generator's state, seeded from ParticleSettings3D::seed; it moves by
    // exactly three draws per thermal bounce and is restored by a rollback.
    std::uint64_t random_state() const { return random_.state(); }

private:
    sim_core::SimulationSettings simulation_settings_;
    ParticleSettings3D particle_settings_;
    std::vector<Particle3D> particles_;
    std::vector<Particle3D> previous_;
    WallImpulse3D wall_impulse_{};
    WallHeat3D wall_heat_{};
    SplitMix64 random_;
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
