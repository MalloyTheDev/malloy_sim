#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec2.hpp>
#include <malloy/particles/particle2d.hpp>
#include <malloy/particles/particle_settings.hpp>
#include <malloy/particles/splitmix64.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::particles
{
// The total normal impulse each wall of the box has delivered since the run
// began (M38): the momentum it took to turn particles around, summed over every
// bounce. A wall is immovable, so this is also the momentum the walls carried
// away, which is what makes it measurable: divided by the elapsed time and the
// length of wall it is the pressure the gas exerts (thermo.hpp).
//
// Every entry is non-negative and counts only the NORMAL momentum change of
// each bounce: m (1 + e) times the speed into a specular wall, and at a thermal
// wall (M39) m times the speed in plus the speed out, since what leaves a
// thermal wall does not depend on what arrived. Tangential impulse, friction's
// or a thermal wall's, pushes along a wall, not into it, so it is not pressure
// and is not tallied.
struct WallImpulse2D
{
    math::Real min_x{0.0};
    math::Real max_x{0.0};
    math::Real min_y{0.0};
    math::Real max_y{0.0};
};

// The total heat each wall has delivered to the gas since the run began (M39):
// at every bounce off a THERMAL wall, the particle's kinetic energy after the
// re-emission minus its kinetic energy before. Positive into the gas, so a hot
// wall's entry grows and a cold wall's falls. A specular wall exchanges no
// heat and its entry stays exactly 0. With restitution 1 and no gravity the
// gas's kinetic energy changes by exactly the sum of these, to rounding. Like
// the impulse tally it is a running sum and is never reset, so a very long run
// loses precision in it before it could overflow.
struct WallHeat2D
{
    math::Real min_x{0.0};
    math::Real max_x{0.0};
    math::Real min_y{0.0};
    math::Real max_y{0.0};
};

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
    //   4. resolve wall contacts, particle by particle in ascending index, each
    //      particle's x walls before its y walls.
    //
    // A THERMAL wall (M39) re-emits a particle instead of reflecting it, drawing
    // three uniforms from the world's generator per bounce, in a fixed order:
    // one for the normal speed, then two for the tangential component
    // (thermo.hpp). That order and count, with the step order above, make the
    // stream a function of the run, so the same seed repeats it (docs/04). No
    // other step draws, and a box without a thermal wall never draws at all.
    //
    // Step 4 clamps a particle that has passed a wall back to the surface,
    // which changes its position without changing its velocity. In a gravity
    // field that MOVES ENERGY: shifting a particle by `delta` against gravity
    // adds `m |g| delta` of potential energy, once per clamping event.
    //
    // That term is real, first order in dt, and comparable to the one the
    // project had already characterized. Free flight loses exactly
    // `(1/2) m |g|^2 dt^2` per step; at dt = 4e-3 an elastic ball bouncing on a
    // floor gains +0.94 from clamping against a free-flight loss of -1.6, so
    // the correction is 59 per cent of the other term and opposite in sign.
    //
    // The clamp distance is bounded by how far the particle overshot in one
    // step, so `0 <= delta < |v| dt` and the whole term vanishes as dt does.
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

    // The impulse each wall has delivered so far. Tallying it reads the bounce
    // without changing it, so every trajectory is exactly what it was before
    // this existed. A step rolled back for non-finite state rolls this back too.
    const WallImpulse2D& wall_impulse() const { return wall_impulse_; }

    // The heat each wall has delivered to the gas so far (M39); rolled back with
    // a rolled-back step, as the impulse tally is.
    const WallHeat2D& wall_heat() const { return wall_heat_; }

    // The state of the world's random number generator, seeded from
    // ParticleSettings::seed. It moves only when a thermal wall draws, by
    // exactly three draws per thermal bounce, and a rolled-back step restores
    // it, so a failed step leaves no trace in the stream.
    std::uint64_t random_state() const { return random_.state(); }

private:
    sim_core::SimulationSettings simulation_settings_;
    ParticleSettings particle_settings_;
    std::vector<Particle2D> particles_;
    // Scratch copy of the pre-step state, used only to roll back a step that
    // produced non-finite values. A member so the common case reuses its
    // allocation rather than allocating once per step.
    std::vector<Particle2D> previous_;
    WallImpulse2D wall_impulse_{};
    WallHeat2D wall_heat_{};
    SplitMix64 random_;
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

// Total angular momentum about the origin (a scalar in 2D):
// sum of m (x_x v_y - x_y v_x).
//
// This domain went without one until issue #20. It runs the same position-only
// correction that malloy_rigid does, and that correction perturbs angular
// momentum by c x (v_b - v_a) without touching a velocity, so the artefact was
// present here too and there was no way to measure it.
//
// Not conserved in the presence of walls, which is not a defect: an immovable
// wall is an external agent, exactly as an immovable body is in malloy_rigid.
math::Real total_angular_momentum(const std::vector<Particle2D>& particles);

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
