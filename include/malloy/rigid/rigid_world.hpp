#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec2.hpp>
#include <malloy/rigid/rigid_body2d.hpp>
#include <malloy/rigid/rigid_settings.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::rigid
{
// A concrete 2D rigid-body simulation: free motion plus impulses.
//
// Concrete and derived from nothing, like every other world here (ADR 0004,
// ADR 0006). It shares only malloy_math, malloy_time and the sim_core
// vocabulary with the other domains, and knows nothing about them.
//
// M11 gave it free motion plus impulses. M14 added contact response: bodies
// with a nonzero radius collide as discs and the resulting impulse generates
// torque, so they tumble rather than merely bouncing.
//
// Still no persistent forces and no force or torque accumulators. A body with
// no contacts translates and spins forever, which remains the cleanest
// invariant here: linear and angular momentum conserved exactly.
class RigidWorld
{
public:
    RigidWorld(sim_core::SimulationSettings simulation_settings,
               std::vector<RigidBody2D> bodies, RigidSettings settings = {});

    const RigidSettings& settings() const { return settings_; }

    // Ok, or InvalidSettings (dt, restitution, gravity) / InvalidState (a body).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step:
    //
    //   1. add gravity * dt to every velocity whose body has finite mass;
    //   2. move the centre of mass by the UPDATED velocity * dt, which keeps
    //      the integration semi-implicit (symplectic) Euler;
    //   3. advance the angle by angular_velocity * dt;
    //   4. resolve body-against-body contacts in ascending pair order;
    //   5. resolve every body against every ground plane, in ascending body
    //      then plane order.
    //
    // Ground comes last on purpose: a body pushed by a pair impulse in step 4
    // gets the floor's answer in the same step rather than the next one.
    //
    // Step 2 integrates EVERY body, including immovable ones. That is
    // deliberate: an immovable body with a velocity you set is a kinematic
    // body, a moving platform that pushes everything and is pushed by nothing.
    // Gravity does not accelerate it, no impulse can slow it, and it carries no
    // momentum in the diagnostics, so it acts as an external agent rather than
    // part of the closed system. Leave its velocity at zero for a wall.
    //
    // The body origin follows from the centre of mass and the new angle, so a
    // body whose origin is offset from its centre of mass orbits correctly
    // rather than rotating about the wrong point.
    //
    // The angle is NOT wrapped here. Canonicalisation is a separate concern
    // with one owner (ADR 0007); forcing it into the state transition would
    // introduce discontinuities and complicate comparison.
    //
    // On validation failure, before or after the update, the state is left
    // unchanged and the failing status is returned. Never throws.
    sim_core::StepResult step();

    // Apply an impulse at a world point. Changes linear velocity by J/m always,
    // and angular velocity by (r x J)/I where r runs from the CENTRE OF MASS to
    // the application point. An impulse through the centre of mass therefore
    // leaves the spin untouched, which is the cleanest test of the torque arm.
    //
    // Returns false and changes nothing when the index is out of range, the
    // body is invalid, or the impulse or point is not finite.
    bool apply_impulse_at(std::size_t index, const math::Vec2& impulse,
                          const math::Vec2& world_point);

    const std::vector<RigidBody2D>& bodies() const { return bodies_; }
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
    RigidSettings settings_;
    std::vector<RigidBody2D> bodies_;
    std::vector<RigidBody2D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Whole-system diagnostics ---
//
// With no forces, no impulses and no contacts, all three are conserved exactly.
// Immovable bodies are skipped: infinite mass times zero velocity would be NaN,
// and an immovable body carries nothing to report.

// Contacts conserve linear momentum exactly: the impulse is equal and opposite,
// and the positional correction moves positions without touching a velocity.
math::Vec2 total_linear_momentum(const std::vector<RigidBody2D>& bodies);

// Angular momentum about the world origin: the spin term I*omega plus the
// orbital term m * (r x v), with r running from the origin to the centre of
// mass. Both terms are needed: a body translating past the origin carries
// angular momentum about it even with zero spin.
//
// The contact IMPULSE conserves this exactly, because it is equal and opposite
// and acts at one shared point, so the two contributions are -(p x J) and
// +(p x J). The positional CORRECTION does not: it moves positions without
// changing velocities, which shifts the orbital term by c x (v_b - v_a), where
// c is the correction. The perturbation is proportional to penetration depth,
// so it shrinks with the timestep rather than accumulating from the impulse.
//
// This is the usual cost of resolving penetration by moving bodies, and it is
// recorded here rather than hidden behind a loose tolerance.
math::Real total_angular_momentum(const std::vector<RigidBody2D>& bodies);

math::Real total_kinetic_energy(const std::vector<RigidBody2D>& bodies);

// Gravitational potential energy in the uniform field: the sum of -m (g . r),
// with r the centre of mass. Zero when gravity is zero. Static bodies are
// skipped, since infinite mass times a position is not a number.
math::Real total_potential_energy(const std::vector<RigidBody2D>& bodies,
                                  const math::Vec2& gravity);

// Kinetic plus potential.
//
// NOT conserved under a constant field. Semi-implicit Euler loses exactly
// (1/2) * (sum of m) * |g|^2 * dt^2 per free-flight step, the same constant
// secular drift derived for the particle domain in M12, because it is the same
// integrator on the same kind of field. Contacts remove more on top of that
// whenever restitution is below 1.
math::Real total_energy(const std::vector<RigidBody2D>& bodies,
                        const math::Vec2& gravity);
} // namespace malloy::rigid
