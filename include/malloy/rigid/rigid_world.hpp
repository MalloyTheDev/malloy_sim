#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec2.hpp>
#include <malloy/rigid/rigid_body2d.hpp>
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
// M11 scope: no forces and no contact response. Nothing accelerates a body
// except an impulse the caller applies, so a body left alone translates and
// spins forever. That is not a placeholder, it is the invariant this milestone
// is built to demonstrate: with no impulse, linear and angular momentum are
// conserved exactly. Persistent forces arrive with ballistics, and contact
// response with rigid-body collision.
class RigidWorld
{
public:
    RigidWorld(sim_core::SimulationSettings simulation_settings,
               std::vector<RigidBody2D> bodies);

    // Ok, or InvalidSettings (dt) / InvalidState (a body).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step:
    //
    //   1. move the centre of mass by velocity * dt;
    //   2. advance the angle by angular_velocity * dt.
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
    std::vector<RigidBody2D> bodies_;
    std::vector<RigidBody2D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Whole-system diagnostics ---
//
// With no forces and no impulses applied, both are conserved exactly.

math::Vec2 total_linear_momentum(const std::vector<RigidBody2D>& bodies);

// Angular momentum about the world origin: the spin term I*omega plus the
// orbital term m * (r x v), with r running from the origin to the centre of
// mass. Both terms are needed: a body translating past the origin carries
// angular momentum about it even with zero spin.
math::Real total_angular_momentum(const std::vector<RigidBody2D>& bodies);

math::Real total_kinetic_energy(const std::vector<RigidBody2D>& bodies);
} // namespace malloy::rigid
