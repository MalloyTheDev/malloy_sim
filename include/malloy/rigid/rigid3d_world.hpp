#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/rigid/rigid_body3d.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::rigid
{
// Torque-free rotation of rigid bodies in three dimensions.
//
// The smallest honest slice of 3D rigid-body dynamics, chosen the same way M19
// chose N-body: this needs no contacts, no collision geometry, no forces and no
// solver, so none of those had to be built speculatively to reach the physics
// that is actually new.
//
// What IS new here cannot be reached from the 2D code by replacing types. In
// two dimensions the inertia is a scalar and the angular velocity lies along a
// fixed axis, so a body's inertia never changes in world space and
// `omega x (I omega)` is identically zero. In three dimensions it is not, and
// the consequence is that a body with three distinct principal moments TUMBLES
// under no torque at all (ADR 0009).
class Rigid3DWorld
{
public:
    Rigid3DWorld(sim_core::SimulationSettings simulation_settings,
                 std::vector<RigidBody3D> bodies);

    // Ok, or InvalidSettings (dt) / InvalidState (a body).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step:
    //
    //   1. integrate the angular velocity with Euler's equations, in the body
    //      frame, where the inertia is diagonal:
    //
    //        Ix wx' = (Iy - Iz) wy wz
    //        Iy wy' = (Iz - Ix) wz wx
    //        Iz wz' = (Ix - Iy) wx wy
    //
    //      Every right-hand side is a DIFFERENCE of principal moments, which is
    //      why a sphere (all three equal) can never precess and why a body with
    //      three distinct moments can tumble.
    //
    //   2. integrate the orientation with the UPDATED angular velocity, which
    //      keeps this semi-implicit Euler like every other world here:
    //
    //        q' = q + (dt/2) * q * (0, omega_body)
    //
    //      The body-frame form puts the quaternion on the left. Using the world
    //      frame would put it on the right, and mixing the two is a silent way
    //      to rotate backwards.
    //
    //   3. RENORMALIZE the orientation. Integrating a unit quaternion does not
    //      keep it unit: the step above is a first-order move along a curve
    //      that lies on the unit sphere, so it leaves tangentially. The size of
    //      that is exact rather than approximate, because the quaternion norm
    //      is multiplicative and the update is a right multiplication:
    //
    //        |q'| = |q| |(1, (dt/2) omega)| = |q| sqrt(1 + |omega|^2 dt^2 / 4)
    //
    //      Unlike the positional correction in the 2D contact code, this costs
    //      nothing physically. Scaling a quaternion does not change the
    //      rotation it represents, so renormalizing moves no body and perturbs
    //      no conserved quantity; it restores the magnitude `rotate` assumes,
    //      which would otherwise scale every rotated vector by |q|^2, and it
    //      stops the factor above from compounding over a long run.
    //
    //   4. move the centre of mass by velocity * dt. Nothing acts on the
    //      translation here, so it is constant, and a body flies straight while
    //      it tumbles.
    //
    // On validation failure, before or after the update, the state is left
    // unchanged and the failing status is returned; this never throws.
    sim_core::StepResult step();

    const std::vector<RigidBody3D>& bodies() const { return bodies_; }
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
    std::vector<RigidBody3D> bodies_;
    std::vector<RigidBody3D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Diagnostics ---

math::Vec3 total_linear_momentum3d(const std::vector<RigidBody3D>& bodies);

// Sum of the world-frame angular momenta.
//
// Conserved by the continuum equations, because nothing applies a torque, but
// NOT by the discrete ones, and the size of the difference is exact rather than
// merely bounded. Writing u = L x omega, so that the update is
// omega' = omega + dt I^-1 u, both L . u and omega . u vanish identically and
// the first-order terms cancel, leaving
//
//   |L'|^2 = |L|^2 + dt^2 |u|^2
//
// per step. That is an equality, it is asserted every step in the tests, and it
// explains the cases that ARE exact: when omega lies along a principal axis, or
// when all three moments are equal, L is parallel to omega, u is zero, and
// nothing drifts at all.
math::Vec3 total_angular_momentum3d(const std::vector<RigidBody3D>& bodies);

// Translational plus rotational.
//
// The translational part is exactly constant, since no force acts. The
// rotational part drifts for the same reason and by the same second-order
// mechanism as the angular momentum above:
//
//   T' = T + dt^2 (I^-1 u) . u / 2
//
// Both corrections are positive, so this is a slow gain rather than a loss.
math::Real total_kinetic_energy3d(const std::vector<RigidBody3D>& bodies);
} // namespace malloy::rigid
