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
    //      frame, where the inertia is diagonal, using FORWARD Euler (the rate
    //      depends on omega alone, so there is no older state to evaluate it
    //      at):
    //
    //        Ix wx' = (Iy - Iz) wy wz
    //        Iy wy' = (Iz - Ix) wz wx
    //        Iz wz' = (Ix - Iy) wx wy
    //
    //      Every right-hand side is a DIFFERENCE of principal moments, which is
    //      why a sphere (all three equal) can never precess and why a body with
    //      three distinct moments can tumble.
    //
    //   2. integrate the orientation with the UPDATED angular velocity. This
    //      coupling is what makes the world semi-implicit in the sense every
    //      other world here uses, velocities before positions, and it is the
    //      only place the word applies: step 1 above is forward Euler.
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
    //      nothing physically. A nonzero quaternion and any scaling of it
    //      represent the same rotation, because q v q^-1 is invariant under
    //      q -> k q, so renormalizing moves no body and perturbs no conserved
    //      quantity.
    //
    //      It is required all the same, because `rotate` uses the optimized
    //      unit-only form and that form is NOT scale invariant: fed a
    //      quaternion of norm k it returns p + k^2 (R p - p), scaling the
    //      displacement instead of rotating. Renormalizing also stops the
    //      factor above from compounding over a long run.
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

// Total angular momentum about the world origin: the SPIN of each body plus
// its ORBITAL term m (r x v), exactly as the 2D `total_angular_momentum` does.
// Dropping either makes a whole class of motion look conserved when it is not.
//
// The orbital term is exactly constant here, since no force acts: for a free
// body d(r x p)/dt = v x mv = 0. That is specific to CONSTANT velocity, and to
// the bit level as well: with v unchanging, r_{n+1} x v = (r_n + dt v) x v =
// r_n x v every step. It stops holding the moment a force changes v within a
// step, so a later milestone that adds forces cannot keep this line and assume
// the orbital term still conserves itself.
//
// --- The drift law, and what it does NOT describe ---
//
// Each body's SPIN angular momentum obeys an exact discrete law. Writing
// u = L x omega, so the update is omega' = omega + dt I^-1 u, both L . u and
// omega . u vanish identically, the first-order terms cancel, and
//
//   |L'|^2 = |L|^2 + dt^2 |u|^2
//
// survives, per step and per body. It is an equality, it is asserted every step
// in the tests, and it comes from the FORWARD-EULER angular velocity update:
// the drift is that method's second-order error, not a property of the
// orientation integration.
//
// It is zero exactly when u = 0, that is when I omega is parallel to omega:
// when the angular velocity lies in an EIGENSPACE of the inertia. For three
// distinct moments that means a principal axis. For a sphere every direction
// qualifies, and for an axisymmetric body (two equal moments) so does every
// direction in the degenerate plane, not only the two principal axes in it.
//
// This function's MAGNITUDE does not inherit that law, and must not be read as
// if it did. It is the magnitude of a vector SUM, and
//
//   |A + B|^2 = |A|^2 + |B|^2 + 2 A . B
//
// so it moves when the world-frame vectors turn relative to each other, which
// they do at O(dt) through the orientation integration. It can and does
// DECREASE over a run in which every per-body |L| is strictly increasing.
// `total_kinetic_energy3d` is different: it is a sum of scalars, so it does
// inherit the law term for term.
math::Vec3 total_angular_momentum3d(const std::vector<RigidBody3D>& bodies);

// Translational plus rotational.
//
// The translational part is exactly constant, since no force acts. The
// rotational part drifts by the same second-order mechanism as the spin angular
// momentum above:
//
//   T' = T + dt^2 (I^-1 u) . u / 2
//
// Both corrections are positive, so this is a slow gain rather than a loss.
//
// Unlike `total_angular_momentum3d`, this figure DOES inherit the per-body law
// exactly, because it sums scalars rather than vectors: the system total moves
// by the sum of the per-body terms and can only rise.
math::Real total_kinetic_energy3d(const std::vector<RigidBody3D>& bodies);
} // namespace malloy::rigid
