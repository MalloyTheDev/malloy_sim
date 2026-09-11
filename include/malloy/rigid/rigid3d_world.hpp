#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec3.hpp>
#include <malloy/rigid/rigid_body3d.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::rigid
{
// Settings for the 3D rigid world: a single constant torque, applied to every
// body each step.
//
// The torque is in the WORLD frame, not the body frame. That is a deliberate
// choice and it is what makes the milestone's headline invariant clean: the
// physical law dL/dt = torque holds in the world frame, so with a constant
// world torque each body's world-frame angular momentum grows along a straight
// line,
//
//   L(t) = L(0) + torque * t,
//
// exactly in the continuum. An external couple fixed in the lab (the kind a
// field produces) is world-frame; a thruster bolted to the body would be
// constant in the BODY frame instead, a different thing that is not reachable
// from this one without knowing the body, and is deferred (M21 keeps to the
// smallest step, as M20 did).
//
// The default is zero torque, which makes a `Rigid3DWorld` built without
// settings identical to the torque-free M20 world, bit for bit.
struct Rigid3DSettings
{
    math::Vec3 torque{};

    // Valid when the torque is SQUARABLE, not merely finite: it enters |u|^2 in
    // the drift law below, so a torque whose square overflows is refused up
    // front, the same bound positions and velocities carry (docs/04, and M18's
    // softening).
    bool is_valid() const { return math::is_squarable(torque); }
};

// Rotation of rigid bodies in three dimensions under a constant applied torque.
//
// M20 was the torque-free case, the smallest honest slice of 3D rigid-body
// dynamics. M21 adds the one thing that was missing to make something act on
// the rotation: a torque. It is still the smallest step that reaches new
// physics, chosen the way M19 and M20 were: no contacts, no collision geometry,
// no force on the translation and no solver, so none of those had to be built
// speculatively.
//
// What torque-free rotation could not show, this can. In M20 nothing pushed on
// a body, so a body spinning about a principal axis stayed there forever. Here
// a torque perpendicular to the spin makes the spin axis move: the gyroscopic
// response, which is the reason a spinning top does not simply fall over. In
// two dimensions there is no such thing, because there is only one rotation
// axis and a torque can only speed the spin up or slow it down.
class Rigid3DWorld
{
public:
    Rigid3DWorld(sim_core::SimulationSettings simulation_settings,
                 std::vector<RigidBody3D> bodies, Rigid3DSettings settings = {});

    // Ok, or InvalidSettings (dt or torque) / InvalidState (a body).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step:
    //
    //   1. integrate the angular velocity with Euler's equations, in the body
    //      frame, where the inertia is diagonal, using FORWARD Euler (the rate
    //      depends on the old angular velocity and orientation, so there is no
    //      older state to evaluate it at):
    //
    //        Ix wx' = (Iy - Iz) wy wz + Tx
    //        Iy wy' = (Iz - Ix) wz wx + Ty
    //        Iz wz' = (Ix - Iy) wx wy + Tz
    //
    //      where (Tx, Ty, Tz) is the applied torque expressed in the BODY
    //      frame. The setting is a WORLD-frame torque, so it is rotated in by
    //      the body's current orientation, T_body = R^-1 T_world, before this
    //      step. Euler's equations are diagonal only in the body frame, which
    //      is why the torque comes to them rather than the other way round.
    //
    //      The gyroscopic part is the DIFFERENCE of principal moments, exactly
    //      as in M20: a sphere still cannot precess, and it is nonzero only in
    //      three dimensions. The torque is what M21 adds to it.
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
    //      it tumbles: the torque acts only on the rotation.
    //
    // On validation failure, before or after the update, the state is left
    // unchanged and the failing status is returned; this never throws.
    sim_core::StepResult step();

    const std::vector<RigidBody3D>& bodies() const { return bodies_; }
    const Rigid3DSettings& settings() const { return settings_; }
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
    Rigid3DSettings settings_;
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
// the orbital term still conserves itself. A torque does not change it: M21
// acts only on the rotation.
//
// --- What the torque does to the SPIN, exactly ---
//
// Each body's SPIN angular momentum obeys an exact discrete law. Writing
// L = I omega and u = L x omega + T_body, so the update is
// omega' = omega + dt I^-1 u, the cross term is orthogonal to both L and omega
// and only the torque survives the first-order part:
//
//   |L'|^2 = |L|^2 + 2 dt (L . T_body) + dt^2 |u|^2
//   T'     = T     +     dt (omega . T_body) + dt^2 (I^-1 u) . u / 2
//
// Both are equalities, both are asserted every step in the tests, and both come
// from the FORWARD-EULER angular step. Note L . T_body = L_world . T_world and
// omega . T_body = omega_world . T_world, since a rotation preserves the dot
// product: the work the torque does is frame independent.
//
// With T = 0 this is the M20 law, |L'|^2 = |L|^2 + dt^2 |u|^2, and the
// first-order term is gone, which is why the drift there is second order and a
// principal-axis spin holds exactly. Under a torque the first-order term is
// present and the WORLD-frame angular momentum grows along the straight line
// L(t) = L(0) + T_world t, which is the exact continuum law the discrete scheme
// tracks to first order (and, for a spin-up about the torque axis, exactly).
//
// This function's MAGNITUDE still does not inherit the per-body law: it is the
// magnitude of a vector SUM, and moves when the world-frame vectors turn
// relative to each other. The world-frame VECTOR is the quantity with the clean
// law, not its magnitude.
math::Vec3 total_angular_momentum3d(const std::vector<RigidBody3D>& bodies);

// Translational plus rotational.
//
// The translational part is exactly constant, since no force acts. The
// rotational part follows the energy law above: a torque does work at rate
// omega . T_body, plus the same second-order forward-Euler term M20 had. With
// no torque it can only rise; with a torque it does whatever the work does.
math::Real total_kinetic_energy3d(const std::vector<RigidBody3D>& bodies);
} // namespace malloy::rigid
