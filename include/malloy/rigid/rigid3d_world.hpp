#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/collide/shapes3d.hpp>
#include <malloy/math/vec3.hpp>
#include <malloy/rigid/rigid_body3d.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::rigid
{
// Settings for the 3D rigid world: a constant torque (M21), a constant applied
// force (M28), and gravity, restitution and ground planes for contacts (M22).
//
// The torque is in the WORLD frame, not the body frame. That is a deliberate
// choice and it is what makes M21's headline invariant clean: the physical law
// dL/dt = torque holds in the world frame, so with a constant world torque each
// body's world-frame angular momentum grows along a straight line,
//
//   L(t) = L(0) + torque * t,
//
// exactly in the continuum. An external couple fixed in the lab (the kind a
// field produces) is world-frame; a thruster bolted to the body would be
// constant in the BODY frame instead, a different thing that is not reachable
// from this one without knowing the body, and is deferred.
//
// Every field defaults to the pre-M22 behaviour: zero torque, zero gravity, no
// ground, so a `Rigid3DWorld` built without settings is the torque-free M20
// world bit for bit, and one given only a torque is the M21 world.
struct Rigid3DSettings
{
    math::Vec3 torque{};

    // A constant applied FORCE, in the world frame: the translational half of a
    // wrench whose rotational half is `torque` above. It is applied as the
    // acceleration F/m it produces, before the position update, so unlike
    // `gravity` it scales with 1/mass. That is the whole difference between a
    // force and an acceleration, and it is why the force reads the body's mass
    // (now computed from geometry, M25 to M27). It acts through the centre of
    // mass, so on its own it makes no torque; an off-centre push is this force
    // plus the couple r x F carried by `torque`, which is what makes the pair a
    // complete wrench. It is a single constant SETTING applied as forcing, NOT
    // the force accumulator rule 5 defers: M28 is that rule's dedicated
    // milestone for a translational force, and adds one force, not an
    // accumulation of them. Zero by default, so a pre-M28 world is unchanged.
    math::Vec3 force{};

    // Bounciness of every contact. 1 is perfectly elastic and conserves kinetic
    // energy across a bounce; 0 is perfectly inelastic, so the body stops
    // separating along the contact normal. The 3D sibling of RigidSettings.
    math::Real restitution{1.0};

    // Coulomb friction coefficient for every contact. 0 is frictionless, which
    // is how M22 behaved and remains the default. The tangential impulse is
    // clamped to `friction` times the normal impulse, so it is self-limiting: a
    // body in mid-air has no normal impulse and so cannot be turned by friction.
    //
    // This is the FIRST 3D contact that imparts spin. A normal impulse on a
    // centred sphere has no lever arm (r x n = 0), but a tangential one does
    // (r x t is nonzero), so friction is where the rotational effective-mass
    // term finally does work. Not capped at 1: a coefficient above 1 is real.
    math::Real friction{0.0};

    // Uniform gravitational ACCELERATION, applied to every body before the
    // position update. An acceleration, not a force (CLAUDE.md rule 5): it does
    // not scale with mass and needs no force accumulator. Defaults to zero.
    //
    // Gravity acts through the centre of mass, so on its own it generates no
    // torque; a body spins under it only through a contact away from the centre
    // of mass, which for a centred sphere means not at all until friction.
    math::Vec3 gravity{};

    // Immovable ground planes: floors, walls, ramps. Each resolves as a body of
    // infinite mass would, but a plane is NOT stored as a body here: it is pure
    // geometry, and the sphere is the only thing that moves, so `RigidBody3D`
    // never has to represent an infinite mass. A plane's normal is its own and
    // never turns, unlike a floor built from spheres. Empty by default.
    std::vector<collide::Plane3> ground{};

    // Valid when the torque is SQUARABLE (it enters |u|^2 in the drift law),
    // the force is SQUARABLE (its square enters the free-flight energy drift, as
    // gravity's does), restitution is in [0, 1] and finite, friction is
    // non-negative and finite, gravity is SQUARABLE, and every ground plane is
    // valid.
    bool is_valid() const;
};

// Rigid bodies in three dimensions: free rotation (M20), a constant applied
// torque (M21), and now gravity and restitution contacts against ground planes
// (M22).
//
// M22 gives a body a collision radius, drops it under a gravity acceleration,
// and bounces it off immovable planes. It is deliberately the 3D echo of M10's
// colliding particles, not M14's rigid contacts: a contact on a CENTRED sphere
// passes through the centre of mass, so it imparts no spin, and the response is
// a pure normal impulse. Spin from a contact needs a tangential component:
// friction, which M23 added, so a sliding sphere now spins up and rolls. Box
// contacts (M30 on the ground, M33 against a box, M34 against a sphere) have
// lever arms of their own and impart spin through the principal moments.
//
// The contacts are impulses and gravity is an acceleration, so the world still
// needs no force accumulator; the applied force and torque (M28, M21) are
// constant settings applied as forcing (CLAUDE.md rule 5).
class Rigid3DWorld
{
public:
    Rigid3DWorld(sim_core::SimulationSettings simulation_settings,
                 std::vector<RigidBody3D> bodies, Rigid3DSettings settings = {});

    // Ok, or InvalidSettings (dt; torque, force or gravity not squarable;
    // restitution outside [0, 1]; friction negative; a ground plane invalid) /
    // InvalidState (a body).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step (M40 order; every constant named here
    // is defined once in rigid3d.cpp with its derivation):
    //
    //   A. DETECT contacts at the START-of-step poses, in a fixed order so a run
    //      repeats (docs/04): every pair of bodies in ascending pair order
    //      (sphere against sphere, M24; box against box, M33; box against
    //      sphere, M34), then every body against every ground plane in
    //      ascending body then plane order (M22, boxes M30). Each pair gives a
    //      MANIFOLD of points in the collide layer's own order (up to 4 for a
    //      box on a plane, 8 for two boxes). The queries are SPECULATIVE
    //      (kMarginFraction): a pair still apart by up to 1e-4 times the
    //      smaller body's size is kept, at a negative penetration (the gap), so
    //      a resting contact cannot drop out on rounding. A body with neither a
    //      radius nor box half-widths does not collide.
    //
    //   B. IMPACT phase, only when the restitution is above 0: on the
    //      start-of-step velocities, a point is an impact when its normal
    //      velocity at the MIDDLE of this step's forcing (v_n + delta/2, delta
    //      the change steps 0 and 1 below will make there) approaches faster
    //      than the threshold 4 dt (|g| + |F| / m_min) (kThresholdFactor), and
    //      it is touching or would close its gap within the step. Those points
    //      are driven by sequential impulses to the Newton target
    //      -e (v_n + delta/2) - delta/2, so after the forcing each leaves at
    //      e times its mid-step approach. Points that were not approaching take
    //      no part, so an impact can never pull energy out of a separating
    //      point (the phase never raises the kinetic energy).
    //
    //   0. add the gravity acceleration to every body's velocity, before the
    //      position update, which keeps this semi-implicit Euler. An
    //      acceleration, so it does not scale with mass and no body is exempt
    //      (every RigidBody3D is movable; the ground is geometry, not a body).
    //      Then add the acceleration F/m of the applied FORCE (M28), which DOES
    //      scale with mass and acts through the centre of mass, so it changes
    //      the velocity but makes no torque. Both are zero by default.
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
    //   C. REST phase, every point: sequential impulses drive each point's
    //      normal velocity to 0 (to -gap/dt for a point still apart, so it may
    //      close its gap exactly but not pass it), with Coulomb friction: a
    //      tangential impulse in the contact plane, solved as a 2x2 block and
    //      kept within a disc of radius friction times that point's ACCUMULATED
    //      normal impulse (M23's cone, per point). Outside the disc the impulse
    //      is the disc's point nearest the sticking one in the metric of the
    //      point's tangent effective mass, which is Coulomb's maximum-
    //      dissipation law: the friction then opposes the slide it leaves
    //      (scaling radially does that only for an isotropic effective mass, a
    //      sphere's, and turned a box sliding diagonally). Impulses are
    //      accumulated and clamped (normal never below 0), WARM-STARTED from the
    //      last step's rest-phase impulses of the matching point (box corner
    //      index on a plane; nearest point in body a's frame within 0.05 of its
    //      size for two boxes; the one point of a sphere pair). Per sweep,
    //      manifolds in detection order, each its normal rows then its friction
    //      rows; each manifold first takes one impulse at its weighted centroid
    //      (the prepass that keeps a symmetric landing exactly symmetric). The
    //      sweeps stop when no point's constraint is violated by more than
    //      1e-8 times the velocity scale (1e-12 in the impact phase), or after
    //      100 sweeps (kMaxSweeps): a deterministic count, set by the state.
    //      Only bodies that received a nonzero impulse are written back, so a
    //      body out of contact is bit for bit what it was without contacts.
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
    //   4. move the centre of mass by velocity * dt, with the solved velocity.
    //      Without gravity and without a contact this is still a body flying
    //      straight while it tumbles.
    //
    //   D. PROJECT positions (nonlinear Gauss-Seidel, 3 passes): each point's
    //      penetration beyond the slop (1e-6 times the pair's size) is removed a
    //      fifth at a time (kProjectionFactor), at most 0.2 times the size per
    //      push, by moving and turning the bodies along the normal (a sphere
    //      is only moved); a LEVEL manifold is first lifted whole at its
    //      centroid. Linear velocities are not touched, and through a turn the
    //      body's WORLD-frame angular momentum is kept, re-expressed in the
    //      turned body frame, so the projection changes no momentum, linear or
    //      angular. It can only restore height a body lost by sinking in.
    //
    // A ground plane is immovable, so only the body changes. Impulses are
    // equal and opposite at one point, so a pair conserves linear momentum.
    //
    // The world's state is the bodies AND the warm-start impulses: a world
    // rebuilt from `bodies()` alone starts without them and so differs for its
    // first steps, deterministically. On validation failure, before or after
    // the update, both are left unchanged and the failing status is returned;
    // this never throws.
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

    // One rest-phase contact point remembered for the next step's warm start
    // (M40), the part of the world's state that is not in `bodies()`.
    struct WarmContact
    {
        std::size_t a{0};
        std::size_t b{0};        // the other body, or the plane index when `ground`
        bool ground{false};
        int key{-1};             // box corner index, 0 for a sphere pair, -1: by position
        math::Vec3 local{};      // the contact point in body a's frame
        math::Real normal_impulse{0.0};
        math::Vec3 tangent_impulse{};
    };

private:
    sim_core::SimulationSettings simulation_settings_;
    Rigid3DSettings settings_;
    std::vector<RigidBody3D> bodies_;
    std::vector<RigidBody3D> previous_;
    // Ordered as the manifolds are detected (pairs, then ground), so matching
    // is one merge per step. Snapshotted and rolled back with the bodies.
    std::vector<WarmContact> warm_;
    std::vector<WarmContact> previous_warm_;
    std::optional<time::FixedStep> step_;
};

// --- Diagnostics ---

math::Vec3 total_linear_momentum3d(const std::vector<RigidBody3D>& bodies);

// Total angular momentum about the world origin: the SPIN of each body plus
// its ORBITAL term m (r x v), exactly as the 2D `total_angular_momentum` does.
// Dropping either makes a whole class of motion look conserved when it is not.
//
// The orbital term m (r x v) is exactly constant ONLY when nothing changes a
// body's velocity: for a free or torque-only body d(r x p)/dt = v x mv = 0, bit
// for bit, because with v unchanging r_{n+1} x v = (r_n + dt v) x v = r_n x v.
// M22 is the milestone the M21 comment warned about: gravity changes v every
// step, so gravity applies a torque about the world origin, m r x g, and a
// CONTACT changes v as well. With gravity set or a bounce in progress this
// total is therefore NOT conserved, and is not meant to be. It is conserved
// again exactly when gravity is zero and no contact fires, which is the M20/M21
// regime.
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

// Translational plus rotational KINETIC energy.
//
// The rotational part follows the energy law above: a torque does work at rate
// omega . T_body, plus the same second-order forward-Euler term M20 had. The
// translational part is constant without gravity or a contact; under gravity it
// trades with the potential below, and a contact removes (1 - e^2) of the normal
// part on each bounce.
math::Real total_kinetic_energy3d(const std::vector<RigidBody3D>& bodies);

// Gravitational potential energy in a uniform field: the sum of -m (g . r), r
// the centre of mass. Zero when gravity is zero. The 3D sibling of the 2D
// rigid `total_potential_energy`. Kinetic plus this is the quantity that would
// be conserved in free flight if the integrator were exact; semi-implicit Euler
// instead sheds the derived constant (1/2)(sum m)|g|^2 dt^2 per free-flight
// step, the same drift M12 and M15 measured in two dimensions.
math::Real total_potential_energy3d(const std::vector<RigidBody3D>& bodies,
                                    const math::Vec3& gravity);

// Potential energy of the constant applied FORCE: the sum of -(F . r), r the
// centre of mass, F the world-frame force setting. Zero when the force is zero.
// The force does work F . dr as a body moves, so this is the potential whose
// decrease matches that work, exactly as the gravity potential above does for
// gravity. Kinetic plus both potentials is the quantity conserved in free
// flight if the integrator were exact; semi-implicit Euler instead sheds
// (1/2)(sum |F|^2 / m) dt^2 per free-flight step, the force's analogue of the
// gravity drift (the per-body acceleration is F/m, so its square carries the
// 1/m the gravity term does not).
math::Real total_force_potential3d(const std::vector<RigidBody3D>& bodies,
                                   const math::Vec3& force);
} // namespace malloy::rigid
