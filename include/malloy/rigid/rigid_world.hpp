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
// torque, so they tumble rather than merely bouncing. M36 added boxes (bodies
// with half extents), which collide with each other by the separating-axis
// test, and M37 put them on halfplane ground. A box against a disc is not
// implemented yet: the pair is skipped (issue #48).
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

    // Ok, or InvalidSettings (dt; restitution outside [0, 1]; gravity not
    // finite; friction negative; a ground halfplane invalid) / InvalidState (a
    // body).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step (M40 order; every constant named here
    // is defined once in rigid.cpp with its derivation, the same values as the
    // 3D world's):
    //
    //   A. DETECT contacts at the START-of-step poses, in a fixed order so a run
    //      repeats (docs/04): every pair of bodies in ascending pair order (two
    //      discs, M14; two boxes, M36; two immovable bodies and a box with a
    //      disc are skipped), then every movable body against every ground
    //      plane in ascending body then plane order (discs M16, boxes M37).
    //      Each pair gives a MANIFOLD of points in the collide layer's order (up
    //      to 2 for a box). The queries are SPECULATIVE (kMarginFraction): a
    //      pair still apart by up to 1e-4 times the smaller body's size is kept,
    //      at a negative penetration (the gap), so a resting contact cannot drop
    //      out on rounding.
    //
    //   B. IMPACT phase, only when the restitution is above 0: on the
    //      start-of-step velocities, a point is an impact when its normal
    //      velocity at the MIDDLE of this step's gravity (v_n + delta/2, delta
    //      what step 1 adds there) approaches faster than the threshold
    //      4 dt |g| (kThresholdFactor), and it is touching or would close its gap
    //      within the step. Those points are driven by sequential impulses to
    //      the Newton target -e (v_n + delta/2) - delta/2, so after the gravity
    //      each leaves at e times its mid-step approach. Points that were not
    //      approaching take no part, so the phase never raises the kinetic
    //      energy.
    //
    //   1. add gravity * dt to every velocity whose body has finite mass;
    //
    //   C. REST phase, every point: sequential impulses drive each point's
    //      normal velocity to 0 (to -gap/dt for a point still apart), with
    //      Coulomb friction clamped to friction times that point's ACCUMULATED
    //      normal impulse (M17, per point). Impulses are accumulated and
    //      clamped (normal never below 0), WARM-STARTED from the last step's
    //      rest-phase impulses of the matching point (box corner index on a
    //      plane; nearest point in body a's frame within 0.05 of its size for
    //      two boxes; the one point of a disc pair). Per sweep, manifolds in
    //      detection order, each its normal rows then its friction rows; each
    //      manifold first takes one impulse at its weighted centroid. The
    //      sweeps stop when no point's constraint is violated by more than
    //      1e-8 of the velocity scale (1e-12 in the impact phase), or after
    //      100 sweeps (kMaxSweeps). Only bodies that received a nonzero impulse
    //      are written back, and an immovable body never is.
    //
    //   2. move the centre of mass by the solved velocity * dt, which keeps
    //      the integration semi-implicit (symplectic) Euler;
    //   3. advance the angle by angular_velocity * dt;
    //
    //   D. PROJECT positions (nonlinear Gauss-Seidel, 3 passes): each point's
    //      penetration beyond the slop (1e-6 times the pair's size) is removed a
    //      fifth at a time (kProjectionFactor), at most 0.2 times the size per
    //      push, by moving and turning the movable bodies along the normal; a
    //      LEVEL manifold is first lifted whole at its centroid. Velocities are
    //      not touched, so this adds no kinetic energy.
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
    // The world's state is the bodies AND the warm-start impulses: a world
    // rebuilt from `bodies()` alone starts without them and so differs for its
    // first steps, deterministically. On validation failure, before or after
    // the update, both are left unchanged and the failing status is returned.
    // Never throws.
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

    // One rest-phase contact point remembered for the next step's warm start
    // (M40), the part of the world's state that is not in `bodies()`.
    struct WarmContact
    {
        std::size_t a{0};
        std::size_t b{0};        // the other body, or the plane index when `ground`
        bool ground{false};
        int key{-1};             // box corner index, 0 for a disc pair, -1: by position
        math::Vec2 local{};      // the contact point in body a's frame
        math::Real normal_impulse{0.0};
        math::Real tangent_impulse{0.0};
    };

    // Read-only, so the whole state can be compared (a failed step restores
    // it with the bodies); there is no way to set it.
    const std::vector<WarmContact>& warm_contacts() const { return warm_; }

private:
    sim_core::SimulationSettings simulation_settings_;
    RigidSettings settings_;
    std::vector<RigidBody2D> bodies_;
    std::vector<RigidBody2D> previous_;
    // Ordered as the manifolds are detected (pairs, then ground), so matching
    // is one merge per step. Snapshotted and rolled back with the bodies.
    std::vector<WarmContact> warm_;
    std::vector<WarmContact> previous_warm_;
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
