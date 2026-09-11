#pragma once

#include <malloy/math/quat.hpp>
#include <malloy/math/vec3.hpp>

namespace malloy::rigid
{
// A rigid body in three dimensions, free to rotate.
//
// This is `RigidBody2D` with its two scalars promoted, exactly as ADR 0007 said
// they would be: the angle becomes a quaternion and the scalar inertia becomes
// three principal moments. Neither was a simplification to be undone; the 2D
// forms are these with two components forced to zero.
//
// It lives in `malloy_rigid` beside the 2D body, on the same reasoning M19 used
// for gravity: rigid-body dynamics is one domain, and the dimension is not a
// domain (ADR 0009).
struct RigidBody3D
{
    math::Vec3 position{};
    math::Vec3 velocity{};

    // Orientation, as a UNIT quaternion. Validity requires it to be one,
    // because every rotation formula assumes it.
    math::Quat orientation{};

    // Angular velocity in the BODY frame, not the world frame.
    //
    // That choice is what makes Euler's equations diagonal below. In the body
    // frame the inertia is constant and diagonal; in the world frame it rotates
    // with the body and the equations need a full tensor at every step.
    math::Vec3 angular_velocity{};

    math::Real mass{1.0};

    // The three PRINCIPAL moments of inertia, in the body frame, about the
    // centre of mass. A general inertia tensor is symmetric and therefore
    // diagonalizable, so every rigid body has such a frame; storing the
    // diagonal is not a restriction, it is a choice of axes.
    //
    // Deliberately NOT a 3x3 matrix. A general tensor is needed once bodies are
    // built from composed shapes and the parallel-axis step moves inertia off
    // the principal axes, which comes with 3D mass properties and 3D contacts.
    // Until then a matrix would be three extra zeros and a speculative
    // abstraction (rule 11).
    //
    // When all three are equal the body is a sphere and cannot precess. When
    // two are equal it is a symmetric top. When all three differ it can tumble,
    // which is what M20 is about.
    math::Vec3 inertia{1.0, 1.0, 1.0};

    // Valid when mass and all three principal moments are strictly positive and
    // finite, the orientation is a unit quaternion, and position, velocity and
    // angular velocity are SQUARABLE rather than merely finite (docs/04).
    bool is_valid() const;
};

// Angular momentum in the WORLD frame: R (I omega_body).
//
// The invariant M20 is built around, and the reason it must be computed in the
// world frame: under torque-free rotation the body-frame vector I omega changes
// as the body tumbles, while the world-frame one does not.
//
// That holds exactly for the continuum equations. The discrete ones gain
// dt^2 |L x omega|^2 in |L|^2 per step, exactly; see `total_angular_momentum3d`
// for the derivation and for the configurations in which it is zero.
math::Vec3 angular_momentum(const RigidBody3D& body);

// Rotational kinetic energy, (1/2) omega . (I omega), computed in the body
// frame where the inertia is diagonal. A scalar, so the frame does not matter,
// and conserved under torque-free rotation on the same terms as the angular
// momentum above: exactly in the continuum, with an exact second-order gain per
// discrete step.
math::Real rotational_energy(const RigidBody3D& body);
} // namespace malloy::rigid
