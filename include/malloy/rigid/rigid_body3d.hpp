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

    // Collision radius of the sphere centred on the body's position, which is
    // its centre of mass. Zero (the default) means the body does not collide,
    // so a body written before M22 keeps its exact behaviour: it tumbles and
    // flies but passes through everything.
    //
    // Centred on the centre of mass, unlike the 2D disc which is centred on the
    // body ORIGIN offset from the centre of mass. That 2D offset existed only
    // to give a contact a moment arm; a sphere needs no such device, and the
    // consequence is that a NORMAL contact on a centred sphere passes through
    // the centre of mass and so imparts no spin. Spin from contact waits for
    // friction, its own milestone.
    math::Real radius{0.0};

    // Valid when mass and all three principal moments are strictly positive and
    // finite, the orientation is a unit quaternion, the radius is non-negative
    // and finite, and position, velocity and angular velocity are SQUARABLE
    // rather than merely finite (docs/04).
    bool is_valid() const;
};

// SPIN angular momentum in the WORLD frame: R (I omega_body).
//
// Spin only. The orbital term m (r x v) about the world origin is NOT here; it
// belongs to the system rather than to the body's rotation, and
// `total_angular_momentum3d` adds it. The name says `spin` because the
// distinction is easy to lose and expensive to lose: dropping the orbital term
// makes a whole class of motion look conserved when it is not, which is the
// same warning the 2D `total_angular_momentum` carries.
//
// The reason it is computed in the WORLD frame is the invariant M20 is built
// around: under torque-free rotation the body-frame vector I omega changes as
// the body tumbles, while the world-frame one does not.
//
// This is also the quantity the discrete drift law is about, and the law is
// PER BODY. The continuum equations conserve it; the discrete ones gain
// dt^2 |L x omega|^2 in |L|^2 per step, exactly. See
// `total_angular_momentum3d` for the derivation and for when it is zero, and
// note there why the magnitude of a SUM of these does not inherit the law.
math::Vec3 spin_angular_momentum(const RigidBody3D& body);

// Rotational kinetic energy, (1/2) omega . (I omega), computed in the body
// frame where the inertia is diagonal. A scalar, so the frame does not matter,
// and conserved under torque-free rotation on the same terms as the angular
// momentum above: exactly in the continuum, with an exact second-order gain per
// discrete step.
math::Real rotational_energy(const RigidBody3D& body);
} // namespace malloy::rigid
