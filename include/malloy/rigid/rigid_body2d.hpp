#pragma once

#include <malloy/math/vec2.hpp>

namespace malloy::rigid
{
// A 2D rigid body: a pose and a mass distribution.
//
// Not a widened nbody::Body2D and not a Particle2D with an angle bolted on. A
// rigid body obeys different equations of motion and exposes different
// invariants, so it owns its own state
// (docs/decisions/0007-rigid-bodies-own-their-state.md).
//
// 2D throughout: orientation is a scalar angle in radians and inertia is a
// scalar, which in 2D is rotation-invariant. There are no quaternions and no
// inertia tensors here; those are 3D concerns deferred to M19 (docs/04).
//
// `local_center_of_mass` is the offset from `position` (the body origin) to the
// centre of mass, expressed in body-local coordinates and therefore rotated by
// `angle` to reach world space. It is deliberately allowed to be nonzero: a
// body whose origin IS its centre of mass hides an implementation that measures
// torque arms from the origin (docs/05, "Prefer asymmetric configurations").
struct RigidBody2D
{
    math::Vec2 position{};        // body origin in world space
    math::Real angle{0.0};        // radians, counter-clockwise positive
    math::Vec2 velocity{};        // linear velocity of the CENTRE OF MASS
    math::Real angular_velocity{0.0};

    math::Real mass{1.0};
    // Moment of inertia about the CENTRE OF MASS, not about the body origin.
    // Keeping one documented reference point is what stops a parallel-axis term
    // being applied twice or not at all.
    math::Real inertia{1.0};

    math::Vec2 local_center_of_mass{};

    // Valid when mass and inertia are strictly positive and every value is
    // finite. Infinite ("static") bodies are not represented: nothing needs
    // them until contact response does, and inventing a representation early
    // would be speculative (CLAUDE.md rule 11).
    bool is_valid() const;
};

// The centre of mass in world space: the body origin plus the local offset
// rotated by the body angle.
math::Vec2 center_of_mass(const RigidBody2D& body);

// Convert between body-local and world coordinates using the body pose.
// Round-tripping a point through both must return it unchanged.
math::Vec2 to_world(const RigidBody2D& body, const math::Vec2& local_point);
math::Vec2 to_local(const RigidBody2D& body, const math::Vec2& world_point);

// Velocity of the material point currently at `world_point`, which is the
// linear velocity of the centre of mass plus the rotational contribution
// omega x r, with r measured from the centre of mass. In 2D that cross product
// is (-omega * r.y, omega * r.x).
math::Vec2 velocity_at(const RigidBody2D& body, const math::Vec2& world_point);
} // namespace malloy::rigid
