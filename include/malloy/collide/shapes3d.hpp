#pragma once

#include <malloy/math/vec3.hpp>

namespace malloy::collide
{
// A sphere in 3D. A zero radius is allowed and means a point.
//
// The 3D sibling of `Circle`, and like it this module knows nothing about
// RigidBody3D or any simulation type: it works on shapes, and pairing a body
// with a shape is the caller's job (docs/03_MODULE_BOUNDARIES.md). The 2D and
// 3D primitives are separate concrete types for the same reason Vec2 and Vec3
// are (ADR 0009): a cross product is a scalar in the plane and a vector in
// space, so one is not the other with a component added.
struct Sphere
{
    math::Vec3 center{};
    math::Real radius{0.0};

    // Valid when the radius is non-negative and every value is finite.
    bool is_valid() const;
};

// A plane: everything on one side of an infinite flat surface, the 3D sibling
// of `Halfplane`. Used for ground and walls, where a sphere is the wrong shape
// and a box is a lie (a floor is not a slab a body pops out the bottom of).
//
// The surface is the set of points p with dot(normal, p) == offset, so the
// signed distance from p to it is dot(normal, p) - offset. `normal` points OUT
// of the solid side, into free space, so that distance is positive for a point
// in free space and negative for one inside the solid. A floor at z = -2 is
// normal (0, 0, 1) with offset -2.
//
// Like the 2D halfplane and unlike every other pair, a sphere against a plane
// has NO degenerate case: the normal is the plane's own and never has to be
// inferred from two centres, so no fallback is ever needed. That is the point
// of the primitive, and it is why M22 uses it before sphere against sphere.
struct Plane3
{
    math::Vec3 normal{0.0, 0.0, 1.0};
    math::Real offset{0.0};

    // Valid when the normal is finite, the offset is finite, and the normal is
    // a UNIT vector. Unit length is required rather than normalized on use, for
    // the reason Halfplane gives: the signed distance is only a distance when
    // the normal is unit, and a caller with arbitrary input (the scenario
    // loader) normalizes once at the boundary and rejects what cannot be.
    bool is_valid() const;
};
} // namespace malloy::collide
