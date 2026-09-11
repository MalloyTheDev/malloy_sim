#pragma once

#include <optional>

#include <malloy/collide/shapes3d.hpp>
#include <malloy/math/vec3.hpp>

namespace malloy::collide
{
// A resolved overlap between two 3D shapes, `a` and `b`, in that order. The 3D
// sibling of `Contact`, with the same conventions:
//
//   normal       unit vector pointing from `a` toward `b`.
//   penetration  how deeply they overlap along `normal`; never negative.
//                Moving `b` by `+normal * penetration`, or `a` by
//                `-normal * penetration`, separates them exactly.
//   point        a representative point in the overlap region, midway between
//                the two surfaces along `normal`.
//
// Shapes that merely touch produce a contact with `penetration == 0` rather
// than no contact, exactly as in 2D.
struct Contact3
{
    math::Vec3 normal{};
    math::Real penetration{0.0};
    math::Vec3 point{};
};

// Cheap "do these touch" tests, no contact data. Touching counts. An invalid
// shape never overlaps anything.
bool overlaps(const Sphere& sphere, const Plane3& plane);
bool overlaps(const Sphere& a, const Sphere& b);

// Sphere against plane. Like circle against halfplane, this has NO fallback:
// the normal is the plane's own, exactly -plane.normal for every configuration,
// including a sphere whose centre lies exactly on the surface. Returns no value
// when they do not overlap or when either shape is invalid; never throws
// (docs/04).
std::optional<Contact3> contact(const Sphere& sphere, const Plane3& plane);

// Sphere against sphere, the 3D sibling of circle against circle. The normal
// points from `a` toward `b`. Coincident centres are the one degenerate case:
// every direction separates them equally, so the normal falls back to a fixed
// +x, exactly as the 2D circle pair does (a documented, deterministic choice
// rather than a divide by zero).
std::optional<Contact3> contact(const Sphere& a, const Sphere& b);
} // namespace malloy::collide
