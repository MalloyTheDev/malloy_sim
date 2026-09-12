#pragma once

#include <optional>
#include <vector>

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
bool overlaps(const Box3& box, const Plane3& plane);
bool overlaps(const Box3& a, const Box3& b);

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

// An oriented box against a plane (M30). Unlike the sphere pairs this returns a
// MANIFOLD, one Contact3 per box corner that is inside the solid (signed
// distance <= 0), because a box resting on a plane touches it at up to four
// corners at once and a single point could not hold it flat. Each entry has:
//
//   normal       -plane.normal, the box's push-out direction (the plane
//                supplies it, so there is no fallback), the same for every
//                corner in the manifold.
//   penetration  how far that corner is below the surface; never negative.
//   point        that corner, moved half its penetration back toward the
//                surface, matching the midway convention of the other pairs.
//
// The corners are enumerated in a fixed order (the sign bits of the local
// axes), so a run repeats (docs/04). Empty when the box is clear of the plane
// or either shape is invalid; never throws.
std::vector<Contact3> contacts(const Box3& box, const Plane3& plane);

// Oriented box against oriented box (M33), by the separating-axis theorem.
// Unlike the box/plane pair this reduces to a SINGLE contact rather than a
// manifold: the goal here is a correct bounce (an exact normal and penetration
// depth), not a stable resting stack, which needs the full contact manifold and
// an iterative solver that rule 12 defers. Fifteen axes are tested (the three
// face normals of each box and the nine pairwise edge cross products); the one
// of least overlap is the separating direction, and its overlap is the
// penetration. `overlaps` shares that same predicate, so the two agree exactly.
//
// The normal points from `a` toward `b`, as for every 3D pair. The contact
// point depends on which axis won: a face axis places it on the deepest vertex
// of the other box (moved half the penetration back toward the surface), an
// edge-edge axis at the midpoint of the closest approach of the two edges. Two
// coincident box centres are the one degenerate case and fall back to the axis
// of least combined radius with a fixed sign, the analogue of the sphere pair's
// fixed +x. Returns no value when they do not overlap or either box is invalid;
// never throws (docs/04).
std::optional<Contact3> contact(const Box3& a, const Box3& b);
} // namespace malloy::collide
