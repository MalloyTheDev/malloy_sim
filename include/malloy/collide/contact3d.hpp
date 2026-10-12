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
//   penetration  how deeply they overlap along `normal`; never negative,
//                except from the speculative queries at the end of this
//                file (M40), where a negative value is the gap.
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

// The contact MANIFOLD of two oriented boxes, the 3D sibling of
// `contacts(Obb2, Obb2)`. When a face axis wins the SAT, the face of the other
// box that points most against the reference face is clipped (Sutherland-
// Hodgman) to the reference face's four side planes and kept only where it lies
// inside the reference box, giving up to eight points: four for a face resting
// flat on a face, fewer for a tilted or overhanging one. When an edge-edge axis
// wins, the contact genuinely is a point, so the manifold is the single contact
// `contact(Box3, Box3)` reports. An edge axis that merely TIES a face axis does
// not win: two faces meeting with a twist about their shared normal have edges
// in the face plane whose cross product is that normal, and the face is the
// contact (issue #22). Each entry has the SAT normal (from `a` toward `b`), its
// own penetration behind the reference face (0 for a point kept within
// rounding in front of it), and the midway point. Clipping can leave a vertex
// on top of another or partway along a straight edge; those are dropped, since
// they would weight the centroid (issue #23).
//
// A rigid body that reduces this to its centroid gets a contact through the
// middle of a face-on-face touch, so two aligned boxes meeting head-on exchange
// momentum along the line of centres with no torque, as M30 does for a box on a
// plane. Points within a rounding-level tolerance in front of the reference
// face are kept, so a bare touch gives the whole face rather than a lopsided
// part of it (issue #23). Never empty when the boxes overlap (if nothing
// survives clipping, it falls back to the single-point contact), so
// `!contacts(a, b).empty()` agrees with `overlaps(a, b)`. Empty when apart or
// either box is invalid; never throws.
std::vector<Contact3> contacts(const Box3& a, const Box3& b);

// An oriented box against a sphere (M34), the last movable pair in 3D. No SAT is
// needed: the nearest point of the box to the sphere's centre is found by
// clamping that centre into the box's local frame, and the two overlap when it
// lies within `radius` of the centre. Touching counts. The order is symmetric:
// `contact(Box3, Sphere)` gives a normal from the box toward the sphere, and
// `contact(Sphere, Box3)` the negation, so both keep the a-toward-b convention.
//
// The normal is the direction from that nearest point to the sphere's centre,
// and the penetration is `radius` minus their distance. The one degenerate case
// is a sphere centre inside the box, where that direction vanishes: it falls
// back to the least-penetrated face axis (its outward direction, ties broken by
// axis order), the deterministic analogue of the sphere pair's fixed normal.
// The contact point is midway between the box surface and the sphere surface
// along the normal. Returns no value when they do not overlap or either shape is
// invalid; never throws (docs/04).
bool overlaps(const Box3& box, const Sphere& sphere);
bool overlaps(const Sphere& sphere, const Box3& box);
std::optional<Contact3> contact(const Box3& box, const Sphere& sphere);
std::optional<Contact3> contact(const Sphere& sphere, const Box3& box);

// --- Speculative queries (M40) ---
//
// The five contact queries the 3D rigid world resolves, with a `margin`: a pair
// is reported while it is separated by up to `margin`, not only once it
// touches. The iterative contact solver uses this to keep a resting contact
// alive across rounding (a corner that drifts 1e-17 clear of the floor would
// otherwise drop out and fall for a step) and to stop an approach exactly at
// the surface. Each entry follows the conventions above with one difference:
//
//   penetration  SIGNED. Negative means the shapes are apart by -penetration
//                (a gap of at most `margin`); the point is still midway
//                between the two surfaces, so it lies in the gap.
//
// With `margin == 0` the sphere and box/plane queries return exactly what the
// plain queries return (those delegate here). The box/box manifold differs from
// `contacts(Box3, Box3)` in two ways: its depths and its single-point fallback
// are signed rather than clamped at 0, and an edge axis that beats a
// PENETRATING face by less than 5 per cent of the face's overlap does not win
// when the face manifold keeps a point (the face preference ODE, Bullet and
// Box2D use against face/edge flicker in resting stacks; see kFacePreference).
// A margin that is negative or not finite is refused like an invalid shape: no
// contact. Never throws.
std::optional<Contact3> contact(const Sphere& sphere, const Plane3& plane,
                                math::Real margin);
std::optional<Contact3> contact(const Sphere& a, const Sphere& b, math::Real margin);
std::vector<Contact3> contacts(const Box3& box, const Plane3& plane, math::Real margin);
std::vector<Contact3> contacts(const Box3& a, const Box3& b, math::Real margin);
std::optional<Contact3> contact(const Box3& box, const Sphere& sphere, math::Real margin);
} // namespace malloy::collide
