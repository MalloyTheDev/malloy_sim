#pragma once

#include <optional>
#include <vector>

#include <malloy/collide/shapes.hpp>
#include <malloy/math/vec2.hpp>

namespace malloy::collide
{
// A resolved overlap between two shapes, `a` and `b`, in that order.
//
// Conventions, which every query below obeys:
//
//   normal       unit vector pointing from `a` toward `b`.
//   penetration  how deeply they overlap along `normal`; never negative.
//                Moving `b` by `+normal * penetration`, or `a` by
//                `-normal * penetration`, separates them exactly.
//   point        a representative point in the overlap region, midway between
//                the two surfaces along `normal`.
//
// Shapes that merely touch produce a contact with `penetration == 0` rather
// than no contact. Reporting the boundary is more useful to a caller than
// silently dropping it, and it is exactly representable, so it is testable.
struct Contact
{
    math::Vec2 normal{};
    math::Real penetration{0.0};
    math::Vec2 point{};
};

// Cheap "do these touch" tests, with no contact data. Touching counts as an
// overlap, matching the Contact convention above. An invalid shape never
// overlaps anything.
bool overlaps(const Circle& a, const Circle& b);
bool overlaps(const Aabb& a, const Aabb& b);
bool overlaps(const Circle& circle, const Aabb& box);
bool overlaps(const Circle& circle, const Halfplane& plane);
bool overlaps(const Obb2& a, const Obb2& b);
bool overlaps(const Obb2& box, const Halfplane& plane);

// Full contact queries. Return no value when the shapes do not overlap or when
// either shape is invalid; these never throw (docs/04).
//
// Where the normal is geometrically undefined the queries fall back to a
// documented, deterministic choice rather than returning NaN:
//
//   circle/circle with coincident centers  -> normal is +x
//   circle inside a box, centered exactly  -> normal is +x
//   box/box overlapping equally on both axes -> the x axis wins
//
// The fallbacks are arbitrary but fixed, so results stay repeatable
// (docs/04 determinism).
std::optional<Contact> contact(const Circle& a, const Circle& b);
std::optional<Contact> contact(const Aabb& a, const Aabb& b);
std::optional<Contact> contact(const Circle& circle, const Aabb& box);

// Circle against halfplane. This one has no fallback, because the normal is
// the plane's own: it is exactly -plane.normal for every configuration,
// including a circle whose centre lies exactly on the line. The convention
// still holds, so it points from the circle toward the solid side.
std::optional<Contact> contact(const Circle& circle, const Halfplane& plane);

// Oriented box against oriented box (M36), by the separating-axis theorem, the
// 2D sibling of the `Box3` pair. Four axes are tested (the two face normals of
// each box; 2D has no edge-edge case, because an edge's separating direction is
// already one of those face normals), and the axis of least overlap is the
// contact normal with its overlap the penetration. `overlaps` runs the same
// predicate, so the two always agree. Like the box/box pair in 3D it reduces to
// a SINGLE contact, placed on the deepest vertex of the other box: enough for a
// bounce, not the manifold a resting stack needs (rule 12). Two coincident
// centres are the degenerate case and fall back to the axis of least combined
// width with a fixed sign, the analogue of the circle pair's +x. Returns no
// value when they do not overlap or either box is invalid; never throws.
std::optional<Contact> contact(const Obb2& a, const Obb2& b);

// An oriented box against a halfplane (M37), the 2D sibling of `Box3` against a
// plane (M30). Like that pair, and unlike the single-point queries, this returns
// a MANIFOLD: one Contact per box corner that lies inside the solid (signed
// distance <= 0), because a box resting on a floor touches it at up to two
// corners at once and a single point could not hold it flat. Each entry has:
//
//   normal       -plane.normal, the box's push-out direction (the plane supplies
//                it, so there is no fallback), the same for every corner.
//   penetration  how far that corner is below the surface; never negative.
//   point        that corner, moved half its penetration back toward the
//                surface, matching the midway convention of the other pairs.
//
// The corners are enumerated in a fixed order (the sign bits of the local axes),
// so a run repeats (docs/04). Empty when the box is clear of the plane or either
// shape is invalid; never throws.
std::vector<Contact> contacts(const Obb2& box, const Halfplane& plane);

// The contact MANIFOLD of two oriented boxes: the points where they actually
// touch, rather than the single representative vertex `contact(Obb2, Obb2)`
// reports. The SAT picks the reference face (the separating face of whichever
// box won the least-overlap axis); the face of the other box that points most
// against it is the incident edge, which is clipped to the reference face's
// extent and kept only where it lies inside the reference box. That gives up to
// two points: two for a face resting on a face, one for a corner. Each has:
//
//   normal       the SAT normal, from `a` toward `b`, the same for every point.
//   penetration  how far that point lies behind the reference face (0 for a
//                point kept within rounding in front of it, issue #23).
//   point        that point, moved half its penetration toward the reference
//                face, the midway convention of every pair.
//
// A rigid body that reduces this to its centroid gets a contact through the
// middle of a face-on-face touch, so two aligned boxes meeting head-on exchange
// momentum along the line of centres with no torque; the single vertex of
// `contact` sits on a corner and would invent one. Points within a
// rounding-level tolerance in front of the reference face are kept, so a bare
// touch gives both ends of the face rather than one (issue #23). Never empty
// when the boxes overlap (if nothing survives clipping, it falls back to the
// single-point contact), so `!contacts(a, b).empty()` agrees with
// `overlaps(a, b)`. Empty when apart or either box is invalid; never throws.
std::vector<Contact> contacts(const Obb2& a, const Obb2& b);
} // namespace malloy::collide
