#pragma once

#include <optional>

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
} // namespace malloy::collide
