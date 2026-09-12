#include <malloy/collide/contact3d.hpp>

#include <algorithm>
#include <cmath>

#include <malloy/collide/shapes3d.hpp>
#include <malloy/math/math.hpp>

namespace malloy::collide
{
bool Sphere::is_valid() const
{
    return radius >= math::Real{0} && math::is_finite(radius) &&
           math::is_finite(center);
}

bool Plane3::is_valid() const
{
    // Unit normal, on purpose: the signed distance below is only a distance
    // when the normal is unit. Same invariant Halfplane keeps in 2D.
    return math::is_finite(normal) && math::is_finite(offset) &&
           std::abs(math::length_squared(normal) - math::Real{1}) <= math::Real{1e-12};
}

bool Aabb3::is_valid() const
{
    return math::is_finite(min) && math::is_finite(max) && min.x <= max.x &&
           min.y <= max.y && min.z <= max.z;
}

bool Box3::is_valid() const
{
    return half_extents.x > math::Real{0} && half_extents.y > math::Real{0} &&
           half_extents.z > math::Real{0} && math::is_finite(half_extents) &&
           math::is_unit(orientation) && math::is_finite(center);
}

namespace
{
// Signed distance from a point to the plane's surface: positive in free space,
// negative inside the solid. The 3D copy of the 2D helper, one component wider.
math::Real signed_distance(const Plane3& plane, const math::Vec3& point)
{
    return math::dot(plane.normal, point) - plane.offset;
}

// Coincident-centre fallback, the 3D copy of the 2D one: arbitrary but fixed,
// so a run repeats.
const math::Vec3 fallback_normal{1.0, 0.0, 0.0};
} // namespace

bool overlaps(const Sphere& sphere, const Plane3& plane)
{
    if (!sphere.is_valid() || !plane.is_valid())
    {
        return false;
    }
    // Touching counts, matching every other query.
    return signed_distance(plane, sphere.center) <= sphere.radius;
}

std::optional<Contact3> contact(const Sphere& sphere, const Plane3& plane)
{
    if (!sphere.is_valid() || !plane.is_valid())
    {
        return std::nullopt;
    }

    const math::Real distance = signed_distance(plane, sphere.center);
    const math::Real depth = sphere.radius - distance;
    if (depth < math::Real{0})
    {
        return std::nullopt;
    }

    Contact3 result;
    // From the sphere toward the plane's solid side, which is the direction the
    // plane's outward normal does not point. No fallback is needed even when
    // the centre lies exactly on the surface: the plane supplies the direction.
    result.normal = -plane.normal;
    result.penetration = depth;
    // Midway between the two surfaces along the normal, as for every other
    // pair: the sphere's surface is `radius` along the normal from its centre,
    // and the plane's surface is `depth` further back.
    result.point =
        sphere.center + result.normal * (sphere.radius - depth / math::Real{2});
    return result;
}

bool overlaps(const Sphere& a, const Sphere& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return false;
    }
    const math::Vec3 delta = b.center - a.center;
    const math::Real distance_squared = math::dot(delta, delta);
    if (!math::is_finite(distance_squared))
    {
        return false; // so far apart the squared distance overflowed
    }
    const math::Real radius_sum = a.radius + b.radius;
    return distance_squared <= radius_sum * radius_sum;
}

std::optional<Contact3> contact(const Sphere& a, const Sphere& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return std::nullopt;
    }

    const math::Vec3 delta = b.center - a.center;
    const math::Real distance_squared = math::dot(delta, delta);
    if (!math::is_finite(distance_squared))
    {
        return std::nullopt;
    }

    const math::Real radius_sum = a.radius + b.radius;
    if (distance_squared > radius_sum * radius_sum)
    {
        return std::nullopt;
    }

    Contact3 result;
    if (distance_squared > math::Real{0})
    {
        const math::Real distance = std::sqrt(distance_squared);
        result.normal = delta / distance; // from a toward b
        result.penetration = radius_sum - distance;
    }
    else
    {
        // Coincident centres: every direction separates them equally, so take
        // the documented one rather than dividing by zero. Same rule as the 2D
        // circle pair.
        result.normal = fallback_normal;
        result.penetration = radius_sum;
    }
    // Midway between the two surfaces along the normal, as for every other pair.
    result.point =
        a.center + result.normal * (a.radius - result.penetration / math::Real{2});
    return result;
}

namespace
{
// The eight corners of an oriented box, in a fixed order (the sign bits of the
// three local axes), so any loop over them repeats a run (docs/04).
void box_corners(const Box3& box, math::Vec3 (&out)[8])
{
    for (int i = 0; i < 8; ++i)
    {
        const math::Real sx = (i & 1) ? math::Real{1} : math::Real{-1};
        const math::Real sy = (i & 2) ? math::Real{1} : math::Real{-1};
        const math::Real sz = (i & 4) ? math::Real{1} : math::Real{-1};
        const math::Vec3 local{sx * box.half_extents.x, sy * box.half_extents.y,
                               sz * box.half_extents.z};
        out[i] = box.center + math::rotate(box.orientation, local);
    }
}
} // namespace

bool overlaps(const Box3& box, const Plane3& plane)
{
    if (!box.is_valid() || !plane.is_valid())
    {
        return false;
    }
    math::Vec3 corners[8];
    box_corners(box, corners);
    for (const math::Vec3& corner : corners)
    {
        if (signed_distance(plane, corner) <= math::Real{0}) // touching counts
        {
            return true;
        }
    }
    return false;
}

std::vector<Contact3> contacts(const Box3& box, const Plane3& plane)
{
    std::vector<Contact3> result;
    if (!box.is_valid() || !plane.is_valid())
    {
        return result;
    }
    math::Vec3 corners[8];
    box_corners(box, corners);
    for (const math::Vec3& corner : corners)
    {
        const math::Real distance = signed_distance(plane, corner);
        if (distance > math::Real{0})
        {
            continue; // this corner is in free space
        }
        const math::Real depth = -distance;
        Contact3 hit;
        hit.normal = -plane.normal; // from the box toward the solid, no fallback
        hit.penetration = depth;
        // Midway between the corner and the surface along the plane normal, the
        // same convention the sphere pairs use.
        hit.point = corner + plane.normal * (depth / math::Real{2});
        result.push_back(hit);
    }
    return result;
}

namespace
{
// An oriented box reduced to the working form the SAT wants: three orthonormal
// world-frame axes `u`, the matching half-extents `e`, and the centre `c`. A
// Box3's orientation columns ARE its world axes (`to_mat3` is built so its
// columns are the rotated basis vectors), so this is just a repackaging.
struct Obb
{
    math::Vec3 u[3];
    math::Real e[3];
    math::Vec3 c;
};

Obb make_obb(const Box3& box)
{
    const math::Mat3 r = math::to_mat3(box.orientation);
    return Obb{{r.col0, r.col1, r.col2},
               {box.half_extents.x, box.half_extents.y, box.half_extents.z},
               box.center};
}

// Half-width of a box's shadow cast onto direction `axis`: each half-extent
// times how much its own axis leans along `axis`. `axis` is unit for every axis
// tested (the face normals are orthonormal and the edge axes are normalised
// before use), so this radius and the centre gap are both true distances and
// directly comparable.
math::Real projected_radius(const Obb& o, const math::Vec3& axis)
{
    return o.e[0] * std::abs(math::dot(o.u[0], axis)) +
           o.e[1] * std::abs(math::dot(o.u[1], axis)) +
           o.e[2] * std::abs(math::dot(o.u[2], axis));
}

// The corner of the box farthest along direction `d`: step from the centre by
// each half-extent, toward `d` on every axis. Used to place a face contact on
// the deepest vertex of the opposing box.
math::Vec3 support(const Obb& o, const math::Vec3& d)
{
    math::Vec3 p = o.c;
    for (int k = 0; k < 3; ++k)
    {
        const math::Real s =
            math::dot(o.u[k], d) >= math::Real{0} ? math::Real{1} : math::Real{-1};
        p = p + o.u[k] * (s * o.e[k]);
    }
    return p;
}

// Centre of the box edge that runs along local axis `along` and sits farthest
// in direction `dir`: pushed to the extreme corner on the OTHER two axes but
// left centred along the edge itself. One of these per box gives the two line
// segments whose closest approach is the edge-edge contact point.
math::Vec3 edge_center(const Obb& o, int along, const math::Vec3& dir)
{
    math::Vec3 p = o.c;
    for (int k = 0; k < 3; ++k)
    {
        if (k == along)
        {
            continue;
        }
        const math::Real s =
            math::dot(o.u[k], dir) >= math::Real{0} ? math::Real{1} : math::Real{-1};
        p = p + o.u[k] * (s * o.e[k]);
    }
    return p;
}

// Midpoint of the shortest segment between two finite line segments, each given
// as a centre, a unit direction and a half-length (Ericson, Real-Time Collision
// Detection section 5.1.9, in centre+half-length form). Callers reach this only
// for edge axes that passed the near-parallel guard below, so the denominator
// 1 - cos^2 is safely positive and there is no divide by zero. The clamp is
// applied once rather than iterated: away from the edge ends the unconstrained
// solution is already interior and this is exact, and a contact that lands near
// an end only shifts the torque lever arm slightly, which a single bounce
// tolerates (the normal and penetration, which set the impulse, stay exact).
math::Vec3 closest_segment_midpoint(const math::Vec3& p1, const math::Vec3& d1,
                                    math::Real h1, const math::Vec3& p2,
                                    const math::Vec3& d2, math::Real h2)
{
    const math::Vec3 r = p1 - p2;
    const math::Real b = math::dot(d1, d2);
    const math::Real c = math::dot(d1, r);
    const math::Real f = math::dot(d2, r);
    const math::Real denom = math::Real{1} - b * b;
    math::Real s = (b * f - c) / denom;
    math::Real u = (f - b * c) / denom;
    s = std::clamp(s, -h1, h1);
    u = std::clamp(u, -h2, h2);
    const math::Vec3 c1 = p1 + d1 * s;
    const math::Vec3 c2 = p2 + d2 * u;
    return (c1 + c2) * math::Real{0.5};
}

// Edge cross products shorter than this (squared) come from near-parallel edges:
// the axis is numerically meaningless and whatever it would report is already
// covered by the six face axes, so it is skipped. Deliberately looser than the
// other squared tolerances here, because it guards the squared sine between two
// unit edge directions, not a distance.
constexpr math::Real kParallelEpsSq = 1e-8;

enum AxisKind
{
    FaceA,
    FaceB,
    Edge
};

// The separating axis of least overlap: the direction, how much the boxes
// overlap along it (the penetration once they are known to touch), which family
// it came from, and for an edge axis which local edge of each box produced it.
struct MinAxis
{
    math::Vec3 axis;
    math::Real overlap;
    AxisKind kind;
    int i;
    int j;
};

// The one SAT core behind BOTH `overlaps` and `contact`, so the two can never
// disagree. Tests all fifteen candidate axes (three face normals per box and
// the nine pairwise edge cross products); returns false the instant one
// separates the boxes (an overlap below zero), otherwise fills `out` with the
// axis of least overlap. Two coincident centres are not special-cased: every
// gap is zero, no axis separates, and the least-overlap axis is simply the one
// of least combined projected radius, a fixed and repeatable choice.
bool box_box_min_axis(const Box3& a, const Box3& b, MinAxis& out)
{
    const Obb box_a = make_obb(a);
    const Obb box_b = make_obb(b);
    const math::Vec3 t = b.center - a.center;

    bool have = false;
    auto test = [&](const math::Vec3& axis, AxisKind kind, int i, int j) -> bool
    {
        const math::Real gap = std::abs(math::dot(t, axis));
        const math::Real overlap =
            projected_radius(box_a, axis) + projected_radius(box_b, axis) - gap;
        if (overlap < math::Real{0})
        {
            return false; // a separating axis exists: the boxes are apart
        }
        if (!have || overlap < out.overlap)
        {
            out = MinAxis{axis, overlap, kind, i, j};
            have = true;
        }
        return true;
    };

    for (int k = 0; k < 3; ++k)
    {
        if (!test(box_a.u[k], FaceA, k, -1))
        {
            return false;
        }
    }
    for (int k = 0; k < 3; ++k)
    {
        if (!test(box_b.u[k], FaceB, -1, k))
        {
            return false;
        }
    }
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            const math::Vec3 raw = math::cross(box_a.u[i], box_b.u[j]);
            const math::Real len_sq = math::length_squared(raw);
            if (len_sq <= kParallelEpsSq)
            {
                continue; // near-parallel edges: the axis is degenerate, skip it
            }
            if (!test(raw / std::sqrt(len_sq), Edge, i, j))
            {
                return false;
            }
        }
    }
    return true;
}
} // namespace

bool overlaps(const Box3& a, const Box3& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return false;
    }
    MinAxis axis;
    return box_box_min_axis(a, b, axis);
}

std::optional<Contact3> contact(const Box3& a, const Box3& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return std::nullopt;
    }
    MinAxis m;
    if (!box_box_min_axis(a, b, m))
    {
        return std::nullopt;
    }

    const Obb box_a = make_obb(a);
    const Obb box_b = make_obb(b);
    const math::Vec3 t = b.center - a.center;

    // The axis is only a line; orient it from `a` toward `b`, like every 3D pair.
    // When the centres coincide the gap is zero, `dot(n, t)` is zero, and the
    // sign is left as computed: a fixed, repeatable choice, the box analogue of
    // the sphere pair's +x fallback.
    math::Vec3 n = m.axis;
    if (math::dot(n, t) < math::Real{0})
    {
        n = -n;
    }
    // A bare touch can round to a hair-negative overlap; never report a negative
    // penetration (docs/04), matching the sphere and box/plane pairs.
    const math::Real depth = std::max(math::Real{0}, m.overlap);

    Contact3 result;
    result.normal = n;
    result.penetration = depth;
    if (m.kind == Edge)
    {
        // The contact rides where the two edges cross: `a`'s edge pushed toward
        // `b` (along +n), `b`'s edge pushed toward `a` (along -n).
        const math::Vec3 p1 = edge_center(box_a, m.i, n);
        const math::Vec3 p2 = edge_center(box_b, m.j, -n);
        result.point = closest_segment_midpoint(p1, box_a.u[m.i], box_a.e[m.i], p2,
                                                 box_b.u[m.j], box_b.e[m.j]);
    }
    else if (m.kind == FaceA)
    {
        // A face of `a` separates: the contact is `b`'s deepest vertex into `a`
        // (its support along -n), nudged half the penetration back to the midway
        // point the other pairs report.
        result.point = support(box_b, -n) + n * (depth / math::Real{2});
    }
    else
    {
        // A face of `b` separates: `a`'s deepest vertex into `b` (support along
        // +n), nudged half the penetration back.
        result.point = support(box_a, n) - n * (depth / math::Real{2});
    }
    return result;
}

namespace
{
// The point of an oriented box nearest to `p`: clamp `p`, expressed in the box's
// local frame, into the half-extents on each axis, then map back to world. The
// building block of the box/sphere test (Ericson, Real-Time Collision Detection
// section 5.1.5).
math::Vec3 closest_point_on_obb(const Obb& o, const math::Vec3& p)
{
    const math::Vec3 d = p - o.c;
    math::Vec3 q = o.c;
    for (int k = 0; k < 3; ++k)
    {
        q = q + o.u[k] * std::clamp(math::dot(d, o.u[k]), -o.e[k], o.e[k]);
    }
    return q;
}
} // namespace

bool overlaps(const Box3& box, const Sphere& sphere)
{
    if (!box.is_valid() || !sphere.is_valid())
    {
        return false;
    }
    const math::Vec3 nearest = closest_point_on_obb(make_obb(box), sphere.center);
    const math::Vec3 delta = sphere.center - nearest;
    const math::Real distance_squared = math::dot(delta, delta);
    if (!math::is_finite(distance_squared))
    {
        return false; // so far away the squared distance overflowed
    }
    // Touching counts, matching every other query.
    return distance_squared <= sphere.radius * sphere.radius;
}

bool overlaps(const Sphere& sphere, const Box3& box)
{
    return overlaps(box, sphere); // the query is symmetric
}

std::optional<Contact3> contact(const Box3& box, const Sphere& sphere)
{
    if (!box.is_valid() || !sphere.is_valid())
    {
        return std::nullopt;
    }
    const Obb o = make_obb(box);
    const math::Vec3 d = sphere.center - o.c;

    // Signed distance along each box axis, and whether the centre lies within
    // the slab on every one of them: that, not the numeric distance, is what
    // separates "centre outside the box" from "centre inside", and it stays
    // exact where a distance-to-zero test would blur a grazing contact into the
    // degenerate one.
    math::Real projection[3];
    bool inside = true;
    for (int k = 0; k < 3; ++k)
    {
        projection[k] = math::dot(d, o.u[k]);
        if (std::abs(projection[k]) > o.e[k])
        {
            inside = false;
        }
    }

    math::Vec3 normal;
    math::Real penetration;
    math::Vec3 box_surface;
    if (!inside)
    {
        // The usual case: the nearest point is on the box surface, and the push
        // direction runs from it to the sphere's centre.
        const math::Vec3 nearest = closest_point_on_obb(o, sphere.center);
        const math::Vec3 delta = sphere.center - nearest;
        const math::Real distance_squared = math::dot(delta, delta);
        if (!math::is_finite(distance_squared) ||
            distance_squared > sphere.radius * sphere.radius)
        {
            return std::nullopt; // clear of the box, or too far to be finite
        }
        const math::Real distance = std::sqrt(distance_squared);
        normal = delta / distance; // from the box toward the sphere
        penetration = sphere.radius - distance;
        box_surface = nearest;
    }
    else
    {
        // The sphere's centre is inside the box (or exactly on a face): every
        // direction to the surface is valid, so push out through the LEAST
        // penetrated face, ties broken by axis order (0 then 1 then 2), the
        // deterministic analogue of the sphere pair's fixed normal.
        int axis = 0;
        math::Real shallowest = o.e[0] - std::abs(projection[0]);
        for (int k = 1; k < 3; ++k)
        {
            const math::Real depth_k = o.e[k] - std::abs(projection[k]);
            if (depth_k < shallowest)
            {
                shallowest = depth_k;
                axis = k;
            }
        }
        const math::Real sign =
            projection[axis] >= math::Real{0} ? math::Real{1} : math::Real{-1};
        normal = o.u[axis] * sign;                 // outward through that face
        penetration = sphere.radius + shallowest;  // its depth plus the radius
        box_surface = sphere.center; // the centre is its own nearest point
    }

    Contact3 result;
    result.normal = normal;
    result.penetration = penetration;
    // Midway between the two surfaces along the normal: the box's at
    // `box_surface`, the sphere's `radius` back along the normal from its centre.
    result.point =
        (box_surface + (sphere.center - normal * sphere.radius)) * math::Real{0.5};
    return result;
}

std::optional<Contact3> contact(const Sphere& sphere, const Box3& box)
{
    // The same contact with the roles swapped, so the normal points from the
    // sphere toward the box instead of from the box toward the sphere.
    std::optional<Contact3> hit = contact(box, sphere);
    if (hit)
    {
        hit->normal = -hit->normal;
    }
    return hit;
}
} // namespace malloy::collide
