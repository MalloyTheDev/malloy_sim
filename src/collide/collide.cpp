#include <malloy/collide/collide.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

#include <malloy/math/math.hpp>

namespace malloy::collide
{
namespace
{
// Used wherever the true normal is geometrically undefined. Fixed rather than
// chosen per call, so degenerate configurations stay repeatable (docs/04).
const math::Vec2 fallback_normal{1.0, 0.0};

// Not in <numbers> because the project targets C++20 without pulling in extras
// it does not otherwise use, and one constant does not justify the include.
const math::Real pi = math::Real{3.14159265358979323846};

math::Real clamp_to(math::Real value, math::Real low, math::Real high)
{
    return std::min(std::max(value, low), high);
}

// Overlap extent of two intervals; negative means they are separated.
math::Real interval_overlap(math::Real a_min, math::Real a_max, math::Real b_min,
                            math::Real b_max)
{
    return std::min(a_max, b_max) - std::max(a_min, b_min);
}

// Signed distance from a point to a halfplane's line: positive in free space,
// negative inside the solid, zero exactly on it. Only a true distance because
// the normal is required to be unit (shapes.hpp).
math::Real signed_distance(const Halfplane& plane, const math::Vec2& point)
{
    return math::dot(plane.normal, point) - plane.offset;
}

// The point on box closest to the given point. Equals it when it is inside.
math::Vec2 closest_point_on(const Aabb& box, const math::Vec2& point)
{
    return math::Vec2{clamp_to(point.x, box.min.x, box.max.x),
                      clamp_to(point.y, box.min.y, box.max.y)};
}
} // namespace

bool Halfplane::is_valid() const
{
    // A unit normal squared is 1. The tolerance is on the SQUARE, so it needs
    // no square root, and it is loose enough for a hand-written direction like
    // (0.6, 0.8) whose components are not exact in binary.
    const math::Real square = math::length_squared(normal);
    return math::is_finite(normal) && math::is_finite(offset) &&
           std::abs(square - math::Real{1}) <= math::Real{1e-12};
}

bool Circle::is_valid() const
{
    return radius >= math::Real{0} && math::is_finite(radius) &&
           math::is_finite(center);
}

bool Aabb::is_valid() const
{
    // NaN fails both the finite check and the ordering check, so it is rejected
    // whichever way it appears.
    return math::is_finite(min) && math::is_finite(max) && min.x <= max.x &&
           min.y <= max.y;
}

bool Obb2::is_valid() const
{
    return half_extents.x > math::Real{0} && half_extents.y > math::Real{0} &&
           math::is_finite(half_extents) && math::is_finite(center) &&
           math::is_finite(orientation);
}

// ----------------------------------------------------------------------------
// Area properties
// ----------------------------------------------------------------------------

math::Real area(const Circle& c)
{
    if (!c.is_valid())
    {
        return math::Real{0};
    }
    return pi * c.radius * c.radius;
}

math::Real area(const Aabb& box)
{
    if (!box.is_valid())
    {
        return math::Real{0};
    }
    return (box.max.x - box.min.x) * (box.max.y - box.min.y);
}

math::Vec2 centroid(const Circle& c)
{
    if (!c.is_valid())
    {
        return math::Vec2{};
    }
    return c.center;
}

math::Vec2 centroid(const Aabb& box)
{
    if (!box.is_valid())
    {
        return math::Vec2{};
    }
    return (box.min + box.max) / math::Real{2};
}

math::Real second_moment_of_area(const Circle& c)
{
    if (!c.is_valid())
    {
        return math::Real{0};
    }
    // pi*R^4/2
    const math::Real r2 = c.radius * c.radius;
    return pi * r2 * r2 / math::Real{2};
}

math::Real second_moment_of_area(const Aabb& box)
{
    if (!box.is_valid())
    {
        return math::Real{0};
    }
    const math::Real w = box.max.x - box.min.x;
    const math::Real h = box.max.y - box.min.y;
    return w * h * (w * w + h * h) / math::Real{12};
}

// ----------------------------------------------------------------------------
// Overlap tests
// ----------------------------------------------------------------------------

bool overlaps(const Circle& a, const Circle& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return false;
    }
    const math::Vec2 delta = b.center - a.center;
    const math::Real distance_squared = math::dot(delta, delta);
    if (!math::is_finite(distance_squared))
    {
        return false; // so far apart the squared distance overflowed
    }
    const math::Real radius_sum = a.radius + b.radius;
    return distance_squared <= radius_sum * radius_sum;
}

bool overlaps(const Aabb& a, const Aabb& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return false;
    }
    return interval_overlap(a.min.x, a.max.x, b.min.x, b.max.x) >= math::Real{0} &&
           interval_overlap(a.min.y, a.max.y, b.min.y, b.max.y) >= math::Real{0};
}

bool overlaps(const Circle& circle, const Aabb& box)
{
    if (!circle.is_valid() || !box.is_valid())
    {
        return false;
    }
    const math::Vec2 delta = closest_point_on(box, circle.center) - circle.center;
    const math::Real distance_squared = math::dot(delta, delta);
    if (!math::is_finite(distance_squared))
    {
        return false;
    }
    return distance_squared <= circle.radius * circle.radius;
}

// ----------------------------------------------------------------------------
// Contact queries
// ----------------------------------------------------------------------------

bool overlaps(const Circle& circle, const Halfplane& plane)
{
    if (!circle.is_valid() || !plane.is_valid())
    {
        return false;
    }
    // Touching counts, matching every other query here.
    return signed_distance(plane, circle.center) <= circle.radius;
}

std::optional<Contact> contact(const Circle& circle, const Halfplane& plane)
{
    if (!circle.is_valid() || !plane.is_valid())
    {
        return std::nullopt;
    }

    const math::Real distance = signed_distance(plane, circle.center);
    const math::Real depth = circle.radius - distance;
    if (depth < math::Real{0})
    {
        return std::nullopt;
    }

    Contact result;
    // From the circle toward the plane's solid side, which is the direction the
    // plane's outward normal does not point. No fallback is needed even when
    // the centre lies exactly on the line: the plane supplies the direction.
    result.normal = -plane.normal;
    result.penetration = depth;
    // Midway between the two surfaces along the normal, as for every other
    // pair: the circle's surface is at distance -radius along the normal from
    // its centre, and the plane's surface is `depth` further back.
    result.point = circle.center + result.normal * (circle.radius - depth / math::Real{2});
    return result;
}

std::optional<Contact> contact(const Circle& a, const Circle& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return std::nullopt;
    }

    const math::Vec2 delta = b.center - a.center;
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

    Contact result;
    if (distance_squared > math::Real{0})
    {
        const math::Real distance = std::sqrt(distance_squared);
        result.normal = delta / distance;
        result.penetration = radius_sum - distance;
    }
    else
    {
        // Coincident centers: every direction separates them equally well, so
        // take the documented one instead of dividing by zero.
        result.normal = fallback_normal;
        result.penetration = radius_sum;
    }
    result.point =
        a.center + result.normal * (a.radius - result.penetration / math::Real{2});
    return result;
}

std::optional<Contact> contact(const Aabb& a, const Aabb& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return std::nullopt;
    }

    const math::Real overlap_x = interval_overlap(a.min.x, a.max.x, b.min.x, b.max.x);
    const math::Real overlap_y = interval_overlap(a.min.y, a.max.y, b.min.y, b.max.y);
    if (overlap_x < math::Real{0} || overlap_y < math::Real{0})
    {
        return std::nullopt;
    }

    const math::Vec2 a_center = (a.min + a.max) / math::Real{2};
    const math::Vec2 b_center = (b.min + b.max) / math::Real{2};

    Contact result;
    // Separate along the axis of least penetration: the shortest way out. A tie
    // goes to x, so an exactly symmetric overlap still gives one fixed answer.
    if (overlap_x <= overlap_y)
    {
        const math::Real direction =
            (b_center.x >= a_center.x) ? math::Real{1} : math::Real{-1};
        result.normal = math::Vec2{direction, math::Real{0}};
        result.penetration = overlap_x;
    }
    else
    {
        const math::Real direction =
            (b_center.y >= a_center.y) ? math::Real{1} : math::Real{-1};
        result.normal = math::Vec2{math::Real{0}, direction};
        result.penetration = overlap_y;
    }

    // Center of the overlapping rectangle.
    result.point = math::Vec2{
        (std::max(a.min.x, b.min.x) + std::min(a.max.x, b.max.x)) / math::Real{2},
        (std::max(a.min.y, b.min.y) + std::min(a.max.y, b.max.y)) / math::Real{2}};
    return result;
}

std::optional<Contact> contact(const Circle& circle, const Aabb& box)
{
    if (!circle.is_valid() || !box.is_valid())
    {
        return std::nullopt;
    }

    const math::Vec2 closest = closest_point_on(box, circle.center);
    const math::Vec2 delta = closest - circle.center; // circle toward box
    const math::Real distance_squared = math::dot(delta, delta);
    if (!math::is_finite(distance_squared))
    {
        return std::nullopt;
    }

    Contact result;
    if (distance_squared > math::Real{0})
    {
        // Center outside the box: the closest face point gives the normal.
        const math::Real distance = std::sqrt(distance_squared);
        if (distance > circle.radius)
        {
            return std::nullopt;
        }
        result.normal = delta / distance;
        result.penetration = circle.radius - distance;
    }
    else
    {
        // Center inside the box, so they always overlap. Leave through the
        // nearest face: the circle must travel that far plus its own radius to
        // clear the box entirely. Comparisons are strict, so a tie keeps the
        // earlier axis and an exactly centered circle gets +x.
        const math::Real to_left = circle.center.x - box.min.x;
        const math::Real to_right = box.max.x - circle.center.x;
        const math::Real to_bottom = circle.center.y - box.min.y;
        const math::Real to_top = box.max.y - circle.center.y;

        math::Real nearest = to_left;
        result.normal = math::Vec2{math::Real{1}, math::Real{0}};
        if (to_right < nearest)
        {
            nearest = to_right;
            result.normal = math::Vec2{math::Real{-1}, math::Real{0}};
        }
        if (to_bottom < nearest)
        {
            nearest = to_bottom;
            result.normal = math::Vec2{math::Real{0}, math::Real{1}};
        }
        if (to_top < nearest)
        {
            nearest = to_top;
            result.normal = math::Vec2{math::Real{0}, math::Real{-1}};
        }
        result.penetration = nearest + circle.radius;
    }

    // The same representative point as circle/circle: midway between the two
    // surfaces along the normal.
    result.point = circle.center +
                   result.normal * (circle.radius - result.penetration / math::Real{2});
    return result;
}

namespace
{
// An oriented 2D box reduced to the working form the SAT wants: two orthonormal
// axes, the matching half-widths, and the centre. The 2D sibling of the `Obb`
// helper behind the box/box test in 3D.
struct Obb2Frame
{
    math::Vec2 u[2];
    math::Real e[2];
    math::Vec2 c;
};

Obb2Frame make_obb2(const Obb2& box)
{
    const math::Real cs = std::cos(box.orientation);
    const math::Real sn = std::sin(box.orientation);
    return Obb2Frame{{math::Vec2{cs, sn}, math::Vec2{-sn, cs}},
                     {box.half_extents.x, box.half_extents.y}, box.center};
}

// Half-width of the box's shadow on `axis` (unit for every axis tested here), so
// this radius and the centre gap are both true distances and comparable.
math::Real projected_radius(const Obb2Frame& o, const math::Vec2& axis)
{
    return o.e[0] * std::abs(math::dot(o.u[0], axis)) +
           o.e[1] * std::abs(math::dot(o.u[1], axis));
}

// The corner of the box farthest along `d`: step from the centre by each
// half-width toward `d`. Places the face contact on the deepest vertex.
math::Vec2 support(const Obb2Frame& o, const math::Vec2& d)
{
    math::Vec2 p = o.c;
    for (int k = 0; k < 2; ++k)
    {
        const math::Real s =
            math::dot(o.u[k], d) >= math::Real{0} ? math::Real{1} : math::Real{-1};
        p = p + o.u[k] * (s * o.e[k]);
    }
    return p;
}

enum Axis2Kind
{
    FaceA2,
    FaceB2
};

struct MinAxis2
{
    math::Vec2 axis;
    math::Real overlap;
    Axis2Kind kind;
};

// The one SAT core behind both `overlaps` and `contact`, so they cannot
// disagree. Tests the four face axes (two per box); returns false the instant
// one separates the boxes, otherwise fills `out` with the axis of least overlap.
// Two coincident centres are not special-cased: every gap is zero, no axis
// separates, and the least-overlap axis is the one of least combined width, a
// fixed repeatable choice.
bool box2_min_axis(const Obb2& a, const Obb2& b, MinAxis2& out)
{
    const Obb2Frame fa = make_obb2(a);
    const Obb2Frame fb = make_obb2(b);
    const math::Vec2 t = b.center - a.center;

    bool have = false;
    auto test = [&](const math::Vec2& axis, Axis2Kind kind) -> bool
    {
        const math::Real gap = std::abs(math::dot(t, axis));
        const math::Real overlap =
            projected_radius(fa, axis) + projected_radius(fb, axis) - gap;
        if (overlap < math::Real{0})
        {
            return false; // a separating axis: the boxes are apart
        }
        if (!have || overlap < out.overlap)
        {
            out = MinAxis2{axis, overlap, kind};
            have = true;
        }
        return true;
    };

    for (int k = 0; k < 2; ++k)
    {
        if (!test(fa.u[k], FaceA2))
        {
            return false;
        }
    }
    for (int k = 0; k < 2; ++k)
    {
        if (!test(fb.u[k], FaceB2))
        {
            return false;
        }
    }
    return true;
}
} // namespace

bool overlaps(const Obb2& a, const Obb2& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return false;
    }
    MinAxis2 axis;
    return box2_min_axis(a, b, axis);
}

std::optional<Contact> contact(const Obb2& a, const Obb2& b)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return std::nullopt;
    }
    MinAxis2 m;
    if (!box2_min_axis(a, b, m))
    {
        return std::nullopt;
    }

    const math::Vec2 t = b.center - a.center;
    // Orient the axis from a toward b, like every pair. Coincident centres leave
    // the sign as computed, a fixed choice (the box analogue of +x).
    math::Vec2 n = m.axis;
    if (math::dot(n, t) < math::Real{0})
    {
        n = -n;
    }
    // A bare touch can round to a hair-negative overlap; never report negative.
    const math::Real depth = std::max(math::Real{0}, m.overlap);

    Contact result;
    result.normal = n;
    result.penetration = depth;
    if (m.kind == FaceA2)
    {
        // A face of a separates: contact on b's deepest vertex into a (its
        // support along -n), nudged half the penetration back to the midpoint.
        result.point = support(make_obb2(b), -n) + n * (depth / math::Real{2});
    }
    else
    {
        // A face of b separates: a's deepest vertex into b (support along +n).
        result.point = support(make_obb2(a), n) - n * (depth / math::Real{2});
    }
    return result;
}

namespace
{
// The four corners of an oriented 2D box, in a fixed order (the sign bits of its
// two local axes), so any loop over them repeats a run (docs/04). The 2D copy of
// `box_corners` in the 3D module.
void box2_corners(const Obb2& box, math::Vec2 (&out)[4])
{
    const Obb2Frame f = make_obb2(box);
    for (int i = 0; i < 4; ++i)
    {
        const math::Real sx = (i & 1) ? math::Real{1} : math::Real{-1};
        const math::Real sy = (i & 2) ? math::Real{1} : math::Real{-1};
        out[i] = f.c + f.u[0] * (sx * f.e[0]) + f.u[1] * (sy * f.e[1]);
    }
}
} // namespace

bool overlaps(const Obb2& box, const Halfplane& plane)
{
    if (!box.is_valid() || !plane.is_valid())
    {
        return false;
    }
    math::Vec2 corners[4];
    box2_corners(box, corners);
    for (const math::Vec2& corner : corners)
    {
        if (signed_distance(plane, corner) <= math::Real{0}) // touching counts
        {
            return true;
        }
    }
    return false;
}

std::vector<Contact> contacts(const Obb2& box, const Halfplane& plane)
{
    std::vector<Contact> result;
    if (!box.is_valid() || !plane.is_valid())
    {
        return result;
    }
    math::Vec2 corners[4];
    box2_corners(box, corners);
    for (const math::Vec2& corner : corners)
    {
        const math::Real distance = signed_distance(plane, corner);
        if (distance > math::Real{0})
        {
            continue; // this corner is in free space
        }
        const math::Real depth = -distance;
        Contact hit;
        hit.normal = -plane.normal; // from the box toward the solid, no fallback
        hit.penetration = depth;
        // Midway between the corner and the surface along the plane normal, the
        // same convention every pair uses.
        hit.point = corner + plane.normal * (depth / math::Real{2});
        result.push_back(hit);
    }
    return result;
}
} // namespace malloy::collide
