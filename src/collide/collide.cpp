#include <malloy/collide/collide.hpp>

#include <algorithm>
#include <cmath>
#include <optional>

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

// The point on box closest to the given point. Equals it when it is inside.
math::Vec2 closest_point_on(const Aabb& box, const math::Vec2& point)
{
    return math::Vec2{clamp_to(point.x, box.min.x, box.max.x),
                      clamp_to(point.y, box.min.y, box.max.y)};
}
} // namespace

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
} // namespace malloy::collide
