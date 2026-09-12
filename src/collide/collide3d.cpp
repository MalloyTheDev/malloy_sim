#include <malloy/collide/contact3d.hpp>
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
} // namespace malloy::collide
