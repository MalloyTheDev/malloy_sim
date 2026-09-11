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

namespace
{
// Signed distance from a point to the plane's surface: positive in free space,
// negative inside the solid. The 3D copy of the 2D helper, one component wider.
math::Real signed_distance(const Plane3& plane, const math::Vec3& point)
{
    return math::dot(plane.normal, point) - plane.offset;
}
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
} // namespace malloy::collide
