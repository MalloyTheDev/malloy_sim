#include <malloy/view3d/mesh.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include <malloy/math/quat.hpp>
#include <malloy/math/vec3.hpp>
#include <malloy/view3d/view3d.hpp>

namespace malloy::view3d
{
namespace
{
constexpr Real pi = 3.14159265358979323846;

// The furthest a grid line may sit from the plane's anchor, in spacings,
// before the spacing is lost to rounding (2^40, far below 2^53).
constexpr Real grid_reach = 1099511627776.0;
} // namespace

std::array<Triangle, 12> box_triangles(const Box& box)
{
    std::array<Triangle, 12> out{};
    for (std::size_t face = 0; face < box_faces.size(); ++face)
    {
        // Faces come in the order -x, +x, -y, +y, -z, +z.
        const Real sign = (face % 2 == 0) ? Real{-1} : Real{1};
        const std::size_t axis = face / 2;
        const Vec3 local{axis == 0 ? sign : Real{0}, axis == 1 ? sign : Real{0},
                         axis == 2 ? sign : Real{0}};
        const Vec3 normal = math::rotate(box.orientation, local);
        const auto corner = [&](std::size_t k) {
            return box.corners[static_cast<std::size_t>(box_faces[face][k])];
        };
        out[2 * face] = Triangle{corner(0), corner(1), corner(2), normal};
        out[2 * face + 1] = Triangle{corner(0), corner(2), corner(3), normal};
    }
    return out;
}

SphereDetail sphere_detail(std::size_t sphere_count)
{
    if (sphere_count <= 64)
    {
        return SphereDetail{16, 24};
    }
    if (sphere_count <= 512)
    {
        return SphereDetail{10, 16};
    }
    if (sphere_count <= 4096)
    {
        return SphereDetail{6, 10};
    }
    return SphereDetail{4, 6};
}

std::vector<Triangle> sphere_triangles(const Vec3& center, Real radius,
                                       const SphereDetail& detail)
{
    const int rings = std::max(detail.rings, 2);
    const int slices = std::max(detail.slices, 3);
    std::vector<Triangle> out;
    out.reserve(static_cast<std::size_t>(slices) * static_cast<std::size_t>(2 * rings - 2));

    // Vertex (i, j): latitude band i from the +z pole, longitude j. The poles
    // are exact and j wraps, so every shared vertex is computed identically and
    // the mesh closes bit for bit.
    const auto vertex = [&](int i, int j) {
        if (i == 0)
        {
            return center + Vec3{0.0, 0.0, radius};
        }
        if (i == rings)
        {
            return center + Vec3{0.0, 0.0, -radius};
        }
        const Real theta = pi * static_cast<Real>(i) / static_cast<Real>(rings);
        const Real phi = Real{2} * pi * static_cast<Real>(j % slices) / static_cast<Real>(slices);
        return center + Vec3{std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
                             std::cos(theta)} *
                            radius;
    };
    const auto push = [&](const Vec3& a, const Vec3& b, const Vec3& c) {
        out.push_back(Triangle{a, b, c, math::normalize(math::cross(b - a, c - a))});
    };
    for (int i = 0; i < rings; ++i)
    {
        for (int j = 0; j < slices; ++j)
        {
            // Seen from outside, north up and east (growing phi) to the right:
            // p00 top-left, p01 top-right, p10 bottom-left, p11 bottom-right,
            // so (p00, p10, p11) and (p00, p11, p01) run counter-clockwise. A
            // pole band keeps only the triangle that is not degenerate.
            const Vec3 p00 = vertex(i, j);
            const Vec3 p01 = vertex(i, j + 1);
            const Vec3 p10 = vertex(i + 1, j);
            const Vec3 p11 = vertex(i + 1, j + 1);
            if (i != 0)
            {
                push(p00, p11, p01);
            }
            if (i != rings - 1)
            {
                push(p00, p10, p11);
            }
        }
    }
    return out;
}

void place_sphere(const std::vector<Triangle>& unit, const Vec3& center, Real radius,
                  std::vector<Triangle>& out)
{
    // No reserve here: reserving exactly what one more sphere needs would
    // reallocate on every call, and a frame of n spheres would cost n^2.
    for (const Triangle& t : unit)
    {
        out.push_back(Triangle{center + t.a * radius, center + t.b * radius,
                               center + t.c * radius, t.normal});
    }
}

std::array<Triangle, 2> patch_triangles(const Patch& patch)
{
    const std::array<Vec3, 4>& c = patch.corners;
    return {Triangle{c[0], c[1], c[2], patch.normal}, Triangle{c[0], c[2], c[3], patch.normal}};
}

Real grid_spacing(Real span, int lines)
{
    if (!math::is_finite(span) || !(span > Real{0}) || lines < 1)
    {
        return Real{0};
    }
    const Real raw = span / static_cast<Real>(lines);
    const Real power = std::pow(Real{10}, std::floor(std::log10(raw)));
    const Real mantissa = raw / power;
    const Real nice = mantissa < Real{1.5}   ? Real{1}
                      : mantissa < Real{3.5} ? Real{2}
                      : mantissa < Real{7.5} ? Real{5}
                                             : Real{10};
    const Real spacing = nice * power;
    return math::is_finite(spacing) && spacing > Real{0} ? spacing : Real{0};
}

std::vector<Segment> patch_grid(const Patch& patch)
{
    std::vector<Segment> lines;
    const Real h = patch.half_size;
    const Real spacing = grid_spacing(Real{2} * h, 20);
    if (!(spacing > Real{0}) || !math::is_finite(patch.offset) ||
        !math::is_finite(patch.center) || !math::is_finite(patch.u) ||
        !math::is_finite(patch.v) || !math::is_finite(patch.normal))
    {
        return lines;
    }
    // Measured from the plane's own anchor, the point nearest the world origin,
    // so a line's place in the world does not depend on where the patch is.
    const Vec3 anchor = patch.normal * patch.offset;
    const Vec3 along[2] = {patch.u, patch.v};
    const Vec3 across[2] = {patch.v, patch.u};
    for (std::size_t k = 0; k < 2; ++k)
    {
        const Real at = math::dot(patch.center - anchor, along[k]);
        const Real first = std::ceil((at - h) / spacing);
        const Real last = std::floor((at + h) / spacing);
        if (!math::is_finite(first) || !math::is_finite(last) ||
            std::abs(first) > grid_reach || std::abs(last) > grid_reach ||
            last - first > Real{1000})
        {
            return {};
        }
        for (Real n = first; n <= last; n += Real{1})
        {
            const Vec3 middle = patch.center + along[k] * (n * spacing - at);
            lines.push_back(Segment{middle - across[k] * h, middle + across[k] * h, 0.0,
                                    patch.source, Role::Ground});
        }
    }
    return lines;
}

std::vector<Segment> arrow_lines(const Arrow& arrow)
{
    const Vec3 shaft = arrow.tip - arrow.origin;
    const Real length = math::length(shaft);
    if (!math::is_finite(length) || !(length > Real{0}))
    {
        return {};
    }
    const Vec3 direction = shaft / length;
    Vec3 u;
    Vec3 v;
    tangents(direction, u, v);
    const Real head = Real{0.2} * length;
    const Vec3 base = arrow.tip - direction * head;
    const Real spread = Real{0.35} * head;
    return {Segment{arrow.origin, arrow.tip, 0.0, 0, arrow.role},
            Segment{arrow.tip, base + u * spread, 0.0, 0, arrow.role},
            Segment{arrow.tip, base - u * spread, 0.0, 0, arrow.role},
            Segment{arrow.tip, base + v * spread, 0.0, 0, arrow.role},
            Segment{arrow.tip, base - v * spread, 0.0, 0, arrow.role}};
}

Real shade(const Vec3& normal, const Vec3& to_light, Real ambient)
{
    return ambient + (Real{1} - ambient) * std::max(Real{0}, math::dot(normal, to_light));
}
} // namespace malloy::view3d
