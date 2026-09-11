#include <malloy/rigid/mass_properties3d.hpp>

#include <vector>

#include <malloy/math/math.hpp>

namespace malloy::rigid
{
namespace
{
// 4/3 pi, the volume of a unit sphere.
constexpr math::Real four_thirds_pi = math::Real{4.18879020478639098461685784437};

bool sphere_is_valid(const SolidSphere& part)
{
    return part.radius > math::Real{0} && math::is_finite(part.radius) &&
           part.density > math::Real{0} && math::is_finite(part.density) &&
           math::is_squarable(part.center);
}

bool box_is_valid(const SolidBox& box)
{
    return box.half_extents.x > math::Real{0} && box.half_extents.y > math::Real{0} &&
           box.half_extents.z > math::Real{0} && math::is_finite(box.half_extents) &&
           box.density > math::Real{0} && math::is_finite(box.density) &&
           math::is_unit(box.orientation) && math::is_squarable(box.center);
}

bool mesh_is_valid(const SolidMesh& mesh)
{
    // A closed mesh needs at least a tetrahedron: four vertices, four faces.
    if (mesh.vertices.size() < 4 || mesh.triangles.size() < 4)
    {
        return false;
    }
    if (!(mesh.density > math::Real{0}) || !math::is_finite(mesh.density))
    {
        return false;
    }
    for (const math::Vec3& vertex : mesh.vertices)
    {
        if (!math::is_squarable(vertex))
        {
            return false;
        }
    }
    for (const MeshTriangle& triangle : mesh.triangles)
    {
        if (triangle.v0 >= mesh.vertices.size() ||
            triangle.v1 >= mesh.vertices.size() ||
            triangle.v2 >= mesh.vertices.size())
        {
            return false;
        }
    }
    return true;
}

// Turn a mass, a centre of mass, and an inertia tensor about that centre into
// principal moments and a principal frame. The one place the diagonalization
// happens, shared by every builder.
MassProperties3D finalize(math::Real mass, const math::Vec3& center_of_mass,
                          const math::Mat3& tensor)
{
    const math::SymmetricEigen eigen = math::eigen_symmetric(tensor);
    MassProperties3D result;
    result.mass = mass;
    result.center_of_mass = center_of_mass;
    result.inertia = eigen.values;
    result.orientation = math::to_quat(eigen.vectors);
    return result;
}

math::Real sphere_mass(const SolidSphere& part)
{
    return part.density * four_thirds_pi * part.radius * part.radius * part.radius;
}
} // namespace

MassProperties3D mass_properties_3d(const std::vector<SolidSphere>& parts)
{
    if (parts.empty())
    {
        return MassProperties3D{};
    }
    for (const SolidSphere& part : parts)
    {
        if (!sphere_is_valid(part))
        {
            return MassProperties3D{};
        }
    }

    math::Real total_mass = 0.0;
    math::Vec3 weighted_center{};
    for (const SolidSphere& part : parts)
    {
        const math::Real mass = sphere_mass(part);
        total_mass += mass;
        weighted_center += part.center * mass;
    }
    if (!(total_mass > math::Real{0}) || !math::is_finite(total_mass))
    {
        return MassProperties3D{};
    }
    const math::Vec3 center_of_mass = weighted_center / total_mass;

    // Inertia tensor about the centre of mass: each sphere's own isotropic
    // tensor, shifted to the common centre by the parallel-axis theorem
    // I += I_own + m (|d|^2 I - d d^T).
    math::Mat3 tensor{};
    for (const SolidSphere& part : parts)
    {
        const math::Real mass = sphere_mass(part);
        const math::Real own = math::Real{0.4} * mass * part.radius * part.radius;
        const math::Vec3 offset = part.center - center_of_mass;
        const math::Real distance_squared = math::dot(offset, offset);
        const math::Mat3 shift =
            (math::diagonal3(math::Vec3{distance_squared, distance_squared,
                                        distance_squared}) -
             math::outer(offset, offset)) *
            mass;
        tensor = tensor + math::diagonal3(math::Vec3{own, own, own}) + shift;
    }
    return finalize(total_mass, center_of_mass, tensor);
}

MassProperties3D mass_properties_3d(const SolidBox& box)
{
    if (!box_is_valid(box))
    {
        return MassProperties3D{};
    }
    const math::Vec3 h = box.half_extents;
    const math::Real mass = box.density * math::Real{8} * h.x * h.y * h.z;
    // Solid box moments about its own axes: (1/12) m (w_j^2 + w_k^2) with w the
    // FULL widths (2h), which is (1/3) m (h_j^2 + h_k^2).
    const math::Real third = math::Real{1} / math::Real{3};
    const math::Vec3 own{third * mass * (h.y * h.y + h.z * h.z),
                         third * mass * (h.x * h.x + h.z * h.z),
                         third * mass * (h.x * h.x + h.y * h.y)};
    // Carry the box-frame diagonal into the lab frame: R diag(own) R^T.
    const math::Mat3 rotation = math::to_mat3(box.orientation);
    const math::Mat3 tensor =
        rotation * math::diagonal3(own) * math::transpose(rotation);
    return finalize(mass, box.center, tensor);
}

MassProperties3D mass_properties_3d(const SolidMesh& mesh)
{
    if (!mesh_is_valid(mesh))
    {
        return MassProperties3D{};
    }

    // Signed sums over the tetrahedra (origin, a, b, c), one per triangle. Each
    // is weighted by six times its signed volume, det = a . (b x c). A closed
    // outward-wound mesh makes these sums equal the integrals over the enclosed
    // solid, whichever side of the origin the mesh lies on: the parts outside
    // the solid cancel between overlapping signed tetrahedra.
    math::Real volume_times_six = 0.0;      // sum det           = 6 V
    math::Vec3 centroid_numerator{};        // sum det (a+b+c)   = 24 V c
    math::Mat3 covariance_times_120{};      // sum det (...)     = 120 C

    for (const MeshTriangle& triangle : mesh.triangles)
    {
        const math::Vec3& a = mesh.vertices[triangle.v0];
        const math::Vec3& b = mesh.vertices[triangle.v1];
        const math::Vec3& c = mesh.vertices[triangle.v2];
        const math::Real det = math::dot(a, math::cross(b, c));
        const math::Vec3 sum = a + b + c;
        volume_times_six += det;
        centroid_numerator += sum * det;
        // The canonical-tetrahedron covariance, mapped to this tetrahedron:
        // integral x x^T over (0,a,b,c) is (det / 120) (a a^T + b b^T + c c^T +
        // (a+b+c)(a+b+c)^T).
        const math::Mat3 moment = math::outer(a, a) + math::outer(b, b) +
                                  math::outer(c, c) + math::outer(sum, sum);
        covariance_times_120 = covariance_times_120 + moment * det;
    }

    const math::Real volume = volume_times_six / math::Real{6};
    if (!(volume > math::Real{0}) || !math::is_finite(volume))
    {
        return MassProperties3D{};
    }
    const math::Vec3 center_of_mass =
        centroid_numerator / (math::Real{24} * volume);

    // Covariance about the origin, C = integral x x^T dV, then the inertia
    // tensor about the origin, I = trace(C) I - C.
    const math::Mat3 covariance =
        covariance_times_120 * (math::Real{1} / math::Real{120});
    const math::Mat3 inertia_origin =
        math::identity3() * math::trace(covariance) - covariance;
    // Parallel-axis to the centre of mass. These are geometric integrals, so
    // the "mass" in the shift is the volume: I_com = I_origin - V (|d|^2 I - d d^T).
    const math::Vec3 offset = center_of_mass;
    const math::Real distance_squared = math::dot(offset, offset);
    const math::Mat3 shift =
        (math::diagonal3(math::Vec3{distance_squared, distance_squared,
                                    distance_squared}) -
         math::outer(offset, offset)) *
        volume;
    const math::Mat3 inertia_com_geometric = inertia_origin - shift;

    // Density scales the mass and every second moment linearly.
    const math::Real mass = mesh.density * volume;
    const math::Mat3 tensor = inertia_com_geometric * mesh.density;
    return finalize(mass, center_of_mass, tensor);
}

MassProperties3D combine(const MassProperties3D& a, const MassProperties3D& b)
{
    // A zero-mass operand is the identity, so a fold can start from nothing.
    if (!(a.mass > math::Real{0}))
    {
        return b;
    }
    if (!(b.mass > math::Real{0}))
    {
        return a;
    }

    const math::Real total_mass = a.mass + b.mass;
    const math::Vec3 center_of_mass =
        (a.center_of_mass * a.mass + b.center_of_mass * b.mass) / total_mass;

    // Each part's tensor about ITS OWN centre, reconstructed from its principal
    // moments and frame, then shifted to the shared centre by parallel-axis.
    const auto shifted = [&](const MassProperties3D& part) {
        const math::Mat3 rotation = math::to_mat3(part.orientation);
        const math::Mat3 own =
            rotation * math::diagonal3(part.inertia) * math::transpose(rotation);
        const math::Vec3 offset = part.center_of_mass - center_of_mass;
        const math::Real distance_squared = math::dot(offset, offset);
        return own + (math::diagonal3(math::Vec3{distance_squared, distance_squared,
                                                 distance_squared}) -
                      math::outer(offset, offset)) *
                         part.mass;
    };
    return finalize(total_mass, center_of_mass, shifted(a) + shifted(b));
}

RigidBody3D rigid_body_from(const MassProperties3D& properties)
{
    RigidBody3D body;
    body.position = properties.center_of_mass;
    body.mass = properties.mass;
    body.inertia = properties.inertia;
    body.orientation = properties.orientation;
    return body;
}
} // namespace malloy::rigid
