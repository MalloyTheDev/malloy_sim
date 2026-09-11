#include <malloy/rigid/mass_properties3d.hpp>

#include <vector>

#include <malloy/math/math.hpp>

namespace malloy::rigid
{
namespace
{
// 4/3 pi, the volume of a unit sphere.
constexpr math::Real four_thirds_pi = math::Real{4.18879020478639098461685784437};

bool part_is_valid(const SolidSphere& part)
{
    return part.radius > math::Real{0} && math::is_finite(part.radius) &&
           part.density > math::Real{0} && math::is_finite(part.density) &&
           math::is_squarable(part.center);
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
        if (!part_is_valid(part))
        {
            return MassProperties3D{};
        }
    }

    // Mass and centre of mass.
    math::Real total_mass = 0.0;
    math::Vec3 weighted_center{};
    for (const SolidSphere& part : parts)
    {
        const math::Real mass =
            part.density * four_thirds_pi * part.radius * part.radius * part.radius;
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
    // I += I_own + m (|d|^2 I - d d^T), with d the offset from the common
    // centre to the part's centre.
    math::Mat3 tensor{};
    for (const SolidSphere& part : parts)
    {
        const math::Real mass =
            part.density * four_thirds_pi * part.radius * part.radius * part.radius;
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

    // Diagonalize: the principal moments and the orientation of the principal
    // frame in the lab frame.
    const math::SymmetricEigen eigen = math::eigen_symmetric(tensor);
    MassProperties3D result;
    result.mass = total_mass;
    result.center_of_mass = center_of_mass;
    result.inertia = eigen.values;
    result.orientation = math::to_quat(eigen.vectors);
    return result;
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
