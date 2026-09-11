#include <malloy/rigid/rigid.hpp>

#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include <malloy/collide/collide.hpp>
#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::rigid
{
namespace
{
// Resolve one sphere body against one immovable ground plane: a normal impulse
// plus a positional correction, and nothing else. This is M22's whole contact
// response, and it is deliberately the translational one.
//
// A centred sphere's contact normal passes through its centre of mass, so the
// arm r from the centre of mass to the contact point is parallel to the normal
// and r x n = 0. A NORMAL impulse therefore imparts no angular velocity, and
// none of the rotational effective-mass machinery is needed or written: the
// effective mass is just 1/m against an immovable plane. Spin from a contact is
// friction's job (a tangential impulse has r x t != 0), a later milestone.
void resolve_ground3d(RigidBody3D& body, const collide::Plane3& plane,
                      const Rigid3DSettings& settings)
{
    if (!(body.radius > math::Real{0})) // zero radius does not collide
    {
        return;
    }
    const std::optional<collide::Contact3> hit =
        collide::contact(collide::Sphere{body.position, body.radius}, plane);
    if (!hit)
    {
        return;
    }

    // Push the sphere out along the normal so it no longer overlaps. The plane
    // is immovable, so the whole correction goes to the body. The Contact3
    // normal points from the sphere INTO the plane, so moving out is -normal.
    if (hit->penetration > math::Real{0})
    {
        body.position -= hit->normal * hit->penetration;
    }

    // Normal impulse, in the same convention as the 2D apply_contact: the
    // contact normal points from the sphere (a) into the plane (b), and the
    // relative velocity is the plane's minus the sphere's, v_b - v_a. The plane
    // is still, so that is -body.velocity, and the closing speed along the
    // normal is -dot(velocity, normal): negative when the sphere is moving into
    // the plane. If it is already separating, do nothing.
    //
    // A centred sphere's spin does not enter this: the spin contribution
    // omega x r is perpendicular to r, and the normal is parallel to r, so it
    // has no component along the normal.
    const math::Real closing = -math::dot(body.velocity, hit->normal);
    if (closing >= math::Real{0})
    {
        return;
    }
    // Effective mass is 1/m (the plane contributes nothing), so the applied
    // velocity change along the normal is exactly -(1 + e) * closing, which
    // reverses the closing speed to -e times itself. Written through the impulse
    // and the inverse mass, as the 2D code is, so it reads as the impulse it is
    // and generalizes to a second movable body later.
    const math::Real inverse_mass = math::Real{1} / body.mass;
    const math::Real magnitude =
        -(math::Real{1} + settings.restitution) * closing / inverse_mass;
    const math::Vec3 impulse = hit->normal * magnitude;
    body.velocity -= impulse * inverse_mass;
}
} // namespace

bool RigidBody3D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           inertia.x > math::Real{0} && inertia.y > math::Real{0} &&
           inertia.z > math::Real{0} && math::is_finite(inertia) &&
           radius >= math::Real{0} && math::is_finite(radius) &&
           math::is_unit(orientation) && math::is_squarable(position) &&
           math::is_squarable(velocity) && math::is_squarable(angular_velocity);
}

bool Rigid3DSettings::is_valid() const
{
    if (!math::is_squarable(torque) || !math::is_squarable(gravity))
    {
        return false;
    }
    if (!(restitution >= math::Real{0}) || !(restitution <= math::Real{1}) ||
        !math::is_finite(restitution))
    {
        return false;
    }
    for (const collide::Plane3& plane : ground)
    {
        if (!plane.is_valid())
        {
            return false;
        }
    }
    return true;
}

math::Vec3 spin_angular_momentum(const RigidBody3D& body)
{
    // I omega in the body frame, where the inertia is diagonal, then rotated
    // into the world frame. The body-frame vector tumbles; this one does not.
    const math::Vec3 body_frame{body.inertia.x * body.angular_velocity.x,
                                body.inertia.y * body.angular_velocity.y,
                                body.inertia.z * body.angular_velocity.z};
    return math::rotate(body.orientation, body_frame);
}

math::Real rotational_energy(const RigidBody3D& body)
{
    const math::Vec3& w = body.angular_velocity;
    // (1/2) omega . (I omega), in the body frame. A scalar, so it does not
    // matter which frame it is computed in, and the body frame is where the
    // inertia is diagonal.
    return math::Real{0.5} * (body.inertia.x * w.x * w.x + body.inertia.y * w.y * w.y +
                              body.inertia.z * w.z * w.z);
}

Rigid3DWorld::Rigid3DWorld(sim_core::SimulationSettings simulation_settings,
                           std::vector<RigidBody3D> bodies,
                           Rigid3DSettings settings)
    : simulation_settings_{simulation_settings},
      settings_{settings},
      bodies_{std::move(bodies)},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus Rigid3DWorld::validate() const
{
    if (!simulation_settings_.is_valid() || // dt > 0, finite
        !settings_.is_valid())              // torque squarable
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const RigidBody3D& body : bodies_)
    {
        if (!body.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

sim_core::StepResult Rigid3DWorld::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    previous_ = bodies_;

    const math::Real dt = simulation_settings_.dt;

    for (RigidBody3D& body : bodies_)
    {
        // (0) gravity, before the position update, which keeps this
        //     semi-implicit Euler. An acceleration, so it does not scale with
        //     mass, and every RigidBody3D is movable (the ground is geometry,
        //     not a body). Zero by default, so a pre-M22 world is unchanged.
        body.velocity += settings_.gravity * dt;

        // (1) Euler's equations, in the body frame, stepped with FORWARD
        //     Euler: the rate is evaluated wholly at the old angular velocity
        //     and orientation, so there is no older state to evaluate it at.
        //     That is where the exact second-order drift in |L|^2 and in the
        //     energy comes from (see rigid3d_world.hpp); the semi-implicit part
        //     of this world is the COUPLING in step (2), not this line.
        //
        //     The applied torque is a WORLD-frame setting, and Euler's
        //     equations are diagonal only in the body frame, so it is rotated
        //     in by the body's CURRENT orientation, before step (2) changes it.
        //     T_body = R^-1 T_world is `rotate` by the conjugate quaternion.
        //
        //     The gyroscopic term is a DIFFERENCE of principal moments, so a
        //     sphere cannot precess and a body with three distinct moments can
        //     tumble; it is identically zero in two dimensions, which is why
        //     none of that physics exists there. The torque is what M21 adds.
        const math::Vec3& w = body.angular_velocity;
        const math::Vec3& I = body.inertia;
        const math::Vec3 torque_body =
            math::rotate(math::conjugate(body.orientation), settings_.torque);
        const math::Vec3 rate{((I.y - I.z) * w.y * w.z + torque_body.x) / I.x,
                              ((I.z - I.x) * w.z * w.x + torque_body.y) / I.y,
                              ((I.x - I.y) * w.x * w.y + torque_body.z) / I.z};
        body.angular_velocity += rate * dt;

        // (2) the orientation, from the UPDATED angular velocity, which keeps
        //     this semi-implicit Euler. Body-frame omega puts the quaternion on
        //     the LEFT of the product.
        const math::Quat spin{math::Real{0}, body.angular_velocity};
        body.orientation =
            body.orientation + (body.orientation * spin) * (dt * math::Real{0.5});

        // (3) renormalize. The step above moves tangentially off the unit
        //     sphere, so the norm grows by exactly
        //     sqrt(1 + |omega|^2 dt^2 / 4) per step: the quaternion norm is
        //     multiplicative and this is a right multiplication.
        //
        //     This is NOT the positional correction's kind of projection and
        //     costs nothing physically. A nonzero quaternion and any scaling of
        //     it represent the same rotation, since q v q^-1 is invariant under
        //     q -> k q, so renormalizing moves no body and perturbs no
        //     conserved quantity. It is required because `rotate` uses the
        //     optimized unit-only form, which is NOT scale invariant: fed a
        //     quaternion of norm k it returns p + k^2 (R p - p), scaling the
        //     displacement rather than performing the rotation.
        body.orientation = math::normalize(body.orientation);

        // (4) translation, with the gravity-updated velocity.
        body.position += body.velocity * dt;
    }

    // (5) resolve contacts: each body against every ground plane, in a fixed
    //     body-then-plane order so a run repeats (docs/04). The plane is
    //     immovable geometry, so only the sphere moves.
    for (RigidBody3D& body : bodies_)
    {
        for (const collide::Plane3& plane : settings_.ground)
        {
            resolve_ground3d(body, plane, settings_);
        }
    }

    for (const RigidBody3D& body : bodies_)
    {
        if (!body.is_valid())
        {
            bodies_.swap(previous_);
            return sim_core::StepResult{sim_core::StepStatus::InvalidState};
        }
    }

    step_->advance();
    return sim_core::StepResult{sim_core::StepStatus::Ok};
}

math::Vec3 total_linear_momentum3d(const std::vector<RigidBody3D>& bodies)
{
    math::Vec3 momentum{};
    for (const RigidBody3D& body : bodies)
    {
        momentum += body.velocity * body.mass;
    }
    return momentum;
}

math::Vec3 total_angular_momentum3d(const std::vector<RigidBody3D>& bodies)
{
    math::Vec3 angular{};
    for (const RigidBody3D& body : bodies)
    {
        // Spin about the centre of mass, plus the orbital term about the world
        // origin. Dropping either makes a whole class of motion look conserved
        // when it is not, which is why the 2D version sums both as well.
        //
        // No infinity guard, unlike the 2D version: there an infinite mass is a
        // legitimate state meaning an immovable body, while a RigidBody3D with
        // an infinite mass or moment is not valid, so a world that steps cannot
        // hold one. `position` is the centre of mass; there is no local offset
        // in three dimensions.
        angular += spin_angular_momentum(body);
        angular += body.mass * math::cross(body.position, body.velocity);
    }
    return angular;
}

math::Real total_kinetic_energy3d(const std::vector<RigidBody3D>& bodies)
{
    math::Real kinetic = 0.0;
    for (const RigidBody3D& body : bodies)
    {
        kinetic += math::Real{0.5} * body.mass * math::dot(body.velocity, body.velocity);
        kinetic += rotational_energy(body);
    }
    return kinetic;
}

math::Real total_potential_energy3d(const std::vector<RigidBody3D>& bodies,
                                    const math::Vec3& gravity)
{
    // Sum of -m (g . r), r the centre of mass. Zero when gravity is zero.
    math::Real potential = 0.0;
    for (const RigidBody3D& body : bodies)
    {
        potential -= body.mass * math::dot(gravity, body.position);
    }
    return potential;
}
} // namespace malloy::rigid
