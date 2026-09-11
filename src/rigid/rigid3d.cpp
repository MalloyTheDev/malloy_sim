#include <malloy/rigid/rigid.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::rigid
{
bool RigidBody3D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           inertia.x > math::Real{0} && inertia.y > math::Real{0} &&
           inertia.z > math::Real{0} && math::is_finite(inertia) &&
           math::is_unit(orientation) && math::is_squarable(position) &&
           math::is_squarable(velocity) && math::is_squarable(angular_velocity);
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
                           std::vector<RigidBody3D> bodies)
    : simulation_settings_{simulation_settings},
      bodies_{std::move(bodies)},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus Rigid3DWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
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
        // (1) Euler's equations, in the body frame, stepped with FORWARD
        //     Euler: the rate is evaluated wholly at the old angular velocity,
        //     because it depends on omega and on nothing else. That is where
        //     the exact second-order drift in |L|^2 and in the energy comes
        //     from (see rigid3d_world.hpp); the semi-implicit part of this
        //     world is the COUPLING in step (2), not this line.
        //
        //     Every right-hand side is a DIFFERENCE of principal moments, so a
        //     sphere cannot precess and a body with three distinct moments can
        //     tumble. This term is identically zero in two dimensions, which is
        //     why none of this physics exists there.
        const math::Vec3& w = body.angular_velocity;
        const math::Vec3& I = body.inertia;
        const math::Vec3 rate{(I.y - I.z) * w.y * w.z / I.x,
                              (I.z - I.x) * w.z * w.x / I.y,
                              (I.x - I.y) * w.x * w.y / I.z};
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

        // (4) translation, on which nothing acts here.
        body.position += body.velocity * dt;
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
} // namespace malloy::rigid
