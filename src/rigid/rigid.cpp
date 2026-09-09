#include <malloy/rigid/rigid.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <malloy/collide/shapes.hpp>
#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::rigid
{
namespace
{
// Rotate a vector by an angle. The one place trigonometry enters this module,
// so a sign error here has exactly one home.
math::Vec2 rotate(const math::Vec2& v, math::Real angle)
{
    const math::Real c = std::cos(angle);
    const math::Real s = std::sin(angle);
    return math::Vec2{v.x * c - v.y * s, v.x * s + v.y * c};
}

// The 2D cross product, a scalar: the z component of a x b.
math::Real cross(const math::Vec2& a, const math::Vec2& b)
{
    return a.x * b.y - a.y * b.x;
}
} // namespace

bool RigidBody2D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           inertia > math::Real{0} && math::is_finite(inertia) &&
           math::is_finite(angle) && math::is_finite(angular_velocity) &&
           math::is_finite(position) && math::is_finite(velocity) &&
           math::is_finite(local_center_of_mass);
}

math::Vec2 center_of_mass(const RigidBody2D& body)
{
    return body.position + rotate(body.local_center_of_mass, body.angle);
}

math::Vec2 to_world(const RigidBody2D& body, const math::Vec2& local_point)
{
    return body.position + rotate(local_point, body.angle);
}

math::Vec2 to_local(const RigidBody2D& body, const math::Vec2& world_point)
{
    return rotate(world_point - body.position, -body.angle);
}

math::Vec2 velocity_at(const RigidBody2D& body, const math::Vec2& world_point)
{
    // r from the centre of mass, not from the body origin. Measuring from the
    // origin is the classic mistake, and it is invisible whenever the two
    // coincide (ADR 0007).
    const math::Vec2 r = world_point - center_of_mass(body);
    const math::Real w = body.angular_velocity;
    return body.velocity + math::Vec2{-w * r.y, w * r.x};
}

// ----------------------------------------------------------------------------
// Mass properties
// ----------------------------------------------------------------------------

namespace
{
MassProperties from_area(math::Real area, math::Real second_moment,
                         const math::Vec2& centroid, math::Real density)
{
    if (!(density > math::Real{0}) || !math::is_finite(density) ||
        !(area > math::Real{0}))
    {
        return MassProperties{};
    }
    return MassProperties{density * area, density * second_moment, centroid};
}
} // namespace

MassProperties mass_properties(const collide::Circle& shape, math::Real density)
{
    return from_area(collide::area(shape), collide::second_moment_of_area(shape),
                     collide::centroid(shape), density);
}

MassProperties mass_properties(const collide::Aabb& shape, math::Real density)
{
    return from_area(collide::area(shape), collide::second_moment_of_area(shape),
                     collide::centroid(shape), density);
}

math::Real shift_inertia(math::Real inertia_about_com, math::Real mass,
                         math::Real distance)
{
    // A negative or non-finite distance is a caller error, and squaring it
    // would quietly turn it into a plausible-looking answer.
    if (!(distance >= math::Real{0}) || !math::is_finite(distance) ||
        !math::is_finite(mass) || !math::is_finite(inertia_about_com))
    {
        return inertia_about_com;
    }
    return inertia_about_com + mass * distance * distance;
}

// ----------------------------------------------------------------------------
// World
// ----------------------------------------------------------------------------

RigidWorld::RigidWorld(sim_core::SimulationSettings simulation_settings,
                       std::vector<RigidBody2D> bodies)
    : simulation_settings_{simulation_settings},
      bodies_{std::move(bodies)},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus RigidWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const RigidBody2D& body : bodies_)
    {
        if (!body.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

bool RigidWorld::apply_impulse_at(std::size_t index, const math::Vec2& impulse,
                                  const math::Vec2& world_point)
{
    if (index >= bodies_.size() || !math::is_finite(impulse) ||
        !math::is_finite(world_point))
    {
        return false;
    }
    RigidBody2D& body = bodies_[index];
    if (!body.is_valid())
    {
        return false;
    }

    // r from the centre of mass. An impulse through the centre of mass gives
    // cross(r, J) == 0 and so leaves the spin untouched.
    const math::Vec2 r = world_point - center_of_mass(body);
    body.velocity += impulse / body.mass;
    body.angular_velocity += cross(r, impulse) / body.inertia;
    return true;
}

sim_core::StepResult RigidWorld::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    previous_ = bodies_;

    const math::Real dt = simulation_settings_.dt;

    for (RigidBody2D& body : bodies_)
    {
        // The centre of mass is what translates; the body origin follows from
        // it. Advancing the origin directly would make an offset body rotate
        // about the wrong point.
        const math::Vec2 com = center_of_mass(body) + body.velocity * dt;
        body.angle += body.angular_velocity * dt;
        body.position = com - rotate(body.local_center_of_mass, body.angle);
    }

    for (const RigidBody2D& body : bodies_)
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

// ----------------------------------------------------------------------------
// Diagnostics
// ----------------------------------------------------------------------------

math::Vec2 total_linear_momentum(const std::vector<RigidBody2D>& bodies)
{
    math::Vec2 momentum{};
    for (const RigidBody2D& body : bodies)
    {
        momentum += body.mass * body.velocity;
    }
    return momentum;
}

math::Real total_angular_momentum(const std::vector<RigidBody2D>& bodies)
{
    math::Real angular = 0.0;
    for (const RigidBody2D& body : bodies)
    {
        // Spin about the centre of mass, plus the orbital term about the world
        // origin. Dropping either makes a whole class of motion look conserved
        // when it is not.
        angular += body.inertia * body.angular_velocity;
        angular += body.mass * cross(center_of_mass(body), body.velocity);
    }
    return angular;
}

math::Real total_kinetic_energy(const std::vector<RigidBody2D>& bodies)
{
    math::Real kinetic = 0.0;
    for (const RigidBody2D& body : bodies)
    {
        kinetic += math::Real{0.5} * body.mass * math::dot(body.velocity, body.velocity);
        kinetic += math::Real{0.5} * body.inertia * body.angular_velocity *
                   body.angular_velocity;
    }
    return kinetic;
}
} // namespace malloy::rigid
