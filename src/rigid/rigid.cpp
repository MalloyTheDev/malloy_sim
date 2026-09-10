#include <malloy/rigid/rigid.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
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

// Resolve one rigid contact: separate the overlap, then apply an impulse along
// the contact normal that also spins both bodies.
//
// This is what distinguishes a rigid contact from a particle contact. The
// impulse acts at a point, not at the centre of mass, so each body gets both a
// linear change J/m and an angular change (r x J)/I, with r running from that
// body's CENTRE OF MASS to the contact point.
//
// The effective mass therefore carries two rotational terms that a particle
// contact does not have:
//
//   k = 1/m_a + 1/m_b + (r_a x n)^2 / I_a + (r_b x n)^2 / I_b
//   j = -(1 + restitution) * v_n / k
//
// where v_n is the relative velocity AT THE CONTACT POINT along the normal,
// which itself includes the omega x r terms. Dropping either rotational term
// makes the pair exchange too much momentum and gain energy.
//
// Immovable bodies fall out for free: their inverse mass and inverse inertia
// are zero, so they contribute nothing to k and receive nothing from j.
// Velocity of the material point at a given ARM from the centre of mass.
//
// velocity_at() takes a world point and subtracts the body's CURRENT centre of
// mass, so calling it before and after the positional correction silently
// measures two different arms. Inside a contact the arm is fixed once, up
// front, and every use goes through here.
math::Vec2 velocity_at_arm(const RigidBody2D& body, const math::Vec2& arm)
{
    const math::Real w = body.angular_velocity;
    return body.velocity + math::Vec2{-w * arm.y, w * arm.x};
}

// Applies one already-resolved contact. Split out of resolve_contact because
// it has a second caller: a body against a ground plane, where the plane stands
// in as a body of infinite mass and inertia. Writing the impulse formula twice
// is exactly the second copy that becomes a third (ADR 0008).
void apply_contact(RigidBody2D& a, RigidBody2D& b, const collide::Contact& contact_data,
                   const RigidSettings& settings)
{
    const math::Real restitution = settings.restitution;
    const std::optional<collide::Contact> hit = contact_data;

    const math::Real inverse_mass_a = a.inverse_mass();
    const math::Real inverse_mass_b = b.inverse_mass();
    const math::Real inverse_inertia_a = a.inverse_inertia();
    const math::Real inverse_inertia_b = b.inverse_inertia();

    // Arms from each centre of mass to the contact point, not from the origins.
    //
    // Measured BEFORE the positional correction below, which is the whole point
    // of the ordering. `hit->point` was computed from the pre-correction
    // positions, so pairing it with a post-correction centre of mass mixes two
    // frames: the correction lifts the centre by the full penetration while the
    // point stays where it was, and the lever arm comes out as R + depth/2
    // instead of the R - depth/2 the contact actually describes. That is the
    // one of the three possible answers that errs outward, and it was an
    // accident of statement order rather than a choice.
    const math::Vec2 arm_a = hit->point - center_of_mass(a);
    const math::Vec2 arm_b = hit->point - center_of_mass(b);

    // Relative velocity AT THE CONTACT POINT, which is what a rigid contact
    // responds to: the centre-of-mass velocities plus each body's omega x r.
    // Same frame, for the same reason.
    const math::Vec2 relative = velocity_at_arm(b, arm_b) - velocity_at_arm(a, arm_a);

    // Separate the overlap in inverse-mass proportion. Positions only, so no
    // momentum changes and an immovable body does not move at all.
    //
    // Still before the impulse, so a pair that is overlapping but already
    // separating is pushed apart even though the early return below skips its
    // impulse.
    const math::Real inverse_sum = inverse_mass_a + inverse_mass_b;
    if (hit->penetration > math::Real{0} && inverse_sum > math::Real{0})
    {
        const math::Vec2 correction = hit->normal * (hit->penetration / inverse_sum);
        a.position -= correction * inverse_mass_a;
        b.position += correction * inverse_mass_b;
    }
    const math::Real closing = math::dot(relative, hit->normal);
    if (closing >= math::Real{0})
    {
        return; // already separating, so a second impulse would add energy
    }

    const math::Real cross_a = cross(arm_a, hit->normal);
    const math::Real cross_b = cross(arm_b, hit->normal);
    const math::Real effective = inverse_mass_a + inverse_mass_b +
                                 cross_a * cross_a * inverse_inertia_a +
                                 cross_b * cross_b * inverse_inertia_b;
    if (!(effective > math::Real{0}))
    {
        return;
    }

    const math::Real magnitude =
        -(math::Real{1} + restitution) * closing / effective;
    const math::Vec2 impulse = hit->normal * magnitude;

    a.velocity -= impulse * inverse_mass_a;
    a.angular_velocity -= cross(arm_a, impulse) * inverse_inertia_a;
    b.velocity += impulse * inverse_mass_b;
    b.angular_velocity += cross(arm_b, impulse) * inverse_inertia_b;

    if (!(settings.friction > math::Real{0}))
    {
        return; // frictionless, which is the pre-M17 behaviour and the default
    }

    // Recomputed AFTER the normal impulse, because friction resists what is
    // left rather than what arrived. The ordering is observable and is fixed in
    // the step contract in rigid_world.hpp.
    const math::Vec2 remaining = velocity_at_arm(b, arm_b) - velocity_at_arm(a, arm_a);
    const math::Vec2 sideways =
        remaining - hit->normal * math::dot(remaining, hit->normal);
    const math::Real sliding = math::length(sideways);
    if (!(sliding > math::Real{0}))
    {
        return; // no sliding, so no direction for friction to act along
    }

    // A unit tangent along the slide. There is no fallback for a zero tangent
    // and there must not be: inventing a direction would push a body that is
    // not sliding, which is the classic way a stack of bodies drifts.
    const math::Vec2 tangent = sideways / sliding;

    const math::Real tangent_cross_a = cross(arm_a, tangent);
    const math::Real tangent_cross_b = cross(arm_b, tangent);
    const math::Real tangent_effective =
        inverse_mass_a + inverse_mass_b +
        tangent_cross_a * tangent_cross_a * inverse_inertia_a +
        tangent_cross_b * tangent_cross_b * inverse_inertia_b;
    if (!(tangent_effective > math::Real{0}))
    {
        return;
    }

    // What it would take to stop the slide outright, then Coulomb's cone. The
    // smaller of the two in magnitude is what actually gets applied, so a
    // contact that CAN stop sliding stops exactly rather than reversing: taking
    // the cone limit unconditionally is how friction ends up adding energy.
    math::Real tangent_magnitude = -sliding / tangent_effective;
    const math::Real limit = settings.friction * magnitude;
    if (tangent_magnitude < -limit)
    {
        tangent_magnitude = -limit;
    }

    const math::Vec2 friction_impulse = tangent * tangent_magnitude;

    a.velocity -= friction_impulse * inverse_mass_a;
    a.angular_velocity -= cross(arm_a, friction_impulse) * inverse_inertia_a;
    b.velocity += friction_impulse * inverse_mass_b;
    b.angular_velocity += cross(arm_b, friction_impulse) * inverse_inertia_b;
}

void resolve_contact(RigidBody2D& a, RigidBody2D& b, const RigidSettings& settings)
{
    if (a.radius <= math::Real{0} || b.radius <= math::Real{0})
    {
        return; // a zero radius means the body does not take part in contacts
    }
    if (a.is_static() && b.is_static())
    {
        return; // nothing to do, and the effective mass would be zero
    }

    // Discs are centred on the body ORIGIN, so a contact normal generally does
    // not pass through the centre of mass. That offset is what creates torque.
    const std::optional<collide::Contact> hit =
        collide::contact(collide::Circle{a.position, a.radius},
                         collide::Circle{b.position, b.radius});
    if (!hit)
    {
        return;
    }
    apply_contact(a, b, *hit, settings);
}

// Resolves one body against one immovable ground plane.
//
// The plane stands in as a body of infinite mass and inertia, which costs
// nothing and keeps a single copy of the impulse formula. Its inverse mass and
// inverse inertia are exactly zero, so it absorbs the impulse without moving,
// takes no share of the positional correction, and its arm cannot matter:
// wherever the stand-in nominally sits, its contribution to the effective mass
// is multiplied by zero. No static early-out is needed either, because a body
// that cannot be moved leaves the effective mass at zero and the guard inside
// apply_contact returns.
void resolve_ground(RigidBody2D& body, const collide::Halfplane& plane,
                    const RigidSettings& settings)
{
    if (body.radius <= math::Real{0})
    {
        return; // a zero radius means the body does not take part in contacts
    }

    const std::optional<collide::Contact> hit =
        collide::contact(collide::Circle{body.position, body.radius}, plane);
    if (!hit)
    {
        return;
    }

    RigidBody2D ground;
    ground.mass = std::numeric_limits<math::Real>::infinity();
    ground.inertia = std::numeric_limits<math::Real>::infinity();
    apply_contact(body, ground, *hit, settings);
}

} // namespace

bool RigidBody2D::is_valid() const
{
    // mass and inertia may be infinite, meaning immovable. NaN still fails,
    // because NaN > 0 is false.
    return mass > math::Real{0} && inertia > math::Real{0} &&
           radius >= math::Real{0} && math::is_finite(radius) &&
           math::is_finite(angle) && math::is_finite(angular_velocity) &&
           math::is_finite(position) && math::is_finite(velocity) &&
           math::is_finite(local_center_of_mass);
}

bool RigidBody2D::has_infinite_mass() const
{
    return !math::is_finite(mass);
}

bool RigidBody2D::has_infinite_inertia() const
{
    return !math::is_finite(inertia);
}

bool RigidBody2D::is_static() const
{
    // AND, not OR. See the header: infinity is per quantity, and a body that
    // is immovable in only one sense is still movable in the other.
    return has_infinite_mass() && has_infinite_inertia();
}

math::Real RigidBody2D::inverse_mass() const
{
    return math::is_finite(mass) ? math::Real{1} / mass : math::Real{0};
}

math::Real RigidBody2D::inverse_inertia() const
{
    return math::is_finite(inertia) ? math::Real{1} / inertia : math::Real{0};
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

bool RigidSettings::is_valid() const
{
    if (!(restitution >= math::Real{0}) || !(restitution <= math::Real{1}) ||
        !math::is_finite(restitution) || !math::is_finite(gravity))
    {
        return false;
    }
    if (!(friction >= math::Real{0}) || !math::is_finite(friction))
    {
        return false;
    }
    for (const collide::Halfplane& plane : ground)
    {
        if (!plane.is_valid())
        {
            return false;
        }
    }
    return true;
}

RigidWorld::RigidWorld(sim_core::SimulationSettings simulation_settings,
                       std::vector<RigidBody2D> bodies, RigidSettings settings)
    : simulation_settings_{simulation_settings},
      settings_{settings},
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
    if (!settings_.is_valid())
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
    // Inverse form rather than division, so an immovable body simply does not
    // move instead of needing a special case.
    body.velocity += impulse * body.inverse_mass();
    body.angular_velocity += cross(r, impulse) * body.inverse_inertia();
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

    // (1) gravity, before the position update, which keeps this semi-implicit
    //     (symplectic) Euler and matches every other world here. It is an
    //     acceleration, so it does NOT scale with mass; a body of infinite
    //     mass is skipped explicitly, because otherwise an immovable wall would
    //     start falling. The test is on the MASS, not on is_static(): a body
    //     that merely cannot spin still falls.
    const math::Vec2 gravity = settings_.gravity;
    for (RigidBody2D& body : bodies_)
    {
        if (!body.has_infinite_mass())
        {
            body.velocity += gravity * dt;
        }
    }

    // (2) integrate the pose with the UPDATED velocities.
    for (RigidBody2D& body : bodies_)
    {
        // The centre of mass is what translates; the body origin follows from
        // it. Advancing the origin directly would make an offset body rotate
        // about the wrong point.
        const math::Vec2 com = center_of_mass(body) + body.velocity * dt;
        body.angle += body.angular_velocity * dt;
        body.position = com - rotate(body.local_center_of_mass, body.angle);
    }

    // (3) resolve contacts in a fixed ascending pair order, so a run repeats.
    //     Floating-point addition is not associative, so the order is part of
    //     the observable behaviour (docs/04).
    const std::size_t count = bodies_.size();
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            resolve_contact(bodies_[i], bodies_[j], settings_);
        }
    }

    // (5) then every body against every ground plane, in ascending body then
    //     plane order. Fixed, so a run repeats (docs/04).
    for (RigidBody2D& body : bodies_)
    {
        for (const collide::Halfplane& plane : settings_.ground)
        {
            resolve_ground(body, plane, settings_);
        }
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
        // Infinite mass with zero velocity is inf * 0 = NaN, and such a body
        // never translates, so it carries no linear momentum to report. This
        // asks about the MASS specifically: a body with infinite inertia but
        // finite mass does translate, and its momentum is real.
        if (body.has_infinite_mass())
        {
            continue;
        }
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
        //
        // The two terms are guarded SEPARATELY, because they depend on
        // different quantities. A body with infinite inertia and finite mass
        // has no spin term (inf * 0 is NaN, and it cannot spin) but a perfectly
        // real orbital one.
        if (!body.has_infinite_inertia())
        {
            angular += body.inertia * body.angular_velocity;
        }
        if (!body.has_infinite_mass())
        {
            angular += body.mass * cross(center_of_mass(body), body.velocity);
        }
    }
    return angular;
}

math::Real total_kinetic_energy(const std::vector<RigidBody2D>& bodies)
{
    math::Real kinetic = 0.0;
    for (const RigidBody2D& body : bodies)
    {
        // Translational and rotational energy are guarded separately, for the
        // same reason as angular momentum: an infinite quantity contributes
        // nothing and would produce NaN, but the other half is still real.
        if (!body.has_infinite_mass())
        {
            kinetic +=
                math::Real{0.5} * body.mass * math::dot(body.velocity, body.velocity);
        }
        if (!body.has_infinite_inertia())
        {
            kinetic += math::Real{0.5} * body.inertia * body.angular_velocity *
                       body.angular_velocity;
        }
    }
    return kinetic;
}

math::Real total_potential_energy(const std::vector<RigidBody2D>& bodies,
                                  const math::Vec2& gravity)
{
    math::Real potential = 0.0;
    for (const RigidBody2D& body : bodies)
    {
        if (body.has_infinite_mass())
        {
            continue; // infinite mass times a position is not a number
        }
        // U = -m (g . r) with r the centre of mass, so with g pointing down a
        // higher body has more.
        potential -= body.mass * math::dot(gravity, center_of_mass(body));
    }
    return potential;
}

math::Real total_energy(const std::vector<RigidBody2D>& bodies,
                        const math::Vec2& gravity)
{
    return total_kinetic_energy(bodies) + total_potential_energy(bodies, gravity);
}
} // namespace malloy::rigid
