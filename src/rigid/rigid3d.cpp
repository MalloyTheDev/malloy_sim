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
// The world-frame inverse inertia applied to a vector: R I^-1 R^T x. The stored
// inertia is the diagonal of principal moments in the body frame, so this
// rotates x into the body frame, divides component-wise, and rotates back. The
// same body-frame bridge Euler's equations use for the torque (M21), now for a
// contact impulse. For an isotropic (spherical) body it reduces to x / I.
math::Vec3 inverse_inertia_world(const RigidBody3D& body, const math::Vec3& x)
{
    const math::Vec3 body_frame = math::rotate(math::conjugate(body.orientation), x);
    const math::Vec3 scaled{body_frame.x / body.inertia.x,
                            body_frame.y / body.inertia.y,
                            body_frame.z / body.inertia.z};
    return math::rotate(body.orientation, scaled);
}

// One side of a contact: a real body with the arm from its centre of mass to
// the contact point, or the immovable plane (body == nullptr), which
// contributes nothing and receives nothing. The plane as a zero-inverse-mass
// participant is the 3D echo of the 2D `resolve_ground` stand-in body: it lets
// one impulse formula serve both a sphere against a plane and a sphere against
// a sphere, so the formula is never written twice (ADR 0008).
struct Participant
{
    RigidBody3D* body;
    math::Vec3 arm;
};

math::Real inverse_mass_of(const Participant& p)
{
    return p.body ? math::Real{1} / p.body->mass : math::Real{0};
}

// Velocity of the participant's material point at the contact: v + omega x r,
// with omega in the world frame. Zero for the immovable plane.
math::Vec3 contact_velocity_of(const Participant& p)
{
    if (!p.body)
    {
        return math::Vec3{};
    }
    const math::Vec3 omega_world =
        math::rotate(p.body->orientation, p.body->angular_velocity);
    return p.body->velocity + math::cross(omega_world, p.arm);
}

// The rotational part of the effective mass along a direction d:
// d . [(I^-1 (r x d)) x r]. Zero for a plane, and zero for any d parallel to
// the arm (so a centred sphere's normal impulse has no rotational term).
math::Real angular_effective_mass(const Participant& p, const math::Vec3& d)
{
    if (!p.body)
    {
        return math::Real{0};
    }
    return math::dot(
        d, math::cross(inverse_inertia_world(*p.body, math::cross(p.arm, d)), p.arm));
}

// Apply an impulse `signed_impulse` to a participant: change its velocity by
// J/m and its angular velocity by I^-1 (r x J). The caller passes +J to one
// side and -J to the other, so the pair exchanges momentum exactly.
void apply_impulse(const Participant& p, const math::Vec3& signed_impulse)
{
    if (!p.body)
    {
        return;
    }
    const math::Real inverse_mass = math::Real{1} / p.body->mass;
    p.body->velocity += signed_impulse * inverse_mass;
    const math::Vec3 angular_world =
        inverse_inertia_world(*p.body, math::cross(p.arm, signed_impulse));
    p.body->angular_velocity +=
        math::rotate(math::conjugate(p.body->orientation), angular_world);
}

// Resolve one contact between two participants: a positional correction, a
// normal impulse with restitution, then a tangential friction impulse clamped
// to the Coulomb cone. `normal` points from a toward b, and the relative
// velocity is taken as b's minus a's, matching the 2D apply_contact.
//
// This is the whole contact response, and it serves both the ground (b is the
// plane) and a sphere pair (both real). A centred sphere's normal impulse has
// no lever arm, so it imparts no spin; friction is tangential, so it does.
void resolve_pair(const Participant& a, const Participant& b,
                  const math::Vec3& normal, math::Real penetration,
                  const Rigid3DSettings& settings)
{
    const math::Real inverse_mass_a = inverse_mass_of(a);
    const math::Real inverse_mass_b = inverse_mass_of(b);
    const math::Real inverse_sum = inverse_mass_a + inverse_mass_b;
    if (!(inverse_sum > math::Real{0}))
    {
        return; // two immovable things never resolve
    }

    // Positional correction, split by inverse mass so the heavier body moves
    // less. The plane, with zero inverse mass, does not move.
    if (penetration > math::Real{0})
    {
        const math::Vec3 correction = normal * (penetration / inverse_sum);
        if (a.body)
        {
            a.body->position -= correction * inverse_mass_a;
        }
        if (b.body)
        {
            b.body->position += correction * inverse_mass_b;
        }
    }

    // Normal impulse. Closing speed is the relative velocity along the normal;
    // negative when the two are approaching. If already separating, do nothing.
    const math::Vec3 relative = contact_velocity_of(b) - contact_velocity_of(a);
    const math::Real closing = math::dot(relative, normal);
    if (closing >= math::Real{0})
    {
        return;
    }
    const math::Real effective = inverse_sum + angular_effective_mass(a, normal) +
                                 angular_effective_mass(b, normal);
    if (!(effective > math::Real{0}))
    {
        return;
    }
    const math::Real normal_magnitude =
        -(math::Real{1} + settings.restitution) * closing / effective;
    const math::Vec3 normal_impulse = normal * normal_magnitude;
    apply_impulse(a, -normal_impulse);
    apply_impulse(b, normal_impulse);

    // Friction: a tangential impulse, from the velocity that REMAINS after the
    // normal impulse (the observable 2D ordering), clamped to the Coulomb cone.
    if (!(settings.friction > math::Real{0}))
    {
        return;
    }
    const math::Vec3 remaining = contact_velocity_of(b) - contact_velocity_of(a);
    const math::Vec3 sideways = remaining - normal * math::dot(remaining, normal);
    const math::Real sliding = math::length(sideways);
    if (!(sliding > math::Real{0}))
    {
        return; // not sliding; invent no tangent (the 2D no-fallback rule)
    }
    const math::Vec3 tangent = sideways / sliding;
    const math::Real tangent_effective = inverse_sum +
                                         angular_effective_mass(a, tangent) +
                                         angular_effective_mass(b, tangent);
    if (!(tangent_effective > math::Real{0}))
    {
        return;
    }
    math::Real tangent_magnitude = sliding / tangent_effective;
    const math::Real limit = settings.friction * normal_magnitude;
    if (tangent_magnitude > limit)
    {
        tangent_magnitude = limit; // Coulomb's cone
    }
    const math::Vec3 friction_impulse = tangent * tangent_magnitude;
    // Friction opposes the slide, which points from a toward b along `tangent`,
    // so a receives +friction_impulse and b receives -friction_impulse (the
    // opposite of the normal impulse's sign).
    apply_impulse(a, friction_impulse);
    apply_impulse(b, -friction_impulse);
}

// A body collides as an oriented box when all three half-extents are positive;
// otherwise it uses the sphere radius. Zero half-extents (the default) is not a
// box, so a pre-M30 body keeps its sphere behaviour.
bool body_is_box(const RigidBody3D& body)
{
    return body.half_extents.x > math::Real{0} &&
           body.half_extents.y > math::Real{0} && body.half_extents.z > math::Real{0};
}

// A body against one immovable ground plane (M22, M23 for a sphere; M30 for a
// box): the plane is the null participant.
void resolve_ground3d(RigidBody3D& body, const collide::Plane3& plane,
                      const Rigid3DSettings& settings)
{
    const Participant ground{nullptr, math::Vec3{}};

    if (body_is_box(body))
    {
        // A box against a plane touches at up to four corners at once. Rather
        // than resolve each corner in sequence, which would inject a spurious
        // torque on a symmetric landing (the corners are handled one at a time,
        // so the body is no longer symmetric by the second), the coplanar
        // manifold is reduced to a SINGLE contact at the CENTROID of the
        // penetrating corners, resolved with the deepest penetration. For a
        // flat face that centroid lies directly below the centre of mass, so
        // the normal impulse passes through it and makes no torque: a flat drop
        // stays flat and its angular momentum is conserved. For a corner or
        // edge landing the centroid is that corner or edge, which is the
        // correct lever arm, so the box tips and tumbles as it should. A full
        // per-corner manifold with an iterative solver, needed to stand a
        // stack, is a later milestone (the same boundary the 2D solver drew).
        const collide::Box3 shape{body.position, body.half_extents, body.orientation};
        const std::vector<collide::Contact3> hits = collide::contacts(shape, plane);
        if (hits.empty())
        {
            return;
        }
        math::Vec3 centroid{};
        math::Real deepest = math::Real{0};
        for (const collide::Contact3& hit : hits)
        {
            centroid += hit.point;
            deepest = std::fmax(deepest, hit.penetration);
        }
        centroid = centroid / static_cast<math::Real>(hits.size());
        const Participant box_body{&body, centroid - body.position};
        resolve_pair(box_body, ground, -plane.normal, deepest, settings);
        return;
    }

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
    const Participant sphere{&body, hit->normal * body.radius};
    resolve_pair(sphere, ground, hit->normal, hit->penetration, settings);
}

// A sphere against another sphere (M24): both participants are real. The normal
// points from a toward b, so each arm is the radius along it toward the
// contact, a's forward and b's back.
void resolve_contact3d(RigidBody3D& a, RigidBody3D& b,
                       const Rigid3DSettings& settings)
{
    // Body against body is spheres only so far. A box collides with ground
    // planes (above) but not yet with another body, which needs a
    // separating-axis test and an edge-edge case; a box body is skipped here.
    if (body_is_box(a) || body_is_box(b))
    {
        return;
    }
    if (!(a.radius > math::Real{0}) || !(b.radius > math::Real{0}))
    {
        return; // either without a shape does not collide
    }
    const std::optional<collide::Contact3> hit = collide::contact(
        collide::Sphere{a.position, a.radius}, collide::Sphere{b.position, b.radius});
    if (!hit)
    {
        return;
    }
    const Participant pa{&a, hit->normal * a.radius};
    const Participant pb{&b, hit->normal * -b.radius};
    resolve_pair(pa, pb, hit->normal, hit->penetration, settings);
}
} // namespace

bool RigidBody3D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           inertia.x > math::Real{0} && inertia.y > math::Real{0} &&
           inertia.z > math::Real{0} && math::is_finite(inertia) &&
           radius >= math::Real{0} && math::is_finite(radius) &&
           half_extents.x >= math::Real{0} && half_extents.y >= math::Real{0} &&
           half_extents.z >= math::Real{0} && math::is_finite(half_extents) &&
           math::is_unit(orientation) && math::is_squarable(position) &&
           math::is_squarable(velocity) && math::is_squarable(angular_velocity);
}

bool Rigid3DSettings::is_valid() const
{
    if (!math::is_squarable(torque) || !math::is_squarable(force) ||
        !math::is_squarable(gravity))
    {
        return false;
    }
    if (!(restitution >= math::Real{0}) || !(restitution <= math::Real{1}) ||
        !math::is_finite(restitution))
    {
        return false;
    }
    if (!(friction >= math::Real{0}) || !math::is_finite(friction))
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
        !settings_.is_valid())              // torque, force, gravity squarable
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

        //     Then the applied FORCE (M28), as the acceleration F/m it
        //     produces. Unlike gravity this scales with 1/mass, which is the
        //     whole difference between a force and an acceleration and is why it
        //     needs the body's mass (computed from geometry, M25 to M27). It
        //     acts through the centre of mass, so it makes no torque; the
        //     rotational half of a wrench is the torque in step (1). Zero by
        //     default, and F = 0 leaves this line adding a zero vector, so a
        //     pre-M28 world is unchanged bit for bit.
        body.velocity += (settings_.force / body.mass) * dt;

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

    // (5) resolve contacts. First every pair of bodies against each other, in a
    //     fixed ascending pair order, then every body against every ground
    //     plane. Fixed order so a run repeats (docs/04): floating-point addition
    //     is not associative, so the order is part of the observable behaviour.
    const std::size_t count = bodies_.size();
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            resolve_contact3d(bodies_[i], bodies_[j], settings_);
        }
    }
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

math::Real total_force_potential3d(const std::vector<RigidBody3D>& bodies,
                                   const math::Vec3& force)
{
    // Sum of -(F . r), r the centre of mass. The applied force is uniform, not
    // mass-scaled like gravity's m g, so there is no mass factor here. Zero when
    // the force is zero.
    math::Real potential = 0.0;
    for (const RigidBody3D& body : bodies)
    {
        potential -= math::dot(force, body.position);
    }
    return potential;
}
} // namespace malloy::rigid
