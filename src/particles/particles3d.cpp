#include <malloy/particles/particles.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include <malloy/collide/collide.hpp>
#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::particles
{
namespace
{
// Reduce a two-component tangential velocity toward zero by at most `max_delta`,
// the Coulomb cap for a wall contact (a wall's tangent plane has two axes in 3D).
// `max_delta` is friction times the normal impulse per unit mass, so friction 0
// leaves the tangent untouched.
void damp_tangential(math::Real& t1, math::Real& t2, math::Real max_delta)
{
    const math::Real speed = std::sqrt(t1 * t1 + t2 * t2);
    if (speed > math::Real{0})
    {
        const math::Real scale = (speed - std::min(speed, max_delta)) / speed;
        t1 *= scale;
        t2 *= scale;
    }
}

// Resolve one particle/particle contact in 3D: separate the overlap, apply an
// equal and opposite impulse along the contact normal, then a friction impulse
// along the tangent clamped to the Coulomb cone. The 3D copy of the 2D resolver,
// using the sphere geometry M22/M24 added.
void resolve_pair(Particle3D& a, Particle3D& b, math::Real restitution,
                  math::Real friction)
{
    // Reject a pair that is clearly apart before building a contact (M38). It
    // is the same squared-distance test collide::contact makes first, on the
    // same arithmetic, so the result is bit-identical; in a gas almost every
    // pair is apart, and skipping the contact keeps the O(N^2) loop cheap.
    // Every pair is still visited, so this is not a broadphase.
    const math::Real ex = b.position.x - a.position.x;
    const math::Real ey = b.position.y - a.position.y;
    const math::Real ez = b.position.z - a.position.z;
    const math::Real reach = a.radius + b.radius;
    if (ex * ex + ey * ey + ez * ez > reach * reach)
    {
        return;
    }
    const std::optional<collide::Contact3> hit =
        collide::contact(collide::Sphere{a.position, a.radius},
                         collide::Sphere{b.position, b.radius});
    if (!hit)
    {
        return;
    }

    const math::Real inverse_a = math::Real{1} / a.mass;
    const math::Real inverse_b = math::Real{1} / b.mass;
    const math::Real inverse_sum = inverse_a + inverse_b;

    // Push apart in inverse-mass proportion; touches only positions, so momentum
    // is untouched.
    if (hit->penetration > math::Real{0})
    {
        const math::Vec3 correction = hit->normal * (hit->penetration / inverse_sum);
        a.position -= correction * inverse_a;
        b.position += correction * inverse_b;
    }

    // Only push back if they are actually closing, or a pair already separated
    // by an earlier contact would gain energy from a second impulse.
    const math::Vec3 relative = b.velocity - a.velocity;
    const math::Real closing = math::dot(relative, hit->normal);
    if (closing >= math::Real{0})
    {
        return;
    }

    // The separating speed becomes restitution times the closing speed; the two
    // velocity changes are equal and opposite in momentum.
    const math::Real impulse = -(math::Real{1} + restitution) * closing / inverse_sum;
    a.velocity -= hit->normal * (impulse * inverse_a);
    b.velocity += hit->normal * (impulse * inverse_b);

    // Coulomb friction: a tangential impulse from the velocity that REMAINS after
    // the normal impulse, clamped to `friction` times the normal impulse. A
    // particle has no orientation, so this only damps the tangential slide
    // between the two and imparts no spin.
    if (!(friction > math::Real{0}))
    {
        return;
    }
    const math::Vec3 remaining = b.velocity - a.velocity;
    const math::Vec3 sideways = remaining - hit->normal * math::dot(remaining, hit->normal);
    const math::Real sliding = math::length(sideways);
    if (!(sliding > math::Real{0}))
    {
        return; // not sliding; invent no tangent (the no-fallback rule)
    }
    const math::Vec3 tangent = sideways / sliding;
    math::Real tangent_impulse = sliding / inverse_sum;
    const math::Real limit = friction * impulse; // `impulse` is the normal magnitude, > 0 here
    if (tangent_impulse > limit)
    {
        tangent_impulse = limit; // Coulomb's cone
    }
    // a is dragged forward along the slide and b back, the opposite of the normal
    // pair's signs, so the tangential relative velocity shrinks.
    a.velocity += tangent * (tangent_impulse * inverse_a);
    b.velocity -= tangent * (tangent_impulse * inverse_b);
}

// Re-emit a particle from a THERMAL wall (M39), the 3D sibling of the 2D
// thermal_bounce: the normal component redrawn as a flux-weighted speed
// pointing `into` the box, then the two tangential components, which go to the
// wall's tangent axes in ascending order (`t1` before `t2`), three draws in all.
// The impulse tally gets m (|v_n,in| + v_n,out) and the heat tally the change in
// kinetic energy, written as total_kinetic_energy3d writes it.
void thermal_bounce(Particle3D& p, math::Real& normal, math::Real& t1, math::Real& t2,
                    math::Real into, math::Real temperature, SplitMix64& random,
                    math::Real& impulse, math::Real& heat)
{
    const math::Real before = math::Real{0.5} * p.mass * math::dot(p.velocity, p.velocity);
    const math::Real speed_in = std::abs(normal);
    const math::Real speed_out = emission_normal_speed(random, temperature, p.mass);
    const math::Vec2 sideways = emission_tangent(random, temperature, p.mass);
    normal = into * speed_out;
    t1 = sideways.x;
    t2 = sideways.y;
    impulse += p.mass * (speed_in + speed_out);
    heat += math::Real{0.5} * p.mass * math::dot(p.velocity, p.velocity) - before;
}

// Keep one particle inside the box. Walls are immovable, so this is a clamp plus
// a damped reflection rather than an exchange of momentum. On a real bounce the
// two tangential axes of the wall are also damped, clamped Coulomb-style to
// `friction` times the normal impulse per unit mass, which is `(1 + restitution)`
// times the incoming normal speed. Friction 0 is exactly the old reflect. Per
// axis, the 3D copy of the 2D resolver with a third axis added, including the
// tally of each bounce's normal impulse against the wall it hit. A wall with a
// temperature above 0 is thermal instead (M39): clamped as before, re-emitted
// by thermal_bounce, with no restitution or friction there; a specular wall's
// statements are unchanged, so a box without a thermal wall never draws.
void resolve_walls(Particle3D& p, const collide::Aabb3& bounds,
                   math::Real restitution, math::Real friction,
                   const WallTemperature3D& walls, SplitMix64& random,
                   WallImpulse3D& impulse, WallHeat3D& heat)
{
    const math::Real kick = p.mass * (math::Real{1} + restitution);
    if (p.position.x - p.radius < bounds.min.x)
    {
        p.position.x = bounds.min.x + p.radius;
        if (p.velocity.x < math::Real{0})
        {
            if (walls.min_x > math::Real{0})
            {
                thermal_bounce(p, p.velocity.x, p.velocity.y, p.velocity.z,
                               math::Real{1}, walls.min_x, random,
                               impulse.min_x, heat.min_x);
            }
            else
            {
                const math::Real normal_speed = -p.velocity.x;
                p.velocity.x = -p.velocity.x * restitution;
                damp_tangential(p.velocity.y, p.velocity.z,
                                friction * (math::Real{1} + restitution) * normal_speed);
                impulse.min_x += kick * normal_speed;
            }
        }
    }
    else if (p.position.x + p.radius > bounds.max.x)
    {
        p.position.x = bounds.max.x - p.radius;
        if (p.velocity.x > math::Real{0})
        {
            if (walls.max_x > math::Real{0})
            {
                thermal_bounce(p, p.velocity.x, p.velocity.y, p.velocity.z,
                               math::Real{-1}, walls.max_x, random,
                               impulse.max_x, heat.max_x);
            }
            else
            {
                const math::Real normal_speed = p.velocity.x;
                p.velocity.x = -p.velocity.x * restitution;
                damp_tangential(p.velocity.y, p.velocity.z,
                                friction * (math::Real{1} + restitution) * normal_speed);
                impulse.max_x += kick * normal_speed;
            }
        }
    }

    if (p.position.y - p.radius < bounds.min.y)
    {
        p.position.y = bounds.min.y + p.radius;
        if (p.velocity.y < math::Real{0})
        {
            if (walls.min_y > math::Real{0})
            {
                thermal_bounce(p, p.velocity.y, p.velocity.x, p.velocity.z,
                               math::Real{1}, walls.min_y, random,
                               impulse.min_y, heat.min_y);
            }
            else
            {
                const math::Real normal_speed = -p.velocity.y;
                p.velocity.y = -p.velocity.y * restitution;
                damp_tangential(p.velocity.x, p.velocity.z,
                                friction * (math::Real{1} + restitution) * normal_speed);
                impulse.min_y += kick * normal_speed;
            }
        }
    }
    else if (p.position.y + p.radius > bounds.max.y)
    {
        p.position.y = bounds.max.y - p.radius;
        if (p.velocity.y > math::Real{0})
        {
            if (walls.max_y > math::Real{0})
            {
                thermal_bounce(p, p.velocity.y, p.velocity.x, p.velocity.z,
                               math::Real{-1}, walls.max_y, random,
                               impulse.max_y, heat.max_y);
            }
            else
            {
                const math::Real normal_speed = p.velocity.y;
                p.velocity.y = -p.velocity.y * restitution;
                damp_tangential(p.velocity.x, p.velocity.z,
                                friction * (math::Real{1} + restitution) * normal_speed);
                impulse.max_y += kick * normal_speed;
            }
        }
    }

    if (p.position.z - p.radius < bounds.min.z)
    {
        p.position.z = bounds.min.z + p.radius;
        if (p.velocity.z < math::Real{0})
        {
            if (walls.min_z > math::Real{0})
            {
                thermal_bounce(p, p.velocity.z, p.velocity.x, p.velocity.y,
                               math::Real{1}, walls.min_z, random,
                               impulse.min_z, heat.min_z);
            }
            else
            {
                const math::Real normal_speed = -p.velocity.z;
                p.velocity.z = -p.velocity.z * restitution;
                damp_tangential(p.velocity.x, p.velocity.y,
                                friction * (math::Real{1} + restitution) * normal_speed);
                impulse.min_z += kick * normal_speed;
            }
        }
    }
    else if (p.position.z + p.radius > bounds.max.z)
    {
        p.position.z = bounds.max.z - p.radius;
        if (p.velocity.z > math::Real{0})
        {
            if (walls.max_z > math::Real{0})
            {
                thermal_bounce(p, p.velocity.z, p.velocity.x, p.velocity.y,
                               math::Real{-1}, walls.max_z, random,
                               impulse.max_z, heat.max_z);
            }
            else
            {
                const math::Real normal_speed = p.velocity.z;
                p.velocity.z = -p.velocity.z * restitution;
                damp_tangential(p.velocity.x, p.velocity.y,
                                friction * (math::Real{1} + restitution) * normal_speed);
                impulse.max_z += kick * normal_speed;
            }
        }
    }
}

// True when every entry of both wall tallies is finite: the 2D world's check
// (issue #50), over six walls. A valid particle can overflow either tally.
bool tallies_finite(const WallImpulse3D& impulse, const WallHeat3D& heat)
{
    return math::is_finite(impulse.min_x) && math::is_finite(impulse.max_x) &&
           math::is_finite(impulse.min_y) && math::is_finite(impulse.max_y) &&
           math::is_finite(impulse.min_z) && math::is_finite(impulse.max_z) &&
           math::is_finite(heat.min_x) && math::is_finite(heat.max_x) &&
           math::is_finite(heat.min_y) && math::is_finite(heat.max_y) &&
           math::is_finite(heat.min_z) && math::is_finite(heat.max_z);
}
} // namespace

bool Particle3D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           radius >= math::Real{0} && math::is_finite(radius) &&
           math::is_squarable(position) && math::is_squarable(velocity);
}

bool WallTemperature3D::is_valid() const
{
    // t >= 0 refuses a NaN as well as a negative; -0.0 is the specular wall.
    const auto usable = [](math::Real t) { return t >= math::Real{0} && math::is_finite(t); };
    return usable(min_x) && usable(max_x) && usable(min_y) && usable(max_y) &&
           usable(min_z) && usable(max_z);
}

bool ParticleSettings3D::is_valid() const
{
    return restitution >= math::Real{0} && restitution <= math::Real{1} &&
           math::is_finite(restitution) && friction >= math::Real{0} &&
           math::is_finite(friction) && math::is_finite(gravity) &&
           bounds.is_valid() && wall_temperature.is_valid();
}

ParticleWorld3D::ParticleWorld3D(sim_core::SimulationSettings simulation_settings,
                                 ParticleSettings3D particle_settings,
                                 std::vector<Particle3D> particles)
    : simulation_settings_{simulation_settings},
      particle_settings_{particle_settings},
      particles_{std::move(particles)},
      random_{particle_settings.seed},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus ParticleWorld3D::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    if (!particle_settings_.is_valid()) // restitution in [0,1], bounds valid
    {
        return sim_core::StepStatus::InvalidSettings;
    }

    const collide::Aabb3& bounds = particle_settings_.bounds;
    const math::Real width = bounds.max.x - bounds.min.x;
    const math::Real height = bounds.max.y - bounds.min.y;
    const math::Real depth = bounds.max.z - bounds.min.z;

    for (const Particle3D& p : particles_)
    {
        if (!p.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
        // A particle wider than the box would be clamped against two opposite
        // walls at once, and the clamps would fight every step.
        if (p.radius * math::Real{2} > width || p.radius * math::Real{2} > height ||
            p.radius * math::Real{2} > depth)
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

sim_core::StepResult ParticleWorld3D::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    previous_ = particles_;
    const WallImpulse3D previous_impulse = wall_impulse_;
    const WallHeat3D previous_heat = wall_heat_;
    const SplitMix64 previous_random = random_;

    const math::Real dt = simulation_settings_.dt;
    const math::Real restitution = particle_settings_.restitution;
    const math::Real friction = particle_settings_.friction;
    const std::size_t count = particles_.size();

    // (1) gravity, then (2) move by the UPDATED velocity: semi-implicit Euler.
    //     With zero gravity step (1) is a no-op.
    const math::Vec3 gravity = particle_settings_.gravity;
    for (std::size_t i = 0; i < count; ++i)
    {
        particles_[i].velocity += gravity * dt;
    }
    for (std::size_t i = 0; i < count; ++i)
    {
        particles_[i].position += particles_[i].velocity * dt;
    }

    // (3) contacts in a fixed ascending pair order, so a run repeats.
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            resolve_pair(particles_[i], particles_[j], restitution, friction);
        }
    }

    // (4) keep everything inside the box.
    for (std::size_t i = 0; i < count; ++i)
    {
        resolve_walls(particles_[i], particle_settings_.bounds, restitution, friction,
                      particle_settings_.wall_temperature, random_, wall_impulse_,
                      wall_heat_);
    }

    // (5) a step that produced non-finite state, in a particle or in a wall
    //     tally (issue #50), must not report success.
    const bool particles_valid =
        std::all_of(particles_.begin(), particles_.end(),
                    [](const Particle3D& p) { return p.is_valid(); });
    if (!particles_valid || !tallies_finite(wall_impulse_, wall_heat_))
    {
        particles_.swap(previous_);
        wall_impulse_ = previous_impulse;
        wall_heat_ = previous_heat;
        random_ = previous_random; // a failed step leaves no trace in the stream
        return sim_core::StepResult{sim_core::StepStatus::InvalidState};
    }

    step_->advance();
    return sim_core::StepResult{sim_core::StepStatus::Ok};
}

math::Vec3 total_momentum3d(const std::vector<Particle3D>& particles)
{
    math::Vec3 momentum{};
    for (const Particle3D& p : particles)
    {
        momentum += p.mass * p.velocity;
    }
    return momentum;
}

math::Real total_kinetic_energy3d(const std::vector<Particle3D>& particles)
{
    math::Real kinetic = 0.0;
    for (const Particle3D& p : particles)
    {
        kinetic += math::Real{0.5} * p.mass * math::dot(p.velocity, p.velocity);
    }
    return kinetic;
}

math::Vec3 total_angular_momentum3d(const std::vector<Particle3D>& particles)
{
    math::Vec3 angular{};
    for (const Particle3D& p : particles)
    {
        angular += p.mass * math::cross(p.position, p.velocity);
    }
    return angular;
}

math::Real total_potential_energy3d(const std::vector<Particle3D>& particles,
                                    const math::Vec3& gravity)
{
    math::Real potential = 0.0;
    for (const Particle3D& p : particles)
    {
        // U = -m (g . r), so with g pointing down a higher particle has more.
        potential -= p.mass * math::dot(gravity, p.position);
    }
    return potential;
}

math::Real total_energy3d(const std::vector<Particle3D>& particles,
                          const math::Vec3& gravity)
{
    return total_kinetic_energy3d(particles) +
           total_potential_energy3d(particles, gravity);
}
} // namespace malloy::particles
