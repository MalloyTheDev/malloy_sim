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
// Reduce one tangential velocity component toward zero by at most `max_delta`,
// the Coulomb cap for a wall contact. `max_delta` is friction times the normal
// impulse per unit mass, so friction 0 leaves the component untouched.
void damp_tangential(math::Real& component, math::Real max_delta)
{
    const math::Real speed = std::abs(component);
    if (speed > math::Real{0})
    {
        component -= (component / speed) * std::min(speed, max_delta);
    }
}

// Resolve one particle/particle contact: separate the overlap, apply an equal
// and opposite impulse along the contact normal, then a friction impulse along
// the tangent, clamped to the Coulomb cone.
void resolve_pair(Particle2D& a, Particle2D& b, math::Real restitution,
                  math::Real friction)
{
    // Reject a pair that is clearly apart before building a contact (M38). It
    // is the same squared-distance test collide::contact makes first, on the
    // same arithmetic, so the result is bit-identical; in a gas almost every
    // pair is apart, and skipping the contact keeps the O(N^2) loop cheap.
    // Every pair is still visited, so this is not a broadphase.
    const math::Real ex = b.position.x - a.position.x;
    const math::Real ey = b.position.y - a.position.y;
    const math::Real reach = a.radius + b.radius;
    if (ex * ex + ey * ey > reach * reach)
    {
        return;
    }
    const std::optional<collide::Contact> hit =
        collide::contact(collide::Circle{a.position, a.radius},
                         collide::Circle{b.position, b.radius});
    if (!hit)
    {
        return;
    }

    const math::Real inverse_a = math::Real{1} / a.mass;
    const math::Real inverse_b = math::Real{1} / b.mass;
    const math::Real inverse_sum = inverse_a + inverse_b;

    // Push the pair apart in inverse-mass proportion, so the lighter particle
    // moves further. This touches only positions, so momentum is untouched.
    if (hit->penetration > math::Real{0})
    {
        const math::Vec2 correction = hit->normal * (hit->penetration / inverse_sum);
        a.position -= correction * inverse_a;
        b.position += correction * inverse_b;
    }

    // Only push back if they are actually closing. Without this check a pair
    // that has already been separated by an earlier contact would get a second
    // impulse and gain energy from nothing.
    const math::Vec2 relative = b.velocity - a.velocity;
    const math::Real closing = math::dot(relative, hit->normal);
    if (closing >= math::Real{0})
    {
        return;
    }

    // j is chosen so the separating speed is restitution times the closing
    // speed. The two velocity changes are equal and opposite in momentum, so
    // the pair conserves it at any restitution.
    const math::Real impulse = -(math::Real{1} + restitution) * closing / inverse_sum;
    a.velocity -= hit->normal * (impulse * inverse_a);
    b.velocity += hit->normal * (impulse * inverse_b);

    // Coulomb friction: a tangential impulse from the velocity that REMAINS
    // after the normal impulse, clamped to `friction` times the normal impulse.
    // A particle has no orientation, so this only damps the tangential slide
    // between the two, imparting no spin (no lever arm, unlike the rigid domain).
    if (!(friction > math::Real{0}))
    {
        return;
    }
    const math::Vec2 remaining = b.velocity - a.velocity;
    const math::Vec2 sideways = remaining - hit->normal * math::dot(remaining, hit->normal);
    const math::Real sliding = math::length(sideways);
    if (!(sliding > math::Real{0}))
    {
        return; // not sliding; invent no tangent (the no-fallback rule)
    }
    const math::Vec2 tangent = sideways / sliding;
    math::Real tangent_impulse = sliding / inverse_sum;
    const math::Real limit = friction * impulse; // `impulse` is the normal magnitude, > 0 here
    if (tangent_impulse > limit)
    {
        tangent_impulse = limit; // Coulomb's cone
    }
    // Friction opposes the slide, which runs from a toward b along `tangent`, so
    // a is dragged forward along it and b back, the opposite of the normal pair.
    a.velocity += tangent * (tangent_impulse * inverse_a);
    b.velocity -= tangent * (tangent_impulse * inverse_b);
}

// Re-emit a particle from a THERMAL wall (M39) it is moving into. The wall
// fully accommodates, so the incoming velocity is forgotten: the normal
// component is redrawn as a flux-weighted speed pointing `into` the box (+1
// from a min wall, -1 from a max one), then the tangential component, in that
// order, three draws in all (thermo.hpp). The impulse tally gets the true normal
// momentum change, m (|v_n,in| + v_n,out), and the heat tally the change in
// kinetic energy, written exactly as total_kinetic_energy writes it so the
// energy books balance to rounding. `normal` and `tangent` are the particle's
// own velocity components for this wall.
void thermal_bounce(Particle2D& p, math::Real& normal, math::Real& tangent,
                    math::Real into, math::Real temperature, SplitMix64& random,
                    math::Real& impulse, math::Real& heat)
{
    const math::Real before = math::Real{0.5} * p.mass * math::dot(p.velocity, p.velocity);
    const math::Real speed_in = std::abs(normal);
    const math::Real speed_out = emission_normal_speed(random, temperature, p.mass);
    const math::Vec2 sideways = emission_tangent(random, temperature, p.mass);
    normal = into * speed_out;
    tangent = sideways.x; // 2D has one tangent axis; the second deviate is unused
    impulse += p.mass * (speed_in + speed_out);
    heat += math::Real{0.5} * p.mass * math::dot(p.velocity, p.velocity) - before;
}

// Keep one particle inside the box. Walls are immovable, so this is a clamp
// plus a damped reflection rather than an exchange of momentum. On a real bounce
// the tangential velocity is also damped, clamped Coulomb-style to `friction`
// times the normal impulse per unit mass, which is `(1 + restitution)` times the
// incoming normal speed. Friction 0 leaves the tangent untouched, so the motion
// is exactly the old reflect. Each bounce also adds its normal impulse,
// m (1 + e) times the incoming normal speed, to the tally for the wall it hit.
// A wall with a temperature above 0 is thermal instead (M39): the particle is
// clamped as before but re-emitted by thermal_bounce, and restitution and
// friction do not apply there. A specular wall's statements are unchanged, so a
// box without a thermal wall behaves exactly as before and never draws.
void resolve_walls(Particle2D& p, const collide::Aabb& bounds,
                   math::Real restitution, math::Real friction,
                   const WallTemperature2D& walls, SplitMix64& random,
                   WallImpulse2D& impulse, WallHeat2D& heat)
{
    const math::Real kick = p.mass * (math::Real{1} + restitution);
    if (p.position.x - p.radius < bounds.min.x)
    {
        p.position.x = bounds.min.x + p.radius;
        if (p.velocity.x < math::Real{0})
        {
            if (walls.min_x > math::Real{0})
            {
                thermal_bounce(p, p.velocity.x, p.velocity.y, math::Real{1}, walls.min_x,
                               random, impulse.min_x, heat.min_x);
            }
            else
            {
                const math::Real normal_speed = -p.velocity.x;
                p.velocity.x = -p.velocity.x * restitution;
                damp_tangential(p.velocity.y,
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
                thermal_bounce(p, p.velocity.x, p.velocity.y, math::Real{-1}, walls.max_x,
                               random, impulse.max_x, heat.max_x);
            }
            else
            {
                const math::Real normal_speed = p.velocity.x;
                p.velocity.x = -p.velocity.x * restitution;
                damp_tangential(p.velocity.y,
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
                thermal_bounce(p, p.velocity.y, p.velocity.x, math::Real{1}, walls.min_y,
                               random, impulse.min_y, heat.min_y);
            }
            else
            {
                const math::Real normal_speed = -p.velocity.y;
                p.velocity.y = -p.velocity.y * restitution;
                damp_tangential(p.velocity.x,
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
                thermal_bounce(p, p.velocity.y, p.velocity.x, math::Real{-1}, walls.max_y,
                               random, impulse.max_y, heat.max_y);
            }
            else
            {
                const math::Real normal_speed = p.velocity.y;
                p.velocity.y = -p.velocity.y * restitution;
                damp_tangential(p.velocity.x,
                                friction * (math::Real{1} + restitution) * normal_speed);
                impulse.max_y += kick * normal_speed;
            }
        }
    }
}

// True when every entry of both wall tallies is finite. The tallies are state a
// step rolls back, and a VALID particle can still overflow them (issue #50):
// mass 1e300 at speed 1e9 makes the impulse m (1 + e) |v| infinite, and at a
// thermal wall at speed 1e5 its kinetic energy, and with it the heat, overflows.
bool tallies_finite(const WallImpulse2D& impulse, const WallHeat2D& heat)
{
    return math::is_finite(impulse.min_x) && math::is_finite(impulse.max_x) &&
           math::is_finite(impulse.min_y) && math::is_finite(impulse.max_y) &&
           math::is_finite(heat.min_x) && math::is_finite(heat.max_x) &&
           math::is_finite(heat.min_y) && math::is_finite(heat.max_y);
}
} // namespace

bool Particle2D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           radius >= math::Real{0} && math::is_finite(radius) &&
           math::is_squarable(position) && math::is_squarable(velocity);
}

bool WallTemperature2D::is_valid() const
{
    // t >= 0 refuses a NaN as well as a negative; -0.0 passes and, not being
    // above 0, is the specular wall.
    const auto usable = [](math::Real t) { return t >= math::Real{0} && math::is_finite(t); };
    return usable(min_x) && usable(max_x) && usable(min_y) && usable(max_y);
}

bool ParticleSettings::is_valid() const
{
    return restitution >= math::Real{0} && restitution <= math::Real{1} &&
           math::is_finite(restitution) && friction >= math::Real{0} &&
           math::is_finite(friction) && math::is_finite(gravity) &&
           bounds.is_valid() && wall_temperature.is_valid();
}

ParticleWorld::ParticleWorld(sim_core::SimulationSettings simulation_settings,
                             ParticleSettings particle_settings,
                             std::vector<Particle2D> particles)
    : simulation_settings_{simulation_settings},
      particle_settings_{particle_settings},
      particles_{std::move(particles)},
      random_{particle_settings.seed},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus ParticleWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    if (!particle_settings_.is_valid()) // restitution in [0,1], bounds valid
    {
        return sim_core::StepStatus::InvalidSettings;
    }

    const collide::Aabb& bounds = particle_settings_.bounds;
    const math::Real width = bounds.max.x - bounds.min.x;
    const math::Real height = bounds.max.y - bounds.min.y;

    for (const Particle2D& p : particles_)
    {
        if (!p.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
        // A particle wider than the box would be clamped against both walls at
        // once, and the two clamps would fight every step. Reject it rather
        // than simulate something meaningless.
        if (p.radius * math::Real{2} > width || p.radius * math::Real{2} > height)
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

sim_core::StepResult ParticleWorld::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    previous_ = particles_;
    const WallImpulse2D previous_impulse = wall_impulse_;
    const WallHeat2D previous_heat = wall_heat_;
    const SplitMix64 previous_random = random_;

    const math::Real dt = simulation_settings_.dt;
    const math::Real restitution = particle_settings_.restitution;
    const math::Real friction = particle_settings_.friction;
    const std::size_t count = particles_.size();

    // (1) gravity first, then (2) move by the UPDATED velocity. Velocities
    //     before positions is what makes this semi-implicit (symplectic) Euler,
    //     matching NBodyWorld. With zero gravity step (1) is a no-op and the
    //     result is bit-identical to before gravity existed.
    const math::Vec2 gravity = particle_settings_.gravity;
    for (std::size_t i = 0; i < count; ++i)
    {
        particles_[i].velocity += gravity * dt;
    }
    for (std::size_t i = 0; i < count; ++i)
    {
        particles_[i].position += particles_[i].velocity * dt;
    }

    // (2) resolve contacts in a fixed ascending pair order, so a run repeats.
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            resolve_pair(particles_[i], particles_[j], restitution, friction);
        }
    }

    // (3) then keep everything inside the box.
    for (std::size_t i = 0; i < count; ++i)
    {
        resolve_walls(particles_[i], particle_settings_.bounds, restitution, friction,
                      particle_settings_.wall_temperature, random_, wall_impulse_,
                      wall_heat_);
    }

    // (4) a step that produced non-finite state must not report success. The
    //     wall tallies are part of that state (issue #50).
    const bool particles_valid =
        std::all_of(particles_.begin(), particles_.end(),
                    [](const Particle2D& p) { return p.is_valid(); });
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

math::Vec2 total_momentum(const std::vector<Particle2D>& particles)
{
    math::Vec2 momentum{};
    for (const Particle2D& p : particles)
    {
        momentum += p.mass * p.velocity;
    }
    return momentum;
}

math::Real total_angular_momentum(const std::vector<Particle2D>& particles)
{
    math::Real angular = 0.0;
    for (const Particle2D& particle : particles)
    {
        angular += particle.mass * (particle.position.x * particle.velocity.y -
                                    particle.position.y * particle.velocity.x);
    }
    return angular;
}

math::Real total_kinetic_energy(const std::vector<Particle2D>& particles)
{
    math::Real kinetic = 0.0;
    for (const Particle2D& p : particles)
    {
        kinetic += math::Real{0.5} * p.mass * math::dot(p.velocity, p.velocity);
    }
    return kinetic;
}

math::Real total_potential_energy(const std::vector<Particle2D>& particles,
                                  const math::Vec2& gravity)
{
    math::Real potential = 0.0;
    for (const Particle2D& p : particles)
    {
        // U = -m (g . r), so with g pointing down a higher particle has more.
        potential -= p.mass * math::dot(gravity, p.position);
    }
    return potential;
}

math::Real total_energy(const std::vector<Particle2D>& particles,
                        const math::Vec2& gravity)
{
    return total_kinetic_energy(particles) + total_potential_energy(particles, gravity);
}
} // namespace malloy::particles
