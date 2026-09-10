#include <malloy/particles/particles.hpp>

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
// Resolve one particle/particle contact: separate the overlap, then apply an
// equal and opposite impulse along the contact normal.
void resolve_pair(Particle2D& a, Particle2D& b, math::Real restitution)
{
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
}

// Keep one particle inside the box. Walls are immovable, so this is a clamp
// plus a damped reflection rather than an exchange of momentum.
void resolve_walls(Particle2D& p, const collide::Aabb& bounds, math::Real restitution)
{
    if (p.position.x - p.radius < bounds.min.x)
    {
        p.position.x = bounds.min.x + p.radius;
        if (p.velocity.x < math::Real{0})
        {
            p.velocity.x = -p.velocity.x * restitution;
        }
    }
    else if (p.position.x + p.radius > bounds.max.x)
    {
        p.position.x = bounds.max.x - p.radius;
        if (p.velocity.x > math::Real{0})
        {
            p.velocity.x = -p.velocity.x * restitution;
        }
    }

    if (p.position.y - p.radius < bounds.min.y)
    {
        p.position.y = bounds.min.y + p.radius;
        if (p.velocity.y < math::Real{0})
        {
            p.velocity.y = -p.velocity.y * restitution;
        }
    }
    else if (p.position.y + p.radius > bounds.max.y)
    {
        p.position.y = bounds.max.y - p.radius;
        if (p.velocity.y > math::Real{0})
        {
            p.velocity.y = -p.velocity.y * restitution;
        }
    }
}
} // namespace

bool Particle2D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           radius >= math::Real{0} && math::is_finite(radius) &&
           math::is_squarable(position) && math::is_squarable(velocity);
}

bool ParticleSettings::is_valid() const
{
    return restitution >= math::Real{0} && restitution <= math::Real{1} &&
           math::is_finite(restitution) && math::is_finite(gravity) &&
           bounds.is_valid();
}

ParticleWorld::ParticleWorld(sim_core::SimulationSettings simulation_settings,
                             ParticleSettings particle_settings,
                             std::vector<Particle2D> particles)
    : simulation_settings_{simulation_settings},
      particle_settings_{particle_settings},
      particles_{std::move(particles)},
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

    const math::Real dt = simulation_settings_.dt;
    const math::Real restitution = particle_settings_.restitution;
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
            resolve_pair(particles_[i], particles_[j], restitution);
        }
    }

    // (3) then keep everything inside the box.
    for (std::size_t i = 0; i < count; ++i)
    {
        resolve_walls(particles_[i], particle_settings_.bounds, restitution);
    }

    // (4) a step that produced non-finite state must not report success.
    for (const Particle2D& p : particles_)
    {
        if (!p.is_valid())
        {
            particles_.swap(previous_);
            return sim_core::StepResult{sim_core::StepStatus::InvalidState};
        }
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
