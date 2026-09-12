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
// Resolve one particle/particle contact in 3D: separate the overlap, then apply
// an equal and opposite impulse along the contact normal. The 3D copy of the 2D
// resolver, using the sphere geometry M22/M24 added.
void resolve_pair(Particle3D& a, Particle3D& b, math::Real restitution)
{
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
}

// Keep one particle inside the box. Walls are immovable, so this is a clamp plus
// a damped reflection rather than an exchange of momentum. Per axis, the 3D copy
// of the 2D resolver with a third axis added.
void resolve_walls(Particle3D& p, const collide::Aabb3& bounds, math::Real restitution)
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

    if (p.position.z - p.radius < bounds.min.z)
    {
        p.position.z = bounds.min.z + p.radius;
        if (p.velocity.z < math::Real{0})
        {
            p.velocity.z = -p.velocity.z * restitution;
        }
    }
    else if (p.position.z + p.radius > bounds.max.z)
    {
        p.position.z = bounds.max.z - p.radius;
        if (p.velocity.z > math::Real{0})
        {
            p.velocity.z = -p.velocity.z * restitution;
        }
    }
}
} // namespace

bool Particle3D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           radius >= math::Real{0} && math::is_finite(radius) &&
           math::is_squarable(position) && math::is_squarable(velocity);
}

bool ParticleSettings3D::is_valid() const
{
    return restitution >= math::Real{0} && restitution <= math::Real{1} &&
           math::is_finite(restitution) && math::is_finite(gravity) &&
           bounds.is_valid();
}

ParticleWorld3D::ParticleWorld3D(sim_core::SimulationSettings simulation_settings,
                                 ParticleSettings3D particle_settings,
                                 std::vector<Particle3D> particles)
    : simulation_settings_{simulation_settings},
      particle_settings_{particle_settings},
      particles_{std::move(particles)},
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

    const math::Real dt = simulation_settings_.dt;
    const math::Real restitution = particle_settings_.restitution;
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
            resolve_pair(particles_[i], particles_[j], restitution);
        }
    }

    // (4) keep everything inside the box.
    for (std::size_t i = 0; i < count; ++i)
    {
        resolve_walls(particles_[i], particle_settings_.bounds, restitution);
    }

    for (const Particle3D& p : particles_)
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
