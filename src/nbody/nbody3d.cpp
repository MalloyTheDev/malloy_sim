#include <malloy/nbody/nbody.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::nbody
{
bool Body3D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           math::is_squarable(position) && math::is_squarable(velocity);
}

NBody3DWorld::NBody3DWorld(sim_core::SimulationSettings simulation_settings,
                           NBodySettings nbody_settings, std::vector<Body3D> bodies)
    : simulation_settings_{simulation_settings},
      nbody_settings_{nbody_settings},
      bodies_{std::move(bodies)},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus NBody3DWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    if (!nbody_settings_.is_valid()) // G, softening, and its square
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const Body3D& body : bodies_)
    {
        if (!body.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

std::vector<math::Vec3> NBody3DWorld::compute_accelerations() const
{
    const std::size_t count = bodies_.size();
    std::vector<math::Vec3> accelerations(count); // cleared to zero

    const math::Real g = nbody_settings_.g;
    const math::Real softening_squared =
        nbody_settings_.softening * nbody_settings_.softening;

    // Ascending pair order, so a run repeats (docs/04). Both members of a pair
    // are written in the same visit, which makes the forces exactly equal and
    // opposite and is what keeps momentum and angular momentum exact.
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            const math::Vec3 delta = bodies_[j].position - bodies_[i].position;
            const math::Real r2 = math::dot(delta, delta) + softening_squared;
            if (!(r2 > math::Real{0}))
            {
                continue; // coincident with no softening: no direction to pull
            }
            const math::Real inv_r = math::Real{1} / std::sqrt(r2);
            const math::Real inv_r3 = inv_r * inv_r * inv_r;

            accelerations[i] += delta * (g * bodies_[j].mass * inv_r3);
            accelerations[j] -= delta * (g * bodies_[i].mass * inv_r3);
        }
    }

    return accelerations;
}

sim_core::StepResult NBody3DWorld::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    previous_ = bodies_;

    const math::Real dt = simulation_settings_.dt;
    const std::vector<math::Vec3> accelerations = compute_accelerations();

    // Velocities before positions: semi-implicit (symplectic) Euler.
    for (std::size_t i = 0; i < bodies_.size(); ++i)
    {
        bodies_[i].velocity += accelerations[i] * dt;
    }
    for (Body3D& body : bodies_)
    {
        body.position += body.velocity * dt;
    }

    for (const Body3D& body : bodies_)
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

math::Vec3 total_momentum(const std::vector<Body3D>& bodies)
{
    math::Vec3 momentum{};
    for (const Body3D& body : bodies)
    {
        momentum += body.velocity * body.mass;
    }
    return momentum;
}

math::Real total_kinetic_energy(const std::vector<Body3D>& bodies)
{
    math::Real kinetic = 0.0;
    for (const Body3D& body : bodies)
    {
        kinetic += math::Real{0.5} * body.mass * math::dot(body.velocity, body.velocity);
    }
    return kinetic;
}

math::Real total_potential_energy(const std::vector<Body3D>& bodies, math::Real g,
                                  math::Real softening)
{
    const math::Real softening_squared = softening * softening;

    math::Real potential = 0.0;
    for (std::size_t i = 0; i < bodies.size(); ++i)
    {
        for (std::size_t j = i + 1; j < bodies.size(); ++j)
        {
            const math::Vec3 delta = bodies[j].position - bodies[i].position;
            const math::Real r = std::sqrt(math::dot(delta, delta) + softening_squared);
            // The same guard as compute_accelerations, so a pair contributing
            // no force contributes no potential energy either (docs/04).
            if (!(r > math::Real{0}))
            {
                continue;
            }
            potential -= g * bodies[i].mass * bodies[j].mass / r;
        }
    }
    return potential;
}

math::Real total_energy(const std::vector<Body3D>& bodies, math::Real g,
                        math::Real softening)
{
    return total_kinetic_energy(bodies) + total_potential_energy(bodies, g, softening);
}

math::Vec3 total_angular_momentum(const std::vector<Body3D>& bodies)
{
    math::Vec3 angular{};
    for (const Body3D& body : bodies)
    {
        angular += math::cross(body.position, body.velocity) * body.mass;
    }
    return angular;
}
} // namespace malloy::nbody
