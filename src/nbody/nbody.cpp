#include <malloy/nbody/nbody.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::nbody
{
bool Body2D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           math::is_finite(position) && math::is_finite(velocity);
}

bool NBodySettings::is_valid() const
{
    return g >= math::Real{0} && softening >= math::Real{0} &&
           math::is_finite(g) && math::is_finite(softening);
}

NBodyWorld::NBodyWorld(sim_core::SimulationSettings simulation_settings,
                       NBodySettings nbody_settings,
                       std::vector<Body2D> bodies)
    : simulation_settings_{simulation_settings},
      nbody_settings_{nbody_settings},
      bodies_{std::move(bodies)}
{
}

sim_core::StepStatus NBodyWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    if (!nbody_settings_.is_valid()) // G >= 0, softening >= 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const Body2D& body : bodies_) // mass > 0, finite position/velocity
    {
        if (!body.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

std::vector<math::Vec2> NBodyWorld::compute_accelerations() const
{
    const std::size_t count = bodies_.size();
    std::vector<math::Vec2> accelerations(count); // (1) cleared to zero

    const math::Real g = nbody_settings_.g;
    const math::Real softening_squared =
        nbody_settings_.softening * nbody_settings_.softening;

    // (2) accumulate softened pairwise gravity (docs/04 formula).
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = 0; j < count; ++j)
        {
            if (i == j)
            {
                continue;
            }
            const math::Vec2 delta = bodies_[j].position - bodies_[i].position;
            const math::Real r2 = math::dot(delta, delta) + softening_squared;
            const math::Real inv_r = math::Real{1} / std::sqrt(r2);
            const math::Real inv_r3 = inv_r * inv_r * inv_r;
            accelerations[i] += g * bodies_[j].mass * delta * inv_r3;
        }
    }
    return accelerations;
}

sim_core::StepResult NBodyWorld::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    const std::vector<math::Vec2> accelerations = compute_accelerations();
    const math::Real dt = simulation_settings_.dt;
    const std::size_t count = bodies_.size();

    // (3) update every velocity using acceleration * dt.
    for (std::size_t i = 0; i < count; ++i)
    {
        bodies_[i].velocity += accelerations[i] * dt;
    }

    // (4) update every position using the now-updated velocity * dt.
    //     Doing velocities first makes this semi-implicit (symplectic) Euler.
    for (std::size_t i = 0; i < count; ++i)
    {
        bodies_[i].position += bodies_[i].velocity * dt;
    }

    ++tick_count_; // (5)

    return sim_core::StepResult{sim_core::StepStatus::Ok}; // (6)
}
} // namespace malloy::nbody
