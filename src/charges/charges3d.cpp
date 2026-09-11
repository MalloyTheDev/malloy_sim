#include <malloy/charges/charges.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::charges
{
bool ChargedParticle3D::is_valid() const
{
    // The charge is deliberately unconstrained apart from being finite. Zero
    // and negative are both meaningful.
    return mass > math::Real{0} && math::is_finite(mass) &&
           math::is_finite(charge) && math::is_squarable(position) &&
           math::is_squarable(velocity);
}

bool Charge3DSettings::is_valid() const
{
    return math::is_finite(k) && math::is_squarable(electric) &&
           math::is_squarable(magnetic) && softening >= math::Real{0} &&
           math::is_finite(softening) && math::is_finite(softening * softening);
}

Charge3DWorld::Charge3DWorld(sim_core::SimulationSettings simulation_settings,
                             Charge3DSettings charge_settings,
                             std::vector<ChargedParticle3D> particles)
    : simulation_settings_{simulation_settings},
      charge_settings_{charge_settings},
      particles_{std::move(particles)},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus Charge3DWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    if (!charge_settings_.is_valid())
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const ChargedParticle3D& particle : particles_)
    {
        if (!particle.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

std::vector<math::Vec3> Charge3DWorld::compute_electric_accelerations() const
{
    const std::size_t count = particles_.size();
    std::vector<math::Vec3> accelerations(count); // cleared to zero

    const math::Real softening_squared =
        charge_settings_.softening * charge_settings_.softening;

    // Ascending pair order, so a run repeats (docs/04). Both members of a pair
    // are written in the same visit, which is what makes the forces exactly
    // equal and opposite.
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            const math::Vec3 delta = particles_[j].position - particles_[i].position;
            const math::Real r2 = math::dot(delta, delta) + softening_squared;
            if (!(r2 > math::Real{0}))
            {
                continue; // coincident with no softening: no direction to push
            }
            const math::Real inv_r = math::Real{1} / std::sqrt(r2);
            const math::Real inv_r3 = inv_r * inv_r * inv_r;

            // F on j from i is +k q_i q_j delta / r^3, pointing from i toward j
            // when the product of the charges is positive: like charges repel,
            // the entire difference from gravity.
            const math::Real scale =
                charge_settings_.k * particles_[i].charge * particles_[j].charge * inv_r3;
            const math::Vec3 force = delta * scale;

            accelerations[j] += force / particles_[j].mass;
            accelerations[i] -= force / particles_[i].mass;
        }
    }

    // The uniform field. Divided by the mass, so this is where a charged
    // particle stops behaving like a body in a gravitational field.
    for (std::size_t i = 0; i < count; ++i)
    {
        accelerations[i] +=
            charge_settings_.electric * (particles_[i].charge / particles_[i].mass);
    }

    return accelerations;
}

sim_core::StepResult Charge3DWorld::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    previous_ = particles_;

    const math::Real dt = simulation_settings_.dt;
    const std::vector<math::Vec3> accelerations = compute_electric_accelerations();

    const math::Real field_magnitude = math::length(charge_settings_.magnetic);

    for (std::size_t i = 0; i < particles_.size(); ++i)
    {
        ChargedParticle3D& particle = particles_[i];

        // (1) the electric kick, an acceleration applied for a time.
        particle.velocity += accelerations[i] * dt;

        // (2) the magnetic turn. Not a kick: rotating the velocity about the
        //     field direction is the exact solution of dv/dt = (q/m)(v x B) for
        //     a uniform B, so the speed is untouched and the component along B
        //     is carried, which together make the motion helical. The rotation
        //     vector is omega dt with omega = -(q/m) B, i.e. an angle
        //     -(q |B| / m) dt about the field axis; from_axis_angle normalizes
        //     the axis, so the field itself is passed as the axis.
        if (field_magnitude > math::Real{0} && particle.charge != math::Real{0})
        {
            const math::Real angle =
                -(particle.charge / particle.mass) * field_magnitude * dt;
            const math::Quat turn = math::from_axis_angle(charge_settings_.magnetic, angle);
            particle.velocity = math::rotate(turn, particle.velocity);
        }
    }

    // (3) positions, from the updated velocities.
    for (ChargedParticle3D& particle : particles_)
    {
        particle.position += particle.velocity * dt;
    }

    for (const ChargedParticle3D& particle : particles_)
    {
        if (!particle.is_valid())
        {
            particles_.swap(previous_);
            return sim_core::StepResult{sim_core::StepStatus::InvalidState};
        }
    }

    step_->advance();
    return sim_core::StepResult{sim_core::StepStatus::Ok};
}

math::Vec3 total_momentum3d(const std::vector<ChargedParticle3D>& particles)
{
    math::Vec3 momentum{};
    for (const ChargedParticle3D& particle : particles)
    {
        momentum += particle.velocity * particle.mass;
    }
    return momentum;
}

math::Real total_kinetic_energy3d(const std::vector<ChargedParticle3D>& particles)
{
    math::Real kinetic = 0.0;
    for (const ChargedParticle3D& particle : particles)
    {
        kinetic += math::Real{0.5} * particle.mass *
                   math::dot(particle.velocity, particle.velocity);
    }
    return kinetic;
}

math::Real total_potential_energy3d(const std::vector<ChargedParticle3D>& particles,
                                    const Charge3DSettings& settings)
{
    const math::Real softening_squared = settings.softening * settings.softening;

    math::Real potential = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        for (std::size_t j = i + 1; j < particles.size(); ++j)
        {
            const math::Vec3 delta = particles[j].position - particles[i].position;
            const math::Real r = std::sqrt(math::dot(delta, delta) + softening_squared);
            // The same guard as the accelerations, so a pair contributing no
            // force contributes no energy either (docs/04).
            if (!(r > math::Real{0}))
            {
                continue;
            }
            potential += settings.k * particles[i].charge * particles[j].charge / r;
        }
    }

    // Uniform field: U = -q (E . r), the same shape as gravity's -m (g . r).
    for (const ChargedParticle3D& particle : particles)
    {
        potential -= particle.charge * math::dot(settings.electric, particle.position);
    }

    return potential;
}

math::Real total_energy3d(const std::vector<ChargedParticle3D>& particles,
                          const Charge3DSettings& settings)
{
    return total_kinetic_energy3d(particles) +
           total_potential_energy3d(particles, settings);
}
} // namespace malloy::charges
