#include <malloy/charges/charges.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::charges
{
namespace
{
// Rotate a vector by an angle, counterclockwise positive. Local, like the same
// helper in malloy_rigid: two callers in two domains is not yet evidence that
// malloy::math should own it (ADR 0008's trigger is three).
math::Vec2 rotate(const math::Vec2& v, math::Real angle)
{
    const math::Real c = std::cos(angle);
    const math::Real s = std::sin(angle);
    return math::Vec2{v.x * c - v.y * s, v.x * s + v.y * c};
}
} // namespace

bool ChargedParticle2D::is_valid() const
{
    // The charge is deliberately unconstrained apart from being finite. Zero
    // and negative are both meaningful.
    return mass > math::Real{0} && math::is_finite(mass) &&
           math::is_finite(charge) && math::is_squarable(position) &&
           math::is_squarable(velocity);
}

bool ChargeSettings::is_valid() const
{
    return math::is_finite(k) && math::is_finite(electric) &&
           math::is_finite(magnetic) && softening >= math::Real{0} &&
           math::is_finite(softening) && math::is_finite(softening * softening);
}

ChargeWorld::ChargeWorld(sim_core::SimulationSettings simulation_settings,
                         ChargeSettings charge_settings,
                         std::vector<ChargedParticle2D> particles)
    : simulation_settings_{simulation_settings},
      charge_settings_{charge_settings},
      particles_{std::move(particles)},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus ChargeWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    if (!charge_settings_.is_valid())
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const ChargedParticle2D& particle : particles_)
    {
        if (!particle.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

std::vector<math::Vec2> ChargeWorld::compute_electric_accelerations() const
{
    const std::size_t count = particles_.size();
    std::vector<math::Vec2> accelerations(count); // cleared to zero

    const math::Real softening_squared =
        charge_settings_.softening * charge_settings_.softening;

    // Ascending pair order, so a run repeats (docs/04). Both members of a pair
    // are written in the same visit, which is what makes the forces exactly
    // equal and opposite.
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            const math::Vec2 delta = particles_[j].position - particles_[i].position;
            const math::Real r2 = math::dot(delta, delta) + softening_squared;
            if (!(r2 > math::Real{0}))
            {
                continue; // coincident with no softening: no direction to push
            }
            const math::Real inv_r = math::Real{1} / std::sqrt(r2);
            const math::Real inv_r3 = inv_r * inv_r * inv_r;

            // F on j from i is +k q_i q_j delta / r^3, pointing from i toward j
            // when the product of the charges is positive. That sign is the
            // entire difference from gravity: like charges repel.
            const math::Real scale =
                charge_settings_.k * particles_[i].charge * particles_[j].charge * inv_r3;
            const math::Vec2 force = delta * scale;

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

sim_core::StepResult ChargeWorld::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    previous_ = particles_;

    const math::Real dt = simulation_settings_.dt;
    const std::vector<math::Vec2> accelerations = compute_electric_accelerations();

    for (std::size_t i = 0; i < particles_.size(); ++i)
    {
        ChargedParticle2D& particle = particles_[i];

        // (1) the electric kick, which is an acceleration applied for a time.
        particle.velocity += accelerations[i] * dt;

        // (2) the magnetic turn, which is not. See the header: applying this as
        //     another kick would inflate the speed by sqrt(1 + (q b dt/m)^2)
        //     every step. Rotating is exact, so the speed is untouched.
        //
        //     The sign: with b positive and a positive charge, q(v x B) is
        //     q*b*(v.y, -v.x), which turns a particle moving along +x toward
        //     -y. That is clockwise, hence the minus.
        if (charge_settings_.magnetic != math::Real{0} &&
            particle.charge != math::Real{0})
        {
            const math::Real angle =
                -(particle.charge * charge_settings_.magnetic / particle.mass) * dt;
            particle.velocity = rotate(particle.velocity, angle);
        }
    }

    // (3) positions, from the updated velocities.
    for (ChargedParticle2D& particle : particles_)
    {
        particle.position += particle.velocity * dt;
    }

    for (const ChargedParticle2D& particle : particles_)
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

math::Vec2 total_momentum(const std::vector<ChargedParticle2D>& particles)
{
    math::Vec2 momentum{};
    for (const ChargedParticle2D& particle : particles)
    {
        momentum += particle.velocity * particle.mass;
    }
    return momentum;
}

math::Real total_kinetic_energy(const std::vector<ChargedParticle2D>& particles)
{
    math::Real kinetic = 0.0;
    for (const ChargedParticle2D& particle : particles)
    {
        kinetic += math::Real{0.5} * particle.mass *
                   math::dot(particle.velocity, particle.velocity);
    }
    return kinetic;
}

math::Real total_potential_energy(const std::vector<ChargedParticle2D>& particles,
                                  const ChargeSettings& settings)
{
    const math::Real softening_squared = settings.softening * settings.softening;

    math::Real potential = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        for (std::size_t j = i + 1; j < particles.size(); ++j)
        {
            const math::Vec2 delta = particles[j].position - particles[i].position;
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
    for (const ChargedParticle2D& particle : particles)
    {
        potential -= particle.charge * math::dot(settings.electric, particle.position);
    }

    return potential;
}

math::Real total_energy(const std::vector<ChargedParticle2D>& particles,
                        const ChargeSettings& settings)
{
    return total_kinetic_energy(particles) + total_potential_energy(particles, settings);
}
} // namespace malloy::charges
