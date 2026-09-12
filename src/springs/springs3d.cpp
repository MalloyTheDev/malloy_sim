#include <malloy/springs/springs.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::springs
{
bool SpringBody3D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           math::is_squarable(position) && math::is_squarable(velocity);
}

void accumulate_spring_forces(const SpringNetwork& network,
                              const std::vector<SpringBody3D>& bodies,
                              std::vector<math::Vec3>& forces)
{
    const std::size_t count = bodies.size();
    if (forces.size() != count)
    {
        return; // caller error: the buffer must already match the bodies
    }

    // Network order, deliberately. Floating-point addition is not associative,
    // so the traversal order is observable (docs/04).
    for (const Spring& spring : network.springs())
    {
        if (!spring.is_valid() || spring.a >= count || spring.b >= count)
        {
            continue;
        }

        const SpringBody3D& a = bodies[spring.a];
        const SpringBody3D& b = bodies[spring.b];

        const math::Vec3 delta = b.position - a.position;
        const math::Real length_squared = math::dot(delta, delta);
        if (!(length_squared > math::Real{0}) || !math::is_finite(length_squared))
        {
            // Coincident endpoints have no axis: contribute nothing rather than
            // inventing a direction.
            continue;
        }

        const math::Real length = std::sqrt(length_squared);
        const math::Vec3 axis = delta / length; // a toward b
        const math::Real extension = length - spring.rest_length;
        // Axial relative velocity, not the full relative speed: a damper resists
        // motion along its own axis, and the sign is what makes it oppose the
        // motion rather than drive it.
        const math::Real axial_speed = math::dot(b.velocity - a.velocity, axis);

        const math::Vec3 force =
            axis * (spring.stiffness * extension + spring.damping * axial_speed);

        // Equal and opposite by construction, so a network cannot change total
        // momentum. ADDED, so a body shared by several springs collects every
        // contribution.
        forces[spring.a] += force;
        forces[spring.b] -= force;
    }
}

SpringWorld3D::SpringWorld3D(sim_core::SimulationSettings simulation_settings,
                             SpringNetwork network, std::vector<SpringBody3D> bodies)
    : simulation_settings_{simulation_settings},
      network_{std::move(network)},
      bodies_{std::move(bodies)},
      forces_(bodies_.size()),
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus SpringWorld3D::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const SpringBody3D& body : bodies_)
    {
        if (!body.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    if (!network_.is_valid_for(bodies_.size()))
    {
        return sim_core::StepStatus::InvalidState;
    }

    const math::Real dt = simulation_settings_.dt;
    for (const Spring& spring : network_.springs())
    {
        const SpringBody3D& a = bodies_[spring.a];
        const SpringBody3D& b = bodies_[spring.b];

        // (1) The separation must be squarable: dot(delta, delta) overflows past
        //     about 1.34e154 even with finite positions, and the force kernel
        //     then treats that like coincident endpoints and silently drops the
        //     spring while the bodies coast apart. Checked every step, because a
        //     diverging network grows into this range during a run.
        const math::Vec3 delta = b.position - a.position;
        if (!math::is_finite(math::dot(delta, delta)))
        {
            return sim_core::StepStatus::InvalidState;
        }

        // (2) The timestep must be inside this spring's symplectic stability
        //     limit, dt < 2 (sqrt(gamma^2 + omega^2) - gamma) / omega^2 with
        //     omega^2 = k/mu and gamma = c/(2 mu). A NECESSARY condition for the
        //     network, not a sufficient one (see the 2D world's derivation).
        const math::Real reduced = (a.mass * b.mass) / (a.mass + b.mass);
        if (!(reduced > math::Real{0}) || !math::is_finite(reduced))
        {
            return sim_core::StepStatus::InvalidState;
        }
        const math::Real omega_squared = spring.stiffness / reduced;
        if (omega_squared > math::Real{0})
        {
            const math::Real gamma = spring.damping / (math::Real{2} * reduced);
            const math::Real limit = math::Real{2} *
                                     (std::sqrt(gamma * gamma + omega_squared) - gamma) /
                                     omega_squared;
            if (!(dt < limit))
            {
                return sim_core::StepStatus::InvalidSettings;
            }
        }
    }

    return sim_core::StepStatus::Ok;
}

sim_core::StepResult SpringWorld3D::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    previous_ = bodies_;

    const math::Real dt = simulation_settings_.dt;
    const std::size_t count = bodies_.size();

    // (1) clear the accumulator.
    forces_.assign(count, math::Vec3{});

    // (2) and (3) evaluate springs in network order, accumulating equal and
    //     opposite endpoint forces.
    accumulate_spring_forces(network_, bodies_, forces_);

    // (4) integrate in stable index order, semi-implicit (symplectic) Euler:
    //     acceleration, then velocity, then position from the UPDATED velocity.
    //     Duplicates the 2D world deliberately (ADR 0008).
    for (std::size_t i = 0; i < count; ++i)
    {
        const math::Vec3 acceleration = forces_[i] / bodies_[i].mass;
        bodies_[i].velocity += acceleration * dt;
    }
    for (std::size_t i = 0; i < count; ++i)
    {
        bodies_[i].position += bodies_[i].velocity * dt;
    }

    for (const SpringBody3D& body : bodies_)
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

math::Vec3 total_momentum3d(const std::vector<SpringBody3D>& bodies)
{
    math::Vec3 momentum{};
    for (const SpringBody3D& body : bodies)
    {
        momentum += body.mass * body.velocity;
    }
    return momentum;
}

math::Real total_kinetic_energy3d(const std::vector<SpringBody3D>& bodies)
{
    math::Real kinetic = 0.0;
    for (const SpringBody3D& body : bodies)
    {
        kinetic += math::Real{0.5} * body.mass * math::dot(body.velocity, body.velocity);
    }
    return kinetic;
}

math::Real total_elastic_energy3d(const SpringNetwork& network,
                                  const std::vector<SpringBody3D>& bodies)
{
    const std::size_t count = bodies.size();
    math::Real elastic = 0.0;
    for (const Spring& spring : network.springs())
    {
        if (!spring.is_valid() || spring.a >= count || spring.b >= count)
        {
            continue;
        }
        const math::Real length =
            math::distance(bodies[spring.b].position, bodies[spring.a].position);
        const math::Real extension = length - spring.rest_length;
        elastic += math::Real{0.5} * spring.stiffness * extension * extension;
    }
    return elastic;
}
} // namespace malloy::springs
