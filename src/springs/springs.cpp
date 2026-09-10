#include <malloy/springs/springs.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::springs
{
bool Spring::is_valid() const
{
    return a != b && rest_length >= math::Real{0} && stiffness >= math::Real{0} &&
           damping >= math::Real{0} && math::is_finite(rest_length) &&
           math::is_finite(stiffness) && math::is_finite(damping);
}

bool SpringBody2D::is_valid() const
{
    return mass > math::Real{0} && math::is_finite(mass) &&
           math::is_finite(position) && math::is_finite(velocity);
}

SpringId SpringNetwork::add(const Spring& spring)
{
    springs_.push_back(spring);
    return springs_.size() - 1;
}

bool SpringNetwork::remove(SpringId id)
{
    if (id >= springs_.size())
    {
        return false;
    }
    // erase rather than swap-and-pop: the relative order of the remaining
    // springs is part of the contract, because it decides accumulation order.
    springs_.erase(springs_.begin() + static_cast<std::ptrdiff_t>(id));
    return true;
}

bool SpringNetwork::is_valid_for(std::size_t body_count) const
{
    for (const Spring& s : springs_)
    {
        if (!s.is_valid() || s.a >= body_count || s.b >= body_count)
        {
            return false;
        }
    }
    return true;
}

void accumulate_spring_forces(const SpringNetwork& network,
                              const std::vector<SpringBody2D>& bodies,
                              std::vector<math::Vec2>& forces)
{
    const std::size_t count = bodies.size();
    if (forces.size() != count)
    {
        return; // caller error: the buffer must already match the bodies
    }

    // Network order, deliberately. See the header: the traversal order is
    // observable because floating-point addition is not associative.
    for (const Spring& spring : network.springs())
    {
        if (!spring.is_valid() || spring.a >= count || spring.b >= count)
        {
            continue;
        }

        const SpringBody2D& a = bodies[spring.a];
        const SpringBody2D& b = bodies[spring.b];

        const math::Vec2 delta = b.position - a.position;
        const math::Real length_squared = math::dot(delta, delta);
        if (!(length_squared > math::Real{0}) || !math::is_finite(length_squared))
        {
            // Coincident endpoints have no axis, so there is no direction in
            // which this spring could act. Contribute nothing rather than
            // inventing one.
            continue;
        }

        const math::Real length = std::sqrt(length_squared);
        const math::Vec2 axis = delta / length;        // a toward b
        const math::Real extension = length - spring.rest_length;
        // Axial relative velocity, not the full relative speed: a damper
        // resists motion along its own axis, and the sign is what makes it
        // oppose rather than drive the motion.
        const math::Real axial_speed = math::dot(b.velocity - a.velocity, axis);

        const math::Vec2 force =
            axis * (spring.stiffness * extension + spring.damping * axial_speed);

        // Equal and opposite by construction, so a spring network cannot change
        // total momentum. ADDED, so a body shared by several springs collects
        // every contribution.
        forces[spring.a] += force;
        forces[spring.b] -= force;
    }
}

SpringWorld::SpringWorld(sim_core::SimulationSettings simulation_settings,
                         SpringNetwork network, std::vector<SpringBody2D> bodies)
    : simulation_settings_{simulation_settings},
      network_{std::move(network)},
      bodies_{std::move(bodies)},
      forces_(bodies_.size()),
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus SpringWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const SpringBody2D& body : bodies_)
    {
        if (!body.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    // A spring pointing at a body that does not exist is a broken world, not
    // something to silently skip at every step.
    if (!network_.is_valid_for(bodies_.size()))
    {
        return sim_core::StepStatus::InvalidState;
    }

    const math::Real dt = simulation_settings_.dt;
    for (const Spring& spring : network_.springs())
    {
        const SpringBody2D& a = bodies_[spring.a];
        const SpringBody2D& b = bodies_[spring.b];

        // (1) The separation must be squarable.
        //
        // The force kernel computes dot(delta, delta), which overflows to
        // infinity once the separation passes sqrt(DBL_MAX), about 1.34e154,
        // even though every position is still perfectly finite. The kernel's
        // guard then treats that exactly like coincident endpoints and skips
        // the spring, so the spring silently ceases to exist and the bodies
        // coast apart forever with every step still reporting Ok.
        //
        // Reported here instead. This is checked every step rather than only at
        // load, because a separation that starts representable can grow into
        // this range during a run, which is precisely what a diverging network
        // does.
        const math::Vec2 delta = b.position - a.position;
        if (!math::is_finite(math::dot(delta, delta)))
        {
            return sim_core::StepStatus::InvalidState;
        }

        // (2) The timestep must be inside this spring's stability limit.
        //
        // Symplectic Euler on a spring pair is stable only for a bounded
        // timestep, and past that bound the amplitude grows geometrically every
        // step. Nothing recovers from it, and until the separation overflows
        // (see above) nothing reports it either.
        //
        // With reduced mass mu, omega^2 = k/mu and gamma = c/(2 mu), the bound
        // is dt < 2 (sqrt(gamma^2 + omega^2) - gamma) / omega^2, which reduces
        // to the familiar dt < 2/omega when the damping is zero.
        //
        // This is a NECESSARY condition, not a sufficient one. A network's
        // stiffness matrix is a sum of positive-semidefinite per-spring terms,
        // so its largest eigenvalue is at least that of any single spring, and
        // a spring that is unstable alone is unstable in company. The converse
        // does not hold: a network of individually safe springs can still be
        // too stiff as a whole.
        const math::Real reduced = (a.mass * b.mass) / (a.mass + b.mass);
        if (!(reduced > math::Real{0}) || !math::is_finite(reduced))
        {
            return sim_core::StepStatus::InvalidState;
        }
        const math::Real omega_squared = spring.stiffness / reduced;
        if (omega_squared > math::Real{0})
        {
            const math::Real gamma = spring.damping / (math::Real{2} * reduced);
            const math::Real limit =
                math::Real{2} *
                (std::sqrt(gamma * gamma + omega_squared) - gamma) / omega_squared;
            if (!(dt < limit))
            {
                return sim_core::StepStatus::InvalidSettings;
            }
        }
    }

    return sim_core::StepStatus::Ok;
}

sim_core::StepResult SpringWorld::step()
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
    forces_.assign(count, math::Vec2{});

    // (2) and (3) evaluate springs in network order, accumulating equal and
    //     opposite endpoint forces.
    accumulate_spring_forces(network_, bodies_, forces_);

    // (4) integrate in stable index order. Semi-implicit (symplectic) Euler:
    //     acceleration, then velocity, then position from the UPDATED velocity,
    //     matching every other world here.
    //
    //     These few lines duplicate ParticleWorld deliberately (ADR 0008). If a
    //     third domain needs the same path, extract a shared integrator.
    for (std::size_t i = 0; i < count; ++i)
    {
        const math::Vec2 acceleration = forces_[i] / bodies_[i].mass;
        bodies_[i].velocity += acceleration * dt;
    }
    for (std::size_t i = 0; i < count; ++i)
    {
        bodies_[i].position += bodies_[i].velocity * dt;
    }

    for (const SpringBody2D& body : bodies_)
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

math::Vec2 total_momentum(const std::vector<SpringBody2D>& bodies)
{
    math::Vec2 momentum{};
    for (const SpringBody2D& body : bodies)
    {
        momentum += body.mass * body.velocity;
    }
    return momentum;
}

math::Real total_kinetic_energy(const std::vector<SpringBody2D>& bodies)
{
    math::Real kinetic = 0.0;
    for (const SpringBody2D& body : bodies)
    {
        kinetic += math::Real{0.5} * body.mass * math::dot(body.velocity, body.velocity);
    }
    return kinetic;
}

math::Real total_elastic_energy(const SpringNetwork& network,
                                const std::vector<SpringBody2D>& bodies)
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
