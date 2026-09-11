#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/nbody/body3d.hpp>
#include <malloy/nbody/nbody_settings.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::nbody
{
// Newtonian gravity in three dimensions.
//
// The same physics and the same integrator as NBodyWorld, one component wider.
// Concrete and derived from nothing, like every other world here (ADR 0006),
// and it shares NBodySettings with the 2D world because G and the softening
// length mean exactly the same thing in either dimension.
class NBody3DWorld
{
public:
    NBody3DWorld(sim_core::SimulationSettings simulation_settings,
                 NBodySettings nbody_settings, std::vector<Body3D> bodies);

    // Ok, or InvalidSettings (dt, G, softening) / InvalidState (a body).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step: accelerations from the current
    // positions, then velocities, then positions. Semi-implicit (symplectic)
    // Euler, the same order and the same guarantees as every other world here
    // (docs/04).
    sim_core::StepResult step();

    std::vector<math::Vec3> compute_accelerations() const;

    const std::vector<Body3D>& bodies() const { return bodies_; }
    std::uint64_t tick_count() const { return step_ ? step_->tick_count() : 0; }

    math::Real elapsed_time() const
    {
        return step_ ? step_->elapsed_time() : math::Real{0};
    }

    const sim_core::SimulationSettings& simulation_settings() const
    {
        return simulation_settings_;
    }
    const NBodySettings& nbody_settings() const { return nbody_settings_; }

private:
    sim_core::SimulationSettings simulation_settings_;
    NBodySettings nbody_settings_;
    std::vector<Body3D> bodies_;
    std::vector<Body3D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Diagnostics ---

math::Vec3 total_momentum(const std::vector<Body3D>& bodies);

math::Real total_kinetic_energy(const std::vector<Body3D>& bodies);

// Uses the same softening as the integrator, so total_energy is a quantity the
// dynamics conserve rather than one they approach (docs/04).
math::Real total_potential_energy(const std::vector<Body3D>& bodies, math::Real g,
                                  math::Real softening);

math::Real total_energy(const std::vector<Body3D>& bodies, math::Real g,
                        math::Real softening);

// Total angular momentum about the origin: the sum of m (r x v).
//
// A VECTOR here, where the 2D world returns a scalar, and that is the first
// place three dimensions stop being two with an extra component rather than
// merely being wider. It is still exactly conserved for pairwise central
// forces, to all orders in dt: the (i,j) and (j,i) torque contributions are
// x_i x x_j and x_j x x_i, which cancel identically, and softening changes
// only the magnitude.
math::Vec3 total_angular_momentum(const std::vector<Body3D>& bodies);
} // namespace malloy::nbody
