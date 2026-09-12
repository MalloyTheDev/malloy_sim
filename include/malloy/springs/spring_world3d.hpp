#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec3.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/springs/spring.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::springs
{
// A translational body in three dimensions: the 3D sibling of SpringBody2D.
// Position and velocity gain a third component; mass is unchanged. No radius
// and no orientation, because a spring needs neither.
//
// `Spring` and `SpringNetwork` are NOT duplicated for 3D: they reference bodies
// by id and carry only scalar rest length, stiffness and damping, so they are
// dimension-agnostic and shared with the 2D world unchanged. Only the body and
// the force are dimensional (ADR 0009: the dimension is not a domain).
struct SpringBody3D
{
    math::Vec3 position{};
    math::Vec3 velocity{};
    math::Real mass{1.0};

    // Valid when mass is strictly positive and position and velocity are
    // SQUARABLE (docs/04): the force squares the separation and the drift law
    // squares the velocity.
    bool is_valid() const;
};

// Accumulate the force each spring exerts on its two endpoints, in three
// dimensions. The 3D overload of the 2D kernel, with the same contract: it ADDS
// into `forces` (already sized to match `bodies`), reads body state only, and
// neither integrates nor mutates a body. The force for one spring from a to b is
//
//   d = x_b - x_a,  L = |d|,  n = d / L,  e = L - rest_length,
//   v_r = (v_b - v_a) . n,   F = (stiffness * e + damping * v_r) * n,
//   F_a = +F,  F_b = -F,
//
// exactly as in 2D (the axis is a Vec3 now, and the cross product plays no part
// because a linear spring acts only along its axis). Coincident endpoints
// (L == 0) contribute nothing, and springs are traversed in network order,
// which is observable because floating-point addition is not associative
// (docs/04).
void accumulate_spring_forces(const SpringNetwork& network,
                              const std::vector<SpringBody3D>& bodies,
                              std::vector<math::Vec3>& forces);

// A concrete 3D spring-network simulation, the 3D sibling of SpringWorld. The
// local composition boundary of ADR 0008, not the engine-wide force
// architecture; see it before generalising anything here.
class SpringWorld3D
{
public:
    SpringWorld3D(sim_core::SimulationSettings simulation_settings,
                  SpringNetwork network, std::vector<SpringBody3D> bodies);

    // Ok, or InvalidSettings (dt, or a timestep past a spring's stability
    // limit) / InvalidState (a body, a spring whose endpoints are invalid or
    // out of range, or a separation too large to square). Identical in meaning
    // to the 2D world's, and derived the same way: the per-spring symplectic
    // stability bound dt < 2 (sqrt(gamma^2 + omega^2) - gamma) / omega^2 with
    // omega^2 = k/mu and gamma = c/(2 mu).
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step, in the fixed order the 2D world uses:
    //
    //   1. clear the per-body force accumulator;
    //   2. evaluate springs in stable network order;
    //   3. accumulate equal and opposite endpoint forces;
    //   4. integrate bodies in stable index order, semi-implicit Euler.
    //
    // Both orders are part of the contract, not incidental (docs/04). On
    // validation failure the state is left unchanged and the failing status is
    // returned; never throws.
    sim_core::StepResult step();

    const std::vector<SpringBody3D>& bodies() const { return bodies_; }
    const SpringNetwork& network() const { return network_; }
    const std::vector<math::Vec3>& forces() const { return forces_; }

    std::uint64_t tick_count() const { return step_ ? step_->tick_count() : 0; }
    math::Real elapsed_time() const
    {
        return step_ ? step_->elapsed_time() : math::Real{0};
    }
    const sim_core::SimulationSettings& simulation_settings() const
    {
        return simulation_settings_;
    }

private:
    sim_core::SimulationSettings simulation_settings_;
    SpringNetwork network_;
    std::vector<SpringBody3D> bodies_;
    std::vector<math::Vec3> forces_;
    std::vector<SpringBody3D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Whole-system diagnostics (3D) ---
//
// Every spring contributes equal and opposite endpoint forces, so a network
// cannot change total momentum: it is conserved exactly, at any stiffness and
// any damping.
math::Vec3 total_momentum3d(const std::vector<SpringBody3D>& bodies);

math::Real total_kinetic_energy3d(const std::vector<SpringBody3D>& bodies);

// Elastic potential stored in the springs, the sum of
// (1/2) stiffness extension^2. Kinetic plus elastic is conserved only for an
// UNDAMPED network, and even then only up to the integrator's drift; damping
// removes it on purpose, which is what makes it testable.
math::Real total_elastic_energy3d(const SpringNetwork& network,
                                  const std::vector<SpringBody3D>& bodies);
} // namespace malloy::springs
