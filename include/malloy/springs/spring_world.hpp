#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <malloy/math/vec2.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/springs/spring.hpp>
#include <malloy/time/fixed_step.hpp>

namespace malloy::springs
{
// A translational body: position, velocity, mass. No radius and no angle,
// because a spring needs neither.
//
// Its own type rather than a reused nbody::Body2D or particles::Particle2D, so
// malloy_springs depends on no other domain (ADR 0006, and the same reasoning
// that gave M11 its own RigidBody2D).
struct SpringBody2D
{
    math::Vec2 position{};
    math::Vec2 velocity{};
    math::Real mass{1.0};

    // Valid when mass is strictly positive and every value is finite.
    bool is_valid() const;
};

// Accumulate the force each spring exerts on its two endpoints.
//
// This is the pure kernel: it reads body state and adds into `forces`, and it
// neither integrates nor mutates a body. Keeping it free of integration is what
// lets a later force producer reuse the same pipeline (ADR 0008).
//
// It ADDS rather than assigns, so a body that is an endpoint of several springs
// receives every contribution. `forces` must already be sized to match
// `bodies`; it is not resized or cleared here, because clearing is the caller's
// first step in the documented contract.
//
// Derivation and sign convention, for one spring from body a to body b:
//
//   d   = x_b - x_a
//   L   = |d|
//   n   = d / L                 unit vector pointing from a toward b
//   e   = L - rest_length       positive when stretched
//   v_r = (v_b - v_a) . n       positive when the endpoints are separating
//
//   F   = (stiffness * e + damping * v_r) * n
//   F_a = +F
//   F_b = -F
//
// Check the signs on both cases. Stretched (e > 0) puts F along +n, so a is
// pulled toward b and b toward a: the spring contracts. Compressed (e < 0) puts
// F along -n and they push apart. Separating (v_r > 0) likewise gives a force
// that opposes the separation. F_a + F_b is exactly zero by construction.
//
// The damping term uses the AXIAL relative velocity, not the full relative
// speed. A damper resists motion along its own axis; motion perpendicular to
// the spring is not the damper's business, and using the full magnitude would
// also lose the sign that makes it oppose rather than drive the motion.
//
// A spring whose endpoints are coincident (L == 0) has no axis, so it
// contributes nothing. Any direction would be equally arbitrary, and picking
// one would invent a force out of a degenerate configuration.
//
// Springs are traversed in network order, and floating-point addition is not
// associative, so that order is part of the observable behaviour (docs/04).
void accumulate_spring_forces(const SpringNetwork& network,
                              const std::vector<SpringBody2D>& bodies,
                              std::vector<math::Vec2>& forces);

// A concrete spring-network simulation.
//
// The M13 composition boundary, not the engine-wide force architecture. See
// ADR 0008 before generalising anything here.
class SpringWorld
{
public:
    SpringWorld(sim_core::SimulationSettings simulation_settings,
                SpringNetwork network, std::vector<SpringBody2D> bodies);

    // Ok, or InvalidSettings (dt, or a timestep past a spring's stability
    // limit) / InvalidState (a body, a spring whose endpoints are invalid or
    // out of range, or a separation too large to square).
    //
    // The last two are easy to violate and impossible to notice otherwise.
    //
    // A spring whose endpoints are more than about 1.34e154 apart has a
    // separation whose SQUARE overflows, and the force kernel cannot tell that
    // from coincident endpoints: it skips the spring, which then silently
    // stops existing while the bodies coast apart forever and every step still
    // reports Ok.
    //
    // Symplectic Euler on a spring pair diverges past a bounded timestep,
    // geometrically, and until the separation reaches the range above nothing
    // reports that either. The bound is
    // dt < 2 (sqrt(gamma^2 + omega^2) - gamma) / omega^2 with omega^2 = k/mu
    // and gamma = c/(2 mu), reducing to dt < 2/omega when undamped.
    //
    // The stability check is per spring, which is NECESSARY but not sufficient
    // for a network: a spring that is unstable alone is unstable in company,
    // but a network of individually safe springs can still be too stiff as a
    // whole.
    sim_core::StepStatus validate() const;

    // Advance by exactly one fixed step, in this fixed order:
    //
    //   1. clear the per-body force accumulator;
    //   2. evaluate springs in stable network order;
    //   3. accumulate equal and opposite endpoint forces;
    //   4. integrate bodies in stable index order.
    //
    // Steps 2 and 3 are one pass through accumulate_spring_forces. Both orders
    // are part of the contract rather than incidental: floating-point addition
    // is not associative, so traversing either list differently changes the
    // result (docs/04).
    //
    // Integration is semi-implicit (symplectic) Euler, matching every other
    // world here: acceleration from the accumulated force, then velocity, then
    // position from the updated velocity.
    //
    // On validation failure, before or after the update, the state is left
    // unchanged and the failing status is returned. Never throws.
    sim_core::StepResult step();

    const std::vector<SpringBody2D>& bodies() const { return bodies_; }
    const SpringNetwork& network() const { return network_; }

    // The forces accumulated during the most recent successful step. Empty
    // before the first one. Exposed so a caller can see what the springs did
    // without re-deriving it.
    const std::vector<math::Vec2>& forces() const { return forces_; }

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
    std::vector<SpringBody2D> bodies_;
    std::vector<math::Vec2> forces_;
    std::vector<SpringBody2D> previous_;
    std::optional<time::FixedStep> step_;
};

// --- Whole-system diagnostics ---
//
// Every spring contributes equal and opposite endpoint forces, so a spring
// network cannot change total momentum. It is conserved exactly, at any
// stiffness and any damping.
math::Vec2 total_momentum(const std::vector<SpringBody2D>& bodies);

math::Real total_kinetic_energy(const std::vector<SpringBody2D>& bodies);

// Elastic potential stored in the springs: the sum of
// (1/2) * stiffness * extension^2.
//
// Kinetic plus elastic is conserved only for an UNDAMPED network, and even then
// only up to the integrator's own drift. Damping removes energy on purpose,
// which is what makes it testable.
math::Real total_elastic_energy(const SpringNetwork& network,
                                const std::vector<SpringBody2D>& bodies);
} // namespace malloy::springs
