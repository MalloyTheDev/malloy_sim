#pragma once

#include <cstdint>
#include <vector>

#include <malloy/math/vec2.hpp>
#include <malloy/nbody/body2d.hpp>
#include <malloy/nbody/nbody_settings.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::nbody
{
// A concrete 2D N-body gravity simulation.
//
// It owns its own fixed-step integration (semi-implicit Euler) and is
// deliberately NOT built on any polymorphic simulation base or engine kernel
// (docs/01_V3_ARCHITECTURE_DECISION.md, ADR 0004).
class NBodyWorld
{
public:
    NBodyWorld(sim_core::SimulationSettings simulation_settings,
               NBodySettings nbody_settings,
               std::vector<Body2D> bodies);

    // Check the settings and every body. Returns Ok, or the first failure
    // category found: InvalidSettings (dt/G/softening) or InvalidState (a body).
    sim_core::StepStatus validate() const;

    // Pairwise gravitational accelerations for the current state, with
    // softening (docs/04). accelerations[i] corresponds to bodies()[i].
    //
    // A pair whose softened separation is zero (coincident bodies with
    // softening == 0) contributes nothing rather than dividing by zero, so the
    // result stays finite for any state that passes validate(). Softening is
    // still the intended way to keep close encounters well behaved.
    std::vector<math::Vec2> compute_accelerations() const;

    // Advance the simulation by exactly one fixed step. On validation failure
    // the state is left unchanged and the failing status is returned; this
    // never throws (docs/04).
    //
    // The state is validated both before and after the update, so a step that
    // produces non-finite values (overflow from an extreme mass or dt, say)
    // returns InvalidState and rolls back instead of reporting success. A
    // successful step therefore always leaves every body finite.
    sim_core::StepResult step();

    const std::vector<Body2D>& bodies() const { return bodies_; }
    std::uint64_t tick_count() const { return tick_count_; }

    const sim_core::SimulationSettings& simulation_settings() const
    {
        return simulation_settings_;
    }
    const NBodySettings& nbody_settings() const { return nbody_settings_; }

private:
    sim_core::SimulationSettings simulation_settings_;
    NBodySettings nbody_settings_;
    std::vector<Body2D> bodies_;
    // Scratch copy of the pre-step state, used only to roll back a step that
    // produced non-finite values. Kept as a member so the common case reuses
    // its allocation instead of allocating once per step.
    std::vector<Body2D> previous_;
    std::uint64_t tick_count_{0};
};
} // namespace malloy::nbody
