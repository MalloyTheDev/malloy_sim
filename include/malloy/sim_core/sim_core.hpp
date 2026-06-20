#pragma once

#include <malloy/math/real.hpp>

namespace malloy::sim_core
{
// Outcome category for a single simulation step.
//
// Deliberately coarse and domain-agnostic: concrete simulations map their
// specific validation failures onto these categories. Keeping this tiny is a
// hard rule (ADR 0004) — no engine kernel vocabulary belongs here.
enum class StepStatus
{
    Ok,
    InvalidSettings, // configuration is unusable (e.g. non-positive timestep)
    InvalidState,    // simulation state is unusable (e.g. non-finite values)
};

// The result of advancing a simulation by one step. A default-constructed
// result means success.
struct StepResult
{
    StepStatus status{StepStatus::Ok};

    bool ok() const;
};

// Simulation-wide settings shared by every simulation. The fixed timestep is
// the one knob they all have in common; domain-specific parameters live in
// their own settings types (e.g. NBodySettings in malloy_nbody).
struct SimulationSettings
{
    math::Real dt{0.0};

    bool is_valid() const;
};
} // namespace malloy::sim_core
