#pragma once

#include <malloy/math/real.hpp>

namespace malloy::nbody
{
// N-body-specific parameters. The fixed timestep is not here: that is generic
// and lives in sim_core::SimulationSettings (docs/03_MODULE_BOUNDARIES.md).
struct NBodySettings
{
    math::Real g{1.0};         // gravitational constant G (>= 0)
    math::Real softening{0.0}; // softening length to bound near-zero distances (>= 0)

    // Valid when G and softening are non-negative and finite.
    bool is_valid() const;
};
} // namespace malloy::nbody
