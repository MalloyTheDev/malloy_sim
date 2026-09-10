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

    // Valid when G and softening are non-negative and finite, and when the
    // SQUARE of the softening is finite too.
    //
    // The squared test is not redundant with the finiteness test. Softening
    // enters the denominator SQUARED, and a finite value above sqrt(DBL_MAX),
    // about 1.34e154, squares to infinity. The consequences are silent and
    // look like success: every separation becomes infinite, so every force is
    // exactly zero and every potential is exactly zero, and a simulation in
    // which nothing happens conserves everything perfectly.
    bool is_valid() const;
};
} // namespace malloy::nbody
