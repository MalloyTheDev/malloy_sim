#pragma once

#include <cstdint>
#include <optional>

#include <malloy/math/real.hpp>

namespace malloy::time
{
// A fixed simulation timestep paired with a running tick counter.
//
// There is no wall clock here: "time" means accumulated *simulation* time
// only. Frame pacing, sleeping and real-time loops are out of scope
// (docs/03_MODULE_BOUNDARIES.md).
class FixedStep
{
public:
    // Returns a FixedStep when dt is a usable timestep (strictly positive and
    // finite), and std::nullopt otherwise. Validation never throws; invalid
    // input is a normal, reportable outcome (docs/04).
    static std::optional<FixedStep> create(math::Real dt)
    {
        if (dt > math::Real{0} && math::is_finite(dt))
        {
            return FixedStep{dt};
        }
        return std::nullopt;
    }

    math::Real dt() const { return dt_; }
    std::uint64_t tick_count() const { return tick_count_; }

    // Accumulated simulation time. Computed as tick_count * dt so it cannot
    // drift the way repeated floating-point addition would.
    math::Real elapsed_time() const
    {
        return static_cast<math::Real>(tick_count_) * dt_;
    }

    // Advance the simulation by exactly one fixed step.
    void advance() { ++tick_count_; }

private:
    explicit FixedStep(math::Real dt) : dt_{dt} {}

    math::Real dt_;
    std::uint64_t tick_count_{0};
};
} // namespace malloy::time
