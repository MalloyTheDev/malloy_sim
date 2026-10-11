#pragma once

#include <cstdint>

#include <malloy/math/real.hpp>

namespace malloy::particles
{
// The random number generator of the particle domain (M39), the project's
// first: thermal walls draw from it, and nothing else does. SplitMix64 (Steele,
// Lea and Flood, OOPSLA 2014), fully specified here rather than taken from
// <random>, whose distributions the C++ standard leaves to the implementation
// (docs/04, ADR 0010). Its integer sequence is the same on every platform;
// doubles that later pass through log, cos or sin need not be.
//
// It is one 64-bit word of state, so a step that is rolled back restores it
// with a copy, and every seed, 0 included, is valid. After n draws the state is
// exactly seed + n * gamma (mod 2^64), which lets a test count draws. It passes
// TestU01's BigCrush, as its authors and Vigna report.
class SplitMix64
{
public:
    static constexpr std::uint64_t gamma = 0x9E3779B97F4A7C15ULL;

    constexpr explicit SplitMix64(std::uint64_t seed) : state_{seed} {}

    // The next 64-bit output.
    constexpr std::uint64_t next()
    {
        state_ += gamma;
        std::uint64_t z = state_;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    // A 64-bit output as a double strictly inside (0, 1): the top 52 bits pick
    // one of 2^52 equal cells and the result is that cell's MIDPOINT,
    // (k + 1/2) 2^-52, computed exactly. The extremes are 2^-53 and
    // 1 - 2^-53, so log(u) is always finite and never 0, a sampled speed is
    // never 0 or infinite, and u and 1 - u have the same distribution. (The
    // 53-bit form ((x >> 11) + 1/2) 2^-53 rounds its top value to exactly 1.)
    static constexpr math::Real to_unit_open(std::uint64_t bits)
    {
        return (static_cast<math::Real>(bits >> 12) + math::Real{0.5}) * 0x1p-52;
    }

    // The next draw as a uniform deviate in (0, 1).
    constexpr math::Real uniform() { return to_unit_open(next()); }

    // The current state, so a test can see whether, and how often, it drew.
    constexpr std::uint64_t state() const { return state_; }

private:
    std::uint64_t state_;
};
} // namespace malloy::particles
