#include <malloy/time/fixed_step.hpp>

#include <test_check.hpp>

#include <cstdint>
#include <iostream>
#include <limits>

int main()
{
    using malloy::math::Real;
    using malloy::time::FixedStep;

    const Real eps = 1e-12;

    // A valid timestep yields a step that starts at tick 0 / elapsed 0.
    {
        const auto step = FixedStep::create(0.5);
        MALLOY_CHECK_TRUE(step.has_value());
        MALLOY_CHECK_NEAR(step->dt(), 0.5, eps);
        MALLOY_CHECK_EQ(step->tick_count(), std::uint64_t{0});
        MALLOY_CHECK_NEAR(step->elapsed_time(), 0.0, eps);
    }

    // Invalid timesteps are rejected without throwing.
    {
        MALLOY_CHECK_FALSE(FixedStep::create(0.0).has_value());
        MALLOY_CHECK_FALSE(FixedStep::create(-1.0).has_value());
        MALLOY_CHECK_FALSE(FixedStep::create(std::numeric_limits<Real>::infinity()).has_value());
        MALLOY_CHECK_FALSE(FixedStep::create(std::numeric_limits<Real>::quiet_NaN()).has_value());
    }

    // advance() moves exactly one fixed step; elapsed_time == ticks * dt.
    {
        auto step = FixedStep::create(0.25);
        MALLOY_CHECK_TRUE(step.has_value());

        for (int i = 0; i < 4; ++i)
        {
            step->advance();
        }

        MALLOY_CHECK_EQ(step->tick_count(), std::uint64_t{4});
        MALLOY_CHECK_NEAR(step->elapsed_time(), 1.0, eps);
    }

    std::cout << "malloy_time_tests passed\n";
    return 0;
}
