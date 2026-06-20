#include <malloy/sim_core/sim_core.hpp>

#include <test_check.hpp>

#include <iostream>
#include <limits>

int main()
{
    using malloy::math::Real;
    using malloy::sim_core::SimulationSettings;
    using malloy::sim_core::StepResult;
    using malloy::sim_core::StepStatus;

    // A default StepResult means success.
    {
        const StepResult result;
        MALLOY_CHECK_TRUE(result.status == StepStatus::Ok);
        MALLOY_CHECK_TRUE(result.ok());
    }

    // A non-Ok status is reported as not-ok.
    {
        const StepResult settings_failure{StepStatus::InvalidSettings};
        MALLOY_CHECK_FALSE(settings_failure.ok());
        MALLOY_CHECK_TRUE(settings_failure.status == StepStatus::InvalidSettings);

        const StepResult state_failure{StepStatus::InvalidState};
        MALLOY_CHECK_FALSE(state_failure.ok());
    }

    // The status categories are distinct.
    {
        MALLOY_CHECK_TRUE(StepStatus::Ok != StepStatus::InvalidSettings);
        MALLOY_CHECK_TRUE(StepStatus::Ok != StepStatus::InvalidState);
        MALLOY_CHECK_TRUE(StepStatus::InvalidSettings != StepStatus::InvalidState);
    }

    // SimulationSettings validates its timestep (positive, finite), no throw.
    {
        MALLOY_CHECK_TRUE(SimulationSettings{0.001}.is_valid());
        MALLOY_CHECK_FALSE(SimulationSettings{0.0}.is_valid());
        MALLOY_CHECK_FALSE(SimulationSettings{-1.0}.is_valid());
        MALLOY_CHECK_FALSE(SimulationSettings{std::numeric_limits<Real>::infinity()}.is_valid());
        MALLOY_CHECK_FALSE(SimulationSettings{std::numeric_limits<Real>::quiet_NaN()}.is_valid());
    }

    std::cout << "malloy_sim_core_tests passed\n";
    return 0;
}
