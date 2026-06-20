#include <malloy/sim_core/sim_core.hpp>

#include <malloy/math/real.hpp>

namespace malloy::sim_core
{
bool StepResult::ok() const
{
    return status == StepStatus::Ok;
}

bool SimulationSettings::is_valid() const
{
    return dt > math::Real{0} && math::is_finite(dt);
}
} // namespace malloy::sim_core
