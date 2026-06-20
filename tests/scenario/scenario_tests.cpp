#include <malloy/scenario/scenario.hpp>

#include <malloy/nbody/nbody.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <test_check.hpp>

#include <cstddef>
#include <iostream>
#include <sstream>

using malloy::nbody::NBodyWorld;
using malloy::scenario::parse_scenario;
using malloy::scenario::ScenarioParseResult;
using malloy::sim_core::StepStatus;

int main()
{
    const double eps = 1e-12;

    // A well-formed scenario parses into the expected settings and bodies.
    {
        std::istringstream input(
            "# two bodies\n"
            "g 2.0\n"
            "softening 0.5\n"
            "dt 0.01\n"
            "steps 500\n"
            "output_every 50\n"
            "body 1.0 0.0 0.0 0.0 0.0\n"
            "body 3.0 4.0 0.0 0.0 2.0   # trailing comment ignored\n");
        const ScenarioParseResult result = parse_scenario(input);
        MALLOY_CHECK_TRUE(result.ok);
        MALLOY_CHECK_NEAR(result.scenario.nbody_settings.g, 2.0, eps);
        MALLOY_CHECK_NEAR(result.scenario.nbody_settings.softening, 0.5, eps);
        MALLOY_CHECK_NEAR(result.scenario.simulation.dt, 0.01, eps);
        MALLOY_CHECK_EQ(result.scenario.steps, 500);
        MALLOY_CHECK_EQ(result.scenario.output_every, 50);
        MALLOY_CHECK_EQ(result.scenario.bodies.size(), std::size_t{2});
        MALLOY_CHECK_NEAR(result.scenario.bodies[1].mass, 3.0, eps);
        MALLOY_CHECK_NEAR(result.scenario.bodies[1].position.x, 4.0, eps);
        MALLOY_CHECK_NEAR(result.scenario.bodies[1].velocity.y, 2.0, eps);
    }

    // Comments and blank lines are ignored; omitted keys take their defaults.
    {
        std::istringstream input(
            "\n# only one body, everything else default\n\nbody 1.0 1.0 0.0 0.0 1.0\n");
        const ScenarioParseResult result = parse_scenario(input);
        MALLOY_CHECK_TRUE(result.ok);
        MALLOY_CHECK_EQ(result.scenario.bodies.size(), std::size_t{1});
        MALLOY_CHECK_NEAR(result.scenario.simulation.dt, 0.001, eps); // default
        MALLOY_CHECK_NEAR(result.scenario.nbody_settings.g, 1.0, eps); // default
    }

    // A parsed scenario builds a world that validates.
    {
        std::istringstream input(
            "g 1.0\nsoftening 0.000001\ndt 0.001\n"
            "body 1.0 0.0 0.0 0.0 0.0\nbody 0.000001 1.0 0.0 0.0 1.0\n");
        const ScenarioParseResult result = parse_scenario(input);
        MALLOY_CHECK_TRUE(result.ok);
        NBodyWorld world{result.scenario.simulation, result.scenario.nbody_settings,
                         result.scenario.bodies};
        MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
    }

    // Syntax errors are reported (never thrown) with a non-empty message.
    {
        std::istringstream bad("dt abc\n"); // not a number
        const ScenarioParseResult result = parse_scenario(bad);
        MALLOY_CHECK_FALSE(result.ok);
        MALLOY_CHECK_FALSE(result.error.empty());
    }
    {
        std::istringstream bad("body 1.0 2.0\n"); // too few fields
        const ScenarioParseResult result = parse_scenario(bad);
        MALLOY_CHECK_FALSE(result.ok);
    }
    {
        std::istringstream bad("gravity 1.0\n"); // unknown key
        const ScenarioParseResult result = parse_scenario(bad);
        MALLOY_CHECK_FALSE(result.ok);
    }

    std::cout << "malloy_scenario_tests passed\n";
    return 0;
}
