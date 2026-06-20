#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include <malloy/nbody/body2d.hpp>
#include <malloy/nbody/nbody_settings.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::scenario
{
// A complete, runnable N-body scenario: the fixed timestep, the N-body
// parameters, the bodies, and how long to run. This is the parsed form of a
// scenario file.
struct Scenario
{
    sim_core::SimulationSettings simulation{0.001};
    nbody::NBodySettings nbody_settings{1.0, 0.0};
    std::vector<nbody::Body2D> bodies;
    int steps{1000};
    int output_every{100};
};

// The outcome of parsing. On failure `error` holds a human-readable message
// (with the offending line number) and `scenario` is unspecified. Parsing
// never throws (docs/04: status/result returns for normal failures).
struct ScenarioParseResult
{
    bool ok{false};
    std::string error;
    Scenario scenario;
};

// Parse a scenario from a text stream. One directive per line; '#' starts a
// comment; blank lines are ignored:
//
//   g <value>             gravitational constant   (default 1.0)
//   softening <value>     softening length         (default 0.0)
//   dt <value>            fixed timestep           (default 0.001)
//   steps <value>         number of steps          (default 1000)
//   output_every <value>  steps between reports    (default 100)
//   body <mass> <px> <py> <vx> <vy>
//
// Only syntax is checked here (bad numbers, wrong field counts, unknown keys).
// Semantic validity (dt > 0, mass > 0, ...) is left to NBodyWorld::validate so
// that rule lives in exactly one place.
ScenarioParseResult parse_scenario(std::istream& input);

// Convenience: open a file and parse it. A missing or unreadable file is
// reported as a parse failure, not an exception.
ScenarioParseResult parse_scenario_file(const std::string& path);
} // namespace malloy::scenario
