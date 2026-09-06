#include <malloy/scenario/scenario.hpp>

#include <cstddef>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>

#include <malloy/math/vec2.hpp>

namespace malloy::scenario
{
namespace
{
ScenarioParseResult make_error(std::size_t line_number, const std::string& message)
{
    return ScenarioParseResult{false,
                               "line " + std::to_string(line_number) + ": " + message,
                               {}};
}

// True when nothing but whitespace remains. A directive that parsed its fields
// and left tokens behind is a malformed line, not a valid one: silently
// dropping the remainder would accept a body line written with 3D fields and
// run a different simulation than the file describes.
bool at_end_of_directive(std::istringstream& tokens)
{
    tokens >> std::ws;
    return tokens.eof();
}
} // namespace

ScenarioParseResult parse_scenario(std::istream& input)
{
    Scenario scenario;
    std::string line;
    std::size_t line_number = 0;

    while (std::getline(input, line))
    {
        ++line_number;

        const std::size_t comment = line.find('#');
        if (comment != std::string::npos)
        {
            line.erase(comment);
        }

        std::istringstream tokens(line);
        std::string key;
        if (!(tokens >> key))
        {
            continue; // blank or comment-only line
        }

        if (key == "body")
        {
            math::Real mass{};
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (!(tokens >> mass >> px >> py >> vx >> vy))
            {
                return make_error(line_number, "body requires: mass px py vx vy");
            }
            scenario.bodies.push_back(
                nbody::Body2D{math::Vec2{px, py}, math::Vec2{vx, vy}, mass});
        }
        else if (key == "dt")
        {
            if (!(tokens >> scenario.simulation.dt))
            {
                return make_error(line_number, "dt requires a number");
            }
        }
        else if (key == "g")
        {
            if (!(tokens >> scenario.nbody_settings.g))
            {
                return make_error(line_number, "g requires a number");
            }
        }
        else if (key == "softening")
        {
            if (!(tokens >> scenario.nbody_settings.softening))
            {
                return make_error(line_number, "softening requires a number");
            }
        }
        else if (key == "steps")
        {
            if (!(tokens >> scenario.steps))
            {
                return make_error(line_number, "steps requires an integer");
            }
            if (scenario.steps < 0)
            {
                return make_error(line_number, "steps must be >= 0");
            }
        }
        else if (key == "output_every")
        {
            if (!(tokens >> scenario.output_every))
            {
                return make_error(line_number, "output_every requires an integer");
            }
            if (scenario.output_every < 0)
            {
                return make_error(line_number, "output_every must be >= 0");
            }
        }
        else
        {
            return make_error(line_number, "unknown key '" + key + "'");
        }

        if (!at_end_of_directive(tokens))
        {
            std::string extra;
            tokens >> extra;
            return make_error(line_number, "unexpected extra token '" + extra + "'");
        }
    }

    // --- #7: getline stops on a read error as well as on end of input, so a
    //     truncated read would otherwise be reported as a successful parse.
    if (input.bad())
    {
        return make_error(line_number, "read error");
    }

    return ScenarioParseResult{true, "", std::move(scenario)};
}

ScenarioParseResult parse_scenario_file(const std::string& path)
{
    // Binary mode on purpose: in text mode the MSVC CRT treats 0x1A as end of
    // file, silently truncating a scenario mid-parse and reporting success.
    // CRLF still parses, because '' is whitespace to operator>>.
    std::ifstream file(path, std::ios::binary);
    if (!file)
    {
        return ScenarioParseResult{false, "cannot open scenario file: " + path, {}};
    }
    return parse_scenario(file);
}
} // namespace malloy::scenario
