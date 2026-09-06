#include <malloy/scenario/scenario.hpp>

#include <cstddef>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>

#include <malloy/collide/shapes.hpp>
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
    // `type` selects which domain keys are legal, so it has to be settled
    // before any of them appear.
    bool saw_domain_key = false;

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

        if (key == "type")
        {
            std::string value;
            if (!(tokens >> value))
            {
                return make_error(line_number, "type requires nbody or particles");
            }
            if (value == "nbody")
            {
                scenario.type = ScenarioType::NBody;
            }
            else if (value == "particles")
            {
                scenario.type = ScenarioType::Particles;
            }
            else
            {
                return make_error(line_number, "unknown type '" + value + "'");
            }
            if (saw_domain_key)
            {
                return make_error(line_number,
                                  "type must come before any domain-specific key");
            }
        }
        else if (key == "body")
        {
            math::Real mass{};
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (scenario.type != ScenarioType::NBody)
            {
                return make_error(line_number, "body belongs to type nbody");
            }
            saw_domain_key = true;
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
            if (scenario.type != ScenarioType::NBody)
            {
                return make_error(line_number, "g belongs to type nbody");
            }
            saw_domain_key = true;
            if (!(tokens >> scenario.nbody_settings.g))
            {
                return make_error(line_number, "g requires a number");
            }
        }
        else if (key == "softening")
        {
            if (scenario.type != ScenarioType::NBody)
            {
                return make_error(line_number, "softening belongs to type nbody");
            }
            saw_domain_key = true;
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
        else if (key == "restitution")
        {
            if (scenario.type != ScenarioType::Particles)
            {
                return make_error(line_number, "restitution belongs to type particles");
            }
            saw_domain_key = true;
            if (!(tokens >> scenario.particle_settings.restitution))
            {
                return make_error(line_number, "restitution requires a number");
            }
        }
        else if (key == "bounds")
        {
            if (scenario.type != ScenarioType::Particles)
            {
                return make_error(line_number, "bounds belongs to type particles");
            }
            saw_domain_key = true;
            math::Real min_x{};
            math::Real min_y{};
            math::Real max_x{};
            math::Real max_y{};
            if (!(tokens >> min_x >> min_y >> max_x >> max_y))
            {
                return make_error(line_number, "bounds requires: minx miny maxx maxy");
            }
            scenario.particle_settings.bounds =
                collide::Aabb{math::Vec2{min_x, min_y}, math::Vec2{max_x, max_y}};
        }
        else if (key == "particle")
        {
            if (scenario.type != ScenarioType::Particles)
            {
                return make_error(line_number, "particle belongs to type particles");
            }
            saw_domain_key = true;
            math::Real mass{};
            math::Real radius{};
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (!(tokens >> mass >> radius >> px >> py >> vx >> vy))
            {
                return make_error(line_number,
                                  "particle requires: mass radius px py vx vy");
            }
            scenario.particle_list.push_back(particles::Particle2D{
                math::Vec2{px, py}, math::Vec2{vx, vy}, mass, radius});
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
    // CRLF still parses, because a carriage return is whitespace to the
    // stream extractor.
    std::ifstream file(path, std::ios::binary);
    if (!file)
    {
        return ScenarioParseResult{false, "cannot open scenario file: " + path, {}};
    }
    return parse_scenario(file);
}
} // namespace malloy::scenario
