#include <malloy/scenario/scenario.hpp>

#include <malloy/nbody/nbody.hpp>
#include <malloy/particles/particles.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <test_check.hpp>

#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>

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

    // --- Issue #4: a run length that cannot be run is rejected, with the line
    //     number, rather than silently producing a zero-step or endless run. ---
    {
        std::istringstream in("dt 0.001\nsteps -5\nbody 1.0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_TRUE(r.error.find("line 2") != std::string::npos);
        MALLOY_CHECK_TRUE(r.error.find("steps") != std::string::npos);
    }
    {
        std::istringstream in("output_every -7\nbody 1.0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_TRUE(r.error.find("line 1") != std::string::npos);
    }
    {
        // Zero remains legal for both: zero steps, and reporting disabled.
        std::istringstream in("steps 0\noutput_every 0\nbody 1.0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_EQ(r.scenario.steps, 0);
        MALLOY_CHECK_EQ(r.scenario.output_every, 0);
    }

    // --- Issue #6: trailing tokens are an error, not silently discarded. The
    //     motivating case is a body line written with 3D fields, which used to
    //     be truncated to the first five and run as a different simulation. ---
    {
        std::istringstream in("body 1.0 0.0 0.0 0.0 0.0 1.0 2.0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_TRUE(r.error.find("line 1") != std::string::npos);
        MALLOY_CHECK_TRUE(r.error.find("1.0") != std::string::npos); // names the token
    }
    {
        std::istringstream in("dt 0.5 GARBAGE\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_TRUE(r.error.find("GARBAGE") != std::string::npos);
    }
    {
        // A trailing comment is not a trailing token: comments are stripped
        // before parsing, so this must still be accepted.
        std::istringstream in("dt 0.5   # the timestep\nbody 1.0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.simulation.dt, 0.5, eps);
    }
    {
        // Trailing whitespace, including CRLF line endings, is not a token.
        std::istringstream in("dt 0.5  \r\nbody 1.0 0 0 0 0  \r\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_EQ(r.scenario.bodies.size(), std::size_t{1});
    }

    // --- Issue #7: a stream that goes bad is a parse failure, not a success
    //     with whatever was read so far. ---
    {
        std::istringstream in("dt 0.001\nbody 1.0 0 0 0 0\n");
        in.setstate(std::ios::badbit);
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_FALSE(r.error.empty());
    }

    // --- Issue #12: every shipped template must parse and build a runnable
    //     world. They are a deliverable (CLAUDE.md rule 16) and were previously
    //     the least-tested files in the repository. ---
#ifdef MALLOY_SCENARIO_DIR
    {
        const char* templates[] = {"two_body.scn", "three_body_triangle.scn",
                                   "bouncing_particles.scn"};
        for (const char* name : templates)
        {
            const std::string path = std::string(MALLOY_SCENARIO_DIR) + "/" + name;
            const ScenarioParseResult r =
                malloy::scenario::parse_scenario_file(path);
            if (!r.ok)
            {
                std::cerr << "template " << name << " failed to parse: " << r.error
                          << '\n';
                return 1;
            }
            MALLOY_CHECK_TRUE(r.scenario.steps > 0);

            // Same dispatch the app performs, so a template is validated by the
            // domain it actually declares.
            if (r.scenario.type == malloy::scenario::ScenarioType::NBody)
            {
                MALLOY_CHECK_TRUE(r.scenario.bodies.size() >= 2);
                NBodyWorld world{r.scenario.simulation, r.scenario.nbody_settings,
                                 r.scenario.bodies};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            else
            {
                MALLOY_CHECK_TRUE(r.scenario.particle_list.size() >= 2);
                malloy::particles::ParticleWorld world{r.scenario.simulation,
                                                      r.scenario.particle_settings,
                                                      r.scenario.particle_list};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                MALLOY_CHECK_TRUE(world.step().ok());
            }
        }
    }
#endif

    // --- Issue #12: parse_scenario_file reports a missing file as a parse
    //     failure rather than throwing (documented in scenario.hpp). ---
    {
        const ScenarioParseResult r =
            malloy::scenario::parse_scenario_file("no_such_scenario_12345.scn");
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_FALSE(r.error.empty());
    }

    // --- Issue: multi-domain dispatch. An absent `type` still means nbody, so
    //     every scenario written before the key existed keeps working. ---
    {
        std::istringstream in("dt 0.5\ng 2.0\nbody 1.0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::NBody);
        MALLOY_CHECK_EQ(r.scenario.bodies.size(), std::size_t{1});
    }
    {
        // Stating it explicitly is identical.
        std::istringstream in("type nbody\ng 2.0\nbody 1.0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::NBody);
    }

    // --- A particles scenario parses into the particle fields. ---
    {
        std::istringstream in("type particles\n"
                              "dt 0.01\n"
                              "steps 50\n"
                              "restitution 0.25\n"
                              "bounds -2 -3 4 5\n"
                              "particle 1.5 0.4 1.0 2.0 3.0 4.0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Particles);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.restitution, 0.25, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.bounds.min.x, -2.0, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.bounds.min.y, -3.0, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.bounds.max.x, 4.0, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.bounds.max.y, 5.0, eps);
        MALLOY_CHECK_EQ(r.scenario.particle_list.size(), std::size_t{1});
        // Every field distinct, so a swapped pair in the constructor shows up.
        MALLOY_CHECK_NEAR(r.scenario.particle_list[0].mass, 1.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_list[0].radius, 0.4, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_list[0].position.x, 1.0, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_list[0].position.y, 2.0, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_list[0].velocity.x, 3.0, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_list[0].velocity.y, 4.0, eps);
    }

    // --- A key from the wrong domain is an error, so a typo in `type` surfaces
    //     immediately instead of silently running the wrong simulation. ---
    {
        std::istringstream in("type particles\nbody 1.0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_TRUE(r.error.find("line 2") != std::string::npos);
    }
    {
        std::istringstream in("type nbody\nparticle 1.0 0.5 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
    }
    {
        std::istringstream in("type nbody\nrestitution 0.5\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
    }
    {
        std::istringstream in("type particles\ng 1.0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
    }

    // --- An unknown type is rejected by name. ---
    {
        std::istringstream in("type fluid\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_TRUE(r.error.find("fluid") != std::string::npos);
    }

    // --- `type` after a domain key is rejected: the keys already read would
    //     have been validated against the wrong domain. ---
    {
        std::istringstream in("body 1.0 0 0 0 0\ntype particles\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_TRUE(r.error.find("line 2") != std::string::npos);
    }

    // --- Field counts on the new keys. ---
    {
        std::istringstream in("type particles\nparticle 1.0 0.5 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok); // one field short
    }
    {
        std::istringstream in("type particles\nbounds -1 -1 1\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
    }
    {
        std::istringstream in("type particles\nbounds -1 -1 1 1 EXTRA\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok); // trailing tokens still rejected
    }

    std::cout << "malloy_scenario_tests passed\n";
    return 0;
}
