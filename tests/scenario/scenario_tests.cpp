#include <malloy/scenario/scenario.hpp>

#include <malloy/charges/charges.hpp>
#include <malloy/nbody/nbody.hpp>
#include <malloy/particles/particles.hpp>
#include <malloy/rigid/rigid.hpp>
#include <malloy/springs/springs.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <test_check.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <cctype>
#include <iterator>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <optional>

using malloy::nbody::NBodyWorld;
using malloy::scenario::parse_scenario;
using malloy::scenario::ScenarioParseResult;
using malloy::math::Real;
using malloy::sim_core::StepStatus;

namespace
{
// One `# check` line from a scenario template, which is issue #14: the
// documented figures are a deliverable, and until now nothing compared them
// against a run.
//
// The form is
//
//     # check step <n> <quantity> <value> tol <t>
//
// with `step 0` meaning the state before any step is taken, matching the
// column the app prints. The tolerance is ABSOLUTE and is required rather than
// defaulted, so a template states the precision it is claiming instead of
// inheriting one.
struct TemplateCheck
{
    int step{0};
    std::string quantity;
    malloy::math::Real value{0.0};
    malloy::math::Real tolerance{0.0};
};

std::vector<TemplateCheck> parse_checks(const std::string& text)
{
    std::vector<TemplateCheck> checks;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line))
    {
        std::istringstream tokens(line);
        std::string hash;
        std::string word;
        if (!(tokens >> hash >> word) || hash != "#" || word != "check")
        {
            continue;
        }
        TemplateCheck check;
        std::string step_word;
        std::string tol_word;
        if (!(tokens >> step_word >> check.step >> check.quantity >> check.value >>
              tol_word >> check.tolerance) ||
            step_word != "step" || tol_word != "tol")
        {
            checks.clear();
            checks.push_back(TemplateCheck{-1, "malformed", 0.0, 0.0});
            return checks;
        }
        checks.push_back(check);
    }
    return checks;
}

// Evaluates one check. `value_of` maps a quantity name to its current value, or
// nothing if this domain has no such quantity: `angular` exists only for nbody
// and rigid, and asking for it elsewhere must FAIL rather than pass silently.
template <typename ValueOf>
bool check_passes(const std::string& name, const TemplateCheck& check, ValueOf value_of)
{
    const std::optional<malloy::math::Real> actual = value_of(check.quantity);
    if (!actual)
    {
        std::cerr << "template " << name << ": unknown or unavailable quantity '"
                  << check.quantity << "' at step " << check.step << '\n';
        return false;
    }
    if (!(std::abs(*actual - check.value) <= check.tolerance))
    {
        std::cerr << "template " << name << ": " << check.quantity << " at step "
                  << check.step << " is " << *actual << ", documented as "
                  << check.value << " (tolerance " << check.tolerance << ")" << '\n';
        return false;
    }
    return true;
}

// Runs a world to `steps`, evaluating checks at the step they name.
//
// A function template over the world type rather than a base class: the five
// worlds share no ancestor (ADR 0006), only the shape of step(). The same
// approach the terminal app uses for drive().
template <typename World, typename ValueOf>
bool run_checked(World& world, const std::string& name, int steps,
                 const std::vector<TemplateCheck>& checks, ValueOf value_of)
{
    for (const TemplateCheck& check : checks)
    {
        if (check.step == 0 && !check_passes(name, check, value_of))
        {
            return false;
        }
    }
    for (int i = 1; i <= steps; ++i)
    {
        if (!world.step().ok())
        {
            std::cerr << "template " << name << " failed at step " << i << '\n';
            return false;
        }
        for (const TemplateCheck& check : checks)
        {
            if (check.step == i && !check_passes(name, check, value_of))
            {
                return false;
            }
        }
    }
    return true;
}
} // namespace

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
        // Enumerated from the directory rather than listed by hand. A hardcoded
        // list means a template added without editing this file is silently
        // untested, which would quietly break the rule that every template must
        // have a test or a documented expected result (CLAUDE.md rule 16).
        std::vector<std::string> templates;
        for (const auto& entry :
             std::filesystem::directory_iterator(MALLOY_SCENARIO_DIR))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".scn")
            {
                templates.push_back(entry.path().filename().string());
            }
        }
        // Sorted, so a failure names the same file on every platform.
        std::sort(templates.begin(), templates.end());
        // The directory must not be empty: an empty enumeration would make this
        // whole block silently vacuous.
        MALLOY_CHECK_TRUE(templates.size() >= 5);

        // Counts stated in prose go stale. These are read back out of the
        // documents and checked against reality, because a wrong number in the
        // file a reader believes is a defect even when the code is correct.
        {
            const auto claimed_before =
                [](const std::string& path,
                   const std::string& marker) -> std::optional<std::size_t> {
                std::ifstream document(path);
                if (!document.good())
                {
                    return std::nullopt;
                }
                const std::string text((std::istreambuf_iterator<char>(document)),
                                       std::istreambuf_iterator<char>());
                const std::size_t at = text.find(marker);
                if (at == std::string::npos)
                {
                    return std::nullopt;
                }
                // Walk back over the digits immediately before the marker.
                std::size_t first = at;
                while (first > 0 &&
                       std::isdigit(static_cast<unsigned char>(text[first - 1])) != 0)
                {
                    --first;
                }
                if (first == at)
                {
                    return std::nullopt;
                }
                return static_cast<std::size_t>(
                    std::stoul(text.substr(first, at - first)));
            };

            // How many DOMAINS the templates cover, counted from the
            // templates themselves rather than from memory. The README's
            // figure for this went stale once already, silently, because only
            // the template count was derived and this one was prose.
            std::set<malloy::scenario::ScenarioType> domains;
            for (const std::string& name : templates)
            {
                const ScenarioParseResult r = malloy::scenario::parse_scenario_file(
                    std::string(MALLOY_SCENARIO_DIR) + "/" + name);
                if (r.ok)
                {
                    domains.insert(r.scenario.type);
                }
            }

            const std::string root = std::string(MALLOY_SCENARIO_DIR) + "/..";
            struct Claim
            {
                std::string path;
                std::string marker;
                std::size_t actual;
            };
            const std::vector<Claim> claims = {
                {root + "/README.md", " templates ship across", templates.size()},
                {root + "/README.md", " domains, and", domains.size()},
#ifdef MALLOY_TEST_EXECUTABLE_COUNT
                {root + "/README.md", " test executables",
                 static_cast<std::size_t>(MALLOY_TEST_EXECUTABLE_COUNT)},
                {root + "/docs/00_START_HERE.md", " test executables",
                 static_cast<std::size_t>(MALLOY_TEST_EXECUTABLE_COUNT)},
#endif
            };

            // The milestone-status line is the other thing that goes stale, and
            // it has done so three times: corrected in the pre-M13 audit,
            // drifted again through M13 and M14, and drifted again through M15.
            // Correcting it a fourth time would just reset the clock, so it is
            // derived instead. CHANGELOG.md is the source of truth, because a
            // milestone is complete exactly when it has an entry there.
            {
                const auto digits_at = [](const std::string& text,
                                          std::size_t at) -> std::string {
                    std::size_t last = at;
                    while (last < text.size() &&
                           std::isdigit(static_cast<unsigned char>(text[last])) != 0)
                    {
                        ++last;
                    }
                    return (last == at) ? std::string() : text.substr(at, last - at);
                };
                const auto read_all = [](const std::string& path) {
                    std::ifstream in(path);
                    return std::string((std::istreambuf_iterator<char>(in)),
                                       std::istreambuf_iterator<char>());
                };

                const std::string changelog = read_all(root + "/CHANGELOG.md");
                MALLOY_CHECK_TRUE(!changelog.empty());

                std::size_t newest = 0;
                for (std::size_t at = changelog.find("## [M");
                     at != std::string::npos; at = changelog.find("## [M", at + 1))
                {
                    const std::string found = digits_at(changelog, at + 5);
                    if (!found.empty())
                    {
                        newest = std::max(
                            newest, static_cast<std::size_t>(std::stoul(found)));
                    }
                }
                MALLOY_CHECK_TRUE(newest >= 15);

                const char* status_docs[] = {"/CLAUDE.md", "/README.md",
                                             "/docs/00_START_HERE.md",
                                             "/docs/07_POST_M5_ROADMAP.md",
                                             "/docs/08_AI_HANDOFF_PROMPT.md"};
                for (const char* relative : status_docs)
                {
                    const std::string path = root + relative;
                    const std::string text = read_all(path);
                    MALLOY_CHECK_TRUE(!text.empty());

                    std::size_t stated = 0;
                    for (std::size_t at = text.find("M1-M"); at != std::string::npos;
                         at = text.find("M1-M", at + 1))
                    {
                        const std::string found = digits_at(text, at + 4);
                        if (found.empty())
                        {
                            continue;
                        }
                        const std::size_t claimed_end =
                            static_cast<std::size_t>(std::stoul(found));
                        // "M1-M5" is the locked original roadmap, a fixed
                        // historical range that is permanently correct. Only a
                        // range beyond it is a claim about the CURRENT state.
                        if (claimed_end <= 5)
                        {
                            continue;
                        }
                        ++stated;
                        if (claimed_end != newest)
                        {
                            std::cerr << path << " claims M1-M" << claimed_end
                                      << " but the newest milestone in CHANGELOG.md is M"
                                      << newest << '\n';
                            return 1;
                        }
                    }
                    // Each of these states the range at least once. One that
                    // stopped would otherwise pass this check vacuously.
                    if (stated == 0)
                    {
                        std::cerr << path
                                  << " no longer states a milestone range, so nothing "
                                     "here checks it" << '\n';
                        return 1;
                    }
                }

                // The ADR range is the same kind of claim, and it went
                // stale the moment ADR 0009 was written. The highest number on
                // disk is the truth. Checked in BOTH documents that state it:
                // fixing only the one that was guarded let the other drift
                // through two milestones before anyone noticed.
                {
                    std::size_t highest = 0;
                    for (const auto& entry : std::filesystem::directory_iterator(
                             root + "/docs/decisions"))
                    {
                        if (!entry.is_regular_file() ||
                            entry.path().extension() != ".md")
                        {
                            continue;
                        }
                        const std::string name = entry.path().filename().string();
                        const std::string leading = digits_at(name, 0);
                        if (!leading.empty())
                        {
                            highest = std::max(
                                highest,
                                static_cast<std::size_t>(std::stoul(leading)));
                        }
                    }
                    MALLOY_CHECK_TRUE(highest >= 9);

                    for (const char* document : {"/docs/00_START_HERE.md",
                                                 "/docs/08_AI_HANDOFF_PROMPT.md"})
                    {
                        const std::string text = read_all(root + document);
                        const std::string marker = "ADRs 0001-";
                        const std::size_t at = text.find(marker);
                        MALLOY_CHECK_TRUE(at != std::string::npos);
                        const std::string claimed_adr =
                            digits_at(text, at + marker.size());
                        MALLOY_CHECK_TRUE(!claimed_adr.empty());
                        if (static_cast<std::size_t>(std::stoul(claimed_adr)) !=
                            highest)
                        {
                            std::cerr << document << " claims ADRs up to "
                                      << claimed_adr << " but " << highest
                                      << " exist" << '\n';
                            return 1;
                        }
                    }
                }
            }

            for (const Claim& claim : claims)
            {
                const std::optional<std::size_t> claimed =
                    claimed_before(claim.path, claim.marker);
                if (!claimed)
                {
                    std::cerr << "no parseable count before \"" << claim.marker
                              << "\" in " << claim.path << '\n';
                    return 1;
                }
                if (*claimed != claim.actual)
                {
                    std::cerr << claim.path << " claims " << *claimed << " for \""
                              << claim.marker << "\" but " << claim.actual
                              << " are present" << '\n';
                    return 1;
                }
            }
        }

        for (const std::string& name : templates)
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
            // Rule 16 requires a template to have a test OR a documented
            // expected result. The parse above is the test half; this is the
            // documentation half, enforced rather than trusted to review.
            std::ifstream source(path);
            std::string text((std::istreambuf_iterator<char>(source)),
                             std::istreambuf_iterator<char>());
            if (text.find("# Expected") == std::string::npos)
            {
                std::cerr << "template " << name
                          << " has no documented expected result\n";
                return 1;
            }

            // Issue #14: prose is not enough. Every template must also carry at
            // least one machine-checkable `# check` line, compared against a
            // real run below. Documented figures have gone stale twice in this
            // repository and were corrected by hand both times.
            const std::vector<TemplateCheck> checks = parse_checks(text);
            if (checks.size() == 1 && checks.front().step < 0)
            {
                std::cerr << "template " << name << " has a malformed # check line\n";
                return 1;
            }
            if (checks.empty())
            {
                std::cerr << "template " << name << " has no # check line\n";
                return 1;
            }
            for (const TemplateCheck& check : checks)
            {
                if (check.step > r.scenario.steps)
                {
                    std::cerr << "template " << name << ": # check names step "
                              << check.step << ", past the run's " << r.scenario.steps
                              << '\n';
                    return 1;
                }
            }

            MALLOY_CHECK_TRUE(r.scenario.steps > 0);

            // Same dispatch the app performs, so a template is validated by the
            // domain it actually declares.
            //
            // Each world is run for the FULL documented number of steps, and
            // every check is evaluated at the step it names. A step that fails
            // validation reports InvalidState, so running to the end is also
            // the NaN check.
            if (r.scenario.type == malloy::scenario::ScenarioType::NBody)
            {
                MALLOY_CHECK_TRUE(r.scenario.bodies.size() >= 2);
                NBodyWorld world{r.scenario.simulation, r.scenario.nbody_settings,
                                 r.scenario.bodies};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto& settings = r.scenario.nbody_settings;
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.bodies();
                    if (q == "energy")
                        return malloy::nbody::total_energy(b, settings.g, settings.softening);
                    if (q == "kinetic") return malloy::nbody::total_kinetic_energy(b);
                    if (q == "momentum_x") return malloy::nbody::total_momentum(b).x;
                    if (q == "momentum_y") return malloy::nbody::total_momentum(b).y;
                    if (q == "angular") return malloy::nbody::total_angular_momentum(b);
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
            }
            else if (r.scenario.type == malloy::scenario::ScenarioType::NBody3D)
            {
                MALLOY_CHECK_TRUE(r.scenario.bodies3d.size() >= 2);
                malloy::nbody::NBody3DWorld world{r.scenario.simulation,
                                                  r.scenario.nbody_settings,
                                                  r.scenario.bodies3d};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto& settings = r.scenario.nbody_settings;
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.bodies();
                    if (q == "energy")
                        return malloy::nbody::total_energy(b, settings.g,
                                                           settings.softening);
                    if (q == "kinetic") return malloy::nbody::total_kinetic_energy(b);
                    if (q == "momentum_x") return malloy::nbody::total_momentum(b).x;
                    if (q == "momentum_y") return malloy::nbody::total_momentum(b).y;
                    if (q == "momentum_z") return malloy::nbody::total_momentum(b).z;
                    // Angular momentum is a vector here, so each component is
                    // nameable. A 2D template asking for one of these, or a 3D
                    // one asking for plain `angular`, fails rather than being
                    // skipped.
                    if (q == "angular_x")
                        return malloy::nbody::total_angular_momentum(b).x;
                    if (q == "angular_y")
                        return malloy::nbody::total_angular_momentum(b).y;
                    if (q == "angular_z")
                        return malloy::nbody::total_angular_momentum(b).z;
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
            }
            else if (r.scenario.type == malloy::scenario::ScenarioType::Rigid3D)
            {
                MALLOY_CHECK_TRUE(!r.scenario.rigid_bodies3d.empty());
                malloy::rigid::Rigid3DWorld world{r.scenario.simulation,
                                                  r.scenario.rigid_bodies3d,
                                                  r.scenario.rigid3d_settings};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.bodies();
                    if (q == "energy") return malloy::rigid::total_kinetic_energy3d(b);
                    if (q == "momentum_x")
                        return malloy::rigid::total_linear_momentum3d(b).x;
                    if (q == "momentum_y")
                        return malloy::rigid::total_linear_momentum3d(b).y;
                    if (q == "momentum_z")
                        return malloy::rigid::total_linear_momentum3d(b).z;
                    if (q == "angular_x")
                        return malloy::rigid::total_angular_momentum3d(b).x;
                    if (q == "angular_y")
                        return malloy::rigid::total_angular_momentum3d(b).y;
                    if (q == "angular_z")
                        return malloy::rigid::total_angular_momentum3d(b).z;
                    // The FIRST body's angular velocity, in its own body frame.
                    // A total would be meaningless here: the interesting claim
                    // in this domain is what one body's spin does, and summing
                    // body-frame vectors across differently-oriented bodies
                    // adds quantities that live in different frames.
                    if (b.empty()) return std::nullopt;
                    if (q == "spin_x") return b.front().angular_velocity.x;
                    if (q == "spin_y") return b.front().angular_velocity.y;
                    if (q == "spin_z") return b.front().angular_velocity.z;
                    if (q == "rotational")
                        return malloy::rigid::rotational_energy(b.front());
                    // The FIRST body's position and velocity, for a contact
                    // scenario where a single sphere's height and speed are the
                    // thing to watch.
                    if (q == "position_x") return b.front().position.x;
                    if (q == "position_y") return b.front().position.y;
                    if (q == "position_z") return b.front().position.z;
                    if (q == "velocity_x") return b.front().velocity.x;
                    if (q == "velocity_y") return b.front().velocity.y;
                    if (q == "velocity_z") return b.front().velocity.z;
                    // Gravitational potential and the total mechanical energy,
                    // which trade off during a fall and a bounce.
                    if (q == "potential")
                        return malloy::rigid::total_potential_energy3d(
                            b, world.settings().gravity);
                    if (q == "total_energy")
                        return malloy::rigid::total_kinetic_energy3d(b) +
                               malloy::rigid::total_potential_energy3d(
                                   b, world.settings().gravity);
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
            }
            else if (r.scenario.type == malloy::scenario::ScenarioType::Springs)
            {
                MALLOY_CHECK_TRUE(r.scenario.spring_bodies.size() >= 2);
                MALLOY_CHECK_TRUE(r.scenario.spring_network.size() >= 1);
                malloy::springs::SpringWorld world{r.scenario.simulation,
                                                   r.scenario.spring_network,
                                                   r.scenario.spring_bodies};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.bodies();
                    const Real kinetic = malloy::springs::total_kinetic_energy(b);
                    const Real elastic =
                        malloy::springs::total_elastic_energy(world.network(), b);
                    if (q == "energy") return kinetic + elastic;
                    if (q == "kinetic") return kinetic;
                    if (q == "elastic") return elastic;
                    if (q == "momentum_x") return malloy::springs::total_momentum(b).x;
                    if (q == "momentum_y") return malloy::springs::total_momentum(b).y;
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
            }
            else if (r.scenario.type == malloy::scenario::ScenarioType::Rigid)
            {
                MALLOY_CHECK_TRUE(r.scenario.rigid_bodies.size() >= 1);
                malloy::rigid::RigidWorld world{r.scenario.simulation,
                                                r.scenario.rigid_bodies,
                                                r.scenario.rigid_settings};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto gravity = r.scenario.rigid_settings.gravity;
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.bodies();
                    if (q == "energy") return malloy::rigid::total_energy(b, gravity);
                    if (q == "kinetic") return malloy::rigid::total_kinetic_energy(b);
                    if (q == "momentum_x") return malloy::rigid::total_linear_momentum(b).x;
                    if (q == "momentum_y") return malloy::rigid::total_linear_momentum(b).y;
                    if (q == "angular") return malloy::rigid::total_angular_momentum(b);
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
            }
            else if (r.scenario.type == malloy::scenario::ScenarioType::Charges)
            {
                MALLOY_CHECK_TRUE(r.scenario.charge_list.size() >= 1);
                malloy::charges::ChargeWorld world{r.scenario.simulation,
                                                   r.scenario.charge_settings,
                                                   r.scenario.charge_list};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto& settings = r.scenario.charge_settings;
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.particles();
                    if (q == "energy") return malloy::charges::total_energy(b, settings);
                    if (q == "kinetic") return malloy::charges::total_kinetic_energy(b);
                    if (q == "momentum_x") return malloy::charges::total_momentum(b).x;
                    if (q == "momentum_y") return malloy::charges::total_momentum(b).y;
                    if (q == "speed") return malloy::math::length(b.front().velocity);
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
            }
            else
            {
                MALLOY_CHECK_TRUE(r.scenario.particle_list.size() >= 2);
                malloy::particles::ParticleWorld world{r.scenario.simulation,
                                                      r.scenario.particle_settings,
                                                      r.scenario.particle_list};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto& settings = r.scenario.particle_settings;
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.particles();
                    if (q == "energy")
                        return malloy::particles::total_energy(b, settings.gravity);
                    if (q == "kinetic") return malloy::particles::total_kinetic_energy(b);
                    if (q == "momentum_x") return malloy::particles::total_momentum(b).x;
                    if (q == "momentum_y") return malloy::particles::total_momentum(b).y;
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
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

    // --- type rigid parses into the rigid fields. Every value is distinct so a
    //     swapped pair in the constructor changes the result. ---
    {
        std::istringstream in("type rigid\n"
                              "dt 0.01\n"
                              "steps 50\n"
                              // 11 fields since M14: radius came third.
                              "rigid_body 2.5 3.75 1.5 0.6 -0.4 1.0 2.0 0.7 "
                              "3.0 4.0 0.9\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Rigid);
        MALLOY_CHECK_EQ(r.scenario.rigid_bodies.size(), std::size_t{1});
        const auto& b = r.scenario.rigid_bodies[0];
        MALLOY_CHECK_NEAR(b.mass, 2.5, eps);
        MALLOY_CHECK_NEAR(b.inertia, 3.75, eps);
        MALLOY_CHECK_NEAR(b.radius, 1.5, eps);
        MALLOY_CHECK_NEAR(b.local_center_of_mass.x, 0.6, eps);
        MALLOY_CHECK_NEAR(b.local_center_of_mass.y, -0.4, eps);
        MALLOY_CHECK_NEAR(b.position.x, 1.0, eps);
        MALLOY_CHECK_NEAR(b.position.y, 2.0, eps);
        MALLOY_CHECK_NEAR(b.angle, 0.7, eps);
        MALLOY_CHECK_NEAR(b.velocity.x, 3.0, eps);
        MALLOY_CHECK_NEAR(b.velocity.y, 4.0, eps);
        MALLOY_CHECK_NEAR(b.angular_velocity, 0.9, eps);
    }

    // --- Domain isolation holds for the third domain too. ---
    {
        std::istringstream in("type rigid\nbody 1.0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type rigid\nparticle 1.0 0.5 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type nbody\nrigid_body 1 1 0 0 0 0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type particles\nrigid_body 1 1 0 0 0 0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // One field short.
        std::istringstream in("type rigid\nrigid_body 1 1 0 0 0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }

    // --- type springs parses bodies and topology, and declaration order is
    //     preserved because it decides accumulation order. ---
    {
        std::istringstream in("type springs\n"
                              "dt 0.01\n"
                              "steps 50\n"
                              "spring_body 1.5 1.0 2.0 3.0 4.0\n"
                              "spring_body 2.5 5.0 6.0 7.0 8.0\n"
                              "spring_body 0.5 9.0 1.5 2.5 3.5\n"
                              "spring 0 1 1.25 30.0 0.5\n"
                              "spring 1 2 2.75 12.0 0.25\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Springs);
        MALLOY_CHECK_EQ(r.scenario.spring_bodies.size(), std::size_t{3});
        MALLOY_CHECK_EQ(r.scenario.spring_network.size(), std::size_t{2});

        const auto& b = r.scenario.spring_bodies[0];
        MALLOY_CHECK_NEAR(b.mass, 1.5, eps);
        MALLOY_CHECK_NEAR(b.position.x, 1.0, eps);
        MALLOY_CHECK_NEAR(b.position.y, 2.0, eps);
        MALLOY_CHECK_NEAR(b.velocity.x, 3.0, eps);
        MALLOY_CHECK_NEAR(b.velocity.y, 4.0, eps);

        // Declared order preserved, and every field distinct.
        const auto& s0 = r.scenario.spring_network.springs()[0];
        MALLOY_CHECK_EQ(s0.a, std::size_t{0});
        MALLOY_CHECK_EQ(s0.b, std::size_t{1});
        MALLOY_CHECK_NEAR(s0.rest_length, 1.25, eps);
        MALLOY_CHECK_NEAR(s0.stiffness, 30.0, eps);
        MALLOY_CHECK_NEAR(s0.damping, 0.5, eps);
        const auto& s1 = r.scenario.spring_network.springs()[1];
        MALLOY_CHECK_NEAR(s1.rest_length, 2.75, eps);
        MALLOY_CHECK_NEAR(s1.stiffness, 12.0, eps);
    }

    // --- Domain isolation holds for the fourth domain. ---
    {
        std::istringstream in("type springs\nbody 1.0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type nbody\nspring 0 1 1 1 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type particles\nspring_body 1 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type rigid\nspring 0 1 1 1 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type springs\nspring 0 1 1 1\n"); // one short
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }

    // --- M14: the rigid format carries a radius, and rigid_static declares an
    //     immovable body without needing the token "inf". ---
    {
        std::istringstream in("type rigid\n"
                              "restitution 0.25\n"
                              "rigid_body 2.5 3.75 0.8 0.6 -0.4 1.0 2.0 0.7 3.0 4.0 0.9\n"
                              "rigid_static 1.5 9.0 8.0 0.3\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_EQ(r.scenario.rigid_bodies.size(), std::size_t{2});
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.restitution, 0.25, eps);

        const auto& moving = r.scenario.rigid_bodies[0];
        MALLOY_CHECK_NEAR(moving.radius, 0.8, eps);
        MALLOY_CHECK_NEAR(moving.mass, 2.5, eps);
        MALLOY_CHECK_FALSE(moving.is_static());

        const auto& wall = r.scenario.rigid_bodies[1];
        MALLOY_CHECK_TRUE(wall.is_static());
        MALLOY_CHECK_NEAR(wall.inverse_mass(), 0.0, 0.0);
        MALLOY_CHECK_NEAR(wall.radius, 1.5, eps);
        MALLOY_CHECK_NEAR(wall.position.x, 9.0, eps);
        MALLOY_CHECK_NEAR(wall.angle, 0.3, eps);
    }
    {
        std::istringstream in("type particles\nrigid_static 1.0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // restitution now serves two domains, but not a third.
        std::istringstream in("type springs\nrestitution 0.5\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }

    // --- M15: gravity reaches the rigid domain. The two components differ in
    //     sign and magnitude, so a transposed pair would show, and the particle
    //     field is checked to stay at zero, so a branch routing the value to
    //     the wrong domain would show too. ---
    {
        std::istringstream in("type rigid\n"
                              "gravity 1.5 -9.81\n"
                              "rigid_body 1 1 0.5 0 0 0 0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.gravity.x, 1.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.gravity.y, -9.81, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.gravity.x, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.gravity.y, 0.0, 0.0);
    }
    {
        // The same key still serves the particle domain it was written for.
        std::istringstream in("type particles\n"
                              "gravity -2.0 0.25\n"
                              "particle 1.0 0.5 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.gravity.x, -2.0, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.gravity.y, 0.25, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.gravity.x, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.gravity.y, 0.0, 0.0);
    }
    {
        // Two components are required, not one.
        std::istringstream in("type rigid\ngravity 1.0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // gravity now serves two domains, but not a third.
        std::istringstream in("type springs\ngravity 0.0 -9.81\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }

    // --- M16: the ground key, normalized on load so the geometry type keeps a
    //     strict unit-normal invariant. The direction below has length 5, and
    //     both components are checked, so a dropped one would show. ---
    {
        std::istringstream in("type rigid\n"
                              "ground 3 4 -2.5\n"
                              "ground 0 -1 -7\n"
                              "rigid_body 1 1 0.5 0 0 0 0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_EQ(r.scenario.rigid_settings.ground.size(), std::size_t{2});

        const auto& first = r.scenario.rigid_settings.ground[0];
        MALLOY_CHECK_NEAR(first.normal.x, 0.6, eps);
        MALLOY_CHECK_NEAR(first.normal.y, 0.8, eps);
        MALLOY_CHECK_NEAR(first.offset, -2.5, eps);
        MALLOY_CHECK_TRUE(first.is_valid());

        // Order is preserved, and an already-unit normal is unchanged.
        const auto& second = r.scenario.rigid_settings.ground[1];
        MALLOY_CHECK_NEAR(second.normal.x, 0.0, eps);
        MALLOY_CHECK_NEAR(second.normal.y, -1.0, eps);
        MALLOY_CHECK_NEAR(second.offset, -7.0, eps);

        MALLOY_CHECK_TRUE(r.scenario.rigid_settings.is_valid());
    }
    {
        // A direction with no direction is refused rather than normalized into
        // one, which would silently invent a floor orientation.
        std::istringstream in("type rigid\nground 0 0 1\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // Three fields are required.
        std::istringstream in("type rigid\nground 0 1\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // ground belongs to the rigid domain and to no other.
        std::istringstream in("type particles\nground 0 1 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type springs\nground 0 1 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // Absent by default, so every rigid scenario written before M16 loads
        // with no ground at all and behaves exactly as it did.
        std::istringstream in("type rigid\nrigid_body 1 1 0.5 0 0 0 0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.rigid_settings.ground.empty());
    }

    // --- M17: the friction key. Above 1 is accepted on purpose; negative is
    //     rejected by the settings rather than by the parser, which is where
    //     the rule lives. ---
    {
        std::istringstream in("type rigid\n"
                              "friction 1.35\n"
                              "rigid_body 1 1 0.5 0 0 0 0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.friction, 1.35, eps);
        MALLOY_CHECK_TRUE(r.scenario.rigid_settings.is_valid());
    }
    {
        // Absent means frictionless, so every rigid scenario written before
        // M17 loads and behaves exactly as it did.
        std::istringstream in("type rigid\nrigid_body 1 1 0.5 0 0 0 0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.friction, 0.0, 0.0);
    }
    {
        // Parsed, but refused by validation, which is reported as a bad world
        // rather than a bad file.
        std::istringstream in("type rigid\nfriction -0.5\nrigid_body 1 1 0.5 0 0 0 0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_FALSE(r.scenario.rigid_settings.is_valid());
    }
    {
        std::istringstream in("type rigid\nfriction\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // friction belongs to the rigid domain and to no other.
        std::istringstream in("type particles\nfriction 0.5\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type springs\nfriction 0.5\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }

    // --- M18: type charges parses into the charge fields. Every value is
    //     distinct and the charges differ in SIGN, so a dropped or transposed
    //     field shows. ---
    {
        std::istringstream in("type charges\n"
                              "coulomb 2.5\n"
                              "efield 0.5 -1.5\n"
                              "bfield 3.25\n"
                              "softening 0.125\n"
                              "charge 2.0  1.5  3.0 -4.0  0.25 -0.75\n"
                              "charge 0.5 -2.5 -6.0  7.0 -0.5   0.125\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Charges);

        MALLOY_CHECK_NEAR(r.scenario.charge_settings.k, 2.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.charge_settings.electric.x, 0.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.charge_settings.electric.y, -1.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.charge_settings.magnetic, 3.25, eps);
        MALLOY_CHECK_NEAR(r.scenario.charge_settings.softening, 0.125, eps);
        MALLOY_CHECK_TRUE(r.scenario.charge_settings.is_valid());

        MALLOY_CHECK_EQ(r.scenario.charge_list.size(), std::size_t{2});
        const auto& first = r.scenario.charge_list[0];
        MALLOY_CHECK_NEAR(first.mass, 2.0, eps);
        MALLOY_CHECK_NEAR(first.charge, 1.5, eps);
        MALLOY_CHECK_NEAR(first.position.x, 3.0, eps);
        MALLOY_CHECK_NEAR(first.position.y, -4.0, eps);
        MALLOY_CHECK_NEAR(first.velocity.x, 0.25, eps);
        MALLOY_CHECK_NEAR(first.velocity.y, -0.75, eps);
        MALLOY_CHECK_TRUE(first.is_valid());

        // The second charge is NEGATIVE, which no other domain's body can be.
        MALLOY_CHECK_NEAR(r.scenario.charge_list[1].charge, -2.5, eps);
        MALLOY_CHECK_TRUE(r.scenario.charge_list[1].is_valid());

        // And nothing leaked into another domain's settings.
        MALLOY_CHECK_NEAR(r.scenario.nbody_settings.softening, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.gravity.x, 0.0, 0.0);
    }
    {
        // softening now serves two domains, and routes to the right one.
        std::istringstream in("type nbody\nsoftening 0.75\nbody 1 0 0 0 0\nbody 1 1 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.nbody_settings.softening, 0.75, eps);
        MALLOY_CHECK_NEAR(r.scenario.charge_settings.softening, 0.0, 0.0);
    }
    {
        // A charge of exactly zero is legal: a neutral particle carried by the
        // fields of others.
        std::istringstream in("type charges\ncharge 1.0 0.0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.charge_list[0].charge, 0.0, 0.0);
    }
    {
        // Six fields are required, not five.
        std::istringstream in("type charges\ncharge 1.0 1.0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type charges\nefield 1.0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // Each key belongs to the charge domain and to no other.
        std::istringstream in("type rigid\nbfield 1.0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type particles\ncharge 1 1 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type springs\ncoulomb 1.0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // And a key belonging to another domain is refused inside charges,
        // which is what makes a typo in `type` surface immediately.
        std::istringstream in("type charges\nrestitution 0.5\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type charges\nground 0 1 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }

    // --- M19: type nbody3d parses into the 3D fields. Every value is distinct
    //     and no component is zero, so a dropped or transposed one shows. The
    //     z components in particular would be invisible in a configuration
    //     that left them at zero. ---
    {
        std::istringstream in("type nbody3d\n"
                              "g 2.5\n"
                              "softening 0.125\n"
                              "body3 3.5   1.0 -2.0  4.0   -0.5  0.25 -0.75\n"
                              "body3 0.25 -6.0  7.5 -8.5    1.5 -2.5   3.5\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::NBody3D);

        // g and softening are shared with the 2D gravity world, because they
        // mean the same thing in either dimension.
        MALLOY_CHECK_NEAR(r.scenario.nbody_settings.g, 2.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.nbody_settings.softening, 0.125, eps);

        MALLOY_CHECK_EQ(r.scenario.bodies3d.size(), std::size_t{2});
        const auto& first = r.scenario.bodies3d[0];
        MALLOY_CHECK_NEAR(first.mass, 3.5, eps);
        MALLOY_CHECK_NEAR(first.position.x, 1.0, eps);
        MALLOY_CHECK_NEAR(first.position.y, -2.0, eps);
        MALLOY_CHECK_NEAR(first.position.z, 4.0, eps);
        MALLOY_CHECK_NEAR(first.velocity.x, -0.5, eps);
        MALLOY_CHECK_NEAR(first.velocity.y, 0.25, eps);
        MALLOY_CHECK_NEAR(first.velocity.z, -0.75, eps);
        MALLOY_CHECK_TRUE(first.is_valid());

        const auto& second = r.scenario.bodies3d[1];
        MALLOY_CHECK_NEAR(second.position.z, -8.5, eps);
        MALLOY_CHECK_NEAR(second.velocity.z, 3.5, eps);

        // Nothing leaked into another domain's state.
        MALLOY_CHECK_TRUE(r.scenario.bodies.empty());
        MALLOY_CHECK_NEAR(r.scenario.charge_settings.softening, 0.0, 0.0);
    }
    {
        // Seven fields are required, not six: a 2D body line is not a 3D one
        // with a component missing.
        std::istringstream in("type nbody3d\nbody3 1.0  0 0 0  0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // And the 2D key is not accepted here, so a scenario cannot silently
        // half-convert.
        std::istringstream in("type nbody3d\nbody 1.0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // body3 belongs to this domain and to no other.
        std::istringstream in("type nbody\nbody3 1.0 0 0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type charges\nbody3 1.0 0 0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // A key belonging to another domain is refused inside nbody3d, which is
        // what makes a typo in `type` surface immediately rather than running
        // the wrong simulation.
        std::istringstream in("type nbody3d\nrestitution 0.5\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type nbody3d\nbfield 1.0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type nbody3d\nground 0 1 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    // --- M20/M22: type rigid3d parses into the 3D rigid fields. Eighteen values,
    //     every one of them distinct and none of them zero except where zero is
    //     the thing being tested, so a dropped or transposed field shows. ---
    {
        std::istringstream in(
            "type rigid3d\n"
            "rigid_body3d 2.5  1.5 2.5 3.5   0.5   1.0 -2.0 4.0   0.0 0.0 1.0 "
            "1.5707963267948966   -0.5 0.25 -0.75   0.125 -0.375 0.625\n"
            "rigid_body3d 0.75  1.0 1.0 1.0   0.25   -6.0 7.5 -8.5   0.0 0.0 0.0 0.0   "
            "1.5 -2.5 3.5   -1.25 2.25 -3.25\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Rigid3D);
        MALLOY_CHECK_EQ(r.scenario.rigid_bodies3d.size(), std::size_t{2});

        const auto& first = r.scenario.rigid_bodies3d[0];
        MALLOY_CHECK_NEAR(first.mass, 2.5, eps);
        MALLOY_CHECK_NEAR(first.inertia.x, 1.5, eps);
        MALLOY_CHECK_NEAR(first.inertia.y, 2.5, eps);
        MALLOY_CHECK_NEAR(first.inertia.z, 3.5, eps);
        MALLOY_CHECK_NEAR(first.radius, 0.5, eps); // M22: the collision radius
        MALLOY_CHECK_NEAR(first.position.x, 1.0, eps);
        MALLOY_CHECK_NEAR(first.position.y, -2.0, eps);
        MALLOY_CHECK_NEAR(first.position.z, 4.0, eps);
        MALLOY_CHECK_NEAR(first.velocity.x, -0.5, eps);
        MALLOY_CHECK_NEAR(first.velocity.y, 0.25, eps);
        MALLOY_CHECK_NEAR(first.velocity.z, -0.75, eps);
        MALLOY_CHECK_NEAR(first.angular_velocity.x, 0.125, eps);
        MALLOY_CHECK_NEAR(first.angular_velocity.y, -0.375, eps);
        MALLOY_CHECK_NEAR(first.angular_velocity.z, 0.625, eps);

        // The orientation was given as an axis and an angle, so what lands in
        // the body is the quarter turn about z that they name, and it is a unit
        // quaternion because from_axis_angle cannot produce anything else.
        MALLOY_CHECK_TRUE(malloy::math::is_unit(first.orientation));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::rotate(first.orientation, malloy::math::Vec3{1.0, 0.0, 0.0}),
            malloy::math::Vec3{0.0, 1.0, 0.0}, 1e-15));
        MALLOY_CHECK_TRUE(first.is_valid());

        // A zero axis is the identity, whatever the angle says, which is how a
        // template writes an unrotated body without having to know quaternions.
        const auto& second = r.scenario.rigid_bodies3d[1];
        MALLOY_CHECK_NEAR(second.radius, 0.25, eps);
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(second.orientation,
                                                     malloy::math::Quat{}, 0.0));
        MALLOY_CHECK_TRUE(second.is_valid());

        // Nothing leaked into the 2D rigid state or any other domain.
        MALLOY_CHECK_TRUE(r.scenario.rigid_bodies.empty());
        MALLOY_CHECK_TRUE(r.scenario.bodies3d.empty());
        MALLOY_CHECK_TRUE(r.scenario.bodies.empty());
    }
    {
        // Eighteen fields are required. Seventeen is not a body with a
        // default somewhere, it is an error.
        std::istringstream in("type rigid3d\n"
                              "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  1 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // Nineteen is an error too: a trailing token is refused rather than
        // ignored, so a line written for some later format does not silently
        // run as this one.
        std::istringstream in(
            "type rigid3d\n"
            "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  1 0 0  9\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // Not a number where a number belongs.
        std::istringstream in(
            "type rigid3d\n"
            "rigid_body3d 1  1 2 three  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // The 2D rigid body line is not accepted here, and the 3D one is not
        // accepted there, so a scenario cannot half-convert.
        std::istringstream in("type rigid3d\nrigid_body 1 1 0.5 0 0 0 0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type rigid\n"
                              "rigid_body3d 1  1 2 3  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type nbody3d\n"
                              "rigid_body3d 1  1 2 3  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // `restitution` and `friction` are accepted for rigid3d (M22, M23),
        // but the 2D `gravity` and `ground` keys are not: 3D uses `gravity3`
        // and `plane3`, so a 2D key here is a mistake, not a silent
        // half-conversion. bfield and spring belong to other domains entirely.
        for (const char* line : {"gravity 0 -9.81\n", "ground 0 1 -2\n",
                                 "bfield 1.0\n", "spring 0 1 1 1 0\n"})
        {
            std::istringstream in(std::string("type rigid3d\n") + line);
            MALLOY_CHECK_FALSE(parse_scenario(in).ok);
        }
    }
    {
        // A body with a non-positive principal moment parses, because the
        // parser checks syntax only, and is then refused by the domain. Each
        // rule lives in exactly one place.
        std::istringstream in("type rigid3d\n"
                              "rigid_body3d 1  1 0 3  0.5  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_FALSE(r.scenario.rigid_bodies3d[0].is_valid());
        malloy::rigid::Rigid3DWorld world{r.scenario.simulation,
                                          r.scenario.rigid_bodies3d};
        MALLOY_CHECK_TRUE(world.validate() == StepStatus::InvalidState);
    }
    // --- M21: the torque key, a world-frame setting for type rigid3d. ---
    {
        std::istringstream in("type rigid3d\n"
                              "torque 0.5 -1.5 2.5\n"
                              "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Rigid3D);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.torque.x, 0.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.torque.y, -1.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.torque.z, 2.5, eps);
        MALLOY_CHECK_TRUE(r.scenario.rigid3d_settings.is_valid());

        std::istringstream none("type rigid3d\n"
                                "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        const ScenarioParseResult r2 = parse_scenario(none);
        MALLOY_CHECK_TRUE(r2.ok);
        MALLOY_CHECK_NEAR(r2.scenario.rigid3d_settings.torque.x, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r2.scenario.rigid3d_settings.torque.y, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r2.scenario.rigid3d_settings.torque.z, 0.0, 0.0);
    }
    {
        std::istringstream in("type rigid3d\ntorque 1 2\n"); // three fields, not two
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type rigid3d\ntorque 1 2 3 4\n"); // trailing token refused
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        std::istringstream in("type rigid3d\ntorque 1 two 3\n"); // not a number
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // torque belongs to rigid3d and to no other domain.
        for (const char* head : {"type rigid\n", "type nbody3d\n", "type springs\n",
                                 "type charges\n"})
        {
            std::istringstream in(std::string(head) + "torque 1 2 3\n");
            MALLOY_CHECK_FALSE(parse_scenario(in).ok);
        }
    }
    {
        // A torque whose square overflows is refused by the settings.
        std::istringstream in("type rigid3d\n"
                              "torque 1e200 0 0\n"
                              "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_FALSE(r.scenario.rigid3d_settings.is_valid());
        malloy::rigid::Rigid3DWorld world{r.scenario.simulation,
                                          r.scenario.rigid_bodies3d,
                                          r.scenario.rigid3d_settings};
        MALLOY_CHECK_TRUE(world.validate() == StepStatus::InvalidSettings);
    }
    // --- M22: the contact keys for type rigid3d: gravity3, plane3, and
    //     restitution. Distinct nonzero values, so a dropped axis shows. ---
    {
        std::istringstream in(
            "type rigid3d\n"
            "gravity3 0.5 -1.5 -9.81\n"
            "restitution 0.8\n"
            "plane3 0 0 2 -3\n"                       // non-unit normal, normalized on load
            "plane3 1 0 0 0\n"
            "rigid_body3d 1  1 2 3  0.5  0 0 5  0 0 1 0  0 0 0  0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Rigid3D);

        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.gravity.x, 0.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.gravity.y, -1.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.gravity.z, -9.81, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.restitution, 0.8, eps);

        MALLOY_CHECK_EQ(r.scenario.rigid3d_settings.ground.size(), std::size_t{2});
        // The first plane's (0,0,2) normal was normalized to (0,0,1), and the
        // offset carried through unchanged.
        const auto& p0 = r.scenario.rigid3d_settings.ground[0];
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(p0.normal,
                                                     malloy::math::Vec3{0.0, 0.0, 1.0}, eps));
        MALLOY_CHECK_NEAR(p0.offset, -3.0, eps);
        MALLOY_CHECK_TRUE(p0.is_valid()); // unit normal after normalization
        MALLOY_CHECK_TRUE(r.scenario.rigid3d_settings.is_valid());

        // The whole thing validates and steps as a world.
        malloy::rigid::Rigid3DWorld world{r.scenario.simulation,
                                          r.scenario.rigid_bodies3d,
                                          r.scenario.rigid3d_settings};
        MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
    }
    {
        // gravity3 needs three components, plane3 needs four, and a trailing
        // token is refused.
        for (const char* line : {"gravity3 0 0\n", "gravity3 0 0 0 0\n",
                                 "plane3 0 0 1\n", "plane3 0 0 1 0 0\n",
                                 "plane3 0 0 z 0\n"})
        {
            std::istringstream in(std::string("type rigid3d\n") + line);
            MALLOY_CHECK_FALSE(parse_scenario(in).ok);
        }
    }
    {
        // A plane3 whose normal cannot be normalized is refused at the boundary,
        // so collide::Plane3 keeps its strict unit-normal invariant.
        std::istringstream in("type rigid3d\nplane3 0 0 0 5\n");
        MALLOY_CHECK_FALSE(parse_scenario(in).ok);
    }
    {
        // gravity3 and plane3 belong to rigid3d and to no other domain, and the
        // 2D gravity and ground keys do not belong to rigid3d.
        for (const char* line : {"gravity3 0 0 -9.81\n", "plane3 0 0 1 0\n"})
        {
            for (const char* head : {"type rigid\n", "type nbody3d\n",
                                     "type particles\n"})
            {
                std::istringstream in(std::string(head) + line);
                MALLOY_CHECK_FALSE(parse_scenario(in).ok);
            }
        }
    }
    {
        // restitution is now accepted for rigid3d and writes only that domain's
        // field, not the 2D rigid or particle one.
        std::istringstream in("type rigid3d\nrestitution 0.3\n"
                              "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.restitution, 0.3, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.restitution, 1.0, eps); // untouched default
    }
    {
        // friction (M23) is accepted for rigid3d and writes only that domain's
        // field, not the 2D rigid one. It requires a value and is not capped.
        std::istringstream in("type rigid3d\nfriction 0.4\n"
                              "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.friction, 0.4, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid_settings.friction, 0.0, eps); // untouched default
        MALLOY_CHECK_TRUE(r.scenario.rigid3d_settings.is_valid());

        std::istringstream missing("type rigid3d\nfriction\n");
        MALLOY_CHECK_FALSE(parse_scenario(missing).ok);
    }
    {
        // An unknown type is still an error, and the message names rigid3d
        // among the ones it could have meant (M22 keys included via rigid3d).
        std::istringstream in("type rigid4d\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_FALSE(r.ok);
        MALLOY_CHECK_TRUE(r.error.find("rigid4d") != std::string::npos);
    }
    {
        // A softening whose square overflows is refused by the settings, the
        // same bound the 2D world has, since both use NBodySettings.
        std::istringstream in("type nbody3d\nsoftening 1e200\nbody3 1 0 0 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_FALSE(r.scenario.nbody_settings.is_valid());
    }

    std::cout << "malloy_scenario_tests passed\n";
    return 0;
}
