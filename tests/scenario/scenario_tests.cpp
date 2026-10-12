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
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <cctype>
#include <iterator>
#include <iostream>
#include <locale>
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
// True when every byte is printable ASCII, 0x20 to 0x7E: no control byte, no
// DEL, nothing from 0x80 up. A parse error must be safe to print to a terminal
// whatever the file contained (issue #28).
bool printable_ascii(const std::string& text)
{
    return std::all_of(text.begin(), text.end(), [](char c) {
        const unsigned char u = static_cast<unsigned char>(c);
        return u >= 0x20 && u < 0x7F;
    });
}

// How a parse error should show one byte of a quoted token: the byte itself
// when it is printable ASCII other than the backslash and the quote, otherwise
// \xNN in lowercase hex. Written with snprintf rather than the parser's own
// table, so the expectation and the implementation cannot share a mistake.
std::string rendered(unsigned char byte)
{
    if (byte >= 0x20 && byte < 0x7F && byte != '\\' && byte != '\'')
    {
        return std::string(1, static_cast<char>(byte));
    }
    char escape[5]{};
    std::snprintf(escape, sizeof escape, "\\x%02x", static_cast<unsigned>(byte));
    return escape;
}

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
                    // The potential energy, of gravity and of the applied force
                    // (M28, issue #47), and the total mechanical energy, which
                    // trade off during a fall, a push and a bounce. The total
                    // is what the app's E column prints.
                    const Real potential =
                        malloy::rigid::total_potential_energy3d(b, world.settings().gravity) +
                        malloy::rigid::total_force_potential3d(b, world.settings().force);
                    if (q == "potential") return potential;
                    if (q == "total_energy")
                        return malloy::rigid::total_kinetic_energy3d(b) + potential;
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
            else if (r.scenario.type == malloy::scenario::ScenarioType::Springs3D)
            {
                MALLOY_CHECK_TRUE(r.scenario.spring_bodies3d.size() >= 2);
                MALLOY_CHECK_TRUE(r.scenario.spring_network3d.size() >= 1);
                malloy::springs::SpringWorld3D world{r.scenario.simulation,
                                                     r.scenario.spring_network3d,
                                                     r.scenario.spring_bodies3d};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.bodies();
                    const Real kinetic = malloy::springs::total_kinetic_energy3d(b);
                    const Real elastic =
                        malloy::springs::total_elastic_energy3d(world.network(), b);
                    if (q == "energy") return kinetic + elastic;
                    if (q == "kinetic") return kinetic;
                    if (q == "elastic") return elastic;
                    if (q == "momentum_x") return malloy::springs::total_momentum3d(b).x;
                    if (q == "momentum_y") return malloy::springs::total_momentum3d(b).y;
                    if (q == "momentum_z") return malloy::springs::total_momentum3d(b).z;
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
            }
            else if (r.scenario.type == malloy::scenario::ScenarioType::Charges3D)
            {
                MALLOY_CHECK_TRUE(r.scenario.charge_list3d.size() >= 1);
                malloy::charges::Charge3DWorld world{r.scenario.simulation,
                                                     r.scenario.charge3d_settings,
                                                     r.scenario.charge_list3d};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto& settings = r.scenario.charge3d_settings;
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.particles();
                    if (q == "energy")
                        return malloy::charges::total_energy3d(b, settings);
                    if (q == "kinetic")
                        return malloy::charges::total_kinetic_energy3d(b);
                    if (q == "momentum_x") return malloy::charges::total_momentum3d(b).x;
                    if (q == "momentum_y") return malloy::charges::total_momentum3d(b).y;
                    if (q == "momentum_z") return malloy::charges::total_momentum3d(b).z;
                    if (q == "speed") return malloy::math::length(b.front().velocity);
                    return std::nullopt;
                };
                if (!run_checked(world, name, r.scenario.steps, checks, value_of))
                {
                    return 1;
                }
            }
            else if (r.scenario.type == malloy::scenario::ScenarioType::Particles3D)
            {
                MALLOY_CHECK_TRUE(r.scenario.particle_list3d.size() >= 2);
                malloy::particles::ParticleWorld3D world{r.scenario.simulation,
                                                         r.scenario.particle3d_settings,
                                                         r.scenario.particle_list3d};
                MALLOY_CHECK_TRUE(world.validate() == StepStatus::Ok);
                const auto& settings = r.scenario.particle3d_settings;
                const auto value_of = [&](const std::string& q) -> std::optional<Real> {
                    const auto& b = world.particles();
                    if (q == "energy")
                        return malloy::particles::total_energy3d(b, settings.gravity);
                    if (q == "kinetic") return malloy::particles::total_kinetic_energy3d(b);
                    if (q == "momentum_x") return malloy::particles::total_momentum3d(b).x;
                    if (q == "momentum_y") return malloy::particles::total_momentum3d(b).y;
                    if (q == "momentum_z") return malloy::particles::total_momentum3d(b).z;
                    // M38: the gas quantities. Pressure is empty before any
                    // time has passed, so a `step 0 pressure` check FAILS.
                    if (q == "temperature") return malloy::particles::temperature3d(b);
                    if (q == "speed_ratio") return malloy::particles::speed_ratio3d(b);
                    if (q == "slow_fraction") return malloy::particles::slow_fraction3d(b);
                    if (q == "pressure")
                        return malloy::particles::pressure3d(
                            world.wall_impulse(), settings.bounds, world.elapsed_time());
                    // M39: the heat each wall has delivered to the gas.
                    const auto& heat = world.wall_heat();
                    if (q == "heat_min_x") return heat.min_x;
                    if (q == "heat_max_x") return heat.max_x;
                    if (q == "heat_min_y") return heat.min_y;
                    if (q == "heat_max_y") return heat.max_y;
                    if (q == "heat_min_z") return heat.min_z;
                    if (q == "heat_max_z") return heat.max_z;
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
                    if (q == "temperature") return malloy::particles::temperature(b);
                    if (q == "speed_ratio") return malloy::particles::speed_ratio(b);
                    if (q == "slow_fraction") return malloy::particles::slow_fraction(b);
                    if (q == "pressure")
                        return malloy::particles::pressure(
                            world.wall_impulse(), settings.bounds, world.elapsed_time());
                    // M39: the heat each wall has delivered to the gas. A 2D
                    // world has no z walls, so heat_min_z FAILS here.
                    const auto& heat = world.wall_heat();
                    if (q == "heat_min_x") return heat.min_x;
                    if (q == "heat_max_x") return heat.max_x;
                    if (q == "heat_min_y") return heat.min_y;
                    if (q == "heat_max_y") return heat.max_y;
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

    // --- M32: a particles3d scenario parses into the 3D particle fields. Every
    //     component distinct, so a dropped or swapped axis in the loader shows. ---
    {
        std::istringstream in("type particles3d\n"
                              "restitution 0.25\n"
                              "gravity3 0.1 0.2 -0.3\n"
                              "bounds3 -2 -3 -4 5 6 7\n"
                              "particle3 1.5 0.4  1.0 2.0 3.0  -4.0 -5.0 -6.0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Particles3D);
        MALLOY_CHECK_NEAR(r.scenario.particle3d_settings.restitution, 0.25, eps);
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            r.scenario.particle3d_settings.gravity, malloy::math::Vec3{0.1, 0.2, -0.3}, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            r.scenario.particle3d_settings.bounds.min, malloy::math::Vec3{-2.0, -3.0, -4.0}, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            r.scenario.particle3d_settings.bounds.max, malloy::math::Vec3{5.0, 6.0, 7.0}, eps));
        MALLOY_CHECK_EQ(r.scenario.particle_list3d.size(), std::size_t{1});
        MALLOY_CHECK_NEAR(r.scenario.particle_list3d[0].mass, 1.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.particle_list3d[0].radius, 0.4, eps);
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            r.scenario.particle_list3d[0].position, malloy::math::Vec3{1.0, 2.0, 3.0}, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            r.scenario.particle_list3d[0].velocity, malloy::math::Vec3{-4.0, -5.0, -6.0}, eps));
    }

    // --- M39: thermal walls and the seed. Every wall name lands in its own
    //     field (distinct temperatures, so a swapped pair shows), in 2D and in
    //     3D; the seed takes plain decimal digits across the whole 64-bit range
    //     and nothing else. ---
    {
        std::istringstream in("type particles\n"
                              "wall_temperature min_x 1.5\n"
                              "wall_temperature max_x 2.5\n"
                              "wall_temperature min_y 3.5\n"
                              "wall_temperature max_y 4.5\n"
                              "seed 18446744073709551615\n"
                              "particle 1 0.1 0 0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        const auto& w = r.scenario.particle_settings.wall_temperature;
        MALLOY_CHECK_TRUE(w.min_x == 1.5 && w.max_x == 2.5 && w.min_y == 3.5 && w.max_y == 4.5);
        MALLOY_CHECK_TRUE(r.scenario.particle_settings.seed == 18446744073709551615ULL);
        // 3D writes the 3D settings, z walls included, and leaves 2D's alone.
        std::istringstream in3("type particles3d\n"
                               "wall_temperature max_z 6.5\n"
                               "wall_temperature min_z 5.5\n"
                               "wall_temperature min_x 1.5\n"
                               "seed 0\n"
                               "seed 42\n"); // a scalar key: the last one wins
        const ScenarioParseResult r3 = parse_scenario(in3);
        MALLOY_CHECK_TRUE(r3.ok);
        const auto& w3 = r3.scenario.particle3d_settings.wall_temperature;
        MALLOY_CHECK_TRUE(w3.min_x == 1.5 && w3.max_x == 0.0 && w3.min_z == 5.5 &&
                          w3.max_z == 6.5);
        MALLOY_CHECK_TRUE(r3.scenario.particle3d_settings.seed == 42);
        MALLOY_CHECK_TRUE(r3.scenario.particle_settings.seed == 0);
        MALLOY_CHECK_TRUE(r3.scenario.particle_settings.wall_temperature.min_x == 0.0);
        // A negative temperature parses, and validate() refuses it, so each rule
        // lives in one place.
        std::istringstream cold("type particles\nwall_temperature max_y -1\n"
                                "particle 1 0.1 0 0 0 0\n");
        const ScenarioParseResult c = parse_scenario(cold);
        MALLOY_CHECK_TRUE(c.ok);
        MALLOY_CHECK_FALSE(c.scenario.particle_settings.is_valid());
    }
    {
        const char* const bad[] = {
            "type nbody\nseed 1\n",                              // wrong domain
            "type rigid\nwall_temperature min_x 1\n",            // wrong domain
            "type particles\nseed\n",                            // missing value
            "type particles\nseed -1\n",                         // a sign (>> would wrap it)
            "type particles\nseed +5\n",                         // a sign
            "type particles\nseed 1.5\n",                        // a fraction
            "type particles\nseed 1e3\n",                        // an exponent
            "type particles\nseed 0x10\n",                       // a hex prefix
            "type particles\nseed 18446744073709551616\n",       // one past 2^64 - 1
            "type particles\nseed 7 8\n",                        // a trailing token
            "type particles\nwall_temperature min_x\n",          // missing temperature
            "type particles\nwall_temperature 1.5\n",            // missing wall
            "type particles\nwall_temperature left 1\n",         // unknown wall
            "type particles\nwall_temperature min_z 1\n",        // a 3D wall in 2D
            "type particles3d\nwall_temperature top 1\n",        // unknown wall in 3D
            "type particles\nwall_temperature min_x 1 2\n",      // a trailing token
            "type particles\nwall_temperature min_x nan\n",      // not a number
            "type particles\nwall_temperature min_x inf\n",      // not finite
            "type particles\nwall_temperature min_x 1e400\n",    // too large for a double
            "type particles\nwall_temperature min_x 1\nwall_temperature min_x 2\n",
            "type particles3d\nwall_temperature max_z 1\nwall_temperature max_z 1\n",
            // Either key marks the file as committed to its domain, so a later
            // `type` line cannot switch it (the seed would land in the wrong
            // settings, or a 2D wall name would have been checked against 2D).
            "type particles\nseed 5\ntype particles3d\n",
            "type particles\nwall_temperature min_x 1\ntype particles3d\n",
        };
        for (const char* text : bad)
        {
            std::istringstream in(text);
            const ScenarioParseResult r = parse_scenario(in);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_FALSE(r.error.empty());
        }
        // The messages name the trouble: the wall list for a wrong name (in the
        // declared dimension), and the duplicate by name on its own line.
        std::istringstream wrong("type particles3d\nwall_temperature top 1\n");
        MALLOY_CHECK_TRUE(parse_scenario(wrong).error ==
                          "line 2: unknown wall 'top' (walls are min_x, max_x, min_y, max_y, "
                          "min_z and max_z)");
        std::istringstream twice("type particles\nwall_temperature min_y 1\n"
                                 "wall_temperature max_y 1\nwall_temperature min_y 2\n");
        MALLOY_CHECK_TRUE(parse_scenario(twice).error ==
                          "line 4: wall_temperature names min_y twice");
        std::istringstream wrapped("type particles\nseed -1\n");
        MALLOY_CHECK_TRUE(parse_scenario(wrapped).error ==
                          "line 2: seed must be a whole number from 0 to "
                          "18446744073709551615, not '-1'");
        // The wall name is resolved before the temperature is read, so a bad
        // name is reported as a bad name even with no temperature after it.
        std::istringstream bare("type particles\nwall_temperature left\n");
        MALLOY_CHECK_TRUE(parse_scenario(bare).error ==
                          "line 2: unknown wall 'left' (walls are min_x, max_x, min_y and "
                          "max_y)");
        std::istringstream huge("type particles\nwall_temperature min_x 1e400\n");
        MALLOY_CHECK_TRUE(parse_scenario(huge).error ==
                          "line 2: wall_temperature requires: wall temperature");
        for (const char* text : {"type particles\nseed 5\ntype particles3d\n",
                                 "type particles\nwall_temperature min_x 1\ntype particles3d\n"})
        {
            std::istringstream in(text);
            MALLOY_CHECK_TRUE(parse_scenario(in).error ==
                              "line 3: type must come before any domain-specific key");
        }
    }

    // --- M38: `lattice` builds a whole gas in one line. It is expanded after
    //     the file is read, so `bounds` may come before or after it, and its
    //     particles follow any listed one by one. ---
    {
        std::istringstream before("type particles\n"
                                  "bounds -2 -1 2 1\n"
                                  "lattice 4 2 1.5 0.2 2.0\n");
        std::istringstream after("type particles\n"
                                 "particle 3.0 0.1 0 0 0 0\n"
                                 "lattice 4 2 1.5 0.2 2.0\n"
                                 "bounds -2 -1 2 1\n");
        const ScenarioParseResult a = parse_scenario(before);
        const ScenarioParseResult b = parse_scenario(after);
        MALLOY_CHECK_TRUE(a.ok);
        MALLOY_CHECK_TRUE(b.ok);
        MALLOY_CHECK_EQ(a.scenario.particle_list.size(), std::size_t{8});
        MALLOY_CHECK_EQ(b.scenario.particle_list.size(), std::size_t{9});
        MALLOY_CHECK_NEAR(b.scenario.particle_list[0].mass, 3.0, eps); // listed one first
        for (std::size_t i = 0; i < 8; ++i)
        {
            const auto& x = a.scenario.particle_list[i];
            const auto& y = b.scenario.particle_list[i + 1];
            MALLOY_CHECK_TRUE(x.position.x == y.position.x && x.position.y == y.position.y);
            MALLOY_CHECK_TRUE(x.velocity.x == y.velocity.x && x.velocity.y == y.velocity.y);
            MALLOY_CHECK_NEAR(x.mass, 1.5, eps);
            MALLOY_CHECK_NEAR(x.radius, 0.2, eps);
            MALLOY_CHECK_NEAR(malloy::math::length(x.velocity), 2.0, eps);
        }
        // Laid out in the FINAL bounds, 1 by 1 cells from (-2, -1).
        MALLOY_CHECK_NEAR(a.scenario.particle_list[0].position.x, -1.5, eps);
        MALLOY_CHECK_NEAR(a.scenario.particle_list[0].position.y, -0.5, eps);
        MALLOY_CHECK_NEAR(a.scenario.particle_list[7].position.x, 1.5, eps);
        MALLOY_CHECK_NEAR(a.scenario.particle_list[7].position.y, 0.5, eps);
    }
    {
        std::istringstream in("type particles3d\n"
                              "lattice3 2 3 4 1.5 0.2 2.0\n"
                              "bounds3 0 0 0 2 3 4\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_EQ(r.scenario.particle_list3d.size(), std::size_t{24});
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(r.scenario.particle_list3d[23].position,
                                                     malloy::math::Vec3{1.5, 2.5, 3.5}, eps));
        // Spacing exactly one diameter is allowed: the particles just touch.
        std::istringstream touching("type particles\nbounds -1 -1 1 1\nlattice 4 4 1 0.25 1\n");
        MALLOY_CHECK_TRUE(parse_scenario(touching).ok);
        // Mass, radius and speed are the world's to judge, not the parser's: a
        // negative radius parses and is then refused by validate().
        std::istringstream negative("type particles\nlattice 2 2 1 -0.1 1\n");
        const ScenarioParseResult n = parse_scenario(negative);
        MALLOY_CHECK_TRUE(n.ok);
        malloy::particles::ParticleWorld world{n.scenario.simulation,
                                              n.scenario.particle_settings,
                                              n.scenario.particle_list};
        MALLOY_CHECK_TRUE(world.validate() == StepStatus::InvalidState);
    }
    // Malformed lattices are refused.
    {
        const char* const bad[] = {
            "type nbody\nlattice 2 2 1 0.1 1\n",                        // wrong domain
            "type particles3d\nlattice 2 2 1 0.1 1\n",                  // 2D key in 3D
            "type particles\nlattice3 2 2 2 1 0.1 1\n",                 // 3D key in 2D
            "type particles\nlattice 2 1 0.1 1\n",                      // too few fields
            "type particles\nlattice 2 2 1 0.1 1 9\n",                  // extra token
            "type particles\nlattice 0 2 1 0.1 1\n",                    // count below 1
            "type particles\nlattice 2 -3 1 0.1 1\n",                   // negative count
            "type particles\nlattice 2.5 2 1 0.1 1\n",                  // not an integer
            "type particles\nlattice 2 2 1 0.1 1\nlattice 2 2 1 0.1 1\n", // a second one
            "type particles\nlattice 1000 1000 1 0 1\n",                // 1e6, over the cap
            "type particles\nlattice 100001 1 1 0 1\n",                 // one count over it
            "type particles3d\nlattice3 100 100 11 1 0 1\n",            // 110000
            "type particles3d\nlattice3 2147483647 2147483647 2147483647 1 0 1\n",
        };
        for (const char* text : bad)
        {
            std::istringstream in(text);
            const ScenarioParseResult r = parse_scenario(in);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_FALSE(r.error.empty());
        }
        // Spacing below a diameter would start the particles overlapping. It is
        // only known once the bounds are final, so it is reported at the end of
        // the file, but against the line holding the lattice. Each case is too
        // tight on ONE axis only, so every axis's check is needed on its own.
        const char* const crowded[] = {
            "type particles\nbounds -1 -2 1 2\nlattice 4 4 1 0.3 1\n",          // x 0.5, y 1
            "type particles\nbounds -2 -1 2 1\nlattice 4 4 1 0.3 1\n",          // x 1, y 0.5
            "type particles3d\nbounds3 -1 -2 -2 1 2 2\nlattice3 4 4 4 1 0.3 1\n", // x only
            "type particles3d\nbounds3 -2 -1 -2 2 1 2\nlattice3 4 4 4 1 0.3 1\n", // y only
            "type particles3d\nbounds3 -2 -2 -1 2 2 1\nlattice3 4 4 4 1 0.3 1\n", // z only
        };
        for (const char* text : crowded)
        {
            std::istringstream in(text);
            const ScenarioParseResult c = parse_scenario(in);
            MALLOY_CHECK_FALSE(c.ok);
            MALLOY_CHECK_TRUE(c.error.find("line 3") != std::string::npos);
            MALLOY_CHECK_TRUE(c.error.find("diameter") != std::string::npos);
        }
        std::istringstream flat("type particles3d\n"
                                "lattice3 2 2 5 1 0.25 1\n"
                                "bounds3 -1 -1 -1 1 1 1\n"); // z spacing 0.4 < 0.5
        const ScenarioParseResult f = parse_scenario(flat);
        MALLOY_CHECK_FALSE(f.ok);
        MALLOY_CHECK_TRUE(f.error.find("line 2") != std::string::npos);
        // The cap is reported on the lattice line itself, before allocating,
        // not left to the end-of-file total (which would name the last line).
        std::istringstream huge("type particles\nlattice 1000 1000 1 0 1\ndt 0.1\n");
        const ScenarioParseResult h = parse_scenario(huge);
        MALLOY_CHECK_FALSE(h.ok);
        MALLOY_CHECK_TRUE(h.error.find("line 2") != std::string::npos);
    }

    // --- Issue #31: entity_count is what the app checks before running. A
    //     settings-only fragment still parses (the parser checks syntax and
    //     settings), but it counts zero, so the app refuses it. ---
    {
        const auto count_of = [](const char* text) {
            std::istringstream in(text);
            const ScenarioParseResult r = parse_scenario(in);
            return r.ok ? malloy::scenario::entity_count(r.scenario) : std::size_t{999999};
        };
        MALLOY_CHECK_EQ(count_of(""), std::size_t{0});
        MALLOY_CHECK_EQ(count_of("# only a comment\n"), std::size_t{0});
        MALLOY_CHECK_EQ(count_of("type charges\nsteps 1\n"), std::size_t{0});
        MALLOY_CHECK_EQ(count_of("type particles\nfriction 0.5\n"), std::size_t{0});
        MALLOY_CHECK_EQ(count_of("body 1 0 0 0 0\nbody 1 1 0 0 0\n"), std::size_t{2});
        MALLOY_CHECK_EQ(count_of("type particles\nbounds -2 -2 2 2\nlattice 3 2 1 0.1 1\n"),
                        std::size_t{6});
        MALLOY_CHECK_EQ(count_of("type springs\nspring_body 1 0 0 0 0\nspring_body 1 1 0 0 0\n"
                                 "spring 0 1 1 1 0\n"),
                        std::size_t{3}); // two bodies and the spring
        MALLOY_CHECK_EQ(count_of("type charges3d\ncharge3 1 1 0 0 0 0 0 0\n"), std::size_t{1});
        // One per remaining list, so a list left out of the sum shows.
        MALLOY_CHECK_EQ(count_of("type charges\ncharge 1 1 0 0 0 0\n"), std::size_t{1});
        MALLOY_CHECK_EQ(count_of("type nbody3d\nbody3 1 0 0 0 0 0 0\n"), std::size_t{1});
        MALLOY_CHECK_EQ(count_of("type particles\nparticle 1 0.1 0 0 0 0\n"), std::size_t{1});
        MALLOY_CHECK_EQ(count_of("type particles3d\nparticle3 1 0.1 0 0 0 0 0 0\n"),
                        std::size_t{1});
        MALLOY_CHECK_EQ(count_of("type rigid\nrigid_static 0.5 0 0 0\n"), std::size_t{1});
        MALLOY_CHECK_EQ(
            count_of("type rigid3d\nrigid_box3d 1 0.5 0.5 0.5 0 0 0 0 0 1 0 0 0 0 0 0 0\n"),
            std::size_t{1});
        MALLOY_CHECK_EQ(count_of("type springs3d\nspring_body3 1 0 0 0 0 0 0\n"
                                 "spring_body3 1 1 0 0 0 0 0\nspring 0 1 1 1 0\n"),
                        std::size_t{3});
    }

    // --- M38: before any time has passed there is no pressure, so a
    //     `# check step 0 pressure` line FAILS rather than comparing a 0. ---
    {
        std::istringstream in("type particles\nbounds -2 -2 2 2\nlattice 2 2 1 0.1 1\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        malloy::particles::ParticleWorld world{r.scenario.simulation,
                                              r.scenario.particle_settings,
                                              r.scenario.particle_list};
        const auto value_of = [&](const std::string&) {
            return malloy::particles::pressure(world.wall_impulse(),
                                               r.scenario.particle_settings.bounds,
                                               world.elapsed_time());
        };
        std::ostringstream sink;
        std::streambuf* const saved = std::cerr.rdbuf(sink.rdbuf());
        const bool passed = check_passes("synthetic", TemplateCheck{0, "pressure", 0.0, 1.0}, value_of);
        std::cerr.rdbuf(saved);
        MALLOY_CHECK_FALSE(passed);
        MALLOY_CHECK_TRUE(sink.str().find("unavailable") != std::string::npos);
        MALLOY_CHECK_TRUE(world.step().ok());
        MALLOY_CHECK_TRUE(value_of("pressure").has_value());
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

    // --- Issue #48: in 2D, box against disc is not built yet, so a box and a
    //     disc used to pass through each other without a word. Until that pair
    //     has its milestone, a `type rigid` file holding a rigid_box and a disc
    //     that collides (rigid_body or rigid_static, radius above 0) is refused,
    //     on the line where the second kind first appears. ---
    {
        const auto parse_text = [](const std::string& text) {
            std::istringstream in(text);
            return parse_scenario(in);
        };
        const std::string box = "rigid_box 2.0 0.5 0.5 -1.2 0.0 0.0 1.0 0.0 0.0\n";
        const std::string body = "rigid_body 1.0 0.125 0.5 0 0 1.2 0.0 0.0 -1.0 0.0 0.0\n";
        const std::string fixed = "rigid_static 0.5 1.2 0.0 0.0\n";
        const std::string body0 = "rigid_body 1.0 0.125 0 0 0 1.2 0.0 0.0 -1.0 0.0 0.0\n";
        const std::string fixed0 = "rigid_static 0 1.2 0.0 0.0\n";

        // The issue's reproduction, then every ordering of each disc kind. The
        // error names the line that brought in the second kind, says why, and
        // points back at the first.
        const struct
        {
            std::string text;
            std::string where;
            std::string names;
        } mixes[] = {
            {"type rigid\n" + box + body, "line 3:", "rigid_box on line 2"},
            {"type rigid\n" + body + box, "line 3:", "rigid_body on line 2"},
            {"type rigid\n" + box + fixed, "line 3:", "rigid_box on line 2"},
            {"type rigid\n" + fixed + box, "line 3:", "rigid_static on line 2"},
            // Comments and blank lines count as lines, and the error is on the
            // FIRST line that mixes the kinds, not a later one.
            {"type rigid\n# a box\n" + box + box + "\n" + body + body, "line 6:",
             "rigid_box on line 3"},
            {"type rigid\n" + fixed + body + "# and now\n" + box + box, "line 5:",
             "rigid_static on line 2"},
        };
        for (const auto& mix : mixes)
        {
            const ScenarioParseResult r = parse_text(mix.text);
            MALLOY_CHECK_FALSE(r.ok);
            if (r.error.rfind(mix.where, 0) != 0 ||
                r.error.find("box against disc is not supported yet") == std::string::npos ||
                r.error.find(mix.names) == std::string::npos)
            {
                std::cerr << "unexpected box/disc error: " << r.error << '\n';
                return 1;
            }
        }
        {
            const ScenarioParseResult r = parse_text("type rigid\n" + box + body);
            MALLOY_CHECK_EQ(r.error, std::string("line 3: box against disc is not supported "
                                                 "yet: this rigid_body would pass through the "
                                                 "rigid_box on line 2 (a disc with a radius "
                                                 "above 0 cannot share a scenario with a box)"));
        }

        // A body of radius 0 collides with nothing, so it may share a file with
        // boxes, before or after them, and so may boxes and discs on their own.
        for (const std::string& text :
             {"type rigid\n" + box + body0 + fixed0, "type rigid\n" + body0 + fixed0 + box,
              "type rigid\n" + box + box, "type rigid\n" + body + fixed})
        {
            const ScenarioParseResult r = parse_text(text);
            if (!r.ok)
            {
                std::cerr << "refused a file that does not mix box and disc: " << r.error
                          << '\n';
                return 1;
            }
        }
        // The refusal is the parser's, so a line's own syntax error still comes
        // first: a trailing token on the mixing line is reported as such.
        {
            const ScenarioParseResult r =
                parse_text("type rigid\n" + box + "rigid_static 0.5 1.2 0.0 0.0 junk\n");
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_TRUE(r.error.find("junk") != std::string::npos);
        }
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
        // Since M35 friction belongs to the particle domains too, in both
        // dimensions, writing the field the declared domain reads.
        std::istringstream in("type particles\nfriction 0.5\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.particle_settings.friction, 0.5, 1e-12);
    }
    {
        std::istringstream in("type particles3d\nfriction 0.25\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.particle3d_settings.friction, 0.25, 1e-12);
    }
    {
        // ...but not springs, which has no friction.
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
        // M30: rigid_box3d builds a box body whose inertia is COMPUTED from the
        // box and mass, not given. Seventeen fields: mass, three half-widths,
        // position, axis+angle, velocity, angular velocity.
        std::istringstream in(
            "type rigid3d\n"
            "rigid_box3d 48.0  3.0 2.0 1.0   4.0 -5.0 6.0   0 0 1 0   1 0 0   0 0 2\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_EQ(r.scenario.rigid_bodies3d.size(), std::size_t{1});
        const auto& body = r.scenario.rigid_bodies3d[0];
        MALLOY_CHECK_NEAR(body.mass, 48.0, 1e-12);
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(body.half_extents,
                                                     malloy::math::Vec3{3.0, 2.0, 1.0}, 0.0));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            body.position, malloy::math::Vec3{4.0, -5.0, 6.0}, 0.0));
        // Inertia is the solid box's, (1/3) m (h_j^2 + h_k^2) with h the half
        // widths: m = 48, so about x it is 16*(2^2 + 1^2) = 80, about y
        // 16*(3^2 + 1^2) = 160, about z 16*(3^2 + 2^2) = 208.
        MALLOY_CHECK_NEAR(body.inertia.x, 80.0, 1e-9);
        MALLOY_CHECK_NEAR(body.inertia.y, 160.0, 1e-9);
        MALLOY_CHECK_NEAR(body.inertia.z, 208.0, 1e-9);
    }
    {
        // Issue #21: the moments must follow the box's OWN axes. The case above
        // has descending half-widths, whose moments happen to ascend, so it
        // could not tell the box-frame moments from the sorted principal ones
        // (the symmetric-configuration trap of docs/05). These half-widths are
        // ascending, (1, 2, 3), so the moments DESCEND: about x it is
        // 16*(2^2 + 3^2) = 208, about y 16*(1^2 + 3^2) = 160, about z
        // 16*(1^2 + 2^2) = 80. Sorted, they came back as (80, 160, 208), the
        // x and z moments swapped. A middle-axis case pins y as well.
        std::istringstream in(
            "type rigid3d\n"
            "rigid_box3d 48.0  1.0 2.0 3.0   0 0 0   0 0 1 0   0 0 0   0 0 0\n"
            "rigid_box3d 48.0  2.0 3.0 1.0   0 0 5   0 0 1 0   0 0 0   0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        const auto& ascending = r.scenario.rigid_bodies3d[0];
        MALLOY_CHECK_NEAR(ascending.inertia.x, 208.0, 1e-9);
        MALLOY_CHECK_NEAR(ascending.inertia.y, 160.0, 1e-9);
        MALLOY_CHECK_NEAR(ascending.inertia.z, 80.0, 1e-9);
        // (2, 3, 1): about x 16*(9 + 1) = 160, about y 16*(4 + 1) = 80,
        // about z 16*(4 + 9) = 208.
        const auto& rotated = r.scenario.rigid_bodies3d[1];
        MALLOY_CHECK_NEAR(rotated.inertia.x, 160.0, 1e-9);
        MALLOY_CHECK_NEAR(rotated.inertia.y, 80.0, 1e-9);
        MALLOY_CHECK_NEAR(rotated.inertia.z, 208.0, 1e-9);
        // The rotational energy a spin about z carries is the box's, with the
        // z moment 80: (1/2) 80 1^2 = 40, not the 104 the swap gave.
        malloy::rigid::RigidBody3D spinning = ascending;
        spinning.angular_velocity = malloy::math::Vec3{0.0, 0.0, 1.0};
        MALLOY_CHECK_NEAR(malloy::rigid::total_kinetic_energy3d({spinning}), 40.0, 1e-9); // at rest, so all rotational
        // A degenerate box still gives zero moments, so validate() rejects it.
        std::istringstream flat("type rigid3d\n"
                                "rigid_box3d 1.0  1.0 0.0 1.0   0 0 0   0 0 1 0   0 0 0   0 0 0\n");
        const ScenarioParseResult f = parse_scenario(flat);
        MALLOY_CHECK_TRUE(f.ok);
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(f.scenario.rigid_bodies3d[0].inertia,
                                                     malloy::math::Vec3{}, 0.0));
    }
    {
        // rigid_box3d belongs to rigid3d only, and needs all seventeen fields.
        std::istringstream wrong_type(
            "type rigid\nrigid_box3d 1 1 1 1 0 0 0 0 0 1 0 0 0 0 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(wrong_type).ok);
        std::istringstream too_few("type rigid3d\nrigid_box3d 1 1 1 1 0 0 0\n");
        MALLOY_CHECK_FALSE(parse_scenario(too_few).ok);
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
    // --- Issue #47: the force key, M28's constant applied force, read exactly
    //     as torque is. Distinct values per axis, so a transposed or dropped
    //     component shows, and the torque is checked to stay untouched. ---
    {
        std::istringstream in("type rigid3d\n"
                              "force 0.5 -1.5 2.5\n"
                              "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.force.x, 0.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.force.y, -1.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.force.z, 2.5, eps);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.torque.x, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.torque.y, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.torque.z, 0.0, 0.0);
        MALLOY_CHECK_TRUE(r.scenario.rigid3d_settings.is_valid());

        // Absent, it is zero: a pre-M28 file is unchanged.
        std::istringstream none("type rigid3d\n"
                                "rigid_body3d 1  1 2 3  0.5  0 0 0  0 0 1 0  0 0 0  1 0 0\n");
        const ScenarioParseResult r2 = parse_scenario(none);
        MALLOY_CHECK_TRUE(r2.ok);
        MALLOY_CHECK_NEAR(r2.scenario.rigid3d_settings.force.x, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r2.scenario.rigid3d_settings.force.y, 0.0, 0.0);
        MALLOY_CHECK_NEAR(r2.scenario.rigid3d_settings.force.z, 0.0, 0.0);
    }
    {
        // The parsed force reaches the world, as the acceleration F/m: one step
        // of a body of mass 4 under (2, -6, 8) changes its velocity by
        // (F/m) dt = (0.5, -1.5, 2) * 0.001, and its momentum by F dt.
        std::istringstream in("type rigid3d\ndt 0.001\n"
                              "force 2 -6 8\n"
                              "rigid_body3d 4  1 1 1  0.5  0 0 0  0 0 1 0  0 0 0  0 0 0\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        malloy::rigid::Rigid3DWorld world{r.scenario.simulation, r.scenario.rigid_bodies3d,
                                          r.scenario.rigid3d_settings};
        MALLOY_CHECK_TRUE(world.step().ok());
        const auto& v = world.bodies().front().velocity;
        MALLOY_CHECK_NEAR(v.x, 0.0005, 1e-18);
        MALLOY_CHECK_NEAR(v.y, -0.0015, 1e-18);
        MALLOY_CHECK_NEAR(v.z, 0.002, 1e-18);
    }
    {
        // Malformed: two fields, four, a word, none. Each is refused on its
        // line, with the message torque's would have.
        const struct
        {
            const char* text;
            const char* message;
        } bad[] = {
            {"type rigid3d\nforce 1 2\n", "line 2: force requires: fx fy fz"},
            {"type rigid3d\nforce 1 two 3\n", "line 2: force requires: fx fy fz"},
            {"type rigid3d\nforce\n", "line 2: force requires: fx fy fz"},
            {"type rigid3d\nforce 1 2 3 4\n", "line 2: unexpected extra token '4'"},
        };
        for (const auto& b : bad)
        {
            std::istringstream in(b.text);
            const ScenarioParseResult r = parse_scenario(in);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_EQ(r.error, std::string(b.message));
        }
    }
    {
        // force belongs to rigid3d and to no other domain, including the
        // default nbody a file gets when its `type` line is missing or late.
        for (const char* head : {"", "type nbody\n", "type rigid\n", "type nbody3d\n",
                                 "type particles3d\n", "type springs3d\n", "type charges3d\n"})
        {
            std::istringstream in(std::string(head) + "force 1 2 3\n");
            const ScenarioParseResult r = parse_scenario(in);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_TRUE(r.error.find("force belongs to type rigid3d") != std::string::npos);
        }
    }
    {
        // A second force line replaces the first, exactly as a second torque
        // line does: both are single world settings, not accumulators.
        std::istringstream in("type rigid3d\nforce 1 2 3\nforce 4 5 6\n"
                              "torque 1 2 3\ntorque 7 8 9\n");
        const ScenarioParseResult r = parse_scenario(in);
        MALLOY_CHECK_TRUE(r.ok);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.force.x, 4.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.force.y, 5.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.force.z, 6.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.torque.x, 7.0, 0.0);
        MALLOY_CHECK_NEAR(r.scenario.rigid3d_settings.torque.z, 9.0, 0.0);
    }
    {
        // A force whose square overflows parses, and is refused by the
        // settings, so that rule lives in one place (as for torque).
        std::istringstream in("type rigid3d\n"
                              "force 0 1e200 0\n"
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

    // --- Security hardening (2026-09-13 review): resource caps on the one
    //     untrusted input (a scenario file) and parser robustness against hostile
    //     content. The parser already rejected non-finite and out-of-range
    //     values; these bound sheer MAGNITUDE (work and size) and confirm no
    //     malformed input crashes, hangs, or reads out of bounds. ---
    {
        // The step cap: exactly at the limit parses, one over is rejected, and a
        // value past INT_MAX is refused by the extractor rather than wrapped.
        {
            // output_every 10000 keeps the run well inside the 1e5-report cap
            // (issue #25).
            std::istringstream at("type nbody\nsteps 100000000\noutput_every 10000\n");
            MALLOY_CHECK_TRUE(parse_scenario(at).ok);
            std::istringstream over("type nbody\nsteps 100000001\n");
            MALLOY_CHECK_FALSE(parse_scenario(over).ok);
            std::istringstream past_int("type nbody\nsteps 99999999999\n");
            MALLOY_CHECK_FALSE(parse_scenario(past_int).ok);
        }

        // The entity cap: exactly at the limit parses, one over is rejected. Built
        // as bare particle lines, which the parser only has to count. `steps 0`
        // keeps the WORK at zero, so this tests the entity cap alone (100000
        // particles for the default 1000 steps is far over the work budget).
        {
            std::ostringstream at;
            at << "type particles\nsteps 0\n";
            for (int i = 0; i < 100000; ++i)
            {
                at << "particle 1 0.01 0 0 0 0\n";
            }
            std::istringstream at_in(at.str());
            MALLOY_CHECK_TRUE(parse_scenario(at_in).ok);

            // Issue #29: the cap is checked as each line is read, so the error
            // names the line that crossed it (line 100003, after the type and
            // steps lines), not the last line of the file, and nothing past it
            // is allocated.
            std::ostringstream over;
            over << "type particles\nsteps 0\n";
            for (int i = 0; i < 100010; ++i)
            {
                over << "particle 1 0.01 0 0 0 0\n";
            }
            std::istringstream over_in(over.str());
            const ScenarioParseResult r = parse_scenario(over_in);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_TRUE(r.error.find("line 100003:") != std::string::npos);
        }

        // Issue #29: a lattice that, with the listed particles, crosses the cap
        // is reported against the lattice's line.
        {
            std::istringstream in("type particles\nsteps 0\nlattice 100 1000 1 0 1\n"
                                  "particle 1 0.01 0 0 0 0\nbounds -1000 -1000 1000 1000\n");
            const ScenarioParseResult r = parse_scenario(in);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_TRUE(r.error.find("line 3:") != std::string::npos);
        }

        // Issue #26: ground planes have their own cap, checked per line.
        {
            std::ostringstream at;
            at << "type rigid\nsteps 1\n";
            for (int i = 0; i < 1000; ++i)
            {
                at << "ground 0 1 -100\n";
            }
            at << "rigid_static 0.5 0 0 0\n";
            std::istringstream at_in(at.str());
            MALLOY_CHECK_TRUE(parse_scenario(at_in).ok);

            std::ostringstream over;
            over << "type rigid3d\nsteps 1\n";
            for (int i = 0; i < 1001; ++i)
            {
                over << "plane3 0 0 1 -100\n";
            }
            std::istringstream over_in(over.str());
            const ScenarioParseResult r = parse_scenario(over_in);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_TRUE(r.error.find("line 1003:") != std::string::npos);
            MALLOY_CHECK_TRUE(r.error.find("plane") != std::string::npos);
        }

        // M40: the contact point cap. A rigid world keeps a step's contact
        // points in memory together, so the parser bounds how many it could
        // have to hold: every colliding pair and every collider against every
        // ground plane, at the most points that contact reports (rigid3d box
        // pair or box-plane 8, anything else 1; rigid box pair 2, box-ground 4,
        // anything else 1), at most 250000. Each file below is built to sit
        // EXACTLY at 250000, so it parses and potential_contact_points reports
        // 250000, and one more body or plane is refused on that line. A weight
        // that drifts either way fails one side or the other. `steps 0` keeps
        // the work far inside its own budget, so the cap is tested alone.
        {
            const auto repeat = [](int count, const std::string& line) {
                std::string lines;
                for (int i = 0; i < count; ++i)
                {
                    lines += line;
                }
                return lines;
            };
            const std::string head3 = "type rigid3d\nsteps 0\noutput_every 0\n";
            const std::string head2 = "type rigid\nsteps 0\noutput_every 0\n";
            const std::string box3 = "rigid_box3d 1 0.5 0.5 0.5 0 0 0 1 0 0 0 0 0 0 0 0 0\n";
            const std::string sphere3 = "rigid_body3d 1 1 1 1 0.5 0 0 0 1 0 0 0 0 0 0 0 0 0\n";
            const std::string ghost3 = "rigid_body3d 1 1 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0\n";
            const std::string plane3 = "plane3 0 0 1 -100\n";
            const std::string box2 = "rigid_box 1 0.5 0.5 0 0 0 0 0 0\n";
            const std::string disc2 = "rigid_body 1 0.125 0.5 0 0 0 0 0 0 0 0\n";
            const std::string static2 = "rigid_static 0.5 0 0 0\n";
            const std::string ghost2 = "rigid_body 1 0.125 0 0 0 0 0 0 0 0 0\n";
            const std::string ground = "ground 0 1 -100\n";
            const std::string refusal = "scenario exceeds the maximum contact point count (250000)";

            // Parses `at`, which must sit exactly at the cap, then `at` plus
            // `extra`, which must be refused on line `line`. Lines follow the
            // one that crosses, so a check made only at the end of the file
            // would name the wrong line.
            const auto at_cap_edge = [&](const std::string& label, const std::string& at,
                                         const std::string& extra, std::size_t line) {
                std::istringstream at_in(at);
                const ScenarioParseResult inside = parse_scenario(at_in);
                std::istringstream over_in(at + extra + "# past the cap\n# and more\n");
                const ScenarioParseResult over = parse_scenario(over_in);
                const std::string expected = "line " + std::to_string(line) + ": " + refusal;
                const std::uint64_t points =
                    inside.ok ? malloy::scenario::potential_contact_points(inside.scenario) : 0;
                if (!inside.ok || points != 250000 || over.ok || over.error != expected)
                {
                    std::cerr << label << ": at the cap ok=" << inside.ok << " ("
                              << inside.error << "), points " << points
                              << "; one over ok=" << over.ok << " (" << over.error << ")\n";
                    return false;
                }
                return true;
            };

            // rigid3d boxes and planes: 8 x (125 x 124 / 2) + 8 x 188 x 125
            // = 62000 + 188000. The 126th box is line 3 + 188 + 126 = 317;
            // a 189th plane after the boxes is line 3 + 188 + 125 + 1 = 317.
            const std::string boxes3 = head3 + repeat(188, plane3) + repeat(125, box3);
            MALLOY_CHECK_TRUE(at_cap_edge("rigid3d boxes, one more box", boxes3, box3, 317));
            MALLOY_CHECK_TRUE(at_cap_edge("rigid3d boxes, one more plane", boxes3, plane3, 317));

            // rigid3d spheres and planes: 625 x 624 / 2 + 88 x 625
            // = 195000 + 55000, one point each.
            const std::string spheres3 = head3 + repeat(88, plane3) + repeat(625, sphere3);
            MALLOY_CHECK_TRUE(at_cap_edge("rigid3d spheres", spheres3, sphere3, 3 + 88 + 626));

            // rigid3d, every kind at once: 176 boxes, 21 spheres and 86 planes,
            // 8 x 15400 + 176 x 21 + 210 + 86 x (8 x 176 + 21)
            // = 123200 + 3696 + 210 + 122894. One more sphere adds 176 + 21 + 86.
            const std::string mixed3 =
                head3 + repeat(86, plane3) + repeat(176, box3) + repeat(21, sphere3);
            MALLOY_CHECK_TRUE(at_cap_edge("rigid3d mixed", mixed3, sphere3, 3 + 86 + 176 + 22));

            // rigid boxes and grounds: 2 x (125 x 124 / 2) + 4 x 469 x 125
            // = 15500 + 234500.
            const std::string boxes2 = head2 + repeat(469, ground) + repeat(125, box2);
            MALLOY_CHECK_TRUE(at_cap_edge("rigid boxes", boxes2, box2, 3 + 469 + 126));

            // rigid discs and grounds, half of them immovable (a static disc
            // still collides with a movable one): 625 x 624 / 2 + 88 x 625.
            const std::string discs2 =
                head2 + repeat(88, ground) + repeat(313, disc2) + repeat(312, static2);
            MALLOY_CHECK_TRUE(at_cap_edge("rigid discs", discs2, static2, 3 + 88 + 626));
            MALLOY_CHECK_TRUE(at_cap_edge("rigid discs, one more ground", discs2, ground,
                                          3 + 88 + 626));

            // One point over is refused: 189 boxes, 2 spheres and 71 planes are
            // 8 x 17766 + 378 + 1 + 71 x 1514 = 250001, crossed by the second
            // sphere on line 3 + 71 + 189 + 2. Without it, 249740, inside.
            {
                const std::string almost =
                    head3 + repeat(71, plane3) + repeat(189, box3) + sphere3;
                std::istringstream inside(almost);
                const ScenarioParseResult in = parse_scenario(inside);
                MALLOY_CHECK_TRUE(in.ok);
                MALLOY_CHECK_TRUE(in.ok &&
                                  malloy::scenario::potential_contact_points(in.scenario) == 249740);
                std::istringstream one_over(almost + sphere3 + "# past the cap\n");
                const ScenarioParseResult r = parse_scenario(one_over);
                MALLOY_CHECK_FALSE(r.ok);
                MALLOY_CHECK_EQ(r.error, "line 265: " + refusal);
            }

            // A body that collides with nothing holds no contact points, so any
            // number of them beside a file at the cap still parses, in either
            // dimension, and leaves the count where it was. That includes a box
            // with a zero half-extent and no radius, which its world does not
            // treat as a box (the degenerate box is left for validate()).
            const std::string flat3 = "rigid_box3d 1 0.5 0.5 0 0 0 0 1 0 0 0 0 0 0 0 0 0\n";
            const std::string flat2 = "rigid_box 1 0.5 0 0 0 0 0 0 0\n";
            for (const std::string& file :
                 {spheres3 + repeat(500, ghost3) + repeat(100, flat3),
                  boxes2 + repeat(500, ghost2) + repeat(100, flat2)})
            {
                std::istringstream in(file);
                const ScenarioParseResult r = parse_scenario(in);
                MALLOY_CHECK_TRUE(r.ok);
                MALLOY_CHECK_TRUE(r.ok &&
                                  malloy::scenario::potential_contact_points(r.scenario) == 250000);
            }

            // Only the rigid domains keep contact points: 1000 particles, which
            // collide, count none.
            std::istringstream gas("type particles\nsteps 0\n" +
                                   repeat(1000, "particle 1 0.01 0 0 0 0\n"));
            const ScenarioParseResult g = parse_scenario(gas);
            MALLOY_CHECK_TRUE(g.ok);
            MALLOY_CHECK_TRUE(g.ok && malloy::scenario::potential_contact_points(g.scenario) == 0);
        }

        // Issue #25: the WORK budget, steps plus reports, times the work a step
        // does, at most 1e11 units (a report's diagnostics cost about a step).
        // Issue #46 weighted it: a step costs a unit per body plus its tests,
        // nbody counts every ORDERED pair, as its world evaluates them, and the
        // rigid domains weight each pair and plane test by kind (box-box 250,
        // box-sphere 25, sphere pair 10, box-plane 40, sphere-plane 4 in 3D;
        // box-box 40, any other pair 5, box-ground 20, disc-ground 4 in 2D).
        //
        // Each case below states its work per step W by hand. With output_every
        // 0 there is one report, so s steps do (s + 1) W: the most steps inside
        // 1e11 is floor(1e11 / W) - 1, and one more is refused. Both sides of
        // that edge are checked, so a weight that drifts either way fails.
        {
            const auto at_budget_edge = [](const std::string& head, const std::string& lines,
                                           long long per_step) {
                const long long inside = 100000000000LL / per_step - 1;
                const auto parse_steps = [&](long long steps) {
                    std::istringstream in(head + "output_every 0\nsteps " +
                                          std::to_string(steps) + "\n" + lines);
                    return parse_scenario(in);
                };
                const ScenarioParseResult in_budget = parse_steps(inside);
                const ScenarioParseResult over_budget = parse_steps(inside + 1);
                if (!in_budget.ok || over_budget.ok ||
                    over_budget.error !=
                        "scenario exceeds the work budget: steps plus reports, times pair "
                        "tests per step, must be at most 1e11 (fewer steps or fewer bodies)")
                {
                    std::cerr << head << "work budget edge at " << inside
                              << " steps: inside ok=" << in_budget.ok << " (" << in_budget.error
                              << "), one over ok=" << over_budget.ok << " ("
                              << over_budget.error << ")\n";
                    return false;
                }
                // Issue #84: the exported step_work is the figure the parser
                // charges, so a caller bounding a run with it (the fuzzers)
                // counts exactly what the parser counts. Every W here is an
                // integer below 2^53, so the comparison is exact.
                const double exported = malloy::scenario::step_work(in_budget.scenario);
                if (exported != static_cast<double>(per_step))
                {
                    std::cerr << head << "step_work gives " << exported << ", the edge needs "
                              << per_step << '\n';
                    return false;
                }
                return true;
            };
            const auto repeat = [](int count, const auto& line_for) {
                std::ostringstream lines;
                for (int i = 0; i < count; ++i)
                {
                    lines << line_for(i);
                }
                return lines.str();
            };
            const auto nbody_line = [](int i) { return "body 1 " + std::to_string(i) + " 0 0 0\n"; };
            const auto body3_line = [](int i) {
                return "body3 1 " + std::to_string(i) + " 0 0 0 0 0\n";
            };
            const auto box3_line = [](int i) {
                return "rigid_box3d 1 1 1 1 " + std::to_string(10 * (i % 10)) + " " +
                       std::to_string(10 * (i / 10)) + " 0 0 0 1 0 0 0 0 0 0 0\n";
            };
            const auto sphere3_line = [](int i) {
                return "rigid_body3d 1 1 1 1 0.5 " + std::to_string(10 * i) +
                       " 50 0 0 0 1 0 0 0 0 0 0 0\n";
            };
            const auto box2_line = [](int i) {
                return "rigid_box 1 0.5 0.5 " + std::to_string(2 * i) + " 0 0 0 0 0\n";
            };
            const auto disc2_line = [](int i) {
                return "rigid_body 1 0.125 0.5 0 0 " + std::to_string(2 * i) + " 9 0 0 0 0\n";
            };
            const auto point2_line = [](int i) {
                return i % 2 ? "rigid_static 0 " + std::to_string(2 * i) + " 9 0\n"
                             : "rigid_body 1 0.125 0 0 0 " + std::to_string(2 * i) + " 9 0 0 0 0\n";
            };

            // nbody: 1000 bodies, 1000 x 999 ordered pairs + 1000 = 1e6 per
            // step, so 99999 steps (1e5 x 1e6 = 1e11) are inside and 100000
            // over. Counted as unordered pairs, as before issue #46, it ran
            // twice what it was charged.
            const std::string bodies = repeat(1000, nbody_line);
            MALLOY_CHECK_TRUE(at_budget_edge("type nbody\n", bodies, 1000000));
            // nbody3d tests each pair once: 499500 + 1000 = 500500.
            MALLOY_CHECK_TRUE(at_budget_edge("type nbody3d\n", repeat(1000, body3_line), 500500));
            // rigid3d, 100 boxes: 250 x 4950 + 100 = 1237600.
            const std::string boxes = repeat(100, box3_line);
            MALLOY_CHECK_TRUE(at_budget_edge("type rigid3d\n", boxes, 1237600));
            // rigid3d, 100 spheres: 10 x 4950 + 100 = 49600.
            MALLOY_CHECK_TRUE(at_budget_edge("type rigid3d\n", repeat(100, sphere3_line), 49600));
            // rigid3d, 10 boxes, 10 spheres and a plane: 250 x 45 + 25 x 100 +
            // 10 x 45 pairs, 40 x 10 + 4 x 10 plane tests, 20 bodies: 14660.
            MALLOY_CHECK_TRUE(at_budget_edge(
                "type rigid3d\nplane3 0 0 1 -100\n",
                repeat(10, box3_line) + repeat(10, sphere3_line), 14660));
            // rigid3d, a box and a sphere under 1000 planes, where the plane
            // tests dominate: 25 + 1000 x (40 + 4) + 2 = 44027.
            MALLOY_CHECK_TRUE(at_budget_edge(
                "type rigid3d\n" + repeat(1000, [](int) { return "plane3 0 0 1 -100\n"; }),
                repeat(1, box3_line) + repeat(1, sphere3_line), 44027));
            // rigid, 100 boxes: 40 x 4950 + 100 = 198100.
            MALLOY_CHECK_TRUE(at_budget_edge("type rigid\n", repeat(100, box2_line), 198100));
            // rigid, 100 discs and 10 grounds: 5 x 4950 + 10 x 4 x 100 + 100 = 28850.
            MALLOY_CHECK_TRUE(at_budget_edge(
                "type rigid\n" + repeat(10, [](int) { return "ground 0 1 -100\n"; }),
                repeat(100, disc2_line), 28850));
            // rigid, 10 boxes, 10 bodies of radius 0 and a ground: 40 x 45 +
            // 5 x (100 + 45) pairs, 20 x 10 + 4 x 10 ground tests, 20 bodies: 2785.
            MALLOY_CHECK_TRUE(at_budget_edge(
                "type rigid\nground 0 1 -100\n", repeat(10, box2_line) + repeat(10, point2_line),
                2785));
            // rigid, a box and a body of radius 0 under 1000 grounds: 5 +
            // 1000 x (20 + 4) + 2 = 24007.
            MALLOY_CHECK_TRUE(at_budget_edge(
                "type rigid\n" + repeat(1000, [](int) { return "ground 0 1 -100\n"; }),
                repeat(1, box2_line) + repeat(1, point2_line), 24007));

            // Issue #46's file: 100 separated boxes for 19799000 steps, 4828
            // bytes, inside the old budget at 5050 pair tests a step and about
            // 2.35 hours in Release. Now 1237600 units a step, far over.
            {
                std::istringstream issue("type rigid3d\nsteps 19799000\noutput_every 100000\n" +
                                         boxes);
                const ScenarioParseResult r = parse_scenario(issue);
                MALLOY_CHECK_FALSE(r.ok);
                MALLOY_CHECK_TRUE(r.error.find("work budget") != std::string::npos);
            }

            // The reports count: 60000 steps alone are 6e10, inside, but a
            // report every step adds 60001 reports, and their diagnostics bring
            // the work to 1.2e11, over. Counted by steps alone, this used to
            // pass.
            std::istringstream reporting("type nbody\noutput_every 1\nsteps 60000\n" + bodies);
            const ScenarioParseResult rr = parse_scenario(reporting);
            MALLOY_CHECK_FALSE(rr.ok);
            MALLOY_CHECK_TRUE(rr.error.find("work budget") != std::string::npos);

            // The 41-byte file the issue reported: 99856 particles for the
            // default 1000 steps is about 5e12 pair tests.
            std::istringstream tiny("type particles\nlattice 316 316 1 0.001 1\n");
            MALLOY_CHECK_FALSE(parse_scenario(tiny).ok);

            // A ground plane costs a test per body per step in the rigid
            // domains: 1000 planes and 100 static discs for 1e6 steps is 4 x
            // 1000 x 100 plane tests plus 5 x 4950 pairs and 100 bodies a
            // step, 4.2e11, over.
            std::ostringstream planes;
            planes << "type rigid\noutput_every 0\nsteps 1000000\n";
            for (int i = 0; i < 1000; ++i)
            {
                planes << "ground 0 1 -100\n";
            }
            for (int i = 0; i < 100; ++i)
            {
                planes << "rigid_static 0.1 " << i << " 0 0\n";
            }
            std::istringstream planes_in(planes.str());
            MALLOY_CHECK_FALSE(parse_scenario(planes_in).ok);

            // Springs are linear: a body and a spring per step. A chain of 1000
            // bodies and 999 springs for 1e6 steps is 2e9, inside; counted as
            // pairs it would be 5e11, far over.
            std::ostringstream chain;
            chain << "type springs\noutput_every 0\nsteps 1000000\n";
            for (int i = 0; i < 1000; ++i)
            {
                chain << "spring_body 1 " << i << " 0 0 0\n";
            }
            for (int i = 0; i < 999; ++i)
            {
                chain << "spring " << i << ' ' << i + 1 << " 1 1 0\n";
            }
            std::istringstream springs(chain.str());
            MALLOY_CHECK_TRUE(parse_scenario(springs).ok);
        }

        // Issue #84: step_work for the domains the edges above do not reach,
        // by the documented formula: N(N-1)/2 + N for the other pairwise
        // domains, N + springs for springs, and nothing for an empty file.
        {
            const auto work_of = [](const std::string& text) {
                std::istringstream in(text);
                const ScenarioParseResult r = parse_scenario(in);
                return r.ok ? malloy::scenario::step_work(r.scenario) : -1.0;
            };
            MALLOY_CHECK_EQ(work_of(""), 0.0);
            MALLOY_CHECK_EQ(work_of("type particles\nparticle 1 0.1 0 0 0 0\n"
                                    "particle 1 0.1 1 0 0 0\nparticle 1 0.1 2 0 0 0\n"),
                            6.0); // 3 pairs + 3 bodies
            MALLOY_CHECK_EQ(work_of("type charges\ncharge 1 1 0 0 0 0\ncharge 1 -1 1 0 0 0\n"
                                    "charge 1 1 2 0 0 0\ncharge 1 1 3 0 0 0\n"),
                            10.0); // 6 pairs + 4 bodies
            MALLOY_CHECK_EQ(work_of("type springs\nspring_body 1 0 0 0 0\nspring_body 1 1 0 0 0\n"
                                    "spring_body 1 2 0 0 0\nspring 0 1 1 1 0\nspring 1 2 1 1 0\n"),
                            5.0); // 3 bodies + 2 springs
        }

        // Issue #25: at most 1e5 reports, which bounds the output volume.
        {
            std::istringstream at("type nbody\nsteps 999990\noutput_every 10\nbody 1 0 0 0 0\n");
            MALLOY_CHECK_TRUE(parse_scenario(at).ok); // 99999 + the step-0 report
            std::istringstream over("type nbody\nsteps 1000000\noutput_every 10\nbody 1 0 0 0 0\n");
            const ScenarioParseResult r = parse_scenario(over);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_TRUE(r.error.find("reports") != std::string::npos);
            std::istringstream none("type nbody\nsteps 100000000\noutput_every 0\nbody 1 0 0 0 0\n");
            MALLOY_CHECK_TRUE(parse_scenario(none).ok); // only the step-0 report
        }

        // Issue #29: no line is read past 4096 characters, so one enormous line
        // cannot allocate without limit. A comment counts: it is still read.
        {
            std::istringstream at("type nbody\n#" + std::string(4095, 'c') +
                                  "\nbody 1 0 0 0 0\n");
            MALLOY_CHECK_TRUE(parse_scenario(at).ok);
            std::istringstream over("type nbody\n#" + std::string(4096, 'c') +
                                    "\nbody 1 0 0 0 0\n");
            const ScenarioParseResult r = parse_scenario(over);
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_TRUE(r.error.find("line 2:") != std::string::npos);
            MALLOY_CHECK_TRUE(r.error.find("4096") != std::string::npos);
            // The last line needs no newline, and a file of one line reads.
            std::istringstream last("type nbody\nbody 1 0 0 0 0");
            MALLOY_CHECK_TRUE(parse_scenario(last).ok);
        }

        // A few inputs that must be rejected outright.
        {
            std::istringstream missing("steps\n");
            MALLOY_CHECK_FALSE(parse_scenario(missing).ok);
            std::istringstream unknown_type("type wobble\n");
            MALLOY_CHECK_FALSE(parse_scenario(unknown_type).ok);
            std::istringstream truncated("type nbody\nbody 1 1\n");
            MALLOY_CHECK_FALSE(parse_scenario(truncated).ok);
        }

        // Robustness battery: each hostile input must RETURN a status rather than
        // crash, hang, or read out of bounds. Reaching the end of the loop is the
        // main assertion; a crash or infinite loop here would fail the test. Some
        // of these parse to a structurally valid scenario carrying values a world
        // would later reject, so their ok flag is deliberately not asserted. What
        // every outcome must have, whatever the verdict, is a message that is
        // safe to print (issue #28): printable ASCII, with no control byte.
        {
            const std::vector<std::string> hostile = {
                "",                                                // empty
                "\n\n\n",                                          // blank lines
                "   \t  ",                                         // whitespace only
                "# just a comment\n",                              // comment only
                "\x1a\x01 control bytes on a line",                // control chars
                "type nbody\nbody 1 nan inf 0 0 0 0\n",            // non-finite tokens
                "type nbody\ndt 1e400\nsteps 1\nbody 1 1 0 0 0 0\n", // overflow to inf
                "type springs\nspring 999999999999 0 1 1 0\n",     // out-of-range index
                "type nbody\ndt -1\nsteps 1\n",                    // non-positive dt
                "type nbody\n" + std::string(200000, '9') + "\n",  // one enormous token
                std::string(100000, 'x') + "\n",                   // enormous unknown key
                "type\ntype\ntype\n",                              // repeated bare key
                "TYPE NBODY\n",                                    // wrong case
                "type nbody\nbody 1 1 0 0 0 0 extra\n",            // trailing token
                // Issue #28: terminal escape sequences in each place a token is
                // repeated, NUL bytes, a no-break space, and the longest line
                // allowed, all ESC.
                "\x1b]0;PWNED-TITLE\x07\x1b[2J\x1b[1;31mFAKE-OK\x1b[0m\n",
                "dt 0.001 \x1b[1A\x1b[2Kscenario-ok\n",
                "type \x1b[2Kparticles\n",
                std::string("ke\0y\0 1\n", 8),
                "dt\xC2\xA0" "0.001\n",
                std::string(4096, '\x1b') + "\n",
                // Issue #32: byte-order marks, UTF-8 and UTF-16.
                "\xEF\xBB\xBF" "dt 0.5\n",
                "\xEF\xBB\xBF\xEF\xBB\xBF" "dt 0.5\n",
                std::string("\xFF\xFE" "t\0y\0p\0e\0\n\0", 12),
            };
            std::size_t returned = 0;
            for (const std::string& s : hostile)
            {
                std::istringstream in(s);
                const ScenarioParseResult r = parse_scenario(in);
                MALLOY_CHECK_TRUE(printable_ascii(r.error));
                ++returned;
            }
            MALLOY_CHECK_TRUE(returned == hostile.size());
        }

        // Issue #28: an error that repeats a token from the file escapes it, so
        // the file cannot write terminal control sequences through the one
        // channel that says it is bad. Every byte value goes through each of the
        // three places a token is repeated, inside the token `k<byte>y`, and the
        // message must be printable ASCII and show the byte exactly as
        // `rendered` says. A byte that ends a line or a token (a newline, a
        // comment, whitespace as the extractor's locale defines it) cannot be
        // inside one and is skipped. Compared with == rather than
        // MALLOY_CHECK_EQ, which would print raw bytes into the log on a failure.
        {
            const std::locale locale; // the global one, which istringstream uses
            for (int b = 0; b < 256; ++b)
            {
                const char c = static_cast<char>(b);
                if (c == '\n' || c == '#' || std::isspace(c, locale))
                {
                    continue;
                }
                const std::string token = std::string("k") + c + "y";
                const std::string shown =
                    "'k" + rendered(static_cast<unsigned char>(b)) + "y'";
                // M39 added two more places a token is repeated: an unknown wall
                // name and a seed that is not a whole number.
                const std::string inputs[5] = {
                    token + " 1\n", "type " + token + "\n", "dt 0.5 " + token + "\n",
                    "type particles\nwall_temperature " + token + " 1\n",
                    "type particles3d\nseed " + token + "\n"};
                const std::string expected[5] = {
                    "line 1: unknown key " + shown, "line 1: unknown type " + shown,
                    "line 1: unexpected extra token " + shown,
                    "line 2: unknown wall " + shown + " (walls are min_x, max_x, min_y and max_y)",
                    "line 2: seed must be a whole number from 0 to 18446744073709551615, not " +
                        shown};
                for (int site = 0; site < 5; ++site)
                {
                    std::istringstream in(inputs[site]);
                    const ScenarioParseResult r = parse_scenario(in);
                    if (r.error != expected[site])
                    {
                        std::cerr << "byte " << b << ", site " << site << '\n';
                    }
                    MALLOY_CHECK_FALSE(r.ok);
                    MALLOY_CHECK_TRUE(printable_ascii(r.error));
                    MALLOY_CHECK_TRUE(r.error == expected[site]);
                }
            }
        }

        // The issue's three reproductions, one for each place a token is
        // repeated: a retitle, clear and fake success line in a key, a cursor-up
        // and erase-line in a trailing token, and an erase-line in a type.
        {
            std::istringstream key("\x1b]0;PWNED-TITLE\x07\x1b[2J\x1b[1;31mFAKE-OK\x1b[0m\n");
            MALLOY_CHECK_TRUE(parse_scenario(key).error ==
                              "line 1: unknown key "
                              "'\\x1b]0;PWNED-TITLE\\x07\\x1b[2J\\x1b[1;31mFAKE-OK\\x1b[0m'");
            std::istringstream extra("dt 0.001 \x1b[1A\x1b[2Kscenario-ok\n");
            MALLOY_CHECK_TRUE(parse_scenario(extra).error ==
                              "line 1: unexpected extra token '\\x1b[1A\\x1b[2Kscenario-ok'");
            std::istringstream type("type \x1b[2Kparticles\n");
            MALLOY_CHECK_TRUE(parse_scenario(type).error ==
                              "line 1: unknown type '\\x1b[2Kparticles'");
        }

        // The echo is cut to the token's first 64 bytes, counted before
        // escaping, with the full length after the closing quote. The token's
        // letter differs from one position to the next, so the test shows WHICH
        // bytes were kept, not only how many.
        {
            const auto error_for = [](const std::string& key) {
                std::istringstream in(key + "\n");
                return parse_scenario(in).error;
            };
            std::string token;
            for (int i = 0; i < 100; ++i)
            {
                token += static_cast<char>('a' + i % 26);
            }
            const std::string first64 = token.substr(0, 64);
            MALLOY_CHECK_TRUE(error_for(token) ==
                              "line 1: unknown key '" + first64 + "'... (100 bytes)");
            // Exactly 64 is shown whole and unmarked; 65 is cut.
            MALLOY_CHECK_TRUE(error_for(first64) == "line 1: unknown key '" + first64 + "'");
            MALLOY_CHECK_TRUE(error_for(token.substr(0, 65)) ==
                              "line 1: unknown key '" + first64 + "'... (65 bytes)");
            // An ESC as the 64th byte is shown whole, as an escape, though its
            // four characters take the echo past 64; the 65th byte is not shown.
            MALLOY_CHECK_TRUE(error_for(std::string(63, 'k') + "\x1b" + "Z") ==
                              "line 1: unknown key '" + std::string(63, 'k') +
                                  "\\x1b'... (65 bytes)");
            // The longest token a line allows, all ESC: the message stays short.
            const std::string longest = error_for(std::string(4096, '\x1b'));
            MALLOY_CHECK_TRUE(printable_ascii(longest));
            MALLOY_CHECK_EQ(longest.size(), std::string("line 1: unknown key '").size() +
                                                std::size_t{64 * 4} +
                                                std::string("'... (4096 bytes)").size());
        }
    }

    // --- Issue #32: a UTF-8 byte-order mark at the very start of the input is
    //     skipped, so a file saved by a Windows editor parses exactly as it
    //     would without one. A UTF-16 file is refused by name. The mark anywhere
    //     else is ordinary text, which the error shows escaped. ---
    {
        const std::string bom = "\xEF\xBB\xBF";
        const auto parse_text = [](const std::string& text) {
            std::istringstream in(text);
            return parse_scenario(in);
        };

        // A value on line 1, right after the mark, and two unequal bodies: line 1
        // is neither dropped nor cut short, and the rest parses as written.
        {
            const ScenarioParseResult r =
                parse_text(bom + "dt 0.25\nsteps 7\nbody 2 1 0 0 1\nbody 3 -1 0.5 0 -1\n");
            MALLOY_CHECK_TRUE(r.ok);
            MALLOY_CHECK_NEAR(r.scenario.simulation.dt, 0.25, eps);
            MALLOY_CHECK_EQ(r.scenario.steps, 7);
            MALLOY_CHECK_EQ(r.scenario.bodies.size(), std::size_t{2});
            MALLOY_CHECK_NEAR(r.scenario.bodies[1].mass, 3.0, eps);
            MALLOY_CHECK_NEAR(r.scenario.bodies[1].position.y, 0.5, eps);
        }

        // The issue's case, the mark before a comment line, which left the mark
        // alone as the first key; with Windows line endings, as Notepad saves.
        {
            const ScenarioParseResult r =
                parse_text(bom + "# a comment\r\ntype rigid\r\nrigid_static 0.5 0 0 0\r\n");
            MALLOY_CHECK_TRUE(r.ok);
            MALLOY_CHECK_TRUE(r.scenario.type == malloy::scenario::ScenarioType::Rigid);
        }

        // Line numbers are unchanged by the mark.
        {
            const ScenarioParseResult r = parse_text(bom + "dt 0.5\nsteps 3\nsteps -1\n");
            MALLOY_CHECK_FALSE(r.ok);
            MALLOY_CHECK_TRUE(r.error.find("line 3:") != std::string::npos);
        }

        // A mark alone is an empty file, and the parser accepts it as one (the
        // app refuses it for having nothing to simulate, issue #31).
        {
            const ScenarioParseResult only = parse_text(bom);
            MALLOY_CHECK_TRUE(only.ok);
            MALLOY_CHECK_EQ(malloy::scenario::entity_count(only.scenario), std::size_t{0});
            MALLOY_CHECK_TRUE(parse_text(bom + "\n").ok);
        }

        // Only a whole mark, only one, and only at the very start is skipped.
        {
            const ScenarioParseResult partial = parse_text("\xEF\xBB" "dt 0.5\n");
            MALLOY_CHECK_FALSE(partial.ok);
            MALLOY_CHECK_TRUE(partial.error == "line 1: unknown key '\\xef\\xbbdt'");

            const ScenarioParseResult twice = parse_text(bom + bom + "dt 0.5\n");
            MALLOY_CHECK_FALSE(twice.ok);
            MALLOY_CHECK_TRUE(twice.error == "line 1: unknown key '\\xef\\xbb\\xbfdt'");

            const ScenarioParseResult later = parse_text("dt 0.5\n" + bom + "steps 3\n");
            MALLOY_CHECK_FALSE(later.ok);
            MALLOY_CHECK_TRUE(later.error == "line 2: unknown key '\\xef\\xbb\\xbfsteps'");

            const ScenarioParseResult indented = parse_text(" " + bom + "dt 0.5\n");
            MALLOY_CHECK_FALSE(indented.ok);
            MALLOY_CHECK_TRUE(indented.error == "line 1: unknown key '\\xef\\xbb\\xbfdt'");
        }

        // UTF-16 in either byte order is refused by name, rather than reported as
        // a key full of NUL bytes.
        {
            const auto utf16 = [](const std::string& ascii, bool little_endian) {
                std::string out = little_endian ? "\xFF\xFE" : "\xFE\xFF";
                for (const char c : ascii)
                {
                    out += little_endian ? c : '\0';
                    out += little_endian ? '\0' : c;
                }
                return out;
            };
            for (const bool little_endian : {true, false})
            {
                const ScenarioParseResult r =
                    parse_text(utf16("type nbody\nbody 1 0 0 0 0\n", little_endian));
                MALLOY_CHECK_FALSE(r.ok);
                MALLOY_CHECK_TRUE(r.error.find("line 1:") != std::string::npos);
                MALLOY_CHECK_TRUE(r.error.find("UTF-8 or ASCII") != std::string::npos);
                MALLOY_CHECK_TRUE(r.error.find("UTF-16") != std::string::npos);
            }

            // The check comes before the length check: a first line with no
            // 0x0A byte in its first 4096 still gets the encoding error.
            const ScenarioParseResult long_line = parse_text("\xFF\xFE" + std::string(5000, 'a'));
            MALLOY_CHECK_FALSE(long_line.ok);
            MALLOY_CHECK_TRUE(long_line.error.find("UTF-8 or ASCII") != std::string::npos);

            // Both bytes must match, and only at the very start.
            const ScenarioParseResult near_miss = parse_text("\xFF\xFF" "dt 0.5\n");
            MALLOY_CHECK_FALSE(near_miss.ok);
            MALLOY_CHECK_TRUE(near_miss.error == "line 1: unknown key '\\xff\\xffdt'");
            const ScenarioParseResult second_line = parse_text("dt 0.5\n\xFF\xFE" "steps 3\n");
            MALLOY_CHECK_FALSE(second_line.ok);
            MALLOY_CHECK_TRUE(second_line.error == "line 2: unknown key '\\xff\\xfesteps'");
        }

#ifdef MALLOY_SCENARIO_DIR
        // The issue's reproduction on a shipped template: with the mark in
        // front, two_body.scn parses exactly as it does without one.
        {
            const std::string path = std::string(MALLOY_SCENARIO_DIR) + "/two_body.scn";
            std::ifstream file(path, std::ios::binary);
            const std::string text((std::istreambuf_iterator<char>(file)),
                                   std::istreambuf_iterator<char>());
            MALLOY_CHECK_FALSE(text.empty());
            const ScenarioParseResult plain = malloy::scenario::parse_scenario_file(path);
            const ScenarioParseResult marked = parse_text(bom + text);
            MALLOY_CHECK_TRUE(plain.ok);
            MALLOY_CHECK_TRUE(marked.ok);
            MALLOY_CHECK_TRUE(marked.scenario.type == plain.scenario.type);
            MALLOY_CHECK_NEAR(marked.scenario.simulation.dt, plain.scenario.simulation.dt, eps);
            MALLOY_CHECK_EQ(marked.scenario.steps, plain.scenario.steps);
            MALLOY_CHECK_EQ(marked.scenario.output_every, plain.scenario.output_every);
            MALLOY_CHECK_NEAR(marked.scenario.nbody_settings.g, plain.scenario.nbody_settings.g, eps);
            MALLOY_CHECK_NEAR(marked.scenario.nbody_settings.softening,
                              plain.scenario.nbody_settings.softening, eps);
            MALLOY_CHECK_EQ(marked.scenario.bodies.size(), plain.scenario.bodies.size());
            MALLOY_CHECK_TRUE(plain.scenario.bodies.size() >= 2);
            for (std::size_t i = 0; i < plain.scenario.bodies.size(); ++i)
            {
                MALLOY_CHECK_NEAR(marked.scenario.bodies[i].mass, plain.scenario.bodies[i].mass, eps);
                MALLOY_CHECK_VEC2_NEAR(marked.scenario.bodies[i].position,
                                       plain.scenario.bodies[i].position, eps);
                MALLOY_CHECK_VEC2_NEAR(marked.scenario.bodies[i].velocity,
                                       plain.scenario.bodies[i].velocity, eps);
            }
        }
#endif
    }

    std::cout << "malloy_scenario_tests passed\n";
    return 0;
}
