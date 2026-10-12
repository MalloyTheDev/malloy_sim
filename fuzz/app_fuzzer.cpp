// libFuzzer entry point for the terminal app's whole path, not just the parser
// (issue #30). The parser fuzzer stops at parse_scenario; everything the app
// does after it, building the declared world, validate(), stepping, and the
// diagnostics and ASCII view it prints, was exercised by nothing, and the
// empty-charges crash (issue #24) lived exactly there. This target runs the real
// app on every input. Build and run it with the recipe in
// fuzz/parse_scenario_fuzzer.cpp; it is `malloy_app_fuzzer` there.
//
// The app is apps/nbody_terminal/main.cpp itself, compiled into this file
// unchanged except for two names: its main() becomes malloy_app_main(), and the
// one call that reads the scenario file, parse_scenario_file(), is answered from
// the fuzz input in memory. That call is only "open the file in binary and
// parse it", and the input is parsed with the same parse_scenario, so the app
// sees exactly what it would have read from disk. Writing a real file per input
// was measured instead: about 5 ms of file I/O each on Windows, which is most of
// the run. Should main.cpp ever stop reading through parse_scenario_file, the
// harness says so and stops, rather than quietly fuzzing nothing.
//
// Every input must finish in milliseconds, so the run is bounded here, on top
// of the parser's own caps. A scenario with more than 64 entities is skipped (it
// still joins the corpus for the parser coverage it found). Any other scenario
// runs, but only for its first steps: as many as fit in 1e6 work units, counted
// by the parser's own step_work (issue #84: this file used to carry a copy of
// the formula, which went stale when issue #46 weighted it), and in 32 reports.
// A step costlier than 1e6 units runs none, only its step-0 report. 1e6 units is
// at most about 10 ms of Release work; it was 1e4 until M40, whose worst-case
// rigid weights charge one pair of boxes 45860 units a step, so 1e4 stepped no
// scenario holding a box at all. Cutting a run
// short rather than skipping it means every shipped template, all of which run
// for thousands of steps, still runs its opening.
//
// The app's output is thrown away through a stream buffer that keeps the
// streams good, so all of its formatting still runs. What it promises, checked
// on every input, is:
//
//   - it returns 0 or 1;
//   - it returns 1 on a file the parser refuses, and on one that defines
//     nothing to simulate (issue #31);
//   - it never fails silently: a return of 1 wrote a message to stderr;
//   - what it writes to stderr is printable ASCII, one message per line, so a
//     hostile file cannot reach the terminal through it (issue #28).

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <ios>
#include <iostream>
#include <optional>
#include <sstream>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

#include <malloy/scenario/scenario.hpp>

// The stand-in for parse_scenario_file, declared before main.cpp is read so
// the renamed call in it has something to bind to. Every header main.cpp uses
// that could see the renames is already included above, so the two macros
// below touch main.cpp alone.
namespace malloy::scenario
{
ScenarioParseResult fuzz_parse_scenario_file(const std::string& path);
}

#define main malloy_app_main
#define parse_scenario_file fuzz_parse_scenario_file
#include "../apps/nbody_terminal/main.cpp"
#undef parse_scenario_file
#undef main

namespace malloy_fuzz
{
// The bounds on one run, described above.
constexpr std::size_t max_entities = 64;
constexpr double max_work = 1e6;
constexpr long long max_reports = 32;

// The parse of the current input, handed to the app by the stand-in, and
// whether the app asked for it.
const malloy::scenario::ScenarioParseResult* current = nullptr;
bool intercepted = false;

[[noreturn]] void fail(const char* why)
{
    std::fprintf(stderr, "app_fuzzer: %s\n", why);
    __builtin_trap();
}

// Stands in for the terminal. Discards what the app writes, but keeps the
// stream good, so every formatting call still runs, and remembers how much was
// written and whether any of it was other than printable ASCII and newlines.
class Sink : public std::streambuf
{
public:
    std::size_t written{0};
    bool printable{true};

protected:
    int_type overflow(int_type c) override
    {
        if (!traits_type::eq_int_type(c, traits_type::eof()))
        {
            note(traits_type::to_char_type(c));
        }
        return traits_type::not_eof(c);
    }

    std::streamsize xsputn(const char* text, std::streamsize count) override
    {
        for (std::streamsize i = 0; i < count; ++i)
        {
            note(text[i]);
        }
        return count;
    }

private:
    void note(char c)
    {
        ++written;
        const auto byte = static_cast<unsigned char>(c);
        if (c != '\n' && (byte < 0x20 || byte > 0x7e))
        {
            printable = false;
        }
    }
};

// How many of its steps a parsed scenario may run here: as many as fit in the
// work and report bounds, never more than it asked for. The work a step does is
// the parser's own figure (issue #84), so the bound cannot drift from the
// budget the parser enforces.
int bounded_steps(const malloy::scenario::Scenario& s)
{
    long long allowed =
        static_cast<long long>(max_work / std::max(1.0, malloy::scenario::step_work(s)));
    if (s.output_every > 0)
    {
        allowed = std::min(allowed,
                           (max_reports - 1) * static_cast<long long>(s.output_every));
    }
    return static_cast<int>(std::min(static_cast<long long>(s.steps), allowed));
}
} // namespace malloy_fuzz

malloy::scenario::ScenarioParseResult
malloy::scenario::fuzz_parse_scenario_file(const std::string&)
{
    malloy_fuzz::intercepted = true;
    return *malloy_fuzz::current;
}

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    std::istringstream input(std::string(reinterpret_cast<const char*>(data), size));
    malloy::scenario::ScenarioParseResult parsed = malloy::scenario::parse_scenario(input);

    std::size_t entities = 0;
    if (parsed.ok)
    {
        entities = malloy::scenario::entity_count(parsed.scenario);
        if (entities > malloy_fuzz::max_entities)
        {
            return 0;
        }
        parsed.scenario.steps = malloy_fuzz::bounded_steps(parsed.scenario);
    }

    malloy_fuzz::current = &parsed;
    malloy_fuzz::intercepted = false;

    // Run the app with its two streams pointed at sinks, then put them back,
    // with the number format the app sets on std::cout, so no input inherits
    // state from the one before.
    malloy_fuzz::Sink out;
    malloy_fuzz::Sink err;
    const std::ios_base::fmtflags flags = std::cout.flags();
    const std::streamsize precision = std::cout.precision();
    std::streambuf* const real_out = std::cout.rdbuf(&out);
    std::streambuf* const real_err = std::cerr.rdbuf(&err);

    char program[] = "malloy_nbody_terminal";
    char path[] = "fuzz-input.scn";
    char* argv[] = {program, path, nullptr};
    const int rc = malloy_app_main(2, argv);

    std::cout.rdbuf(real_out);
    std::cerr.rdbuf(real_err);
    std::cout.flags(flags);
    std::cout.precision(precision);

    if (!malloy_fuzz::intercepted)
    {
        malloy_fuzz::fail("main.cpp no longer reads its scenario through "
                          "parse_scenario_file, so this harness no longer feeds it; "
                          "update fuzz/app_fuzzer.cpp");
    }
    if (rc != 0 && rc != 1)
    {
        malloy_fuzz::fail("the app returned something other than 0 or 1");
    }
    if (!parsed.ok && rc != 1)
    {
        malloy_fuzz::fail("the app ran a file the parser refused");
    }
    if (parsed.ok && entities == 0 && rc != 1)
    {
        malloy_fuzz::fail("the app ran a scenario that defines nothing (issue #31)");
    }
    if (rc != 0 && err.written == 0)
    {
        malloy_fuzz::fail("the app failed without a message on stderr");
    }
    if (!err.printable)
    {
        malloy_fuzz::fail("the app wrote a byte outside printable ASCII to stderr");
    }
    return 0;
}
