// libFuzzer entry point for the scenario parser, the one untrusted input this
// project takes (a scenario file named on the command line). Its sibling,
// fuzz/app_fuzzer.cpp, runs the whole app behind the parser (issue #30). Neither
// is built by the default MSVC presets: they need clang++ with libFuzzer and
// the sanitizers. Build and run them with (the first command is one line,
// wrapped here):
//
//   cmake -S . -B out/fuzz -G Ninja -DCMAKE_BUILD_TYPE=Release
//         -DCMAKE_CXX_COMPILER=clang++ -DMALLOY_BUILD_FUZZERS=ON
//         -DMALLOY_BUILD_TESTS=OFF -DMALLOY_BUILD_APPS=OFF
//   cmake --build out/fuzz
//   cmake -E make_directory out/fuzz/parse_corpus out/fuzz/app_corpus
//   out/fuzz/malloy_parse_scenario_fuzzer out/fuzz/parse_corpus scenarios
//         -max_total_time=60 -artifact_prefix=out/fuzz/
//   out/fuzz/malloy_app_fuzzer out/fuzz/app_corpus scenarios
//         -max_total_time=60 -timeout=10 -artifact_prefix=out/fuzz/
//
// The FIRST directory is the corpus libFuzzer writes new inputs into; the ones
// after it only seed it. Never put `scenarios` first, or the fuzzer writes into
// the shipped templates. A reproducer for anything found is written as
// out/fuzz/crash-<hash> (without -artifact_prefix, into the current directory).
//
// CI runs this recipe on Linux on every push and pull request to main, for 30 s
// per target (.github/workflows/ci.yml). On Ubuntu, libFuzzer and the sanitizer
// runtimes come in a package of their own, libclang-rt-<version>-dev.
//
// On Windows, configure and build from a shell with the LLVM bin folder and
// Ninja on PATH (clang++ finds the MSVC headers and libraries itself). ASan's
// runtime there is a DLL, so before running a fuzzer also put its folder on
// PATH: <clang resource dir>\lib\windows, which for clang 22 installed in the
// default place is C:\Program Files\LLVM\lib\clang\22\lib\windows (PowerShell:
// $env:PATH = "C:\Program Files\LLVM\lib\clang\22\lib\windows;" + $env:PATH).
//
// MALLOY_BUILD_FUZZERS builds EVERY target with coverage and with ASan and
// UBSan (issue #27), so libFuzzer is steered by the parser's own branches and a
// crash, an out-of-bounds read or write, or undefined behaviour anywhere in it
// stops the run with a reproducer. UBSan is built not to recover, so undefined
// behaviour stops the run rather than scrolling past; to see where it happened,
// replay the reproducer (pass its path instead of the corpus directories) with
// UBSAN_OPTIONS=print_stacktrace=1 in the environment. The parser must never do
// more than return a status (docs/04: never throw), so the only accepted
// outcomes are a clean parse or a reported error, and the checks below hold
// whatever the input.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <string>

#include <malloy/scenario/scenario.hpp>

namespace
{
// A violated invariant is a real defect, so say which one and trip the fuzzer,
// which then saves the input as a reproducer.
[[noreturn]] void fail(const char* why)
{
    std::fprintf(stderr, "parse_scenario_fuzzer: %s\n", why);
    __builtin_trap();
}
} // namespace

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    std::string text(reinterpret_cast<const char*>(data), size);
    std::istringstream input(text);

    const malloy::scenario::ScenarioParseResult result =
        malloy::scenario::parse_scenario(input);

    if (result.ok)
    {
        const malloy::scenario::Scenario& s = result.scenario;
        if (s.steps < 0 || s.output_every < 0)
        {
            fail("a successful parse yielded a negative count");
        }

        // Every resource cap scenario.hpp documents holds for a file the parser
        // accepted (issue #66). The counts are recomputed here from the
        // documented formulas, not taken from the parser, so a slip in the
        // parser's own arithmetic shows up as a disagreement. The weighted work
        // and the contact point count come from the parser's exported functions
        // (issue #84), each held against an independent lower bound.
        if (s.steps > 100000000)
        {
            fail("a successful parse exceeds the step cap (1e8)");
        }
        const std::size_t entities = malloy::scenario::entity_count(s);
        if (entities > 100000)
        {
            fail("a successful parse exceeds the entity cap (1e5)");
        }
        const std::size_t planes2 = s.rigid_settings.ground.size();
        const std::size_t planes3 = s.rigid3d_settings.ground.size();
        if (planes2 > 1000 || planes3 > 1000)
        {
            fail("a successful parse exceeds the ground plane cap (1000)");
        }
        const long long reports =
            s.output_every > 0 ? static_cast<long long>(s.steps) / s.output_every + 1 : 1;
        if (reports > 100000)
        {
            fail("a successful parse exceeds the report cap (1e5)");
        }
        // The work budget, through the parser's own step_work (issue #84), so
        // this checks the cap the parser actually enforces rather than a copy
        // of its formula, which went stale once when issue #46 weighted it.
        const double per_step = malloy::scenario::step_work(s);
        const double work =
            (static_cast<double>(s.steps) + static_cast<double>(reports)) * per_step;
        if (work > 1e11)
        {
            fail("a successful parse exceeds the work budget (1e11)");
        }
        // And an independent lower bound on it, from the documented formula
        // with every weight at 1: N + springs for the linear springs domains,
        // and N(N-1)/2 pairs plus N per body and N per ground plane otherwise,
        // where N counts every entity but the springs. Every weight the parser
        // applies is at least 1, so a step_work below this (a weight under 1,
        // a dropped term, a NaN) is a slip in the parser's own arithmetic.
        const std::size_t springs =
            s.spring_network.springs().size() + s.spring_network3d.springs().size();
        const double n = static_cast<double>(entities - springs);
        const bool linear = s.type == malloy::scenario::ScenarioType::Springs ||
                            s.type == malloy::scenario::ScenarioType::Springs3D;
        const double unweighted =
            linear ? n + static_cast<double>(springs)
                   : n * (n - 1.0) / 2.0 + n * static_cast<double>(planes2 + planes3) + n;
        if (!(per_step >= unweighted))
        {
            fail("the parser's work per step is below the unweighted pair count");
        }

        // The contact point cap (M40), through the parser's own count, and an
        // independent lower bound on that count: every pair of rigid bodies
        // that collide (a box, or a radius above 0), and every such body
        // against every ground plane, holds at least one point.
        const std::uint64_t points = malloy::scenario::potential_contact_points(s);
        if (points > 250000)
        {
            fail("a successful parse exceeds the contact point cap (250000)");
        }
        std::uint64_t colliders = 0;
        for (const malloy::rigid::RigidBody2D& body : s.rigid_bodies)
        {
            if ((body.half_extents.x > 0.0 && body.half_extents.y > 0.0) || body.radius > 0.0)
            {
                ++colliders;
            }
        }
        for (const malloy::rigid::RigidBody3D& body : s.rigid_bodies3d)
        {
            if ((body.half_extents.x > 0.0 && body.half_extents.y > 0.0 &&
                 body.half_extents.z > 0.0) ||
                body.radius > 0.0)
            {
                ++colliders;
            }
        }
        const std::uint64_t collider_pairs = colliders > 0 ? colliders * (colliders - 1) / 2 : 0;
        if (points < collider_pairs + colliders * (planes2 + planes3))
        {
            fail("the parser's contact point count is below one point per colliding pair "
                 "and plane");
        }

        // No line, counted in bytes up to its newline, may pass 4096.
        std::size_t line_length = 0;
        for (std::size_t i = 0; i < size; ++i)
        {
            line_length = (data[i] == '\n') ? 0 : line_length + 1;
            if (line_length > 4096)
            {
                fail("a successful parse accepted a line longer than 4096 bytes");
            }
        }
    }
    else
    {
        // Every refusal says why, and says it in printable ASCII, however hostile
        // the file was: the message goes straight to a terminal (issue #28).
        if (result.error.empty())
        {
            fail("a failed parse has an empty error message");
        }
        for (const char c : result.error)
        {
            const auto byte = static_cast<unsigned char>(c);
            if (byte < 0x20 || byte > 0x7e)
            {
                fail("a parse error message holds a byte outside printable ASCII");
            }
        }
    }
    return 0;
}
