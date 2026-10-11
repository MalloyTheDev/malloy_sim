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
        if (result.scenario.steps < 0 || result.scenario.output_every < 0)
        {
            fail("a successful parse yielded a negative count");
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
