// libFuzzer entry point for the scenario parser, the one untrusted input this
// project takes (a scenario file named on the command line). It is NOT built by
// the default MSVC presets; it needs a fuzzer- and sanitizer-capable toolchain
// (clang). Turn it on with:
//
//   cmake -S . -B out/fuzz -G Ninja \
//       -DCMAKE_CXX_COMPILER=clang++ -DMALLOY_BUILD_FUZZERS=ON \
//       -DMALLOY_BUILD_TESTS=OFF -DMALLOY_BUILD_APPS=OFF
//   cmake --build out/fuzz
//   out/fuzz/malloy_parse_scenario_fuzzer -max_total_time=60
//
// The CMake target adds -fsanitize=fuzzer,address,undefined, so a crash, a
// read/write out of bounds, or undefined behaviour reachable from any byte
// sequence is reported with a reproducer. The parser must never do more than
// return a status (docs/04: never throw), so the only accepted outcomes are a
// clean parse or a reported error.

#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>

#include <malloy/scenario/scenario.hpp>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    std::string text(reinterpret_cast<const char*>(data), size);
    std::istringstream input(text);

    const malloy::scenario::ScenarioParseResult result =
        malloy::scenario::parse_scenario(input);

    // Invariants the parser promises regardless of input. A violation is a real
    // defect, so trip the sanitizer rather than return quietly.
    if (result.ok)
    {
        if (result.scenario.steps < 0 || result.scenario.output_every < 0)
        {
            __builtin_trap(); // a successful parse never yields a negative count
        }
    }
    return 0;
}
