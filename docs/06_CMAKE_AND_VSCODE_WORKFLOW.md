# 06 - CMake and VS Code Workflow

## Required tools

System:

- Visual Studio Build Tools 2026 or Visual Studio Community 2026
- Workload: Desktop development with C++
- Git for Windows

VS Code extensions:

- CMake Tools
- C/C++

Optional, only for the fuzz targets: LLVM (`clang++` with libFuzzer and the
sanitizers) and Ninja.

## CMake strategy

Use CMake 3.21+ for the project itself. The shipped presets additionally
require CMake 4.2+, because the `Visual Studio 18 2026` generator was added in
4.2; on an older CMake, configure without a preset or add a
`CMakeUserPresets.json` naming a generator your version has.

What the build uses, all in force at M39:

- one root `CMakeLists.txt` (still the only one; no nested CMake files);
- `CMakePresets.json`, with `windows-msvc-debug` and `windows-msvc-release`
  configure, build and test presets;
- out-of-source build prevention (an in-source configure is a fatal error);
- `target_compile_features(... cxx_std_20)` on `malloy_project_options`;
- `CMAKE_EXPORT_COMPILE_COMMANDS ON`;
- on MSVC, `/W4 /WX /permissive- /fp:precise` on `malloy_project_warnings`;
- on GCC and Clang, `-Wall -Wextra -Wpedantic -Werror -ffp-contract=off`
  (warnings are errors there too since issue #66);
- four options: `MALLOY_BUILD_APPS` and `MALLOY_BUILD_TESTS` (both ON),
  `MALLOY_BUILD_FUZZERS` (OFF; see Fuzzing below), and `MALLOY_BUILD_VIEWER`
  (OFF in `CMakeLists.txt`, ON in both presets; M41). The viewer's configure
  step fetches raylib 6.0 with `FetchContent`, from its release archive and
  checked against a pinned SHA-256; it is declared SYSTEM and EXCLUDE_FROM_ALL,
  so its warnings never reach `/W4 /WX`, and built static, without audio and
  with Control Flow Guard. `FETCHCONTENT_SOURCE_DIR_RAYLIB` pointed at a
  verified copy builds offline. The viewer and the fuzzers cannot be on in the
  same build directory.

Avoid package managers, install/export rules, complex CMake functions, source
globbing, and dependency fetching until a milestone needs them. Two planned
milestones will: 3D rendering brings the first graphics dependency
([#80](https://github.com/MalloyTheDev/malloy_sim/issues/80)), and packaging the libraries for other projects brings
install and export rules ([#82](https://github.com/MalloyTheDev/malloy_sim/issues/82)).

## Build commands

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

The same three commands with `windows-msvc-release` give the optimized build.
CTest runs 46 tests: the 15 test executables, fifteen checks that run the
viewer (refusals, flags, and the headless report lines against the terminal
app's at three frame rates), plus sixteen checks that run the
terminal app itself and require both its exit code and its output (four inputs
it must refuse, one scenario path outside the ANSI code page that it must open,
one shipped template of every scenario type, and the force template's energy
column). The app checks are registered after the
test-executable count is taken, so they do not change it (`docs/05`).

CI (`.github/workflows/ci.yml`, issue #66) runs these commands on every push
and pull request to `main`: the Debug and Release presets with every CTest test
on the `windows-2025-vs2026` runner, and a 30-second libFuzzer run of each fuzz
target under ASan and UBSan on `ubuntu-26.04`. It has read-only permissions,
needs no secrets, and pins every action to a commit SHA. Branch protection on
`main` is not set.

## MSVC hardening

Issue #37. Every MSVC target gets these in every configuration:

- `/guard:cf` at compile and link (Control Flow Guard) and `/CETCOMPAT` at
  link (CET shadow-stack compatible), on `malloy_project_options`, which every
  library links PUBLIC, so every object file and executable is covered.
- `/WX` beside `/W4`, on `malloy_project_warnings`. The tree builds
  warning-free, so a new warning fails the build instead of accumulating.

Neither mitigation changes a result: all shipped templates produce
byte-identical output with and without them, at the same speed.

Two further hardenings were measured and are deliberately off, because the
Release pair loop is what the simulator spends its time on. Release build,
100 disks for 40000 steps (compute bound), median of 5 runs:

| Build | Time | Cost |
|---|---|---|
| none of the four | 581 ms | |
| `/guard:cf` + `/CETCOMPAT` (shipped) | 577 ms | 0% |
| `_MSVC_STL_HARDENING=1` | 699 ms | +20% |
| `/sdl` | 822 ms | +41% |
| all four | 964 ms | +66% |

- `_MSVC_STL_HARDENING=1` turns a standard-library precondition violation
  (`front()` of an empty vector, an index past the end) into a fail-fast.
  Debug builds already make the same checks through iterator debugging, in
  every test run; the parser bounds its input, and CI runs the fuzz targets
  below under ASan and UBSan on every push.
- `/sdl` would add only runtime checks: every warning it promotes (C4146,
  C4996 and the rest) is already an error under `/W4 /WX`.

To turn either on for one build, without editing anything:

```powershell
cmake --preset windows-msvc-release --fresh "-DCMAKE_CXX_FLAGS_INIT=/D_MSVC_STL_HARDENING=1"
cmake --preset windows-msvc-release --fresh "-DCMAKE_CXX_FLAGS_INIT=/sdl"
```

To turn one on permanently, add the define or the flag to the `if(MSVC)`
hardening block at the top of `CMakeLists.txt`; its comment gives the line.

To see the mitigations in a built executable (LLVM's `llvm-readobj`):

```powershell
llvm-readobj --file-headers --coff-load-config --coff-debug-directory out/build/windows-msvc-release/Release/malloy_nbody_terminal.exe
```

Look for `IMAGE_DLL_CHARACTERISTICS_GUARD_CF`, a nonzero
`GuardCFFunctionCount`, `CF_FUNCTION_TABLE_PRESENT`, and
`IMAGE_DLL_CHARACTERISTICS_EX_CET_COMPAT` (the CET bit is in the debug
directory, not the file header).

## Unicode paths

Issue #33. The terminal app embeds `apps/nbody_terminal/nbody_terminal.manifest`,
which makes UTF-8 the process code page, so a scenario path outside the system
code page (Greek, CJK, emoji) opens. CMake merges it into the linker's own
manifest because it is listed among the app's sources. It takes effect on
Windows 10 1903 and later; older systems ignore it and keep the system code
page. The app prints the path in its run title as UTF-8, so a console whose
output code page is not 65001 shows non-ASCII characters there garbled; the
run itself is unaffected. The `malloy_app_opens_unicode_path` test runs the app
on such a path.

## Fuzzing

The scenario file is the one untrusted input, so two libFuzzer targets in
`fuzz/` exercise it (issues #27 and #30):

- `malloy_parse_scenario_fuzzer` (`fuzz/parse_scenario_fuzzer.cpp`) parses each
  input and requires either a clean parse with no negative step or report count,
  or a non-empty error in printable ASCII.
- `malloy_app_fuzzer` (`fuzz/app_fuzzer.cpp`) compiles
  `apps/nbody_terminal/main.cpp` into itself, answers its one file read from the
  input in memory, and runs the whole app: refusals, dispatch, validation,
  stepping, diagnostics and ASCII frames. Each run is cut to its first steps
  (1e4 pair tests, 32 reports) and a scenario with more than 64 entities is
  skipped, so every input finishes in milliseconds. It requires a return of 0 or
  1, a message on stderr for every 1, and only printable ASCII on stderr.

How they are built (`CMakeLists.txt`, the `MALLOY_BUILD_FUZZERS` block):

- The option is OFF by default and the MSVC presets never set it.
- It needs `clang++`, the GNU-style driver. Configuring it with MSVC or
  `clang-cl` stops with an error.
- It instruments EVERY target, not only the two harnesses, because libFuzzer
  steers by coverage counters compiled into the code and the sanitizers see
  only code built with them: `-g -O1 -fno-omit-frame-pointer
  -fsanitize=fuzzer-no-link,address,undefined -fno-sanitize-recover=undefined`
  on `malloy_project_options`. Only the two fuzzers link libFuzzer itself
  (`-fsanitize=fuzzer`). UBSan does not recover, so undefined behaviour stops
  the run with a reproducer.
- On Windows (clang++ targeting the MSVC ABI) every target uses the static
  release CRT and the MSVC STL's ASan container annotations are off, which
  libFuzzer's prebuilt runtime requires. The cost: ASan cannot see a read past
  `size()` but inside the capacity.
- The fuzzers are not registered with CTest, so `ctest` never runs them and the
  test count is unchanged.

The recipe, from the header of `fuzz/parse_scenario_fuzzer.cpp` (the first
command is one line):

```text
cmake -S . -B out/fuzz -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++ -DMALLOY_BUILD_FUZZERS=ON -DMALLOY_BUILD_TESTS=OFF -DMALLOY_BUILD_APPS=OFF
cmake --build out/fuzz
cmake -E make_directory out/fuzz/parse_corpus out/fuzz/app_corpus
out/fuzz/malloy_parse_scenario_fuzzer out/fuzz/parse_corpus scenarios -max_total_time=60 -artifact_prefix=out/fuzz/
out/fuzz/malloy_app_fuzzer out/fuzz/app_corpus scenarios -max_total_time=60 -timeout=10 -artifact_prefix=out/fuzz/
```

- The FIRST directory is the corpus libFuzzer writes into; the ones after it
  only seed it. Never put `scenarios` first, or the fuzzer writes into the
  shipped templates.
- On Windows, configure and build from a shell with the LLVM `bin` folder and
  Ninja on `PATH`, and before running a fuzzer also put ASan's runtime folder,
  `<clang resource dir>\lib\windows`, on `PATH`.
- A finding is written as `out/fuzz/crash-<hash>`. Replay it by passing its
  path instead of the corpus directories, with `UBSAN_OPTIONS=print_stacktrace=1`
  set to see where undefined behaviour happened.

CI runs both fuzzers briefly on every push. A fuzzer binary left in
`out/fuzz` from an earlier build does not cover parsing added since, so rebuild
from the current tree before a local run. The parser fuzzer asserts every
documented resource cap on a successful parse (issue #66), though its copy of
the work formula predates the weighted budget of issue #46
([#84](https://github.com/MalloyTheDev/malloy_sim/issues/84)).

## VS Code workflow

1. Open the repository root.
2. CMake Tools should detect `CMakePresets.json`.
3. Select `windows-msvc-debug`.
4. Configure.
5. Build.
6. Run CTest.
7. Run the terminal target, `malloy_nbody_terminal`. With no argument it runs
   the built-in scenarios; with one positional path it runs that scenario file
   (there are no flags and no CLI parser).

## Build directory

The presets use:

```text
out/build/windows-msvc-debug     (app: Debug/malloy_nbody_terminal.exe)
out/build/windows-msvc-release   (app: Release/malloy_nbody_terminal.exe)
```

The fuzz recipe above uses `out/fuzz`. All of `out/` is ignored by git.

Do not commit personal `.vscode` configuration or `CMakeUserPresets.json`.
