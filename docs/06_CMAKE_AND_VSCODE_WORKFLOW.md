# 06 - CMake and VS Code Workflow

## Required tools

System:

- Visual Studio Build Tools 2026 or Visual Studio Community 2026
- Workload: Desktop development with C++
- Git for Windows

VS Code extensions:

- CMake Tools
- C/C++

## CMake strategy

Use CMake 3.21+ for the project itself. The shipped presets additionally
require CMake 4.2+, because the `Visual Studio 18 2026` generator was added in
4.2; on an older CMake, configure without a preset or add one naming a
generator your version has, one root `CMakeLists.txt` through M5, `CMakePresets.json` from M1, out-of-source build prevention, `target_compile_features(... cxx_std_20)`, `CMAKE_EXPORT_COMPILE_COMMANDS ON`, `/W4 /WX` on MSVC, and explicit `/fp:precise` on MSVC.

Avoid package managers, install/export rules, complex CMake functions, source globbing, nested CMake files before M5, and dependency fetching before M5.

## Build commands

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

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
  every test run; the parser is fuzzed under ASan and UBSan and bounds its
  input.
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

## VS Code workflow

1. Open the repository root.
2. CMake Tools should detect `CMakePresets.json`.
3. Select `windows-msvc-debug`.
4. Configure.
5. Build.
6. Run CTest.
7. Run the terminal target.

## Build directory

The preset uses:

```text
out/build/windows-msvc-debug
```

Do not commit personal `.vscode` configuration or `CMakeUserPresets.json`.
