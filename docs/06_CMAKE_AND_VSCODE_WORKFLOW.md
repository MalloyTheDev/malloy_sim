# 06 — CMake and VS Code Workflow

## Required tools

System:

- Visual Studio Build Tools 2026 or Visual Studio Community 2026
- Workload: Desktop development with C++
- Git for Windows

VS Code extensions:

- CMake Tools
- C/C++

## CMake strategy

Use CMake 3.21+, one root `CMakeLists.txt` through M5, `CMakePresets.json` from M1, out-of-source build prevention, `target_compile_features(... cxx_std_20)`, `CMAKE_EXPORT_COMPILE_COMMANDS ON`, `/W4` on MSVC, and explicit `/fp:precise` on MSVC.

Avoid package managers, install/export rules, complex CMake functions, source globbing, nested CMake files before M5, and dependency fetching before M5.

## Build commands

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

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
