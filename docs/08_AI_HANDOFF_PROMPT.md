# 08 — AI Handoff Prompt

Use this prompt when asking Claude or another AI assistant to continue work on MalloySim.

```text
You are helping me build MalloySim, a simulation-first C++20 project.

Before doing anything, read:

- README.md
- CLAUDE.md
- docs/00_START_HERE.md
- docs/01_V3_ARCHITECTURE_DECISION.md
- docs/02_MILESTONE_ROADMAP_M1_M5.md
- docs/03_MODULE_BOUNDARIES.md
- docs/04_NUMERIC_AND_PHYSICS_CONVENTIONS.md
- docs/05_TESTING_STRATEGY.md
- docs/06_CMAKE_AND_VSCODE_WORKFLOW.md

Current milestone:

M1: CMake skeleton and smoke test.

Do not jump to M2 unless M1 is green and committed.

Do not add rendering, ECS, collision, rigid bodies, 3D, quantum, package managers, Catch2, plugins, scripting, editor systems, or config files.

Use C++20, CMake, MSVC, VS Code CMake Tools, CTest, and tiny custom check macros.

The immediate build commands are:

cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug

Your task is to work only on the requested milestone.

If I ask for something outside the current milestone, tell me what milestone it belongs to and suggest the smallest correct next step.
```
