# 00 — Start Here

This repository is the start of **MalloySim C++**, a simulation-first C++20 project.

Read order:

1. `README.md`
2. `CLAUDE.md`
3. `docs/01_V3_ARCHITECTURE_DECISION.md`
4. `docs/02_MILESTONE_ROADMAP_M1_M5.md`
5. `docs/03_MODULE_BOUNDARIES.md`
6. `docs/04_NUMERIC_AND_PHYSICS_CONVENTIONS.md`
7. `docs/05_TESTING_STRATEGY.md`
8. `docs/06_CMAKE_AND_VSCODE_WORKFLOW.md`
9. `docs/07_POST_M5_ROADMAP.md`
10. `docs/08_AI_HANDOFF_PROMPT.md`

## Immediate goal

Do not write physics first.

The immediate goal is **M1**:

```text
CMake configures.
MSVC builds.
CTest runs.
The terminal app prints one line.
Commit green.
```

## M1 commands

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

## Commit M1 only after the commands pass

```powershell
git add .
git commit -m "M1: CMake skeleton and smoke test"
```
