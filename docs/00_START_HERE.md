# 00 — Start Here

**MalloySim C++** is a simulation-first C++20 project. The locked **M1-M5
roadmap is complete**, plus the first post-M5 milestone (M6: N-body
diagnostics): the math, time, sim-core, and N-body libraries are built and the
terminal N-body demo runs. See `README.md` and `CHANGELOG.md`.

## Read order

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

## Build and test

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

All five test executables should pass, and
`out/build/windows-msvc-debug/Debug/malloy_nbody_terminal.exe` should run to
completion.

## What's next

The locked M1-M5 roadmap is done. Post-M5 work is gated and taken one milestone
at a time -- see `docs/07_POST_M5_ROADMAP.md`. Do not expand scope immediately
after M5 (`docs/09_MISTAKES_TO_AVOID.md`, #10).
