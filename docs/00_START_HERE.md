# 00 - Start Here

**MalloySim C++** is an all-in-one science physics simulation in C++20, built
simulation-first and terminal-first. It grows one finished physics domain at a
time, each shipping tested scenario templates, with no engine kernel and no
simulation base class (`docs/decisions/0006-multi-domain-dispatch.md`).

M1-M9 are complete: the math, time, sim-core, N-body, scenario, ASCII, and
collide libraries are built and the terminal N-body demo runs. Gravity is the
first finished domain. See `README.md` and `CHANGELOG.md`.

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
11. `docs/09_MISTAKES_TO_AVOID.md`
12. `docs/10_CHANGELOG_TEMPLATE.md`
13. `docs/decisions/` (ADRs 0001-0006)

## Build and test

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

All eight test executables should pass, and
`out/build/windows-msvc-debug/Debug/malloy_nbody_terminal.exe` should run to
completion.

## What's next

M1-M9 are done and gravity is the first finished domain. The active track is
classical mechanics depth: collision geometry shipped in M9, and rigid bodies,
ballistics, and springs remain. Work stays gated one milestone at a time; see
`docs/07_POST_M5_ROADMAP.md`.

The all-in-one goal is not a license to build ahead. Breadth is earned by
finishing domains, not by scaffolding for them
(`docs/09_MISTAKES_TO_AVOID.md`, #10).
