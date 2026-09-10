# 00 - Start Here

**MalloySim C++** is an all-in-one science physics simulation in C++20, built
simulation-first and terminal-first. It grows one finished physics domain at a
time, each shipping tested scenario templates, with no engine kernel and no
simulation base class (`docs/decisions/0006-multi-domain-dispatch.md`).

M1-M18 are complete: the math, time, sim-core, N-body, scenario, ASCII, collide,
particles, rigid, springs and charges libraries are built, and the terminal app
runs all five domains. Gravity, colliding particles, rigid bodies, spring
networks and charged particles are the finished domains. See `README.md` and `CHANGELOG.md`.

## Where this is going

MalloySim is 2D today and is intended to become a 3D simulator. That has not
started and is gated behind its own milestone, but it is the destination, and
it is why some of the current code looks the way it does: scalar angles and
scalar inertia are 2D specializations of quaternions and inertia tensors, not
an alternative to them. See
`docs/decisions/0009-three-dimensions-are-the-destination.md`.

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
13. `docs/decisions/` (ADRs 0001-0009)

## Build and test

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

All 13 test executables should pass, and
`out/build/windows-msvc-debug/Debug/malloy_nbody_terminal.exe` should run to
completion.

## What's next

M1-M18 are done. The classical mechanics track is complete, including the
rigid-body contact response that M11 promised and never numbered, which became
M14, the uniform gravity field for that domain, which became M15, and the
halfplane ground of M16 and the Coulomb friction of M17, and M18 added
electromagnetism beside it. Five
finished domains sit behind multi-domain dispatch. Work stays gated one milestone at a time; see
`docs/07_POST_M5_ROADMAP.md`.

The all-in-one goal is not a license to build ahead. Breadth is earned by
finishing domains, not by scaffolding for them
(`docs/09_MISTAKES_TO_AVOID.md`, #10).
