# 08 - AI Handoff Prompt

Use this prompt when asking Claude or another AI assistant to continue work on MalloySim.

```text
You are helping me build MalloySim, an all-in-one science physics simulation in
C++20. It is terminal-first, and it grows one finished physics domain at a time.

Before doing anything, read:

- README.md
- CLAUDE.md
- CHANGELOG.md
- docs/00_START_HERE.md
- docs/01_V3_ARCHITECTURE_DECISION.md
- docs/02_MILESTONE_ROADMAP_M1_M5.md
- docs/03_MODULE_BOUNDARIES.md
- docs/04_NUMERIC_AND_PHYSICS_CONVENTIONS.md
- docs/05_TESTING_STRATEGY.md
- docs/06_CMAKE_AND_VSCODE_WORKFLOW.md
- docs/07_POST_M5_ROADMAP.md
- docs/09_MISTAKES_TO_AVOID.md
- docs/decisions/ (ADRs 0001-0007)

Current status:

M1-M12 are complete. Gravity, colliding particles (with ballistics as uniform
gravity in that domain) and 2D rigid bodies are the three finished domains, and
the scenario format dispatches between them with a type key. The active track is
classical mechanics depth (docs/07_POST_M5_ROADMAP.md); springs remain, and
rigid-body contact response is promised but was never given a milestone number.

Do not start any further milestone unless I explicitly ask, and then work only
on that one milestone, one at a time. The all-in-one goal is not a license to
build ahead: breadth is earned by finishing domains, not by scaffolding for
them.

Architecture rule that matters most: every physics domain is its own library
with its own concrete world type. There is no simulation base class, no virtual
step(), no registry, and no engine kernel. Adding a domain must leave
malloy_sim_core unchanged (ADR 0006).

Do not add rendering, ECS, collision, rigid bodies, 3D, quantum, package
managers, Catch2, plugins, scripting, or editor systems unless the milestone I
asked for explicitly requires it.

A domain is not finished until it has validation, an invariant checked by tests,
malformed-input tests, and at least one scenario template in scenarios/.

Use C++20, CMake, MSVC, VS Code CMake Tools, CTest, and tiny custom check macros.

The immediate build commands are:

cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug

Your task is to work only on the requested milestone.

If I ask for something outside the current milestone, tell me what milestone it
belongs to and suggest the smallest correct next step.
```
