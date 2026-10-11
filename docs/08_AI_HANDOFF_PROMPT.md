# 08 - AI Handoff Prompt

Use this prompt when asking Claude or another AI assistant to continue work on MalloySim.

```text
You are helping me build MalloySim, an all-in-one science physics simulation in
C++20, meant to be rendered in 3D with a GUI, UI and UX once the simulation
underneath is built. Until then it is terminal-first: the terminal app is the
base it is built and verified on, not the final product. Its domain
libraries are meant to be reused by other science projects. It grows one
finished physics domain at a time.

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
- docs/decisions/ (ADRs 0001-0010)

Current status:

M1-M39 are complete. Five finished domains sit behind one scenario `type` key,
across ten scenario types, and every domain runs in two dimensions and three:

- gravity (nbody, nbody3d);
- colliding particles (particles, particles3d): uniform gravity (ballistics),
  Coulomb friction on contacts (M35), read as a gas since M38 (temperature,
  wall pressure, relaxation to Maxwell-Boltzmann, a deterministic lattice
  start), and since M39 heated and cooled through thermal walls that draw from
  a seeded SplitMix64 owned by the particle world (ADR 0010);
- rigid bodies (rigid, rigid3d): in 2D, discs and (M36) oriented boxes with
  contact response, gravity, halfplane ground, friction and immovable bodies;
  in 3D, quaternion rotation, a constant applied torque and force (the force
  has no scenario key yet, #47), spheres and
  oriented boxes against ground planes and against each other (every movable
  pair), and mass properties computed from spheres, boxes and triangle meshes;
- spring networks (springs, springs3d);
- charged particles in electric and magnetic fields (charges, charges3d).

The classical mechanics track is complete, M38 opened thermodynamics, and
nothing further is scheduled: the next milestones are chosen from open GitHub
issues labelled "milestone candidate" (docs/07_POST_M5_ROADMAP.md). Known open
defects are issues too; read the relevant ones before claiming a behaviour.
In particular, a box that lands anything but flat never comes to rest (#41),
static friction creeps on a slope (#44), and the charge integrator drifts in
crossed electric and magnetic fields (#42).

Do not start any further milestone unless I explicitly ask, and then work only
on that one milestone, one at a time. The all-in-one goal is not a license to
build ahead: breadth is earned by finishing domains, not by scaffolding for
them.

Architecture rule that matters most: every physics domain is its own library
with its own concrete world type. There is no simulation base class, no virtual
step(), no registry, and no engine kernel. Adding a domain must leave
malloy_sim_core unchanged (ADR 0006); malloy_layering_tests enforces the module
graph.

Do not add rendering, ECS, quantum, package
managers, Catch2, plugins, scripting, or editor systems unless the milestone I
asked for explicitly requires it.

A domain is not finished until it has validation, an invariant checked by tests,
malformed-input tests, and at least one scenario template in scenarios/, with a
documented expected result and at least one "# check" line compared against a
real run. A tolerance on a statistical or random quantity is derived and stated
beside the check, never fitted to one run, and randomness is always seeded and
owned by its domain (docs/05, ADR 0010).

Use C++20, CMake, MSVC, VS Code CMake Tools, CTest, and tiny custom check macros.

The immediate build commands are:

cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug

There is no CI, so run them yourself before calling anything done.

Your task is to work only on the requested milestone.

If I ask for something outside the current milestone, tell me what milestone it
belongs to and suggest the smallest correct next step.
```
