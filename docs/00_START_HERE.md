# 00 - Start Here

**MalloySim C++** is an all-in-one science physics simulation in C++20, built
simulation-first and terminal-first. It grows one finished physics domain at a
time, each shipping tested scenario templates, with no engine kernel and no
simulation base class (`docs/decisions/0006-multi-domain-dispatch.md`).

M1-M27 are complete: the math, time, sim-core, N-body, scenario, ASCII, collide,
particles, rigid, springs and charges libraries are built, and the terminal app
runs all seven scenario types. Gravity, colliding particles, rigid bodies,
spring networks and charged particles are the finished domains, and two of them
now run in three dimensions as well as two. See `README.md` and `CHANGELOG.md`.

## Where this is going

MalloySim is intended to become a 3D simulator, and it is under way. M19
shipped `Vec3` and a 3D gravity world; M20 shipped `Quat` and torque-free
rigid-body rotation, the first piece of 3D with no 2D form at all (a body with
three distinct principal moments tumbles under no torque, and in two dimensions
it cannot); M21 applied a constant torque to that rotation, so a gyroscope
precesses instead of toppling; M22 gave a body a collision radius and bounced
it off a ground plane under gravity, the project's first 3D contact; M23 added
friction, so a sliding sphere spins up and rolls; M24 made two spheres collide
and exchange momentum; M25 computes a body's mass and inertia from its geometry;
M26 does the same for boxes and assembles compound bodies from primitives; and
M27 does it for an arbitrary closed triangle mesh. Everything else is still 2D,
which is most of the project, and each remaining
piece is its own milestone.

That destination is why some of the current code looks the way it does: the
scalar angle and scalar inertia in `RigidBody2D` are 2D specializations of the
quaternion and the principal moments in `RigidBody3D` rather than an
alternative to them. See
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

All 14 test executables should pass, and
`out/build/windows-msvc-debug/Debug/malloy_nbody_terminal.exe` should run to
completion.

## What's next

M1-M27 are done. The classical mechanics track is complete, including the
rigid-body contact response that M11 promised and never numbered, which became
M14, the uniform gravity field for that domain, which became M15, and the
halfplane ground of M16 and the Coulomb friction of M17, and M18 added
electromagnetism beside it. M19 took the first step off the plane, with Vec3
and a 3D gravity world, M20 added quaternions and torque-free 3D rotation, M21
put a constant torque on it, M22 bounced a sphere off a ground plane, M23
added friction, M24 collided two spheres, M25 computes mass and inertia from
geometry, M26 adds box mass properties and `combine` for compound bodies, and
M27 computes them for an arbitrary triangle mesh.
Five finished domains sit behind multi-domain dispatch, across seven
scenario types.
Work stays gated one milestone at a time; see
`docs/07_POST_M5_ROADMAP.md`.

The all-in-one goal is not a license to build ahead. Breadth is earned by
finishing domains, not by scaffolding for them
(`docs/09_MISTAKES_TO_AVOID.md`, #10).
