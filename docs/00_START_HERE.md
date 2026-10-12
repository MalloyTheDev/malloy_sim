# 00 - Start Here

**MalloySim C++** is an all-in-one science physics simulation in C++20, meant
to be seen: its destination is a simulator rendered in 3D, with a GUI, UI and
UX. It is built simulation-first: the terminal app is the base it is built
and verified on until the rendering milestone, not the final product, and its
domain libraries are meant to be the foundation of other science projects. It grows one finished physics domain at a time, each shipping tested scenario templates, with no engine kernel and no
simulation base class (`docs/decisions/0006-multi-domain-dispatch.md`).

M1-M40 are complete: the math, time, sim-core, N-body, scenario, ASCII, collide,
particles, rigid, springs and charges libraries are built, and the terminal app
runs all ten scenario types. Gravity, colliding particles, rigid bodies,
spring networks and charged particles are the finished domains, and all five of
them run in three dimensions as well as two. Thirty-three templates in `scenarios/`
cover all ten types. See `README.md` and `CHANGELOG.md`.

## Where this is going

MalloySim is intended to become a 3D simulator, and it is well under way.
M19 shipped `Vec3` and a 3D gravity world; M20 `Quat` and torque-free
rigid-body rotation, the first piece of 3D with no 2D form at all (a body with
three distinct principal moments tumbles under no torque); M21 a constant
torque, so a gyroscope precesses; M22 to M24 a sphere bouncing on a ground
plane, rolling with friction, and colliding with another sphere; M25 to M27 a
body's mass and inertia computed from its geometry, for spheres, boxes and any
closed triangle mesh; M28 a constant applied force; M29, M31 and M32 the charge,
spring and particle domains in three dimensions; and M30, M33 and M34 the
oriented box, resting on a plane and colliding with boxes and spheres. Every
domain now runs in space as well as the plane, and any two 3D primitives
collide. M36 and M37 brought the box back to the plane: 2D box against box, and
a 2D box on the ground.

M40 gave both rigid worlds a per-point iterative contact solver, so a box that
lands tilted comes to rest and short stacks stand (ADR 0011). What remains of
the 3D arc is listed in `docs/07_POST_M5_ROADMAP.md`: immovable 3D bodies and
rendering, each its own milestone.

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
13. `docs/decisions/` (ADRs 0001-0011)

`CHANGELOG.md` is the record of what each milestone actually shipped, and the
open GitHub issues are the record of known defects and milestone candidates.

## Build and test

```powershell
cmake --preset windows-msvc-debug
cmake --build --preset windows-msvc-debug
ctest --preset windows-msvc-debug
```

All 15 test executables should pass, and
`out/build/windows-msvc-debug/Debug/malloy_nbody_terminal.exe` should run to
completion. CTest runs 30 tests in all: the executables plus sixteen checks
that run the app itself. CI runs the Debug and Release builds with every test,
and a short fuzz run, on every push and pull request to `main`
(`.github/workflows/ci.yml`).

To run one template:

```powershell
out/build/windows-msvc-debug/Debug/malloy_nbody_terminal.exe scenarios/heat_conduction.scn
```

## What's next

M1-M40 are done. The classical mechanics track is complete, including the
rigid-body contact response that M11 promised and never numbered, which became
M14, and every domain has its 3D form. M38 then opened thermodynamics, reading a
box of colliding particles as a gas: its temperature, the pressure on its walls,
and its relaxation to Maxwell-Boltzmann. M39 gave its walls a temperature, so
heat flows through the gas, using the project's first random numbers: seeded,
fully specified and owned by the particle domain (ADR 0010). M40 returned to
the rigid domain with a per-point iterative contact solver (ADR 0011). Five
finished domains sit behind multi-domain dispatch, across ten scenario types.

Nothing further is scheduled. The candidates for the next milestone are open
issues labelled `milestone candidate`, listed in
`docs/07_POST_M5_ROADMAP.md`: 3D rendering and then a GUI, UI and UX (the
destination), packaging the libraries for other projects, immovable 3D bodies, box against disc in 2D, a piston, accommodation
coefficients, a temperature profile, a particle broadphase, gravity for
springs, and a shared `# check` evaluator. Work stays gated one milestone at a
time, started only when asked.

The all-in-one goal is not a license to build ahead. Breadth is earned by
finishing domains, not by scaffolding for them
(`docs/09_MISTAKES_TO_AVOID.md`, #10).
