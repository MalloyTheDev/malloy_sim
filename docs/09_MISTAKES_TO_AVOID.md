# 09 - Mistakes To Avoid

1. Turning `malloy_sim_core` into an engine kernel.
2. Adding rendering before terminal N-body works.
3. Adding ECS before `std::vector<Body2D>` actually fails.
4. Using real astronomical SI values in the first demo.
5. Using `float` for positions/velocities.
6. Enabling fast-math for physics.
7. Letting `main.cpp` contain physics formulas.
8. Creating empty folders and fake architecture.
9. Making CMake clever too early.
10. Expanding scope immediately after M5.
11. Writing a symmetric test, which a wrong implementation can pass
    (`docs/05`, "Prefer asymmetric configurations").
12. Fitting a statistical tolerance to the run it checks, or reseeding until a
    check passes (`docs/05`, ADR 0010).
13. Reaching for `<random>`, `rand`, a time-based seed or a global generator
    (ADR 0010).
14. Trusting a documented figure that no `# check` line compares against a run
    (`CLAUDE.md`, Template library).

## Scope creep watchlist

Do not add any of the following until a dedicated post-M5 milestone explicitly
calls for it:

Raylib, SDL, GLFW, SFML, rendering, input, camera, quantum, ECS, serialization,
CLI parser, logging framework, plugins, scripting, editor, asset manager,
threading, package manager, or install/export rules.

Several items that were on this list have since had their milestone and are no
longer scope creep: config files (the scenario format, M7), collision geometry
(M9), rigid bodies (M11), and 3D, which began with `Vec3` and 3D gravity (M19)
and has since given every domain a 3D form, quaternions, the inertia tensor as a
construction step, and 3D contact geometry (M20 to M34). What remains of 3D,
rendering, a per-point contact solver and immovable 3D bodies among it, is still
gated, and each piece needs its own milestone (`CLAUDE.md` rule 6,
`docs/07_POST_M5_ROADMAP.md`).

The rule is the "until" clause, not the list. A list of names goes stale as
milestones land; the discipline of requiring a milestone first does not.
