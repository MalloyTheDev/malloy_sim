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

## Scope creep watchlist

Do not add any of the following until a dedicated post-M5 milestone explicitly
calls for it:

Raylib, SDL, GLFW, SFML, rendering, input, camera, quantum, ECS, serialization,
config files, CLI parser, logging framework, plugins, scripting, editor, asset
manager, threading, package manager, or install/export rules.

Several items that were on this list have since had their milestone and are no
longer scope creep: collision geometry (M9), rigid bodies (M11), and the first
3D work, `Vec3` and 3D gravity (M19). The rest of 3D, meaning quaternions,
inertia tensors, 3D contact geometry and 3D versions of the other four domains,
is still gated and each piece needs its own milestone (`CLAUDE.md` rule 6).

The rule is the "until" clause, not the list. A list of names goes stale as
milestones land; the discipline of requiring a milestone first does not.
