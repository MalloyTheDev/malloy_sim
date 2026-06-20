# 09 — Mistakes To Avoid

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

Do not add Raylib, SDL, GLFW, SFML, rendering, input, camera, collision, rigid bodies, 3D, quantum, ECS, serialization, config files, CLI parser, logging framework, plugins, scripting, editor, asset manager, threading, package manager, or install/export rules until a dedicated post-M5 milestone explicitly calls for them.
