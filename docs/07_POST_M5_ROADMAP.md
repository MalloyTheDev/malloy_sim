# 07 — Post-M5 Roadmap

M1-M5 are complete, so this is now the live roadmap for what may come next.
Treat each entry as its own milestone: start it only when explicitly asked,
build it one milestone at a time, and do not expand scope immediately after M5
(`docs/09_MISTAKES_TO_AVOID.md`, #10).

Recommended order:

```text
M6  — stabilize N-body diagnostics
M7  — minimal scenario/config loading
M8  — simple 2D debug visualization
M9  — collision primitives
M10 — 2D rigid body basics
M11 — ballistics/projectiles
M12 — vehicle experiments
M13 — broader 2D sandbox integration
M14 — 3D math
M15 — 3D simulation experiments
M16 — quantum side project/module
```

Raylib is the likely first visualization choice later because it gets pixels on screen quickly without turning this into a graphics project.

Quantum should start later as a separate repo or isolated prototype. Only merge it as a module if it develops clean boundaries and does not pollute classical simulation types.
