# ADR 0012 - The 3D Viewer Is A Separate Front End Over Shared, Pure Layers

## Status

Accepted, recorded at M41. Supersedes, in part, ADR 0002: terminal-first stays
as the build order, and rendering now exists.

## Context

MalloySim is meant to be seen. Its destination is a simulator rendered in 3D,
with a GUI, UI and UX on top (`CLAUDE.md`, "Where it is going"). Terminal-first
was the build order, not the goal, and by M40 the physics was ready to show:

- every domain has a 3D form;
- every 3D pair collides;
- boxes come to rest and short stacks stand.

Rendering brings what ADR 0002 deferred: windowing, input, timing, a camera, a
third-party dependency and platform questions. The risks are specific:

- the graphics code reaches the physics;
- a second front end copies the first one's setup and drifts;
- a dependency's warnings or flags break the hardened build;
- a window that nothing in CI can check.

## Decision

**Three layers, and only the outermost one knows about graphics:**

- **The run layer** (`include/malloy/scenario/run.hpp`, in `malloy_scenario`).
  It holds everything a front end needs that is not presentation:
  - building each scenario type's world from a parsed `Scenario`
    (`make_<type>_world`);
  - reading back what a front end reports (a plain readout struct per type);
  - the report line the terminal prints.

  Each front end keeps its own `switch` over the scenario type. There is no
  simulation base class, no type-erased world and no `std::variant` of worlds
  (rule 12, ADR 0006). The terminal app now calls this layer, and its output
  stayed byte-identical.
- **The render model** (`malloy_view3d`), a pure static library and the 3D
  sibling of `malloy_ascii`.
  - It turns each 3D world's state into final world-space primitives and
    triangle meshes, fits a camera to them, and paces fixed steps against a
    wall clock.
  - Every number a viewer could get wrong is computed and tested here: box
    corners, winding, sphere meshes, plane patches, the camera fit, and frame
    pacing.
  - It never links or includes a graphics library.
- **The viewer** (`malloy_viewer`, `apps/viewer`), a thin executable that only
  draws what the render model describes. It is the only target that includes
  `raylib.h`.

**The graphics library is raylib 6.0.**

- **Source:** `https://github.com/raysan5/raylib/archive/refs/tags/6.0.tar.gz`,
  fetched by CMake `FetchContent` and pinned by its SHA-256,
  `2b3ee1e2120c7a0796b33062c7e9a694dd8a8caa56a96319ac8c8ecf54a90d0b`.
  - That value was checked by hand against the independently published Conan
    and vcpkg values.
  - The archive's embedded commit is the tag's,
    `dbc56a87da87d973a9c5baa4e7438a9d20121d28`.
  - A changed archive fails the configure step.
- **Isolation:** it is declared `SYSTEM` and `EXCLUDE_FROM_ALL`, so its headers
  and warnings never reach the project's `/W4 /WX`.
- **Build settings:** it is built static, with audio off and with Control Flow
  Guard added.
- **The option:** `MALLOY_BUILD_VIEWER` is OFF in `CMakeLists.txt`, so a bare
  configure, the fuzzers and future library consumers (#82) never fetch
  anything. It is ON in the presets and CI, which build the viewer.

**The layering test enforces the boundary:**

- no file in `include/` or `src/` includes a graphics header;
- a library links only `malloy::` targets;
- no domain links or includes the run layer or the render model;
- every library sits on a layer and links only lower ones.

**The viewer never changes the physics.**

- **One thread.** It steps the world only through `step()`, on one thread,
  between frames, at the scenario's fixed `dt`. So "reads state between steps"
  holds by construction.
- **Bounded runs.** It stops at the scenario's `steps` and holds the last
  frame. A reset rebuilds the world from the parsed `Scenario`, never from
  `bodies()`, because warm-start impulses and random streams are world state
  (ADRs 0010, 0011).
- **Same caps.** An untrusted file is bounded exactly as in the terminal app,
  by the same parser caps.
- **Scope.** It draws the five 3D scenario types; a 2D type is refused with a
  pointer to the terminal app.

**The window is checked without a window where possible:**

- a headless mode runs the same loop with a synthetic clock and prints the
  terminal's report line. CTest requires that line to equal the terminal app's,
  and to be identical at 20, 60 and 144 frames a second;
- a pixel self-test runs under Xvfb and Mesa on the Ubuntu CI job;
- a written manual check-list covers what only a person at a display can judge.

## Rejected alternatives

- **sokol, SDL3, GLFW plus a GL loader.**
  - sokol has no release tags, so it can only be pinned by commit.
  - Beyond immediate mode, sokol and SDL3 need shaders and an offline shader
    toolchain.
  - SDL3 is large.
  - GLFW plus a loader means writing shaders, buffers, matrices and a font by
    hand.
  - raylib's 3D view needs a few hundred lines on top of a pure model.
- **bgfx, Magnum, VTK, Polyscope.** Heavy, or they own the main loop.
- **raylib's prebuilt binaries.** Release only, no Control Flow Guard, and built
  for AVX2.
- **Vendoring the source or a package manager.** Vendoring adds about 15 MB of
  third-party source and large diffs on every upgrade. A package manager breaks
  rule 8.
- **A `std::variant` of worlds, or a world interface.** The variant duplicates
  the type key and invites engine vocabulary (ADR 0006, rule 12). Drawing and
  formatting need per-type code anyway.
- **A separate simulation thread.** Real time is easy on one thread: the
  heaviest template's physics costs 0.4 per cent of a 60 Hz frame. A second
  thread would turn "reads between steps" into a synchronization argument.
- **Drawing the 2D types in M41.** That is a follow-up on the same drawers.
- **Looping a run past `steps`.** It needs the runtime work budget (#85).

## Consequences

- **The physics is now visible.** The terminal app remains the development base
  and the reference the viewer is checked against.
- **A configure that enables the viewer downloads 52.5 MB.** That covers the
  presets and CI. It needs network access, or a verified local copy through
  `FETCHCONTENT_SOURCE_DIR_RAYLIB`.
- **Coverage differs by runner.** CI builds the viewer on Windows, where the
  hosted runner has no OpenGL 3.3, so only the headless checks run there. It
  draws real frames only on Ubuntu, under software rendering.
- **The run layer is the shared contract for later front ends.** The GUI (#81)
  builds on it and on the render model.
