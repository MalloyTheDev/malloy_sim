# 03 - Module Boundaries

> All modules below are implemented (M2-M8). These boundaries are in force
> in the shipped code; keep them when extending the project.

> New physics domains follow the same shape: one library, one concrete world
> type, one settings type, one test executable, depending only on
> `malloy_math` and `malloy_sim_core`. No domain depends on another and none
> derives from a shared base (`docs/decisions/0006-multi-domain-dispatch.md`).

## `malloy_math`

Responsible for `Real`, `Vec2`, scalar/vector helpers, approximate equality, and finite checks.

Not responsible for mass, force, velocity semantics, acceleration semantics, time, gravity, collision, or rendering.

## `malloy_time`

Responsible for fixed timestep representation, tick count, elapsed simulation time, and positive dt validation.

Currently exercised only by its own test executable. No library or app links it: `NBodyWorld` keeps its own tick count and `SimulationSettings` holds its own `dt`.

Not responsible for wall-clock time, frame pacing, sleeping, render interpolation, or physics formulas.

## `malloy_sim_core`

Responsible for tiny shared simulation vocabulary: `SimulationSettings`, `StepResult`, and `StepStatus`.

Not responsible for polymorphic simulation interfaces, generic world base classes, entity systems, schedulers, event buses, service locators, plugin systems, scene graphs, or engine kernels.

## `malloy_nbody`

Responsible for `Body2D`, `NBodySettings`, `NBodyWorld`, pairwise gravity, softening, semi-implicit Euler integration, deterministic update order, and N-body-specific validation.

Not responsible for printing, rendering, app loop policy, config files, input, or GUI.

## `malloy_scenario`

Responsible for parsing a scenario/config text file into `SimulationSettings`, `NBodySettings`, bodies, and run length, and for reporting syntax errors with a line number.

Not responsible for semantic validation (that stays in `NBodyWorld`), physics, printing, rendering, or CLI argument handling.

## `malloy_ascii`

Responsible for turning 2D points into a character grid: viewport fitting (`fit_viewport`) and framed grid rendering (`render`).

Not responsible for physics, body types, simulation state, terminal control sequences, color, animation, or any graphics API.

## `malloy_nbody_terminal`

Responsible for hardcoded demo setup, reading a single positional scenario path from argv and choosing between a file and the built-in scenarios, calling library APIs, fixed number of steps, formatted terminal output, view framing policy (which viewport each frame shows), and returning nonzero on validation or step failure.

Not responsible for reusable physics logic.
