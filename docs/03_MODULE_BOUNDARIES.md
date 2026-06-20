# 03 — Module Boundaries

> All four modules below are implemented (M2-M4). These boundaries are in force
> in the shipped code; keep them when extending the project.

## `malloy_math`

Responsible for `Real`, `Vec2`, scalar/vector helpers, approximate equality, and finite checks.

Not responsible for mass, force, velocity semantics, acceleration semantics, time, gravity, collision, or rendering.

## `malloy_time`

Responsible for fixed timestep representation, tick count, elapsed simulation time, and positive dt validation.

Not responsible for wall-clock time, frame pacing, sleeping, render interpolation, or physics formulas.

## `malloy_sim_core`

Responsible for tiny shared simulation vocabulary: `SimulationSettings`, `StepResult`, and `StepStatus`.

Not responsible for polymorphic simulation interfaces, generic world base classes, entity systems, schedulers, event buses, service locators, plugin systems, scene graphs, or engine kernels.

## `malloy_nbody`

Responsible for `Body2D`, `NBodySettings`, `NBodyWorld`, pairwise gravity, softening, semi-implicit Euler integration, deterministic update order, and N-body-specific validation.

Not responsible for printing, rendering, app loop policy, config files, input, or GUI.

## `malloy_nbody_terminal`

Responsible for hardcoded demo setup, calling library APIs, fixed number of steps, formatted terminal output, and returning nonzero on validation failure.

Not responsible for reusable physics logic.
