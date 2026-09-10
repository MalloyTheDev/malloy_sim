# ADR 0009 - Three Dimensions Are The Destination

## Status

Accepted, recorded at M16. Nothing 3D has been started, and none is permitted
before its own milestone (`CLAUDE.md` rule 6).

## Context

Every mention of 3D across the documentation was a prohibition. `CLAUDE.md`
listed "a 3D project yet" in a block of things the project is not, immediately
below "an editor project" and "a plugin framework". `README.md` listed a bare
"3D" in a section headed "Out of scope (gated)", between "vehicles" and
"quantum" and above "ECS", "plugins", "editor" and "scripting". The only place
that framed it as future work was one line in Track 4 of the roadmap.

That put two very different things in one list. Plugins, an editor and a
scripting layer are rejected because adopting them would change what the
project is. 3D is not rejected at all: it is where the project is going. A
reader had no way to tell those apart, and the natural conclusion from the
README was that MalloySim will always be a 2D simulator.

The distinction matters for more than tone. Decisions have already been made
whose justification is that a 3D form exists and is being deliberately
postponed. ADR 0007 chose a scalar angle and a scalar inertia for
`RigidBody2D`, and it did so on the grounds that quaternions and inertia
tensors are the same quantities in three dimensions rather than a different
model. Someone who believes 3D is out of scope reads that as an arbitrary
simplification instead of a specialization, and may then "simplify" further in
ways that do not generalize.

## Decision

**MalloySim is intended to become a 3D simulator.** Three dimensions are the
destination, not a door left open.

This changes no code and relaxes no gate. 3D arrives as its own milestone,
after the 2D mechanics it rests on are finished, and nothing 3D is added
speculatively before then. What it changes is how the present code should be
read, and what a future 2D decision has to be checked against.

## What carries over, and what does not

Recorded so that a 2D decision can be checked against it rather than guessed at.

**Carries over essentially unchanged.** The fixed timestep and `malloy_time`.
The `sim_core` vocabulary: `SimulationSettings`, `StepStatus`, `StepResult`,
and the rule that a world validates and returns status rather than throwing.
Semi-implicit Euler and the velocities-before-positions ordering. The
multi-domain dispatch of ADR 0006, since a `type` key selects a concrete world
whatever its dimension. The determinism policy in `docs/04`. The testing
approach, including mutation testing and asymmetric configurations.

**Specializes, meaning the 2D form is the 3D one with a dimension removed.**
`Vec2` becomes a three-component vector. The scalar angle becomes a quaternion.
The scalar inertia becomes a 3x3 tensor. Scalar angular velocity and scalar
torque become vectors. The scalar `cross(a, b)` used in `malloy_rigid` becomes
the vector cross product. The impulse formulas keep their structure: the
effective mass stays `1/mA + 1/mB` plus rotational terms, with
`(r x n)^2 / I` becoming `n . (I^-1 (r x n)) x r`. `shift_inertia`, the
parallel-axis step, becomes the tensor form of the same theorem.

**Does not carry over.** `malloy_ascii` projects points that are already in the
view plane; drawing a 3D scene is a different problem, not a wider one. The
collision primitives are 2D area shapes: `Circle` becomes a sphere and the
`Halfplane` of M16 generalizes cleanly, but `second_moment_of_area` becomes a
volume integral and the area properties are replaced rather than extended. The
scenario format's fixed field lists all grow.

**Is genuinely new, and is the reason 3D is a milestone rather than a
widening.** In 2D, angular velocity points along a fixed axis and inertia is a
scalar, so a body's inertia never changes in world space and there is no
gyroscopic term. In 3D the inertia tensor rotates with the body, and Euler's
equations carry `omega x (I omega)`. A spinning body precesses and can tumble
about its intermediate axis with no torque applied at all. That is physics the
2D code does not contain in any form, and it cannot be reached by replacing
types.

## Deliberately not decided here

Left to the milestone, so this ADR records a direction and not a design:

- whether 3D types are separate (`Vec3`, `RigidBody3D`) or the existing ones
  become templates on dimension. `docs/02` banned `Vec3` and templates for
  M1-M5 and that ban has served; it is not evidence either way now.
- whether 3D domains are new concrete worlds beside the 2D ones, dispatched by
  the same `type` key, or replacements for them. ADR 0006 constrains only that
  neither may be reached through a shared base class.
- whether 2D remains supported once 3D exists. It probably should, since the 2D
  templates are a shipped deliverable, but that is a decision with a cost and
  it has not been made.

## Consequences

The documentation now separates what is planned and gated from what is
rejected, in `README.md` and `CLAUDE.md`, and 3D appears in the first as the
project's stated destination.

A 2D decision can be checked against the lists above. If a proposed
simplification has no 3D form, that is a reason to reconsider it, not a reason
to defer thinking about it. If it specializes cleanly, it can be made without
worrying about the eventual transition.

Nothing here licenses building ahead. Rules 6 and 7 stand, the milestone
sequence is unchanged, and a scalar angle remains the correct representation
for a 2D rigid body today.
