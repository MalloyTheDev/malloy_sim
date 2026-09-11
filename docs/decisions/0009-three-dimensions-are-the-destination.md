# ADR 0009 - Three Dimensions Are The Destination

## Status

Accepted, recorded at M16, amended at M19 and M20.

As of M20 the first two pieces of 3D have shipped: `math::Vec3` and
`nbody::NBody3DWorld`, then `math::Quat` and `rigid::Rigid3DWorld`. Everything
else is still 2D, and each remaining piece needs its own milestone
(`CLAUDE.md` rule 6). The amendments below record how M19 answered the three
questions this ADR deliberately left open, and what M20 found when it built the
part this ADR called genuinely new.

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

## Amendment at M19: the three open questions, answered

M19 shipped `math::Vec3` and `nbody::NBody3DWorld`, and building them settled
all three questions below. Each is answered by what the work actually required
rather than by choosing in the abstract, which is why they were left open.

**Separate types, not templates.** `Vec3` is its own concrete type beside
`Vec2`. The decisive reason is that the two are not one algebra with a
different component count: `cross` returns a SCALAR in two dimensions and a
VECTOR in three, so a template would need specializing for the one operation
that matters most to rigid-body dynamics, and `perp` and rotation by a scalar
angle have no 3D form at all. A template would have unified the component-wise
arithmetic, which is the easy half, and then needed specializing for every part
that is interesting. What is shared is spelled identically on purpose, so code
that reads one reads the other.

**3D worlds sit inside their domain's library, not beside it.** `Body3D` and
`NBody3DWorld` live in `malloy_nbody`. Gravity is one domain, and two and three
dimensions are the same physics with a different component count, not two
domains that share a name. Splitting them would put the same equations in two
libraries and require keeping them in step by hand. The unit ADR 0006 makes a
library is the domain; the dimension is not one. The dispatch switch grew by
one branch, `type nbody3d`, exactly as a new domain would have cost it.

**2D stays supported.** It follows from the answer above rather than being a
separate decision: the 2D types are untouched, both worlds are selected by the
same `type` key, and the twelve shipped templates still run. The strongest
evidence the 3D code is right is a test that embeds a 2D configuration at
z = 0 and runs both worlds: their positions agree bit for bit, and the third
component stays exactly zero.

What M19 did NOT do, deliberately: quaternions, inertia tensors, 3D contact
geometry, and 3D versions of the other four domains. N-body is the only domain
with no contacts and no orientation, so it needed `Vec3` and nothing else,
which let the rest stay deferred rather than be built speculatively (rule 11).

## Amendment at M20: the genuinely new part, built

The section above singled out one thing as the reason 3D is a milestone rather
than a widening: that in 3D a body precesses and can tumble about its
intermediate axis with no torque applied, and that this is physics the 2D code
does not contain in any form. M20 built it, and that prediction was right in
every respect. `Rigid3DWorld` needs no contacts, no collision geometry, no
forces and no solver, so it reached the new physics without building anything
ahead of its milestone, the same way M19 reached `Vec3` through N-body.

**One prediction above needs correcting.** This ADR said the scalar inertia
"becomes a 3x3 tensor". `RigidBody3D` stores three PRINCIPAL moments instead,
and that is a decision rather than an omission. A symmetric tensor is always
diagonalizable, so every rigid body has a frame in which the inertia is
diagonal; storing the diagonal is a choice of axes, not a restriction. A
general tensor becomes necessary only once bodies are built from composed
shapes and the parallel-axis step moves inertia off the principal axes, which
arrives with 3D mass properties and 3D contacts. Written now it would be three
extra zeros and a speculative abstraction (rule 11). The tensor form of
`shift_inertia` is deferred with it.

Two predictions held exactly as written. The scalar angle became a quaternion,
with the 2D form untouched and still shipping. And "carries over" was accurate
about the integrator: `Rigid3DWorld` is semi-implicit Euler with
velocities-before-positions, validates and returns status rather than throwing,
and is selected by the same `type` key through the same dispatch switch, which
grew by one branch.

**What this did NOT settle.** Whether `malloy_rigid` eventually holds a shared
integrator once a third domain needs the same translational path, which ADR
0008 rule 17 asks to be reassessed then. `Rigid3DWorld` duplicates the
`position += velocity * dt` line and nothing more, which is not yet evidence.

## Deliberately not decided here (resolved above at M19)

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
