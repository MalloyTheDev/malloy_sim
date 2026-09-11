# ADR 0009 - Three Dimensions Are The Destination

## Status

Accepted, recorded at M16, amended at M19, M20, M22, M23, M24 and M25.

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

## Amendment at M22: the collision prediction, confirmed

The section above put contacts in the SPECIALIZES tier, not the genuinely-new
one: "`Circle` becomes a sphere and the `Halfplane` of M16 generalizes cleanly,
but `second_moment_of_area` becomes a volume integral and the area properties
are replaced rather than extended." M22 built the first 3D contact and every
clause held.

`collide::Sphere` and `collide::Plane3` are direct specializations of `Circle`
and `Halfplane`, in new headers inside the same library; the 2D primitives are
untouched. The contact query and its `Contact3` carry the same normal / never-
negative-depth / midpoint convention with no conceptual change, and the plane,
like the halfplane, has no degenerate case. No area or volume property was
needed at all: a centred sphere's normal impulse never reads the inertia, so
the milestone did not have to build the volume integral the prediction warned
about.

The one thing worth recording that the ADR did not foresee: a centred sphere's
normal contact imparts no spin, because the contact point is on the line
through the centre of mass, so `r x n = 0`. That is why M22 is the 3D echo of
M10 (translational contacts) rather than M14 (rotational contacts), and why the
rotational effective-mass term this ADR sanctioned (`n . (I^-1 (r x n)) x r`,
line 66) is not yet exercised: it is identically zero here and first bites with
friction, whose tangential impulse has `r x t != 0`.

## Amendment at M23: the impulse formula, put to work

The SPECIALIZES section predicted that the impulse formulas keep their
structure, with the rotational effective-mass term `(r x n)^2 / I` becoming
`n . (I^-1 (r x n)) x r` (line 66). M22 built the normal impulse but could not
exercise that term: a centred sphere's normal contact has `r x n = 0`, so the
term is identically zero there. M23's friction is the first contact with a
nonzero arm, `r x t`, and it uses exactly that formula, in exactly that shape.

Building it confirmed the prediction and added one detail the ADR did not
spell out: the `I^-1` in that formula is the WORLD-frame inverse inertia,
`R I^-1_body R^T`, computed by the same body-frame bridge Euler's equations use
for the torque (M21). For a sphere the inertia is isotropic and the rotation is
a no-op, but the formula is written in the general form and a test with a
non-isotropic body pins it, so the machinery is correct before the general
inertia tensor it will eventually pair with exists.

The physics it reaches is the 3D echo of M17's rolling: a sliding sphere rolls
without slipping at `5/7` of its sliding speed, independent of the coefficient
and of gravity, the same shape as the 2D `2/3` disc ratio. So 3D contact
response has now specialized from the 2D code exactly as this ADR said it would,
across both the normal impulse (M22) and the tangential one (M23).

## Amendment at M24: one impulse core for both contact types

M22 and M23 resolved a sphere against an immovable plane. M24 added a sphere
against a second movable sphere, and in doing so unified the two: the contact
response is now one two-body core, with the plane expressed as a participant of
zero inverse mass. This is the 3D echo of the 2D `resolve_ground` stand-in body
(ADR 0008), and it is the answer to whether the impulse formula would be written
twice. It is not. Refactoring the sphere-plane path onto the shared core left
every M22/M23 result bit for bit unchanged.

The new physics a second movable body brings is conservation: an immovable plane
is a momentum sink, so a sphere-plane contact says nothing about momentum, while
two real bodies exchange an equal and opposite impulse and so conserve total
linear momentum exactly. That is the sphere-sphere headline, and it is the 3D
form of what M10's colliding particles first showed in the plane.

Still deferred, and still each its own milestone: a translational force on a 3D
body, a general inertia tensor, and 3D versions of the particle, spring and
charge domains. Oriented boxes and the SAT question remain untouched in both
dimensions.

## Amendment at M25: the inertia tensor, as a construction step

The SPECIALIZES section said the scalar inertia "becomes a 3x3 tensor" and that
`shift_inertia` "becomes the tensor form of the same theorem." M25 built both:
`mass_properties_3d` composes solid spheres, shifts each by the 3D parallel-axis
theorem m (|d|^2 I - d d^T), and forms the compound body's inertia tensor, with
`math::Mat3` and a symmetric eigensolver new in `malloy_math` to hold and
diagonalize it.

The refinement the M20 amendment already flagged holds all the way through: the
tensor is where inertia is COMPUTED, not where it is STORED. Diagonalizing it
gives principal moments and an orientation, which is what `RigidBody3D` carries,
so the body's dynamics never touch a 3x3. A symmetric tensor always
diagonalizes, so this is not a shortcut; it is the reason the principal-moments
choice was correct from M20 on.

What this unlocks is the vision the project is heading toward: a body's physical
properties FOLLOWING from its geometry rather than being typed in. M25 does it
for compound spheres; other shapes (boxes, meshes, via volume integrals) are
later milestones, each its own.

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
