# ADR 0008 - SpringWorld Is A Local Composition Boundary

## Status

Accepted, implemented in M13.

## Context

M13 introduces springs, which are the first force model with a topology: a
spring connects two specific bodies, and a body may be an endpoint of several.
That makes it the first many-to-one force pipeline in the project.

Springs could have gone into `ParticleWorld`, which already owns bodies,
collision, walls, and uniform gravity. That was rejected. Gravity is
environmental configuration, one setting that applies to everything. A spring
network is interaction topology, which is new state. Putting it into
`ParticleWorld` would make the particle domain own a specific force model, and
the same argument would then admit cloth, rods, distance constraints, cables and
breakable joints, until `ParticleWorld` became the general physics engine rather
than one domain.

## Decision

`malloy_springs` owns `Spring`, `SpringNetwork`, `SpringBody2D`, and
`accumulate_spring_forces`, which is a pure kernel: it reads body state and
writes into a force buffer, and it neither integrates nor mutates any body.

Beside it, `SpringWorld` composes that kernel with the minimal translational
integration needed to make the domain runnable and testable. Its step contract
is fixed and documented:

```text
1. clear the per-body force accumulator
2. evaluate springs in stable network order
3. accumulate equal and opposite endpoint forces
4. integrate bodies in stable index order
```

`ParticleWorld` is not modified. No generic force-provider API is introduced.

## SpringWorld is not the engine-wide force architecture

`SpringWorld` is the M13 domain-level composition boundary. It exists to make
spring dynamics independently runnable and testable without modifying existing
world APIs. Generalized external-force composition across domains is
intentionally deferred until additional force-producing systems provide
sufficient requirements for a stable shared abstraction.

This paragraph is here so that a later milestone does not conclude that
`SpringWorld` already solves generic force composition. It does not. It is
deliberately local.

The eventual architecture may well converge on a shared pipeline:

```text
gravity  ---+
springs  ---+
drag     ---+--> force buffer --> body integration
buoyancy ---+
thrusters --+
contacts ---+
```

but one force producer is not enough evidence to design that interface
correctly. Designing it now would mean choosing between
`step(dt, forces)`, `apply_forces(forces)`, `accumulator()`,
`ExternalForceProvider`, `ForceGenerator` and similar on the strength of a
single example. Two or three genuinely different producers, such as springs plus
drag plus buoyancy, give empirical grounds for the choice.

## Duplicated integration is intentional, with a trigger

`SpringWorld` repeats roughly five transparent lines that `ParticleWorld` also
has:

```text
acceleration = force / mass
velocity += acceleration * dt
position += velocity * dt
```

That duplication is accepted for M13. Five obvious lines carry less
architectural risk than a premature abstraction.

The refactoring trigger is explicit: **if a third domain independently needs the
same translational integration path, reassess extracting a shared translational
integrator or a generic force-composition mechanism.** Two copies can be
coincidence. Three copies are evidence of an abstraction.

## Deferred

Rigid-body spring attachment points and the torque they generate stay deferred.
A rigid-body spring needs local anchors transformed into world space, and each
endpoint then generates a moment arm and a torque, which would turn M13 into
force accumulation plus torque accumulation plus attachment transforms plus
angular integration plus spring topology plus damping. That is too many concepts
to validate in one milestone. M13 connects body centres only.
