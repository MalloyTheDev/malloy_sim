# ADR 0005 — Normalized Demo Units

## Status

Accepted.

## Decision

The first N-body demo uses normalized units:

```text
G = 1.0
sun_mass = 1.0
planet_mass = 0.000001
radius = 1.0
planet_velocity = 1.0
dt = 0.001
```

## Rationale

Normalized values make terminal output easy to read and make orbital behavior easy to eyeball.
