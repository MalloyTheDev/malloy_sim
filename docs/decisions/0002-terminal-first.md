# ADR 0002 - Terminal First

## Status

Accepted.

## Decision

The first real simulation milestone is a terminal-based 2D N-body demo.

No graphics library is used until the milestone that introduces rendering.
That milestone had not arrived as of M8; see `docs/07_POST_M5_ROADMAP.md`.

## Rationale

Rendering introduces windowing, input, timing, camera, dependency, and platform questions. MalloySim is simulation-first.
