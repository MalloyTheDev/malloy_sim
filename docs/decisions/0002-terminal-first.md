# ADR 0002 - Terminal First

## Status

Accepted.

## Decision

The first real simulation milestone is a terminal-based 2D N-body demo.

No graphics library is used until the milestone that introduces rendering.
That milestone had not arrived as of M8; see `docs/07_POST_M5_ROADMAP.md`.

Note, 2026-10-11: still true at M39. No graphics library is used; the 3D
domains are drawn in the terminal as an orthographic flattening onto the xy
plane of the ASCII view, and rendering remains its own milestone.

## Rationale

Rendering introduces windowing, input, timing, camera, dependency, and platform questions. MalloySim is simulation-first.
