# 10 - Changelog Template

This is the template/reference. The live changelog is `CHANGELOG.md` at the
repo root; keep it updated with this format.

Two kinds of entry exist. Work between milestones (fixes, hardening, review
passes) collects under `[Unreleased]`, citing the issue it closes. A milestone
gets its own entry, with its title in parentheses after the date. Since M28 a
milestone entry has followed the shape below (M20 to M27 close to it): what was
added, a section in plain prose on why the milestone is shaped the way it is
(its heading names the point, for example "The last movable pair in 3D"), the
invariants its tests pin, the mutations applied and what caught them, and,
since M38, the documents it changed. A milestone may add a section of its own
where it has more to say (M39's `Calibration`).

```markdown
# Changelog

## [Unreleased]

### Security

### Added

### Fixed

- Issue #N: what was wrong, how it showed, and what now prevents it.

## [Mn] - YYYY-MM-DD  (short title)

### Added

- New types, settings, scenario keys, check quantities, templates (and the new
  template count).

### <A heading that states the milestone's point>

Why it is built this way, what it deliberately does not do, and which rule or
ADR that follows from.

### Invariants

What the tests pin, in closed form where possible.

### Mutation testing

How many mutations were applied, which were caught, and for each survivor
whether it is a true equivalent or a gap (and how the gap was closed).

### Docs

Which documents, ADRs and headers changed.

## [M1] - YYYY-MM-DD

### Added

- Initial CMake skeleton.
- `malloy_nbody_terminal` app stub.
- `malloy_smoke_tests` CTest target.
- `CMakePresets.json` for Windows MSVC debug/release.
- Project planning docs.

### Tests

- Smoke test passes through CTest.

### Architecture Notes

- Locked C++20 + CMake + VS Code + MSVC baseline.
- No physics code yet.
```

Earlier entries (M1 to M19) mostly use `Tests` and `Architecture Notes`
sections, as in the M1 example; they are left as written.
