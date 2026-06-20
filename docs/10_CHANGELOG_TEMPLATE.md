# 10 — Changelog Template

This is the template/reference. The live changelog is `CHANGELOG.md` at the
repo root; keep it updated with this format.

```markdown
# Changelog

## [Unreleased]

### Added

### Changed

### Fixed

### Removed

### Tests

### Architecture Notes

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
