# ADR 0001 — C++20 + CMake + MSVC

## Status

Accepted.

## Decision

MalloySim starts with C++20, CMake, MSVC on Windows, and VS Code + CMake Tools.

## Rationale

C++20 is mature enough across MSVC/GCC/Clang and is appropriate for a long-term simulation/engine project.

CMake provides a cross-platform build foundation without locking the project to Visual Studio.

MSVC is the first compiler because the primary development machine is Windows.
