# ADR 0013 - The Roadmap After M41 Is Gated: Publication, Correctness, Library

## Status

**Proposed. Not accepted.** This ADR records the owner-reviewed roadmap for the
work after M41. It stays Proposed, in the repository or not, until the owner
explicitly accepts it, and it authorizes nothing. In particular it does not
authorize:

- a push or a merge;
- an issue;
- an extraction;
- the start of any milestone or fix.

**M41 is merged** (pull request #87, merge commit `2910de8`, 2026-10-12), so
its description of M41 as merged holds.

**Still the owner's to decide:**

- the package name and the version policy (section 6);
- the open items under "Consequences".

## Context

M41 finished the first rendering milestone. The next question is no longer
"what can be added" but "in what order, and what belongs here". Three inputs
bear on it:

1. **An expansion audit and three research reports,** reconciled against the
   repository on 2026-10-12. The reports covered foundations, disciplines, and
   platform and consumers. Together they proposed:
   - about ten foundation libraries;
   - a dozen new sciences;
   - packaging, language bindings and coupling.
2. **A master research specification** (revision 2), which states the
   mathematics each proposal would need.
3. **Measurements on the local M41 build** (`35bdde6`, whose code is the
   merged `main`'s), run on 2026-10-12:
   - **Rigid 3D rotation.** Free rotation can gain energy without bound,
     because the angular update is forward Euler: +16.6% over 1000 s at dt
     1e-3, and a step failure at 627.6 s at dt 1/60.
   - **Springs.** The spring timestep check is per spring, so a network it
     accepts can diverge silently.
   - **Fast pairs** pass through each other between steps.

The pressures to resist are the ones the rules already name:

- building foundations before anything consumes them (rules 10 and 11);
- an engine kernel or a coupling runtime (rule 12, ADR 0006);
- numbering work before it starts;
- publishing a versioned library whose validation and template numbers are
  about to change.

## Decision

### 1. Three gates, in order

| Gate | Passes when |
|---|---|
| **Publication** | M41 is reviewed; its Linux CI job has run green on a branch and pull request; the owner approves the merge. Manual input and HiDPI checks stay explicitly open until done. **Passed for M41 on 2026-10-12:** CI was green on pull request #87 and again on `main` after the owner merged it (`2910de8`), and the manual input and HiDPI checks were done (`CHANGELOG.md`, [M41], Verified after merge). |
| **Correctness** | The measured numerical defects are fixed, each as a separately reviewed change with regression tests and updated scientific documentation: the rigid 3D rotation update, the spring network timestep check, #86 (physical inertia), #67 and #63 (checks that cannot detect what they document), and #57 (box and sphere contact point). Fast-pair tunneling gets a scoped decision (a documented envelope, a margin, or continuous detection). |
| **Library** | #82 meets the acceptance criteria in section 5, including the minimal external consumer test: a separate CMake project that installs the packages, finds them with `find_package`, links only `malloy::nbody`, and reproduces a verified N-body test case directly through the installed public API, without the scenario dispatcher or the repository source tree. |

**After the gates,** work is scheduled by dependency, not in a fixed sequence:

- **A real application integration.** One of the owner's own projects
  consumes the installed packages and reports what they lack. This is distinct
  from the library gate's minimal consumer test.
- **The shared spatial index,** a design candidate for #79 and #85 (section 3).
- **Small extensions of existing domains,** each when a demonstrated dependency
  or consumer requires it.
- **Continuum heat,** the preferred first new-domain candidate, subject to
  separate approval.

**The GUI (#81) is not scheduled by this ADR.** It is an application milestone
that needs its own authorization.

### 2. Roadmap status is not authorization

This ADR records priority and classification, not approval. Every milestone,
fix or extraction still needs:

- its own scope report (section 4);
- the owner's explicit go before it starts;
- its own acceptance after it is done.

| Capability | Roadmap status |
|---|---|
| Continuum heat | Preferred first new-domain candidate, pending scope approval |
| Spatial index | Shared-design candidate for #79 and #85 |
| Existing-domain extensions (#74, orbital elements, the piston #75, joints) | Candidates, scheduled by demonstrated dependency or consumer |
| Waves | Later candidate |
| SPH | Deferred until neighbour-search needs are established |
| FEM | Deferred; rotational objectivity must be validated for finite-rotation use |
| GUI (#81) | Independently authorized application milestone |
| Units, numerics, linear algebra, fields, statistics, random streams | Not candidates until a consumer triggers them |
| Ecology, epidemiology, planetary systems | Potential sibling projects, not implicitly approved |
| Quantum | Separate (rule 7) |

### 3. Every capability is classified before it is scheduled

**Foundation:** reusable mathematics, numerics, spatial, or testing and
validation support, with no physics of its own.

- **The rule:** a foundation library is extracted only when a second real
  consumer demonstrates the need. Until then the code lives inside its first
  consumer (rules 10 and 11, ADR 0008, ADR 0010).
- **The spatial index.** #79 and #85 are two PLANNED consumers, not two
  implemented ones. A shared implementation is justified when both
  demonstrate compatible requirements in their scope reports. The extraction
  trigger has not been met yet.
- **Shared testing and validation infrastructure:** #64, one `# check`
  evaluator for the scenario tests and the front ends. It is not a scientific
  application.
- **Waiting for a trigger:** units, numerics, linear algebra, fields,
  statistics and random streams.
- **Not candidates:** signals, optimization and graphs, until a consumer
  exists.

**Scientific domain:** an independent library with its own world, settings,
validation, invariants and templates (rule 16), or an extension or
measurement of an existing domain (the M12 and M38 pattern).

- **Existing-domain extension candidates:** #74, orbital elements, the piston
  (#75) and joints.
- **Continuum heat, the preferred first new-domain candidate.** It is preferred
  because its scheme has an exact discrete reference.
  - Under forward-time centred-space differencing, a discrete sine mode of
    wavenumber k is multiplied every step by
    G(k) = 1 - 4 r sin^2(k dx / 2), with r = alpha dt / dx^2.
  - That checks the implemented recurrence exactly, to rounding. It is not an
    exact reproduction of continuous diffusion, which needs its own
    convergence study.
- **Waves, a later candidate.** For the standard centred one-dimensional scheme
  at Courant number 1, a consistently initialized wave travelling one way on
  a periodic (or effectively unbounded) grid is reproduced exactly: it shifts
  one cell per step. Boundaries and other initial conditions need their own
  tests.
- **SPH:** deferred until neighbour-search needs are established.
- **FEM:** deferred. A future FEM proposal must show rotational objectivity
  (no stress under a rigid rotation) within its stated validity regime, and
  validate it wherever finite rotations are in scope.

**Application or integration:** a consumer of the libraries.

- **In this repository:** the viewer, the GUI and packaging.
- **Separate repositories:** language bindings.
- **Example applications, never a runtime (ADR 0006):** couplers between
  domains. Reading each domain's accumulated heat, work or impulse tallies is
  diagnostics only. Physically conservative coupling would need separately
  approved exchange timing, state and control interfaces, and conservation
  checks, none of which the current domain APIs provide.
- **Potential sibling projects** on the packaged libraries: ecology,
  epidemiology, economics, planetary generators and other population or agent
  models. They are valuable and mathematically checkable, but they have
  different ownership and lifecycles. None is implicitly approved.
- **Quantum** stays separate (rule 7).

### 4. Evidence before code, and verification after it

**Before implementation:** a scope report, reviewed before any code. It covers:

- the need, the named consumer, and what is excluded;
- the evidence: commit, paths, equations and unknowns, each with its
  evidence class (run locally, reported, derived, inferred);
- the continuous model and its discretization, with signs, units and frames;
- the predicted guarantees and limits: stability, accuracy, edge cases and
  invariants;
- the files, APIs, dependencies and licences affected;
- the test plan, with an exact or analytic oracle and a mutation or
  frozen-system check;
- a numerical budget;
- its ADR impact, and the alternatives rejected;
- the review gate.

**After implementation, before acceptance:** a verification record. It covers:

- the exact commands and their results;
- the mutation outcomes;
- every changed template number, with its reason;
- the limitations found.

**Throughout:**

- A reported number is not adopted until it is reproduced.
- Verification is not validation.
- Statistical tolerances are derived, never fitted (`docs/05`).

**The numerical claims this ADR relies on, by evidence class:**

- **Proved** (derivations, under their stated assumptions):
  - forward Euler on the torque-free Euler equations never decreases |L|:
    |L_{n+1}|^2 - |L_n|^2 = h^2 |L_n x w_n|^2 >= 0;
  - the two undamped spring-network bounds are sufficient upper estimates of
    the local maximum stiffness eigenvalue;
  - a splitting of the free rigid body into exact sub-flows conserves the
    world-frame angular momentum in exact arithmetic.
- **Measured** (local runs, 2026-10-12, on the cases tested only):
  - +16.6% energy over 1000 s at dt 1e-3, and a step failure at 627.6 s at dt
    1/60, for a free asymmetric top;
  - a spring network accepted by the per-spring check that diverges while
    every step reports Ok;
  - two particles on a collision course passing through each other in one
    step;
  - the splitting's energy error staying bounded and second order, and the
    body-frame midpoint moving the world-frame angular momentum.
- **Unverified proposals:**
  - a stability condition for damped spring networks;
  - the choice of rotation integrator;
  - the tunneling policy;
  - the package name and the linkage and version policy.

### 5. The library gate's acceptance criteria (#82)

These criteria are recorded here so that the gate depends on no document
outside the repository. #82's scope report may tighten them, not loosen them.

1. **Target isolation.** The link closure of `malloy::nbody` contains no
   other domain, no scenario layer, no raylib and no viewer code.
2. **Relocatable headers.** Include paths are given as build-tree and install
   paths; no source or build path appears in the installed usage
   requirements.
3. **Versioned discovery.** A package configuration file, a version file and
   exported targets, with the supported compilers and C++ version stated.
4. **A clean consumer.** A separate project builds against the installed
   prefix only: no `add_subdirectory` and no source-tree headers. It links
   only `malloy::nbody`. It constructs and steps an N-body world through the
   installed public API, and reproduces a verified N-body test case with its
   reference results. It does not depend on the scenario dispatcher, which
   criterion 1 excludes, or on the repository source tree.
5. **Relocation.** After the installed prefix is moved, the consumer still
   configures, builds and runs.
6. **Offline by default.** Configuring and consuming a library never
   downloads raylib and never needs an OpenGL context.
7. **Library linkage policy.** The installed package implements the
   owner-approved static/shared linkage policy, with correct transitive
   dependencies in Debug and Release. Static-only distribution is the initial
   recommendation, pending approval.
8. **Compiler flags.** Consumers on MSVC, GCC and Clang do not inherit the
   project's warnings-as-errors or hardening flags, and the floating-point
   contraction policy for inlined headers is documented.
9. **Documentation and tests.** Install instructions, licence notices, a
   consumer example and version-bump rules exist, and Debug and Release
   consumer tests run in ordinary CI.

### 6. Names, numbers and versions

**No milestone number is assigned before the work starts** (docs/07, Tracks 2
to 4).

**The CMake package name is pending the owner's decision.** `malloy` is
recommended, because it matches:

- the `malloy::` namespace;
- the `<malloy/...>` include root;
- the tested packaging prototype.

`MalloySim` is the alternative. Exactly one name is chosen, with no alias, and
#82 records it when it starts.

**The version policy is a separate pending decision.** The recommendation is:

- pre-1.0 (0.x) versions;
- static libraries only at first;
- compiler flags kept private;
- consumer compatibility tests.

### 7. What this ADR does not change

ADR 0006 stands: every physics domain stays its own library with its own
concrete world, settings and tests, selected by the scenario `type` key with a
plain switch. Nothing here adds a simulation base class, an engine kernel, a
scheduler, a coupling runtime or a speculative foundation library (rules 10 to
13), and no milestone number is chosen.

## Rejected alternatives

- **Foundation libraries first** (the audit's catalog of about ten). It would
  create empty or single-consumer libraries, which rules 10 and 11 forbid.
- **Packaging (#82) before correctness.** It would publish an API whose
  validation (#86) and template numbers (the rotation and spring fixes) are
  about to change.
- **A coupled-domain runtime or engine kernel.** That is rule 12 and ADR 0006.
  Couplers are example applications. Reading each domain's own tallies serves
  diagnostics; conservative coupling would need its own approved exchange
  timing, interfaces and conservation checks.
- **Ecology and planetary systems inside this repository.** They are better
  served as potential sibling projects on the packaged libraries.
- **A numbered multi-milestone plan.** Numbers written ahead have collided with
  what shipped four times (docs/07).

## Consequences

**Template numbers will change.** The correctness gate re-derives the late-step
rotation checks in `intermediate_axis.scn` and `gyroscope.scn` (the latter also
for #86), and any other template whose numbers move, each with its reason. The
spring check may refuse new inputs, though no shipped template. Each change is
recorded in `CHANGELOG.md` with its reason, and terminal and viewer parity are
re-checked after each.

**#82 starts later** than the audit proposed, and on a stabler API.

**docs/07 carries the gates and the classification** so that a reader sees the
order without reading the research.

**Open items:**

- **Spring network stability.** The two undamped bounds provide sufficient
  upper estimates of the local maximum stiffness eigenvalue under the stated
  assumptions. They do not establish stability against negative-stiffness
  modes, nonlinear state changes, or arbitrary damping. The damped-network
  stability condition remains unproved. A timestep gate may use the smaller
  of the two undamped bounds.
- **The choice of rotation integrator** is the scope report's. The candidates
  are splitting (preferred) and Simo and Wong's energy-momentum method.
- **The tunneling policy** is the scope report's. If a warning-only envelope
  is chosen, the limitation is stated wherever the domain is documented.
