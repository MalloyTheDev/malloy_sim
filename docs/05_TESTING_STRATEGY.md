# 05 - Testing Strategy

## Final decision

MalloySim uses:

```text
CTest + tiny custom CHECK macros
```

Decided for M1-M5 and still in force at M39: Catch2 and GoogleTest are not
added unless explicitly asked for (`CLAUDE.md` rule 9, ADR 0003).

## Why

- zero external dependency friction
- no FetchContent/network failures
- no package manager
- CTest is already available
- simple tests have stayed enough, through ten scenario types
- per-module test executables catch LINK-time dependency leaks

The link step alone has a limit: every module puts the whole `include/` tree on
its consumers' include path, so a header-only leak (using another module's
header without linking it) still compiles. Since issue #13,
`malloy_layering_tests` closes that gap for the libraries: it reads
`CMakeLists.txt` and every library's sources, and fails if a library includes
a module it does not link, if any physics domain links or includes another
(ADR 0006), or if `malloy_sim_core` depends on anything but `malloy_math`. It
does not check the test executables themselves.

## Test executable layout

All fourteen exist and pass (one per module, plus `malloy_quat_tests` for the
quaternion half of `malloy_math` and `malloy_layering_tests`, which checks the
build graph against the module boundaries):

```text
malloy_smoke_tests
malloy_math_tests
malloy_quat_tests
malloy_time_tests
malloy_sim_core_tests
malloy_nbody_tests
malloy_scenario_tests
malloy_ascii_tests
malloy_collide_tests
malloy_particles_tests
malloy_rigid_tests
malloy_springs_tests
malloy_charges_tests
malloy_layering_tests
```

Each domain's one executable covers both of its dimensions: the 3D worlds live
in their domain's library (ADR 0009), so their tests do too.

Each test executable should return 0 on success, return nonzero on failure, print useful failure messages, and link only the module it tests.

CTest also runs five app checks, registered after the executables are counted
(`CMakeLists.txt`): the real terminal app is run on an empty file, a
comment-only file, an empty `charges` and an empty `charges3d` scenario, each of
which it must refuse with "defines nothing to simulate", and on a scenario
whose path lies outside the ANSI code page, which it must open. They match the
output only: CTest ignores the exit code of a test with
`PASS_REGULAR_EXPRESSION`, and no CTest yet runs a template through the app
([#65](https://github.com/MalloyTheDev/malloy_sim/issues/65)).

## Test helpers

Use `MALLOY_CHECK_TRUE`, `MALLOY_CHECK_FALSE`, `MALLOY_CHECK_EQ`,
`MALLOY_CHECK_NEAR`, and `MALLOY_CHECK_VEC2_NEAR` (all in `tests/test_check.hpp`).


## Prefer asymmetric configurations

A test whose setup is symmetric can pass against a wrong implementation, because
the symmetry cancels the error. This has already happened twice in this project,
both times found by mutation testing rather than by review:

- `distance(a, b)` was only ever tested with `b` at the origin, where `a - b`
  and `a + b` are identical. An implementation using the wrong operator passed.
- `approx_equal(Vec2)` had no case where x was the sole difference, so an
  implementation checking only y passed everything.

The same shape of mistake hides sign errors, swapped operands, reversed normals,
and reversed impulses. When writing a physics test, deliberately break symmetry:

- unequal masses, unequal radii, unequal extents;
- positions and velocities that are off-axis and off-origin;
- values distinct in every field, so a swapped pair changes the result;
- for anything rotational: an off-centre application point, a body whose origin
  is not its centre of mass, and a nonzero initial angular velocity.

A useful check on a new test is to ask what wrong implementation would still
pass it, and then either add a case that would not, or confirm by mutation.

## Mutation testing

Every milestone since M9 has been checked by mutation: a deliberate defect is put into the code
(a dropped term, a flipped sign, a swapped operand, a skipped validation), the
whole suite is run, and a test must fail. The defect is then reverted. Each
milestone's `CHANGELOG.md` entry lists its mutations and what caught them; M39
ran thirty-seven, all caught, three at compile time.

A mutation that survives is one of two things, and the entry says which. Either
it is a true equivalent (the program's behaviour cannot differ), or it is a
coverage gap, which is closed by a new assertion that is itself confirmed by
re-applying the mutation. The mutations are applied by hand or by a local
harness; no mutation tool is part of the repository, and nothing runs them
automatically.

## Scenario templates and `# check` lines

The shipped templates in `scenarios/` are a deliverable, and
`malloy_scenario_tests` treats them as one. For every template it requires:

- that it parses;
- a `# Expected` section, the documented result;
- at least one machine-checkable line of the form

  ```text
  # check step <n> <quantity> <value> tol <t>
  ```

  where `step 0` is the state before any step and the tolerance is ABSOLUTE and
  required (`CLAUDE.md`, Template library);

and then runs the world the template declares for its FULL length, evaluating
each check at the step it names. A step that fails validation fails the run, so
running to the end is also the NaN check. A check that names a quantity its
domain does not have FAILS rather than being skipped, and so does one that
names a step past the end of the run.

The same test executable guards the documents that state counts: the template,
domain and test-executable counts in `README.md`, the test-executable count in
`docs/00_START_HERE.md`, the `M1-M<n>` milestone range in `CLAUDE.md`,
`README.md`, `docs/00`, `docs/07` and `docs/08` (checked against the newest
`CHANGELOG.md` entry), and the ADR range in `docs/00` and `docs/08` (checked
against the newest file in `docs/decisions/`).

Known limits of the harness, all open:

- A `#check` line with no space after the `#`, and a check at a negative step
  (unless it is the only one), are skipped without a word; and in the 3D rigid
  domain `energy` is the kinetic energy only and there is no `kinetic`
  ([#63](https://github.com/MalloyTheDev/malloy_sim/issues/63)).
- Several templates' checks would still pass if their headline behaviour broke,
  for example a box that falls through the floor in `box_drop.scn`
  ([#67](https://github.com/MalloyTheDev/malloy_sim/issues/67)).
- The checks run only for files in `scenarios/`, through a copy of the app's
  dispatch in the test; a user's own scenario's checks run nowhere
  ([#64](https://github.com/MalloyTheDev/malloy_sim/issues/64), a milestone
  candidate).

## Statistical checks (M38)

A gas is chaotic, but a run is deterministic (docs/04), so a check on a
statistical quantity returns the same number every time on the same binary and
platform. That makes it tempting to set the tolerance from the observed value.
Do not: a tolerance fitted to one run proves nothing about the physics.

Instead, every tolerance on a statistical quantity is DERIVED and the derivation
is written beside the check:

- three standard deviations of the quantity's own scatter, from its sampling
  distribution where one is known (binomial for a fraction, Beta for an energy
  share), or measured as the run-to-run spread when the timestep is perturbed
  slightly (which reshuffles every collision) for a time average like pressure;
- plus an explicit budget for model error when the comparison is with an
  approximate theory (a first-order wall correction, a finite step);
- and it must EXCLUDE the obvious wrong answer. Say what it excludes: the ideal
  gas pressure, the unrelaxed start, the starting energy ratio.

Exact quantities keep tight tolerances even in a statistical template: a
conserved temperature, or a statistic of the deterministic start.

Since M39 some quantities are random, drawn from a seeded generator. A
template states its seed, and its tolerances come from the spread measured
across many seeds (40 to 60 for the M39 templates), never from the one seed it
ships with. The generator itself is pinned by known answers rather than
statistics, at compile time and at run time, and the draw order by exact
comparison with the public samplers.

M39 also re-ran every statistical check on the real code across 1000 to 5000
seeds other than its own. Where a quantity turned out heavy-tailed, so that
three standard deviations passed fewer than 99.73 per cent of seeds, the
tolerance became the measured 99.73 per cent quantile instead. A shipped seed
should be typical; M39's are within 1.3 standard deviations of their ensembles
on every check.

A bound can replace a standard deviation when the error has one. The point
particle ideal-gas test computes its own tolerance from the run (whole bounces
counted, plus the clamp lag), and also asserts that tolerance is tight enough to
mean something.

## Fuzzing

The scenario file is the one untrusted input. Besides the hostile-input battery
in `malloy_scenario_tests` (control bytes, non-finite and overflowing numbers,
enormous tokens and lines, the resource caps at their boundaries, all 256 byte
values through every place a token is echoed), two libFuzzer targets run the
parser and the whole app under ASan and UBSan. They are built only with clang
and run only by hand; `docs/06` has the recipe.

## When to migrate to Catch2

Only if explicitly asked (`CLAUDE.md` rule 9). The original trigger, test count,
diagnostics, parameterization or CI reporting becoming painful, has not been
reached.
