# ADR 0010 - Randomness Is Seeded And Owned By Its Domain

## Status

Accepted, recorded at M39.

## Context

Until M39 nothing in MalloySim was random. Every run was a pure function of
its input, and the determinism definition (`docs/04`) said so: same binary,
platform, initial state, fixed timestep and update order give the same output.
M38 kept it that way on purpose. Its gas starts from a `lattice`, not a random
draw, and `docs/04` added that a future need for randomness would use "a fully
specified generator inside its domain, never the std distributions".

M39 brought that need. A thermal wall forgets how a particle arrived and
re-emits it with a velocity drawn from the distribution of a wall at its
temperature. There is no deterministic substitute: the emission must be
random, or the wall imposes a pattern rather than a temperature.

Randomness threatens two things the project relies on. The first is
repeatability: a run nobody can repeat is a run nobody can debug, and the
scenario tests compare a real run with fixed `# check` values. The second is
honest tests: a statistical check whose tolerance was fitted to one random run
proves nothing, and a check can be made to pass by reseeding until it does.

## Decision

Randomness is seeded, fully specified, owned by the domain that needs it, and
consumed in a fixed order.

**The generator is SplitMix64** (Steele, Lea and Flood, OOPSLA 2014), written
out in `include/malloy/particles/splitmix64.hpp`:

- **Pinned by known answers.** Its integer output is fixed by a known-answer
  test at compile time and at run time; the values match an independent
  implementation, Java's `SplittableRandom`.
- **Uniform deviates.** A deviate is the midpoint of one of 2^52 cells,
  `((x >> 12) + 1/2) 2^-52`, strictly inside (0, 1), so a logarithm of it is
  always finite and never zero.
- **Small state.** One 64-bit word, so rolling back a failed step is one copy.
- **Any seed works.** Every value, 0 included, is a valid seed.
- **Draws can be counted.** After n draws the state is exactly
  `seed + n * gamma`, so a test can count draws.

**What is not used.** No `<random>` engine or distribution, no `rand`, no
`random_device`, no time-based seed, and no global or static generator. The std
distributions are rejected outright: the standard specifies their interfaces,
not their algorithms, so the same seed would give different runs on different
standard libraries.

**One stream per world.**

- The stream is a member of the world, seeded from a setting (`seed` in the
  scenario file), copied with the world, and restored when a step is rolled
  back.
- A failed step therefore leaves no trace in the stream, and the stream's
  position is a function of the successful steps alone.

**Draws are part of the observable behaviour**, like the spring order in ADR
0008:

- **Where:** a draw happens only where a feature needs one (in M39, a bounce
  off a thermal wall).
- **Order:** draws come in a fixed order, by particle index, then by wall axis.
- **Count:** each event takes a fixed number of draws (three per thermal
  bounce), with no rejection sampling and no cached spare deviate.
- **Consequence:** a box without a thermal wall never draws, so every scenario
  written before M39 runs bit-identically.
- **Changing any of these** (a sampler formula, the order, the count) changes
  every thermal template's output, and is a deliberate act.

**The generator lives in the domain.** It sits in `malloy_particles`, not in
`malloy_math` and not in `malloy_sim_core` (rule 13, enforced by the layering
test). When a second domain needs randomness, that is the trigger to move it to
`malloy_math`, as ADR 0008's trigger governs a shared integrator. Until then,
one consumer is not a pattern.

**A second consumer must not share a seed.** A future random feature in the
same world, a random initial state for example, must draw from a stream seeded
differently, or it replays the walls' numbers and correlates the two.

**The determinism definition gains the seed.** The same binary, platform,
initial state, fixed timestep, update order AND seed give the same output. The
integer stream is the same on every platform. The doubles derived from it pass
through `log`, `cos` and `sin`, which are not correctly rounded, so
cross-platform bitwise repeatability is still not claimed (`docs/04`).

**Statistical checks are derived, never fitted** (`docs/05`). A tolerance on a
random quantity comes from its sampling distribution, or from the spread
measured across many seeds, plus an explicit model budget. It is stated beside
the check, and it must exclude the obvious wrong answer. A template states its
seed. Another seed then passes each check about 997 times in 1000 rather than
by construction.

## Rejected alternatives

- **A random-number service in sim_core.** That would grow sim_core with a
  facility one domain uses (rule 13), and invite every domain to reach for it.
- **`std::mt19937_64`.** The engine is specified, but its 2.5 KB of state would
  be copied every step for rollback. It fails BigCrush's linear-complexity
  tests, `std::generate_canonical` can return 1.0, and the distributions on top
  of it are unspecified.
- **xoshiro256\*\* or PCG32.** Both are good generators. Neither is needed:
  xoshiro's authors seed it with SplitMix64 anyway, and PCG32 needs two calls
  per double.
- **A counter-based generator** (Philox, or a hash of seed, tick, particle and
  wall). It needs no rollback, but it needs a multi-key mixing design and more
  code than the problem calls for.
- **A polymorphic wall model** (a `Wall` base class with specular and thermal
  subclasses). Rule 12. A temperature per wall, with 0 meaning specular, is a
  setting, not a hierarchy.

## Consequences

- **A run with a thermal wall depends on its seed.** The app prints the seed in
  the run header whenever some wall is thermal.
- **The tests can see the stream.** `random_state()` shows whether, and how
  often, a world drew, so tests check no-draw, repeatability, rollback and the
  exact draw order against the public samplers.
- **Randomness opens the next pieces of thermodynamics:**
  - accommodation coefficients;
  - random initial states (on their own seed);
  - eventually DSMC-style methods.

  Each is its own milestone.
