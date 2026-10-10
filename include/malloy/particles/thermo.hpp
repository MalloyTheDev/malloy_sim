#pragma once

#include <optional>
#include <vector>

#include <malloy/collide/shapes.hpp>
#include <malloy/collide/shapes3d.hpp>
#include <malloy/math/real.hpp>
#include <malloy/particles/particle2d.hpp>
#include <malloy/particles/particle3d.hpp>
#include <malloy/particles/particle3d_world.hpp>
#include <malloy/particles/particle_world.hpp>

// Kinetic theory over the particle lists (M38): what a box of colliding
// particles looks like as a GAS. Nothing here changes a trajectory. These are
// pure functions of the particle state (and, for pressure, of the wall tally),
// so the same dynamics that M10 and M32 already validated are simply read as
// temperature, pressure and a velocity distribution. Like ballistics in M12,
// this is measurement on an existing domain, not a new one: no new world,
// settings type or scenario `type`.
//
// Units are the project's normalized ones with Boltzmann's constant k_B = 1, so
// temperature is an energy (docs/04).
namespace malloy::particles
{
// Temperature: the mean kinetic energy per degree of freedom, 2 KE / (d N) with
// d = 2 here (d = 3 for the 3D overload). Measured in the lab frame, so a gas
// drifting as a whole counts its drift as heat; a box gas has no net drift once
// the walls have absorbed it. 0 for an empty list.
math::Real temperature(const std::vector<Particle2D>& particles);
math::Real temperature3d(const std::vector<Particle3D>& particles);

// The ratio <v^4> / <v^2>^2 over all particles, a one-number shape of the
// speed distribution. It is exactly 1 when every particle has the same speed
// (a `lattice` start) and is 2 for a 2D Maxwell-Boltzmann gas (5/3 in 3D),
// whatever the temperature, so its climb from 1 is relaxation made visible.
// Unweighted by mass, so those targets hold for a gas of equal masses. 0 for an
// empty list or one where nothing moves.
math::Real speed_ratio(const std::vector<Particle2D>& particles);
math::Real speed_ratio3d(const std::vector<Particle3D>& particles);

// The fraction of particles whose v^2 is below the mean v^2. For a 2D
// Maxwell-Boltzmann gas v^2 is exponentially distributed, so the fraction is
// 1 - 1/e = 0.632; in 3D it is 0.608 (a chi-squared with three degrees of
// freedom, below its mean). Ill-conditioned when every speed is equal, because
// then each comparison is decided by rounding, so read it after relaxation and
// use speed_ratio for the start. 0 for an empty list.
math::Real slow_fraction(const std::vector<Particle2D>& particles);
math::Real slow_fraction3d(const std::vector<Particle3D>& particles);

// Pressure on the walls: the total normal impulse they have delivered, divided
// by the elapsed time and by the total length of wall, the box's perimeter
// 2 (W + H). It is a time average from the start of the run, which is what a
// gauge on the wall reads. For point particles of an ideal gas it obeys
// P A = N T (A = W H) to within a bound: the fixed step's wall clamp lands a
// particle ON the wall and lowers each bounce rate by at most v dt / L, and a
// finite run counts whole bounces (docs/04). Finite radii raise it (excluded
// area).
//
// Empty when `elapsed` is not a positive finite time: before any time has
// passed no impulse has been measured, so no pressure exists rather than a 0.
std::optional<math::Real> pressure(const WallImpulse2D& impulse,
                                   const collide::Aabb& bounds,
                                   math::Real elapsed);

// The 3D sibling: impulse over elapsed time and total wall AREA,
// 2 (W H + H D + D W). An ideal gas of point particles obeys P V = N T.
std::optional<math::Real> pressure3d(const WallImpulse3D& impulse,
                                     const collide::Aabb3& bounds,
                                     math::Real elapsed);

// A deterministic gas start: nx by ny particles of one mass and radius, centred
// in the cells of a uniform grid that tiles `bounds`, every one at the same
// `speed`. Directions come in exactly opposite pairs (particles 2k and 2k+1),
// with pair k at angle k times the golden angle, so neighbouring pairs point
// unrelated ways, the directions cover the circle evenly, and the total momentum
// is exactly zero for an even count (an odd count leaves its last particle
// unpaired, carrying one particle's momentum). No random numbers are involved,
// so the start, and therefore the run, is repeatable (docs/04).
//
// Preconditions, enforced by the scenario parser rather than here: nx and ny are
// at least 1 and the grid spacing is at least 2 * radius, so no two particles
// start overlapping. Mass, radius and speed are left to validate().
std::vector<Particle2D> lattice_gas(const collide::Aabb& bounds, int nx, int ny,
                                    math::Real mass, math::Real radius,
                                    math::Real speed);

// The 3D sibling on an nx by ny by nz grid. Pair k points along azimuth k times
// the golden angle and height z = 1 - 2 frac((k + 1/2) sqrt 2), two independent
// equidistributed sequences, so the directions cover the sphere evenly; pairs
// are exactly opposite, as in 2D.
std::vector<Particle3D> lattice_gas3d(const collide::Aabb3& bounds, int nx, int ny,
                                      int nz, math::Real mass, math::Real radius,
                                      math::Real speed);
} // namespace malloy::particles
