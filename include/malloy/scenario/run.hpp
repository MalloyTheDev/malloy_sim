#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include <malloy/charges/charges.hpp>
#include <malloy/math/math.hpp>
#include <malloy/nbody/nbody.hpp>
#include <malloy/particles/particles.hpp>
#include <malloy/rigid/rigid.hpp>
#include <malloy/scenario/scenario.hpp>
#include <malloy/springs/springs.hpp>

// The run layer (M41, issue #80): what every front end needs to RUN a parsed
// scenario, so that none of it lives only in one app. It builds the world a
// scenario declares, reads back the quantities a front end reports for it, and
// formats them as the report line the terminal app prints, so the 3D viewer's
// headless mode prints the identical line.
//
// Ten concrete functions of each kind, one per scenario type. There is no world
// base class, no variant of worlds and no registry (rules 11 and 12, ADR 0006):
// each front end keeps its own switch over ScenarioType and calls the one
// function for its case, exactly as the terminal app always has. The readouts
// are sums and selections of each domain's own diagnostics; the run layer has
// no physics formula of its own.
namespace malloy::scenario
{
// --- Building a world ---
//
// Each reads only the fields belonging to its own type, in the order its
// world's constructor takes them; called for a scenario of another type it
// builds an empty world of its own type. The world is NOT validated: the caller
// calls validate() and reports a refusal in its own way.
//
// To start a run again, build the world again from the same Scenario. A world
// rebuilt from bodies() alone is not the same run: the rigid warm-start
// impulses (ADR 0011) and the particle generator (ADR 0010) are world state.
nbody::NBodyWorld make_nbody_world(const Scenario& scenario);
nbody::NBody3DWorld make_nbody3d_world(const Scenario& scenario);
particles::ParticleWorld make_particle_world(const Scenario& scenario);
particles::ParticleWorld3D make_particle3d_world(const Scenario& scenario);
rigid::RigidWorld make_rigid_world(const Scenario& scenario);
rigid::Rigid3DWorld make_rigid3d_world(const Scenario& scenario);
springs::SpringWorld make_spring_world(const Scenario& scenario);
springs::SpringWorld3D make_spring3d_world(const Scenario& scenario);
charges::ChargeWorld make_charge_world(const Scenario& scenario);
charges::Charge3DWorld make_charge3d_world(const Scenario& scenario);

// Whether a particle box has a thermal wall (a wall temperature above 0, M39):
// the case in which a run depends on its seed, so a front end shows it.
bool has_thermal_wall(const particles::ParticleSettings& settings);
bool has_thermal_wall(const particles::ParticleSettings3D& settings);

// --- Reading back ---
//
// One plain struct per world type, holding the quantities a front end reports
// for it, computed exactly as the terminal app has always printed them. A
// quantity that does not exist in the current state is EMPTY rather than 0:
// the pressure before any time has passed (M38), the speed of the first
// particle of an empty world (issue #24), the separation of a world with fewer
// than two bodies, and the heat of a wall that is not thermal (M39).

struct NBodyReadout
{
    std::optional<math::Real> separation; // of the first two bodies
    math::Real energy{};                  // kinetic plus softened potential
    math::Real angular_momentum{};
};

struct NBody3DReadout
{
    math::Real energy{};
    math::Vec3 angular_momentum{};
    math::Vec3 momentum{};
};

// The heat each wall has delivered to the gas, present for a thermal wall only.
struct WallHeatReadout2D
{
    std::optional<math::Real> min_x;
    std::optional<math::Real> max_x;
    std::optional<math::Real> min_y;
    std::optional<math::Real> max_y;
};

struct WallHeatReadout3D
{
    std::optional<math::Real> min_x;
    std::optional<math::Real> max_x;
    std::optional<math::Real> min_y;
    std::optional<math::Real> max_y;
    std::optional<math::Real> min_z;
    std::optional<math::Real> max_z;
};

struct ParticleReadout
{
    math::Real energy{}; // kinetic plus gravitational potential
    math::Vec2 momentum{};
    math::Real temperature{};
    std::optional<math::Real> pressure;
    WallHeatReadout2D heat;
};

struct Particle3DReadout
{
    math::Real energy{};
    math::Vec3 momentum{};
    math::Real temperature{};
    std::optional<math::Real> pressure;
    WallHeatReadout3D heat;
};

struct RigidReadout
{
    math::Real energy{}; // kinetic plus gravitational potential
    math::Real angular_momentum{};
    math::Vec2 momentum{};
};

struct Rigid3DReadout
{
    // Kinetic plus the potentials of gravity and of the applied force (M28,
    // issue #47), the quantity conserved in free flight, summed in that order.
    math::Real energy{};
    math::Vec3 angular_momentum{};
    math::Vec3 momentum{};
};

struct SpringReadout
{
    math::Real kinetic{};
    math::Real elastic{};
    math::Real energy{}; // kinetic plus elastic
    math::Vec2 momentum{};
};

struct Spring3DReadout
{
    math::Real kinetic{};
    math::Real elastic{};
    math::Real energy{};
    math::Vec3 momentum{};
};

struct ChargeReadout
{
    math::Real energy{};
    // The first particle's speed, the column a magnetic field must leave alone.
    std::optional<math::Real> first_speed;
    math::Vec2 momentum{};
};

struct Charge3DReadout
{
    math::Real energy{};
    std::optional<math::Real> first_speed;
    math::Vec3 momentum{};
};

// One overload per concrete world: overload resolution on ten unrelated types,
// not a virtual call.
NBodyReadout readout(const nbody::NBodyWorld& world);
NBody3DReadout readout(const nbody::NBody3DWorld& world);
ParticleReadout readout(const particles::ParticleWorld& world);
Particle3DReadout readout(const particles::ParticleWorld3D& world);
RigidReadout readout(const rigid::RigidWorld& world);
Rigid3DReadout readout(const rigid::Rigid3DWorld& world);
SpringReadout readout(const springs::SpringWorld& world);
Spring3DReadout readout(const springs::SpringWorld3D& world);
ChargeReadout readout(const charges::ChargeWorld& world);
Charge3DReadout readout(const charges::Charge3DWorld& world);

// --- The report line ---
//
// The line the terminal app prints for one report, without its newline:
// `step` and the step index right-aligned in 6 columns, then each quantity
// after its label, right-aligned in 16 columns (15 for the two springs types)
// in scientific notation with 8 digits after the point. The nbody separation
// is the one fixed-point column, in 11 columns, also with 8 digits. A 3D vector
// is printed as its length (|L|, |p|); a 2D one as its two components. An empty
// pressure or first speed prints "-" in its column; an empty separation or
// wall heat prints no column at all, so a run without them prints what it
// printed before they existed. For example, scenarios/two_body.scn at step 0:
//
//   step      0   sep  1.00000000   E_total  -5.00000000e-07   L_total   1.00000000e-06
//
// Formatted in the classic "C" locale whatever the global locale is, so every
// front end prints the same bytes.
std::string report_line(std::int64_t step, const NBodyReadout& r);
std::string report_line(std::int64_t step, const NBody3DReadout& r);
std::string report_line(std::int64_t step, const ParticleReadout& r);
std::string report_line(std::int64_t step, const Particle3DReadout& r);
std::string report_line(std::int64_t step, const RigidReadout& r);
std::string report_line(std::int64_t step, const Rigid3DReadout& r);
std::string report_line(std::int64_t step, const SpringReadout& r);
std::string report_line(std::int64_t step, const Spring3DReadout& r);
std::string report_line(std::int64_t step, const ChargeReadout& r);
std::string report_line(std::int64_t step, const Charge3DReadout& r);
} // namespace malloy::scenario
