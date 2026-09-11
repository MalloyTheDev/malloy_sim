#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include <malloy/nbody/body2d.hpp>
#include <malloy/nbody/body3d.hpp>
#include <malloy/nbody/nbody_settings.hpp>
#include <malloy/particles/particle2d.hpp>
#include <malloy/charges/charge_settings.hpp>
#include <malloy/charges/charged_particle2d.hpp>
#include <malloy/particles/particle_settings.hpp>
#include <malloy/math/real.hpp>
#include <malloy/rigid/rigid_body2d.hpp>
#include <malloy/rigid/rigid_body3d.hpp>
#include <malloy/rigid/rigid3d_world.hpp>
#include <malloy/rigid/rigid_settings.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <malloy/springs/springs.hpp>

namespace malloy::scenario
{
// A complete, runnable N-body scenario: the fixed timestep, the N-body
// parameters, the bodies, and how long to run. This is the parsed form of a
// scenario file.
// Which domain a scenario describes. Selected by the `type` key, which
// defaults to NBody when absent so every scenario written before this key
// existed keeps working unchanged.
enum class ScenarioType
{
    NBody,
    Particles,
    Rigid,
    Springs,
    Charges,
    NBody3D,
    Rigid3D,
};

// A complete, runnable scenario for one domain.
//
// Only the fields belonging to `type` are meaningful; the rest keep their
// defaults. This is a plain tagged struct on purpose (ADR 0006): a base class
// with a virtual load() is exactly what this project does not build, and the
// cost of a few unused fields is far lower than the abstraction it replaces.
struct Scenario
{
    ScenarioType type{ScenarioType::NBody};

    sim_core::SimulationSettings simulation{0.001};
    int steps{1000};
    int output_every{100};

    // type == NBody
    nbody::NBodySettings nbody_settings{1.0, 0.0};
    std::vector<nbody::Body2D> bodies;

    // type == Particles
    particles::ParticleSettings particle_settings{};
    std::vector<particles::Particle2D> particle_list;

    // type == Rigid
    std::vector<rigid::RigidBody2D> rigid_bodies;
    rigid::RigidSettings rigid_settings{};

    // type == NBody3D
    std::vector<nbody::Body3D> bodies3d;

    // type == Rigid3D
    std::vector<rigid::RigidBody3D> rigid_bodies3d;
    rigid::Rigid3DSettings rigid3d_settings{};

    // type == Charges
    std::vector<charges::ChargedParticle2D> charge_list;
    charges::ChargeSettings charge_settings{};

    // type == Springs
    std::vector<springs::SpringBody2D> spring_bodies;
    springs::SpringNetwork spring_network;
};

// The outcome of parsing. On failure `error` holds a human-readable message
// (with the offending line number) and `scenario` is unspecified. Parsing
// never throws (docs/04: status/result returns for normal failures).
struct ScenarioParseResult
{
    bool ok{false};
    std::string error;
    Scenario scenario;
};

// Parse a scenario from a text stream. One directive per line; '#' starts a
// comment; blank lines are ignored:
//
// One comment form is given meaning elsewhere: a line starting `# check` is a
// machine-checkable assertion, evaluated by the scenario tests against a real
// run (CLAUDE.md, Template library). The parser ignores it like any other
// comment.
//
// Common to every domain:
//
//   type <nbody|particles|rigid|springs|charges|nbody3d|rigid3d>
//                           which domain           (default nbody)
//   dt <value>              fixed timestep          (default 0.001)
//   steps <value>           number of steps         (default 1000)
//   output_every <value>    steps between reports   (default 100)
//
// type nbody:
//
//   g <value>             gravitational constant   (default 1.0)
//   softening <value>     softening length         (default 0.0)
//   body <mass> <px> <py> <vx> <vy>
//
// type nbody3d:
//
//   Newtonian gravity in three dimensions. Takes the same g and softening keys
//   as type nbody, because they mean exactly the same thing in either
//   dimension, and differs only in the body line.
//
//   body3 <mass> <px> <py> <pz> <vx> <vy> <vz>
//
// type rigid3d:
//
//   Torque-free rotation of rigid bodies in three dimensions. Takes no
//   restitution, gravity, friction or ground keys, because M20 has no contacts
//   and no forces: a body here tumbles and flies straight, and that is all.
//
//   rigid_body3d <mass> <Ix> <Iy> <Iz>
//                <px> <py> <pz>
//                <axisx> <axisy> <axisz> <angle>
//                <vx> <vy> <vz>
//                <wx> <wy> <wz>
//
//     Seventeen fields, written above in the five groups they form, though the
//     line itself is a single line like every other.
//
//     Ix, Iy and Iz are the PRINCIPAL moments of inertia about the centre of
//     mass, so the body frame is already the one in which the inertia is
//     diagonal. Three equal moments are a sphere, two a symmetric top, three
//     distinct ones a body that can tumble.
//
//     The orientation is given as an axis and an angle in radians rather than
//     as a quaternion, because a hand-written quaternion is almost never a unit
//     one and an invalid body is a worse error message than a converted one.
//     The axis is normalized on load and need not be unit; a zero axis means no
//     rotation, whatever the angle says.
//
//     wx, wy and wz are the angular velocity in the BODY frame, not the world
//     frame. That is the frame in which Euler's equations are diagonal, and for
//     a body loaded unrotated the two coincide.
//
//   torque <tx> <ty> <tz>   a constant WORLD-frame torque on every body
//                           (default 0 0 0, the torque-free M20 case)
//
//     World frame, not body frame: an external couple fixed in the lab. Under
//     it each body's world-frame angular momentum grows along the straight line
//     L(t) = L(0) + torque t. There is one torque for the whole world, the way
//     there is one gravity in the 2D rigid domain.
//
// type particles:
//
//   restitution <value>   bounciness in [0, 1]     (default 1.0)
//   gravity <gx> <gy>     uniform acceleration     (default 0 0)
//   bounds <minx> <miny> <maxx> <maxy>             (default -1 -1 1 1)
//   particle <mass> <radius> <px> <py> <vx> <vy>
//
// type rigid:
//
//   restitution <value>   bounciness in [0, 1]     (default 1.0)
//   gravity <gx> <gy>     uniform acceleration     (default 0 0)
//   friction <mu>         Coulomb coefficient      (default 0, frictionless)
//
//     Not capped at 1: a coefficient above 1 is physically real. Only a
//     negative or non-finite one is refused.
//
//   ground <nx> <ny> <offset>   an immovable infinite plane
//
//     The plane is the line dot((nx, ny), p) == offset, and (nx, ny) points
//     OUT of the solid side, into free space. A floor at y = -2 is
//     `ground 0 1 -2`; a wall at x = 6 that keeps bodies to its left is
//     `ground -1 0 -6`. The direction need not be unit length here: it is
//     normalized once on load, and only a zero or non-finite one is refused.
//     May appear more than once.
//
//   rigid_body <mass> <inertia> <radius> <comx> <comy> <px> <py> <angle> <vx> <vy> <omega>
//
//     mass and inertia are about the centre of mass; radius is the collision
//     disc centred on the body ORIGIN, and zero means the body does not
//     collide; comx/comy is the local offset from the body origin to the
//     centre of mass; angle is in radians. An infinite mass and inertia mean
//     an immovable body.
//
//   rigid_static <radius> <px> <py> <angle>
//
//     An immovable body: infinite mass and inertia. A separate key because
//     MSVC's stream extraction does not accept the token "inf".
//
// type charges:
//
//   coulomb <k>           Coulomb constant         (default 1.0)
//   efield <ex> <ey>      uniform electric field   (default 0 0)
//   bfield <b>            uniform magnetic field, out of plane (default 0)
//   softening <value>     enters the denominator SQUARED       (default 0)
//                         rejected above about 1.34e154, where the square
//                         would overflow
//   charge <mass> <q> <px> <py> <vx> <vy>
//
//     q is SIGNED, and may be zero: a neutral particle is carried by the
//     fields of others while contributing nothing to them.
//
//     bfield is a scalar rather than a vector because in two dimensions only
//     the out-of-plane component of B produces an in-plane force.
//
// type springs:
//
//   spring_body <mass> <px> <py> <vx> <vy>
//   spring <a> <b> <rest_length> <stiffness> <damping>
//
//     a and b are indices into the spring_body list, in the order they appear.
//     Springs are evaluated in the order they are declared, which is part of
//     the observable behaviour (ADR 0008).
//
// A key belonging to another domain is a parse error, so a typo in `type`
// surfaces immediately rather than silently running the wrong simulation.
//
// Only syntax is checked here: bad numbers, too few or too many fields on a
// line, unknown keys, and run lengths that cannot be run at all (a negative
// `steps` or `output_every`). Trailing tokens are an error rather than being
// ignored, so a body line written with 3D fields is rejected instead of
// silently running a different simulation.
// Semantic validity (dt > 0, mass > 0, ...) is left to the target domain's own
// World::validate, so each rule lives in exactly one place.
ScenarioParseResult parse_scenario(std::istream& input);

// Convenience: open a file and parse it. A missing or unreadable file is
// reported as a parse failure, not an exception.
ScenarioParseResult parse_scenario_file(const std::string& path);
} // namespace malloy::scenario
