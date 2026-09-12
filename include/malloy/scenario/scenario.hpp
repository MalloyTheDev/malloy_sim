#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include <malloy/nbody/body2d.hpp>
#include <malloy/nbody/body3d.hpp>
#include <malloy/nbody/nbody_settings.hpp>
#include <malloy/particles/particle2d.hpp>
#include <malloy/particles/particle3d.hpp>
#include <malloy/particles/particle3d_settings.hpp>
#include <malloy/charges/charge3d_settings.hpp>
#include <malloy/charges/charge_settings.hpp>
#include <malloy/charges/charged_particle2d.hpp>
#include <malloy/charges/charged_particle3d.hpp>
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
    Charges3D,
    Springs3D,
    Particles3D,
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

    // type == Particles3D
    particles::ParticleSettings3D particle3d_settings{};
    std::vector<particles::Particle3D> particle_list3d;

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

    // type == Charges3D
    std::vector<charges::ChargedParticle3D> charge_list3d;
    charges::Charge3DSettings charge3d_settings{};

    // type == Springs
    std::vector<springs::SpringBody2D> spring_bodies;
    springs::SpringNetwork spring_network;

    // type == Springs3D
    std::vector<springs::SpringBody3D> spring_bodies3d;
    springs::SpringNetwork spring_network3d;
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
//   type <nbody|particles|rigid|springs|charges|nbody3d|rigid3d|charges3d|
//         springs3d|particles3d>  which domain     (default nbody)
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
//   Rigid bodies in three dimensions: free rotation (M20), a constant torque
//   (M21), and gravity and restitution contacts against ground planes (M22).
//   Friction is not here; a normal contact on a centred sphere imparts no spin,
//   so this domain's contacts are translational until friction gets its own
//   milestone.
//
//   rigid_body3d <mass> <Ix> <Iy> <Iz> <radius>
//                <px> <py> <pz>
//                <axisx> <axisy> <axisz> <angle>
//                <vx> <vy> <vz>
//                <wx> <wy> <wz>
//
//     Eighteen fields, written above in the groups they form, though the line
//     itself is a single line like every other.
//
//     radius is the collision sphere centred on the body's centre of mass.
//     Zero (a body written before M22) means the body does not collide.
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
//   rigid_box3d <mass> <hx> <hy> <hz>
//               <px> <py> <pz>
//               <axisx> <axisy> <axisz> <angle>
//               <vx> <vy> <vz>
//               <wx> <wy> <wz>
//
//     A body that collides as an ORIENTED BOX rather than a sphere (M30), of
//     half-widths hx, hy, hz along its own axes. Its inertia is COMPUTED from
//     the box and mass (M26), not given, so shape and inertia always agree;
//     there is no radius and no Ix Iy Iz. The pose and velocity fields are
//     exactly as for rigid_body3d. A box collides with ground planes only so
//     far; box against box is a later milestone.
//
//   torque <tx> <ty> <tz>   a constant WORLD-frame torque on every body
//                           (default 0 0 0, the torque-free M20 case)
//
//     World frame, not body frame: an external couple fixed in the lab. Under
//     it each body's world-frame angular momentum grows along the straight line
//     L(t) = L(0) + torque t. There is one torque for the whole world, the way
//     there is one gravity in the 2D rigid domain.
//
//   gravity3 <gx> <gy> <gz>   a uniform gravitational ACCELERATION
//                             (default 0 0 0)
//
//     An acceleration, not a force, applied before the position update. The 3D
//     key is separate from the 2D `gravity` because it takes three components.
//
//   restitution <value>       bounciness of every contact, in [0, 1]
//                             (default 1, perfectly elastic)
//
//   friction <mu>             Coulomb coefficient (default 0, frictionless)
//
//     The tangential impulse is clamped to `mu` times the normal impulse, so a
//     sphere sliding on a plane spins up and rolls without slipping at 5/7 of
//     its sliding speed. Not capped at 1. Shared spelling with the 2D `rigid`
//     key, which writes the 2D settings; here it writes the 3D ones.
//
//   plane3 <nx> <ny> <nz> <offset>   an immovable ground plane
//
//     The plane is dot((nx, ny, nz), p) == offset, with the normal pointing OUT
//     of the solid, into free space, exactly like the 2D `ground` halfplane. A
//     floor at z = -2 is `plane3 0 0 1 -2`. The normal need not be unit: it is
//     normalized once on load, and only a zero or non-finite one is refused.
//     May appear more than once, for a corner or a box of walls.
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
// type charges3d:
//
//   Charged particles in three dimensions. The Coulomb term and the electric
//   field mean the same as in 2D; the magnetic field becomes a full VECTOR,
//   because in three dimensions B has a direction and q(v x B) turns the
//   velocity about it while leaving the component along it alone, which is what
//   makes the motion helical (docs/decisions/0009).
//
//   coulomb <k>                   Coulomb constant             (default 1.0)
//   efield3 <ex> <ey> <ez>        uniform electric field       (default 0 0 0)
//   bfield3 <bx> <by> <bz>        uniform magnetic field       (default 0 0 0)
//   softening <value>             enters the denominator SQUARED (default 0)
//   charge3 <mass> <q> <px> <py> <pz> <vx> <vy> <vz>
//
//     q is SIGNED and may be zero, exactly as in the 2D charges domain.
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
// type springs3d:
//
//   Spring networks in three dimensions. The `spring` key is dimension-agnostic
//   and is reused unchanged (a, b, rest_length, stiffness, damping); only the
//   body line gains a third component, so 3D structures (a chain, a cloth, a
//   lattice) are built the same way as 2D ones.
//
//   spring_body3 <mass> <px> <py> <pz> <vx> <vy> <vz>
//   spring <a> <b> <rest_length> <stiffness> <damping>
//
// type particles3d:
//
//   Colliding spheres in a 3D box under gravity. The 3D form of type particles;
//   the restitution key is shared, and the box and gravity take three
//   components.
//
//   restitution <value>              bounciness in [0, 1]        (default 1.0)
//   gravity3 <gx> <gy> <gz>          uniform acceleration        (default 0 0 0)
//   bounds3 <minx> <miny> <minz> <maxx> <maxy> <maxz>            (default the
//                                    unit box -1 -1 -1 to 1 1 1)
//   particle3 <mass> <radius> <px> <py> <pz> <vx> <vy> <vz>
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
