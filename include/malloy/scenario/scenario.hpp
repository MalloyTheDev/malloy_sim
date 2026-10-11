#pragma once

#include <cstddef>
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
//
// A message from parse_scenario is plain printable ASCII, safe to print to a
// terminal whatever the file contained (issue #28). Where it repeats a token
// from the file, the token is quoted, every byte outside printable ASCII and
// every backslash and quote in it is written as \xNN (two lowercase hex
// digits), and only its first 64 bytes are shown, followed by its full length:
// `unknown key '\x1b[2Kbody'`, or `unknown key '<first 64 bytes>'... (300
// bytes)` for a longer one. The path in parse_scenario_file's open error is
// the caller's, not the file's, and is repeated as given.
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
// Encoding: a scenario file is UTF-8 or ASCII; every key is ASCII, and other
// UTF-8 text may appear only in comments. A UTF-8 byte-order mark at the very
// start of the input is skipped, since Windows editors write one. A UTF-16
// byte-order mark (FF FE or FE FF) there is refused with an error that names
// the encoding. Anywhere else the mark's bytes are ordinary text, so a second
// mark, or one at the start of a later line, is an error that shows its bytes
// (issue #32).
//
// Common to every domain:
//
//   type <nbody|particles|rigid|springs|charges|nbody3d|rigid3d|charges3d|
//         springs3d|particles3d>  which domain     (default nbody)
//   dt <value>              fixed timestep          (default 0.001)
//   steps <value>           number of steps         (default 1000)
//   output_every <value>    steps between reports   (default 100)
//
//     A report is printed at step 0 and then every `output_every` steps; 0
//     prints only the step-0 report, which is how a long run stays quiet.
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
//   (M21) and a constant force (M28), gravity and restitution contacts against
//   ground planes (M22),
//   Coulomb friction (M23, the `friction` key below), and collisions between
//   bodies: spheres (M24), oriented boxes (M33) and a box against a sphere
//   (M34). A normal contact on a centred sphere imparts no spin; friction and
//   off-centre box contacts do.
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
//     exactly as for rigid_body3d. A box collides with ground planes, with
//     other boxes and with spheres.
//
//   torque <tx> <ty> <tz>   a constant WORLD-frame torque on every body
//                           (default 0 0 0, the torque-free M20 case)
//
//     World frame, not body frame: an external couple fixed in the lab. Under
//     it each body's world-frame angular momentum grows along the straight line
//     L(t) = L(0) + torque t. There is one torque for the whole world, the way
//     there is one gravity in the 2D rigid domain.
//
//   force <fx> <fy> <fz>    a constant WORLD-frame force on every body
//                           (default 0 0 0, no force)
//
//     The translational half of the wrench whose rotational half is `torque`
//     (M28), read exactly as `torque` is: one force for the whole world, every
//     body feeling all of it, and a later line replacing an earlier one. It
//     acts through each centre of mass, so it turns nothing, and it is applied
//     as the acceleration F/m, so unlike `gravity3` a heavier body speeds up
//     less: each body's momentum grows along the straight line p(t) = p(0) +
//     force t, whatever its mass. Its potential, the sum of -(force . r), is
//     part of the energy the app prints.
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
//   friction <mu>         Coulomb coefficient      (default 0, frictionless)
//   bounds <minx> <miny> <maxx> <maxy>             (default -1 -1 1 1)
//   particle <mass> <radius> <px> <py> <vx> <vy>
//   lattice <nx> <ny> <mass> <radius> <speed>
//
//     A whole gas in one line (M38): nx by ny identical particles centred in
//     the cells of a grid that tiles `bounds`, all at `speed`, in exactly
//     opposite pairs, so the net momentum is zero for an even count; an odd
//     count leaves one particle unpaired (particles::lattice_gas). It is
//     expanded after the file is read, so `bounds` may come before or after it,
//     and its particles follow any `particle` lines. At most one lattice line;
//     counts must be at least 1, the total is capped at 100000 BEFORE anything
//     is allocated, and a grid spacing below one particle diameter is refused
//     because the particles would start overlapping.
//
//   wall_temperature <wall> <T>      a thermal wall           (default 0)
//
//     wall is min_x, max_x, min_y or max_y, the names of the impulse and heat
//     tallies. A wall at T > 0 is THERMAL (M39): it forgets how a particle
//     arrived and re-emits it with a velocity drawn from the distribution a
//     wall at temperature T emits, so heat flows between the wall and the gas.
//     Restitution and friction do not apply at a thermal wall. T = 0, the
//     default, is the ordinary specular wall. Each wall may be named once. A
//     negative T parses and is refused by validate(); a value that is not a
//     finite number does not parse.
//
//   seed <n>                         random seed              (default 0)
//
//     Seeds the world's random number generator (M39), a whole number from 0
//     to 18446744073709551615 in plain decimal digits, with no sign. Only a
//     thermal wall draws random numbers, so without one the seed changes
//     nothing. The same seed gives the same run (docs/04), and every template
//     with a thermal wall states its seed.
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
//   rigid_box <mass> <hx> <hy> <px> <py> <angle> <vx> <vy> <omega>
//
//     A body that collides as an ORIENTED BOX rather than a disc (M36), of
//     half-widths hx and hy along its own axes, centred on the body origin,
//     which is also its centre of mass. Its inertia is COMPUTED from the box
//     and mass, m ((2 hx)^2 + (2 hy)^2) / 12, not given, so shape and inertia
//     always agree. A box collides with other boxes (M36) and rests on the
//     ground planes (M37). Box against a disc is not built yet, so a file that
//     declares both a `rigid_box` and a disc that collides (a `rigid_body` or
//     `rigid_static` with a radius above 0) is refused, naming the line where
//     the second kind first appears, rather than letting the two pass through
//     each other unseen (issue #48). A body of radius 0 collides with nothing
//     and may share a file with boxes.
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
//   the restitution and friction keys are shared, and the box and gravity take
//   three components.
//
//   restitution <value>              bounciness in [0, 1]        (default 1.0)
//   friction <mu>                    Coulomb coefficient         (default 0)
//   gravity3 <gx> <gy> <gz>          uniform acceleration        (default 0 0 0)
//   bounds3 <minx> <miny> <minz> <maxx> <maxy> <maxz>            (default the
//                                    unit box -1 -1 -1 to 1 1 1)
//   particle3 <mass> <radius> <px> <py> <pz> <vx> <vy> <vz>
//   lattice3 <nx> <ny> <nz> <mass> <radius> <speed>
//
//     The 3D lattice: an nx by ny by nz grid, with the same rules as `lattice`
//     (particles::lattice_gas3d).
//
//   wall_temperature <wall> <T>      a thermal wall           (default 0)
//   seed <n>                         random seed              (default 0)
//
//     As for type particles (M39), with six walls: min_x, max_x, min_y, max_y,
//     min_z and max_z.
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
//
// A scenario file is the one untrusted input, so the parser also bounds what a
// file may ask for (the 2026-09-13 security review, and issues #25, #26, #29):
//
//   steps                at most 1e8
//   entities             at most 1e5 bodies, particles, charges, spring bodies
//                        and springs in total, checked as each line is read
//   ground planes        at most 1000 `ground` or `plane3` lines
//   work                 steps plus reports, times the work a step does (a
//                        report's diagnostics cost about a step), at most 1e11
//                        units of about 10 ns of Release time (under 20 minutes
//                        there; up to about three hours in Debug). A step
//                        is N units for its N bodies, N + springs for springs,
//                        plus its pair tests: N(N-1) for nbody, which evaluates
//                        every ordered pair, N(N-1)/2 for the other pairwise
//                        domains, and in the rigid domains each pair and each
//                        body-plane test weighted by kind, as measured (issue
//                        #46): rigid3d box-box 250, box-sphere 25, sphere pair
//                        10, box-plane 40, sphere-plane 4; rigid box-box 40,
//                        any other pair 5, box-ground 20, disc-ground 4. A body
//                        is a box when all its half-extents are positive
//   reports              at most 1e5 (steps / output_every + 1); output_every
//                        0 prints only the step-0 report
//   line length          at most 4096 characters, comments included; the
//                        reader counts bytes, so a leading byte-order mark
//                        leaves line 1 room for 4093
//
// Each refusal is an error naming the line that crossed the limit, except the
// work and report limits, whose inputs can sit on different lines.
ScenarioParseResult parse_scenario(std::istream& input);

// Convenience: open a file and parse it. A missing or unreadable file is
// reported as a parse failure, not an exception.
ScenarioParseResult parse_scenario_file(const std::string& path);

// How many entities a scenario defines: bodies, particles, charges, rigid
// bodies, spring bodies and springs, summed over every domain. The parser caps
// it, and the app refuses to run a scenario where it is zero (issue #31): an
// empty, comment-only or wrong file parses as a default nbody scenario with
// nothing in it, and running that would report success on a broken input. The
// parser itself still accepts such a fragment, because it checks syntax and
// settings, not whether there is something to simulate.
std::size_t entity_count(const Scenario& scenario);
} // namespace malloy::scenario
