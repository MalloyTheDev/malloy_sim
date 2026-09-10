#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include <malloy/nbody/body2d.hpp>
#include <malloy/nbody/nbody_settings.hpp>
#include <malloy/particles/particle2d.hpp>
#include <malloy/particles/particle_settings.hpp>
#include <malloy/rigid/rigid_body2d.hpp>
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
// Common to every domain:
//
//   type <nbody|particles|rigid|springs>  which domain  (default nbody)
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
// type particles:
//
//   restitution <value>   bounciness in [0, 1]     (default 1.0)
//   gravity <gx> <gy>     uniform acceleration     (default 0 0)
//   bounds <minx> <miny> <maxx> <maxy>             (default -1 -1 1 1)
//   particle <mass> <radius> <px> <py> <vx> <vy>
//
// type rigid:
//
//   rigid_body <mass> <inertia> <comx> <comy> <px> <py> <angle> <vx> <vy> <omega>
//
//     mass and inertia are about the centre of mass; comx/comy is the local
//     offset from the body origin to the centre of mass; angle is in radians.
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
// Semantic validity (dt > 0, mass > 0, ...) is left to NBodyWorld::validate so
// that rule lives in exactly one place.
ScenarioParseResult parse_scenario(std::istream& input);

// Convenience: open a file and parse it. A missing or unreadable file is
// reported as a parse failure, not an exception.
ScenarioParseResult parse_scenario_file(const std::string& path);
} // namespace malloy::scenario
