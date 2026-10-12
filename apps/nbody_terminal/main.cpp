#include <malloy/ascii/grid.hpp>
#include <malloy/math/math.hpp>
#include <malloy/charges/charges.hpp>
#include <malloy/nbody/nbody.hpp>
#include <malloy/particles/particles.hpp>
#include <malloy/rigid/rigid.hpp>
#include <malloy/springs/springs.hpp>
#include <malloy/scenario/run.hpp>
#include <malloy/scenario/scenario.hpp>
#include <malloy/sim_core/sim_core.hpp>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

// The terminal app (M5, M8): runs a scenario of any of the ten `type`s and
// reports that domain's diagnostics together with an ASCII view of the body
// positions (a 3D run is drawn as its projection onto the xy plane). The
// scenario is loaded from a file given as the single argument, or, with no
// argument, the built-in N-body scenarios run. The app stays dumb: it owns no
// physics, no rendering math, and has no CLI parser (docs/03, rules 14-15).
//
// Since M41 it builds each world, reads each world's diagnostics and formats
// each report line through the run layer (malloy/scenario/run.hpp), which the
// 3D viewer shares, so nothing a second front end needs lives only here. What
// stays is the terminal's own part: the run headers, the ASCII view and its
// framing, and the loop that steps to the end.

using malloy::ascii::fit_viewport;
using malloy::ascii::render;
using malloy::ascii::Viewport;
using malloy::charges::Charge3DWorld;
using malloy::charges::ChargeWorld;
using malloy::math::Real;
using malloy::math::Vec2;
using malloy::nbody::Body2D;
using malloy::nbody::NBodySettings;
using malloy::nbody::NBodyWorld;
using malloy::particles::ParticleSettings;
using malloy::particles::ParticleSettings3D;
using malloy::particles::ParticleWorld;
using malloy::particles::ParticleWorld3D;
using malloy::rigid::RigidBody2D;
using malloy::rigid::RigidWorld;
using malloy::scenario::parse_scenario_file;
using malloy::scenario::ScenarioParseResult;
using malloy::scenario::ScenarioType;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;
using malloy::springs::SpringWorld;
using malloy::springs::SpringWorld3D;

namespace
{
// Size of the ASCII debug view. Roughly 2:1 because terminal cells are about
// twice as tall as they are wide, so a square region of world space reads as
// square on screen.
constexpr int view_width = 61;
constexpr int view_height = 25;

// The positions of any bodies that have one, as the plain point list the ASCII
// view works with.
//
// At M10 this was two copies and a comment saying the duplication was three
// lines and the abstraction would be a template or a base class, so the copies
// won. It is now three copies, which is the threshold this project uses as
// evidence rather than coincidence (ADR 0008). A template over "has a position"
// is not a base class and invents no vocabulary, so rule 12 is untouched.
//
// malloy_rigid keeps its own extractor below, because it emits two points per
// body rather than one and is genuinely a different function.
template <typename Body>
std::vector<Vec2> positions_of(const std::vector<Body>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size());
    for (const Body& body : bodies)
    {
        points.push_back(body.position);
    }
    return points;
}

// Draws the current body positions, first expanding `view` so that every body
// is inside it. The view only ever grows: refitting it from scratch each frame
// would rescale the picture every report and make a stationary body appear to
// drift, while growing keeps the framing steady and still guarantees no body is
// silently clipped. The extents are printed so a change of scale is visible.
void print_view(Viewport& view, const std::vector<Vec2>& points)
{
    const Viewport frame = fit_viewport(points);

    view.min_x = std::min(view.min_x, frame.min_x);
    view.max_x = std::max(view.max_x, frame.max_x);
    view.min_y = std::min(view.min_y, frame.min_y);
    view.max_y = std::max(view.max_y, frame.max_y);

    std::cout << render(points, view, view_width, view_height) << "  view x ["
              << view.min_x << ", " << view.max_x << "]  y [" << view.min_y << ", "
              << view.max_y << "]\n";
}

// Prints one report: the shared report line for the world's current state,
// then the ASCII view of it.
template <typename World, typename Points>
void print_report(std::int64_t step_index, const World& world, Viewport& view,
                  const Points& points_of)
{
    std::cout << malloy::scenario::report_line(step_index, malloy::scenario::readout(world))
              << '\n';
    print_view(view, points_of(world));
}

// Drives one world to completion, reporting every output_every steps. Returns 0
// on success, 1 on step failure.
//
// One copy rather than four. It is a template over the world type, not a base
// class: RigidWorld, ParticleWorld, NBodyWorld and SpringWorld share nothing but
// a step() returning StepResult, and giving them a common base is exactly what
// CLAUDE.md rule 12 forbids. Duck typing here costs nothing and invents no
// vocabulary.
//
// The int64 counter is the reason this is worth factoring at all. With
// steps == INT_MAX an int counter reaches INT_MAX, passes the loop test, and
// overflows on ++, which is undefined behaviour and in practice never
// terminates. That fix was replicated in all four runners but explained in only
// one, so three copies carried a correctness-critical detail with no note
// saying why it mattered.
template <typename World, typename Report>
int drive(const char* title, World& world, const Report& report, int steps,
          int output_every)
{
    report(0);
    for (std::int64_t step = 1; step <= steps; ++step)
    {
        if (!world.step().ok())
        {
            std::cerr << title << ": step " << step << " failed\n";
            return 1;
        }
        if (output_every > 0 && step % output_every == 0)
        {
            report(step);
        }
    }
    return 0;
}

// Validates the world, then drives it to completion with the shared report line
// and this app's ASCII view. Returns 1 on an invalid configuration or a failed
// step.
template <typename World, typename Points>
int validate_and_drive(const char* title, World& world, const Points& points_of,
                       int steps, int output_every)
{
    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    // One view for the whole run, seeded from the initial state and grown as
    // needed, so successive frames share a scale and can be compared.
    Viewport view = fit_viewport(points_of(world));

    const auto report = [&world, &view, &points_of](std::int64_t step_index) {
        print_report(step_index, world, view, points_of);
    };

    return drive(title, world, report, steps, output_every);
}

// The orthographic flattening of 3D positions onto the xy plane for the debug
// view: plot (x, y), drop z.
//
// An orthographic projection along z, and nothing cleverer. malloy_ascii plots
// points that are already in the view plane, and drawing a 3D scene properly is
// a different problem rather than a wider one (ADR 0009); the 3D viewer
// (malloy_viewer, M41) draws these worlds properly, and this flattening is only
// the terminal base's view.
//
// The consequence is worth knowing while reading a frame: motion along z is
// invisible here, so an orbit tilted out of the xy plane looks like an ellipse
// rather than a circle, two bodies at different depths can overlap, and a body
// tumbling in place does not move on screen at all (the view plots CENTRES; the
// rotation is in the diagnostics columns). A helical charge orbit reads as a
// circle from this angle.
//
// One template over "has a 3D position", for the five 3D body types, which
// share that field name and nothing else.
template <typename Body>
std::vector<Vec2> projected(const std::vector<Body>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size());
    for (const Body& body : bodies)
    {
        points.push_back(Vec2{body.position.x, body.position.y});
    }
    return points;
}

// The points the ASCII view shows for a rigid body: its centre of mass and its
// body origin. Both are real state rather than an invented marker, and when the
// two are offset the origin visibly orbits the centre, which is how rotation
// becomes legible in a character grid.
std::vector<Vec2> rigid_points(const std::vector<RigidBody2D>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size() * 2);
    for (const RigidBody2D& body : bodies)
    {
        points.push_back(malloy::rigid::center_of_mass(body));
        points.push_back(body.position);
    }
    return points;
}

// Runs one particle scenario to completion. Deliberately a separate function
// from run_scenario: ParticleWorld and NBodyWorld share no base class, and
// ADR 0006 says a dispatch switch is the whole mechanism.
int run_particles(const char* title, ParticleWorld world, int steps, int output_every)
{
    const ParticleSettings& settings = world.particle_settings();
    std::cout << "\n== " << title << " ==  particles=" << world.particles().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps
              << "  restitution=" << settings.restitution;
    if (malloy::scenario::has_thermal_wall(settings))
    {
        std::cout << "  seed=" << settings.seed; // a thermal run depends on it
    }
    std::cout << '\n';

    const auto points = [](const ParticleWorld& w) { return positions_of(w.particles()); };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one 3D colliding-particle scenario. A tenth concrete runner and the
// tenth branch of the dispatch switch: ADR 0006 once more, still no base class.
int run_particles3d(const char* title, ParticleWorld3D world, int steps,
                    int output_every)
{
    const ParticleSettings3D& settings = world.particle_settings();
    std::cout << "\n== " << title << " ==  particles=" << world.particles().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps
              << "  restitution=" << settings.restitution << "  (3D)";
    if (malloy::scenario::has_thermal_wall(settings))
    {
        std::cout << "  seed=" << settings.seed; // a thermal run depends on it
    }
    std::cout << '\n';

    const auto points = [](const ParticleWorld3D& w) { return projected(w.particles()); };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one rigid-body scenario. A third concrete runner: RigidWorld shares no
// base class with the other two, and the dispatch switch is the whole
// mechanism (ADR 0006).
int run_rigid(const char* title, RigidWorld world, int steps, int output_every)
{
    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps << '\n';

    const auto points = [](const RigidWorld& w) { return rigid_points(w.bodies()); };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one 3D N-body scenario. A sixth concrete runner and a sixth branch of
// the dispatch switch, with no base class: ADR 0006 holds in three dimensions
// exactly as it does in two.
int run_nbody3d(const char* title, malloy::nbody::NBody3DWorld world, int steps,
                int output_every)
{
    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps
              << "  (3D)" << '\n';

    const auto points = [](const malloy::nbody::NBody3DWorld& w) {
        return projected(w.bodies());
    };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one 3D rigid-body scenario. A seventh concrete runner and a seventh
// branch of the dispatch switch: ADR 0006 again, and still no base class.
int run_rigid3d(const char* title, malloy::rigid::Rigid3DWorld world, int steps,
                int output_every)
{
    std::cout << "\n\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps
              << "  (3D)" << '\n';

    const auto points = [](const malloy::rigid::Rigid3DWorld& w) {
        return projected(w.bodies());
    };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one charged-particle scenario. A fifth concrete runner, and the fifth
// branch of the dispatch switch: no base class, exactly as ADR 0006 intends.
int run_charges(const char* title, ChargeWorld world, int steps, int output_every)
{
    std::cout << "\n== " << title << " ==  charges=" << world.particles().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps << '\n';

    const auto points = [](const ChargeWorld& w) { return positions_of(w.particles()); };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one 3D charged-particle scenario. An eighth concrete runner and an
// eighth branch of the dispatch switch: ADR 0006 once more, still no base class.
int run_charges3d(const char* title, Charge3DWorld world, int steps, int output_every)
{
    std::cout << "\n== " << title << " ==  charges=" << world.particles().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps
              << "  (3D)" << '\n';

    const auto points = [](const Charge3DWorld& w) { return projected(w.particles()); };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one spring scenario. A fourth concrete runner: SpringWorld shares no
// base class with the other three (ADR 0006, ADR 0008).
int run_springs(const char* title, SpringWorld world, int steps, int output_every)
{
    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  springs=" << world.network().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps << '\n';

    const auto points = [](const SpringWorld& w) { return positions_of(w.bodies()); };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one 3D spring scenario. A ninth concrete runner and a ninth branch of
// the dispatch switch: ADR 0006 once more, still no base class.
int run_springs3d(const char* title, SpringWorld3D world, int steps, int output_every)
{
    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  springs=" << world.network().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps
              << "  (3D)" << '\n';

    const auto points = [](const SpringWorld3D& w) { return projected(w.bodies()); };
    return validate_and_drive(title, world, points, steps, output_every);
}

// Runs one N-body scenario to completion, printing system diagnostics every
// output_every steps. Returns 0 on success, 1 on validation/step failure.
int run_scenario(const char* title, NBodyWorld world, int steps, int output_every)
{
    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << world.simulation_settings().dt << "  steps=" << steps << '\n';

    const auto points = [](const NBodyWorld& w) { return positions_of(w.bodies()); };
    return validate_and_drive(title, world, points, steps, output_every);
}

// The built-in scenarios, run when no scenario file is given.
int run_builtin_scenarios()
{
    const int rc = run_scenario(
        "two-body orbit (normalized units)",
        NBodyWorld{SimulationSettings{0.001}, NBodySettings{1.0, 0.000001},
                   {Body2D{Vec2{0.0, 0.0}, Vec2{0.0, 0.0}, 1.0},
                    Body2D{Vec2{1.0, 0.0}, Vec2{0.0, 1.0}, 0.000001}}},
        10000, 1000);
    if (rc != 0)
    {
        return rc;
    }

    const Real omega = 0.75984; // circular angular velocity for the unit triangle
    return run_scenario(
        "equilateral three-body (Lagrange)",
        NBodyWorld{SimulationSettings{0.001}, NBodySettings{1.0, 0.001},
                   {Body2D{Vec2{1.0, 0.0}, Vec2{0.0, omega}, 1.0},
                    Body2D{Vec2{-0.5, 0.8660254},
                           Vec2{-0.8660254 * omega, -0.5 * omega}, 1.0},
                    Body2D{Vec2{-0.5, -0.8660254},
                           Vec2{0.8660254 * omega, -0.5 * omega}, 1.0}}},
        10000, 1000);
}
} // namespace

int main(int argc, char** argv)
{
    std::cout << "MalloySim nbody_terminal\n";
    std::cout << std::fixed << std::setprecision(8);

    // With a path argument, load and run that scenario file; otherwise run the
    // built-in scenarios. A single positional path, with no flags and no CLI
    // parser.
    if (argc >= 2)
    {
        const ScenarioParseResult parsed = parse_scenario_file(argv[1]);
        if (!parsed.ok)
        {
            std::cerr << "scenario error: " << parsed.error << '\n';
            return 1;
        }
        // Refuse a file that defines nothing to simulate, rather than "run" an
        // empty default world and report success on a broken input (issue #31).
        if (malloy::scenario::entity_count(parsed.scenario) == 0)
        {
            std::cerr << "scenario error: the file defines nothing to simulate "
                         "(no bodies, particles, charges or springs)\n";
            return 1;
        }
        const auto& s = parsed.scenario;
        // The whole multi-domain mechanism: one switch, one concrete runner per
        // domain, no base class (ADR 0006). Each world is built by the run
        // layer's function for its type, the one the viewer calls too.
        namespace run = malloy::scenario;
        switch (s.type)
        {
        case ScenarioType::NBody:
            return run_scenario(argv[1], run::make_nbody_world(s), s.steps,
                                s.output_every);
        case ScenarioType::Particles:
            return run_particles(argv[1], run::make_particle_world(s), s.steps,
                                 s.output_every);
        case ScenarioType::Rigid:
            return run_rigid(argv[1], run::make_rigid_world(s), s.steps, s.output_every);
        case ScenarioType::Springs:
            return run_springs(argv[1], run::make_spring_world(s), s.steps,
                               s.output_every);
        case ScenarioType::Charges:
            return run_charges(argv[1], run::make_charge_world(s), s.steps,
                               s.output_every);
        case ScenarioType::NBody3D:
            return run_nbody3d(argv[1], run::make_nbody3d_world(s), s.steps,
                               s.output_every);
        case ScenarioType::Rigid3D:
            return run_rigid3d(argv[1], run::make_rigid3d_world(s), s.steps,
                               s.output_every);
        case ScenarioType::Charges3D:
            return run_charges3d(argv[1], run::make_charge3d_world(s), s.steps,
                                 s.output_every);
        case ScenarioType::Springs3D:
            return run_springs3d(argv[1], run::make_spring3d_world(s), s.steps,
                                 s.output_every);
        case ScenarioType::Particles3D:
            return run_particles3d(argv[1], run::make_particle3d_world(s), s.steps,
                                   s.output_every);
        }
        return 1;
    }

    return run_builtin_scenarios();
}
