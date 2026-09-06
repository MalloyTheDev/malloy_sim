#include <malloy/ascii/grid.hpp>
#include <malloy/math/math.hpp>
#include <malloy/nbody/nbody.hpp>
#include <malloy/scenario/scenario.hpp>
#include <malloy/sim_core/sim_core.hpp>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

// M8 terminal demo: runs an N-body scenario and reports the conserved system
// diagnostics (separation, total energy, total angular momentum) together with
// an ASCII view of the body positions. The scenario is either loaded from a
// file given as the single argument, or one of the built-in scenarios when no
// argument is given. The app stays dumb -- it owns no physics, no rendering
// math, and has no CLI parser (docs/03, rules 14-15).

using malloy::ascii::fit_viewport;
using malloy::ascii::render;
using malloy::ascii::Viewport;
using malloy::math::distance;
using malloy::math::Real;
using malloy::math::Vec2;
using malloy::nbody::Body2D;
using malloy::nbody::NBodySettings;
using malloy::nbody::NBodyWorld;
using malloy::nbody::total_angular_momentum;
using malloy::nbody::total_energy;
using malloy::scenario::parse_scenario_file;
using malloy::scenario::ScenarioParseResult;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;

namespace
{
// Size of the ASCII debug view. Roughly 2:1 because terminal cells are about
// twice as tall as they are wide, so a square region of world space reads as
// square on screen.
constexpr int view_width = 61;
constexpr int view_height = 25;

// The body positions, as the plain point list the ASCII view works with.
std::vector<Vec2> positions_of(const std::vector<Body2D>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size());
    for (const Body2D& body : bodies)
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
void print_view(Viewport& view, const std::vector<Body2D>& bodies)
{
    const std::vector<Vec2> points = positions_of(bodies);
    const Viewport frame = fit_viewport(points);

    view.min_x = std::min(view.min_x, frame.min_x);
    view.max_x = std::max(view.max_x, frame.max_x);
    view.min_y = std::min(view.min_y, frame.min_y);
    view.max_y = std::max(view.max_y, frame.max_y);

    std::cout << render(points, view, view_width, view_height) << "  view x ["
              << view.min_x << ", " << view.max_x << "]  y [" << view.min_y << ", "
              << view.max_y << "]\n";
}

// Runs one scenario to completion, printing system diagnostics every
// output_every steps. Returns 0 on success, 1 on validation/step failure.
int run_scenario(const char* title, const SimulationSettings& sim,
                 const NBodySettings& nbody, std::vector<Body2D> bodies, int steps,
                 int output_every)
{
    NBodyWorld world{sim, nbody, std::move(bodies)};

    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << sim.dt << "  steps=" << steps << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    // One view for the whole run, seeded from the initial state and grown as
    // needed, so successive frames share a scale and can be compared.
    Viewport view = fit_viewport(positions_of(world.bodies()));

    const auto report = [&world, &nbody, &view](std::int64_t step_index) {
        const auto& bodies_now = world.bodies();
        std::cout << "step " << std::setw(6) << step_index;
        if (bodies_now.size() >= 2)
        {
            std::cout << "   sep " << std::setw(11)
                      << distance(bodies_now[0].position, bodies_now[1].position);
        }
        std::cout << "   E_total " << std::setw(13)
                  << total_energy(bodies_now, nbody.g, nbody.softening)
                  << "   L_total " << std::setw(13) << total_angular_momentum(bodies_now)
                  << '\n';
        print_view(view, bodies_now);
    };

    report(0);
    // int64 counter: with steps == INT_MAX an int counter reaches INT_MAX,
    // passes the test, and overflows on ++, which is undefined behavior and
    // in practice loops forever.
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

// The built-in scenarios, run when no scenario file is given.
int run_builtin_scenarios()
{
    const int rc = run_scenario(
        "two-body orbit (normalized units)", SimulationSettings{0.001},
        NBodySettings{1.0, 0.000001},
        {Body2D{Vec2{0.0, 0.0}, Vec2{0.0, 0.0}, 1.0},
         Body2D{Vec2{1.0, 0.0}, Vec2{0.0, 1.0}, 0.000001}},
        10000, 1000);
    if (rc != 0)
    {
        return rc;
    }

    const Real omega = 0.75984; // circular angular velocity for the unit triangle
    return run_scenario(
        "equilateral three-body (Lagrange)", SimulationSettings{0.001},
        NBodySettings{1.0, 0.001},
        {Body2D{Vec2{1.0, 0.0}, Vec2{0.0, omega}, 1.0},
         Body2D{Vec2{-0.5, 0.8660254}, Vec2{-0.8660254 * omega, -0.5 * omega}, 1.0},
         Body2D{Vec2{-0.5, -0.8660254}, Vec2{0.8660254 * omega, -0.5 * omega}, 1.0}},
        10000, 1000);
}
} // namespace

int main(int argc, char** argv)
{
    std::cout << "MalloySim nbody_terminal\n";
    std::cout << std::fixed << std::setprecision(8);

    // With a path argument, load and run that scenario file; otherwise run the
    // built-in scenarios. A single positional path -- no flags, no CLI parser.
    if (argc >= 2)
    {
        const ScenarioParseResult parsed = parse_scenario_file(argv[1]);
        if (!parsed.ok)
        {
            std::cerr << "scenario error: " << parsed.error << '\n';
            return 1;
        }
        const auto& s = parsed.scenario;
        return run_scenario(argv[1], s.simulation, s.nbody_settings, s.bodies, s.steps,
                            s.output_every);
    }

    return run_builtin_scenarios();
}
