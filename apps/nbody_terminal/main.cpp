#include <malloy/math/math.hpp>
#include <malloy/nbody/nbody.hpp>
#include <malloy/sim_core/sim_core.hpp>

#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

// M6 terminal demo: runs hardcoded N-body scenarios and reports the conserved
// system diagnostics (separation, total energy, total angular momentum) every
// few thousand steps. The app stays dumb -- it owns no physics and has no CLI;
// it just sets up bodies and drives NBodyWorld (docs/03, rules 14-15).

using malloy::math::distance;
using malloy::math::Real;
using malloy::math::Vec2;
using malloy::nbody::Body2D;
using malloy::nbody::NBodySettings;
using malloy::nbody::NBodyWorld;
using malloy::nbody::total_angular_momentum;
using malloy::nbody::total_energy;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;

namespace
{
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
        std::cerr << title << ": invalid initial configuration\n";
        return 1;
    }

    const auto report = [&world, &nbody](int step_index) {
        const auto& bodies_now = world.bodies();
        const Real separation = distance(bodies_now[0].position, bodies_now[1].position);
        const Real energy = total_energy(bodies_now, nbody.g, nbody.softening);
        const Real angular = total_angular_momentum(bodies_now);
        std::cout << "step " << std::setw(6) << step_index
                  << "   sep " << std::setw(11) << separation
                  << "   E_total " << std::setw(13) << energy
                  << "   L_total " << std::setw(13) << angular << '\n';
    };

    report(0);
    for (int step = 1; step <= steps; ++step)
    {
        if (!world.step().ok())
        {
            std::cerr << title << ": step " << step << " failed\n";
            return 1;
        }
        if (step % output_every == 0)
        {
            report(step);
        }
    }
    return 0;
}
} // namespace

int main()
{
    std::cout << "MalloySim nbody_terminal\n";
    std::cout << std::fixed << std::setprecision(8);

    // Scenario 1: normalized sun + planet, a near-circular orbit (the M5 demo).
    int rc = run_scenario(
        "two-body orbit (normalized units)", SimulationSettings{0.001},
        NBodySettings{1.0, 0.000001},
        {Body2D{Vec2{0.0, 0.0}, Vec2{0.0, 0.0}, 1.0},
         Body2D{Vec2{1.0, 0.0}, Vec2{0.0, 1.0}, 0.000001}},
        10000, 1000);
    if (rc != 0)
    {
        return rc;
    }

    // Scenario 2: three equal masses in the rotating equilateral-triangle
    // (Lagrange) configuration. Exercises N>2 and shows total energy and
    // angular momentum staying put with comparable masses.
    const Real omega = 0.75984; // circular angular velocity for the unit triangle
    rc = run_scenario(
        "equilateral three-body (Lagrange)", SimulationSettings{0.001},
        NBodySettings{1.0, 0.001},
        {Body2D{Vec2{1.0, 0.0}, Vec2{0.0, omega}, 1.0},
         Body2D{Vec2{-0.5, 0.8660254}, Vec2{-0.8660254 * omega, -0.5 * omega}, 1.0},
         Body2D{Vec2{-0.5, -0.8660254}, Vec2{0.8660254 * omega, -0.5 * omega}, 1.0}},
        10000, 1000);
    return rc;
}
