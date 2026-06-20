#include <malloy/math/math.hpp>
#include <malloy/nbody/nbody.hpp>
#include <malloy/sim_core/sim_core.hpp>

#include <iomanip>
#include <iostream>
#include <vector>

// M5 terminal demo: a hardcoded, normalized two-body orbit driven by
// NBodyWorld. This app is intentionally dumb -- it owns no physics, just the
// demo setup, the fixed-step loop, and formatted output (docs/03, rules 14-15).
int main()
{
    using malloy::math::distance;
    using malloy::math::Real;
    using malloy::math::Vec2;
    using malloy::nbody::Body2D;
    using malloy::nbody::NBodySettings;
    using malloy::nbody::NBodyWorld;
    using malloy::nbody::specific_orbital_energy;
    using malloy::sim_core::SimulationSettings;
    using malloy::sim_core::StepStatus;

    // M5 constants: normalized units (docs/02_MILESTONE_ROADMAP_M1_M5.md, ADR 0005).
    constexpr Real g = 1.0;
    constexpr Real softening = 0.000001;
    constexpr Real dt = 0.001;
    constexpr int steps = 10000;
    constexpr int output_every = 1000;

    std::vector<Body2D> bodies = {
        Body2D{Vec2{0.0, 0.0}, Vec2{0.0, 0.0}, 1.0},      // sun
        Body2D{Vec2{1.0, 0.0}, Vec2{0.0, 1.0}, 0.000001}, // planet
    };

    NBodyWorld world{SimulationSettings{dt}, NBodySettings{g, softening},
                     std::move(bodies)};

    std::cout << "MalloySim nbody_terminal\n";
    std::cout << "two-body orbit (normalized units): G=" << g << " dt=" << dt
              << " steps=" << steps << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << "Invalid initial configuration\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(8);

    const auto print_state = [&world, g](int step_index) {
        const Body2D& sun = world.bodies()[0];
        const Body2D& planet = world.bodies()[1];
        const Real radius = distance(planet.position, sun.position);
        const Real energy = specific_orbital_energy(sun, planet, g);
        std::cout << "step " << std::setw(6) << step_index
                  << "   radius " << std::setw(12) << radius
                  << "   specific_energy " << std::setw(14) << energy << '\n';
    };

    print_state(0);

    for (int step = 1; step <= steps; ++step)
    {
        const auto result = world.step();
        if (!result.ok())
        {
            std::cerr << "Step " << step << " failed\n";
            return 1;
        }
        if (step % output_every == 0)
        {
            print_state(step);
        }
    }

    return 0;
}
