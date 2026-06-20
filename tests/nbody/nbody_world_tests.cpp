#include <malloy/nbody/nbody.hpp>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <test_check.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

using malloy::math::distance;
using malloy::math::is_finite;
using malloy::math::Real;
using malloy::math::Vec2;
using malloy::nbody::Body2D;
using malloy::nbody::NBodySettings;
using malloy::nbody::NBodyWorld;
using malloy::nbody::center_of_mass;
using malloy::nbody::total_angular_momentum;
using malloy::nbody::total_energy;
using malloy::nbody::total_kinetic_energy;
using malloy::nbody::total_momentum;
using malloy::nbody::total_potential_energy;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;

namespace
{
// A valid two-body baseline: a normalized "sun" and a light "planet" placed on
// a near-circular orbit (v = sqrt(G*M/r) = 1 at r = 1). Matches the M5 demo.
NBodyWorld make_orbit_world()
{
    const SimulationSettings sim{0.001};
    const NBodySettings nbody{1.0, 1e-6};
    std::vector<Body2D> bodies = {
        Body2D{Vec2{0.0, 0.0}, Vec2{0.0, 0.0}, 1.0},  // sun
        Body2D{Vec2{1.0, 0.0}, Vec2{0.0, 1.0}, 1e-6}, // planet
    };
    return NBodyWorld{sim, nbody, std::move(bodies)};
}

std::vector<Body2D> two_unit_bodies()
{
    return {Body2D{Vec2{0.0, 0.0}, Vec2{}, 1.0},
            Body2D{Vec2{1.0, 0.0}, Vec2{}, 1.0}};
}
} // namespace

int main()
{
    const Real eps = 1e-12;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    // --- Settings validation: dt <= 0 is InvalidSettings, no state change. ---
    {
        NBodyWorld zero_dt{SimulationSettings{0.0}, NBodySettings{1.0, 1e-6}, two_unit_bodies()};
        MALLOY_CHECK_TRUE(zero_dt.step().status == StepStatus::InvalidSettings);
        MALLOY_CHECK_EQ(zero_dt.tick_count(), std::uint64_t{0});

        NBodyWorld neg_dt{SimulationSettings{-0.001}, NBodySettings{1.0, 1e-6}, two_unit_bodies()};
        MALLOY_CHECK_TRUE(neg_dt.step().status == StepStatus::InvalidSettings);
    }

    // --- G < 0 is InvalidSettings. ---
    {
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{-1.0, 1e-6}, two_unit_bodies()};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
    }

    // --- softening < 0 is InvalidSettings. ---
    {
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, -1e-6}, two_unit_bodies()};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
    }

    // --- mass <= 0 is InvalidState. ---
    {
        std::vector<Body2D> zero_mass = {Body2D{Vec2{}, Vec2{}, 0.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, 1e-6}, zero_mass};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);

        std::vector<Body2D> neg_mass = {Body2D{Vec2{}, Vec2{}, -1.0}};
        NBodyWorld w2{SimulationSettings{0.001}, NBodySettings{1.0, 1e-6}, neg_mass};
        MALLOY_CHECK_TRUE(w2.step().status == StepStatus::InvalidState);
    }

    // --- non-finite position is InvalidState. ---
    {
        std::vector<Body2D> b = {Body2D{Vec2{inf, 0.0}, Vec2{}, 1.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, 1e-6}, b};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
    }

    // --- non-finite velocity is InvalidState. ---
    {
        std::vector<Body2D> b = {Body2D{Vec2{}, Vec2{nan, 0.0}, 1.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, 1e-6}, b};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
    }

    // --- A single body feels zero acceleration. ---
    {
        std::vector<Body2D> b = {Body2D{Vec2{3.0, -2.0}, Vec2{}, 5.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, 1e-6}, b};
        const auto acc = w.compute_accelerations();
        MALLOY_CHECK_EQ(acc.size(), std::size_t{1});
        MALLOY_CHECK_VEC2_NEAR(acc[0], Vec2(0.0, 0.0), eps);
    }

    // --- Coincident bodies do not produce NaN when softening > 0. ---
    {
        std::vector<Body2D> b = {Body2D{Vec2{0.0, 0.0}, Vec2{}, 1.0},
                                 Body2D{Vec2{0.0, 0.0}, Vec2{}, 1.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, 1e-3}, b};
        const auto acc = w.compute_accelerations();
        MALLOY_CHECK_TRUE(is_finite(acc[0]));
        MALLOY_CHECK_TRUE(is_finite(acc[1]));
    }

    // --- Pairwise force symmetry: m_a * a_a == -(m_b * a_b). ---
    {
        std::vector<Body2D> b = {Body2D{Vec2{0.0, 0.0}, Vec2{}, 2.0},
                                 Body2D{Vec2{1.5, 0.5}, Vec2{}, 3.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, 1e-6}, b};
        const auto acc = w.compute_accelerations();
        const Vec2 force_a = 2.0 * acc[0];
        const Vec2 force_b = 3.0 * acc[1];
        MALLOY_CHECK_VEC2_NEAR(force_a, -force_b, eps);
    }

    // --- Determinism: identical worlds stepped identically end identical. ---
    {
        NBodyWorld a = make_orbit_world();
        NBodyWorld b = make_orbit_world();
        for (int i = 0; i < 100; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        MALLOY_CHECK_EQ(a.tick_count(), std::uint64_t{100});

        const auto& bodies_a = a.bodies();
        const auto& bodies_b = b.bodies();
        MALLOY_CHECK_EQ(bodies_a.size(), bodies_b.size());
        for (std::size_t i = 0; i < bodies_a.size(); ++i)
        {
            MALLOY_CHECK_VEC2_NEAR(bodies_a[i].position, bodies_b[i].position, eps);
            MALLOY_CHECK_VEC2_NEAR(bodies_a[i].velocity, bodies_b[i].velocity, eps);
        }
    }

    // --- Near-circular orbit: the radius stays close to 1 over a short run. ---
    {
        NBodyWorld w = make_orbit_world();
        const Real tolerance = 0.02;
        Real max_deviation = 0.0;
        for (int i = 0; i < 1000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
            const auto& bodies = w.bodies();
            const Real radius = distance(bodies[1].position, bodies[0].position);
            const Real deviation = std::abs(radius - 1.0);
            if (deviation > max_deviation)
            {
                max_deviation = deviation;
            }
        }
        MALLOY_CHECK_TRUE(max_deviation < tolerance);
    }

    // --- Specific orbital energy: known value, then approximate conservation. ---
    {
        const Body2D sun{Vec2{0.0, 0.0}, Vec2{0.0, 0.0}, 1.0};
        const Body2D planet{Vec2{1.0, 0.0}, Vec2{0.0, 1.0}, 1e-6};
        MALLOY_CHECK_NEAR(malloy::nbody::specific_orbital_energy(sun, planet, 1.0),
                          -0.500001, 1e-9);

        NBodyWorld w = make_orbit_world();
        const Real initial_energy =
            malloy::nbody::specific_orbital_energy(w.bodies()[0], w.bodies()[1], 1.0);
        for (int i = 0; i < 1000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const Real final_energy =
            malloy::nbody::specific_orbital_energy(w.bodies()[0], w.bodies()[1], 1.0);
        MALLOY_CHECK_NEAR(final_energy, initial_energy, 0.01);
    }

    // --- System diagnostics: known values for a hand-computed configuration. ---
    {
        const std::vector<Body2D> b = {Body2D{Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, 2.0},
                                       Body2D{Vec2{4.0, 0.0}, Vec2{0.0, 2.0}, 3.0}};
        MALLOY_CHECK_NEAR(total_kinetic_energy(b), 7.0, 1e-9);
        MALLOY_CHECK_NEAR(total_potential_energy(b, 1.0, 0.0), -1.5, 1e-9);
        MALLOY_CHECK_NEAR(total_energy(b, 1.0, 0.0), 5.5, 1e-9);
        MALLOY_CHECK_VEC2_NEAR(total_momentum(b), Vec2(2.0, 6.0), 1e-9);
        MALLOY_CHECK_VEC2_NEAR(center_of_mass(b), Vec2(2.4, 0.0), 1e-9);
        MALLOY_CHECK_NEAR(total_angular_momentum(b), 24.0, 1e-9);
    }

    // --- N>2 conservation: a rotating equilateral triangle (Lagrange config)
    //     keeps its conserved quantities bounded across a run; linear momentum
    //     is held to machine precision. ---
    {
        const Real omega = 0.75984;
        const std::vector<Body2D> triangle = {
            Body2D{Vec2{1.0, 0.0}, Vec2{0.0, omega}, 1.0},
            Body2D{Vec2{-0.5, 0.8660254}, Vec2{-0.8660254 * omega, -0.5 * omega}, 1.0},
            Body2D{Vec2{-0.5, -0.8660254}, Vec2{0.8660254 * omega, -0.5 * omega}, 1.0}};
        NBodyWorld world{SimulationSettings{0.001}, NBodySettings{1.0, 0.001}, triangle};

        const Real energy0 = total_energy(world.bodies(), 1.0, 0.001);
        const Vec2 momentum0 = total_momentum(world.bodies());
        const Real angular0 = total_angular_momentum(world.bodies());

        for (int i = 0; i < 1000; ++i)
        {
            MALLOY_CHECK_TRUE(world.step().ok());
        }

        MALLOY_CHECK_VEC2_NEAR(total_momentum(world.bodies()), momentum0, 1e-9);
        MALLOY_CHECK_NEAR(total_energy(world.bodies(), 1.0, 0.001), energy0, 0.02);
        MALLOY_CHECK_NEAR(total_angular_momentum(world.bodies()), angular0, 0.02);
    }

    std::cout << "malloy_nbody_tests passed\n";
    return 0;
}
