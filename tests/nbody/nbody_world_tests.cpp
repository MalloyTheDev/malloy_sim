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

    // --- Issue #2: coincident bodies with softening == 0 stay finite. The
    //     inverse-cube law is undefined there, so the pair contributes nothing
    //     rather than producing NaN. ---
    {
        std::vector<Body2D> b = {Body2D{Vec2{0.0, 0.0}, Vec2{}, 1.0},
                                 Body2D{Vec2{0.0, 0.0}, Vec2{}, 1.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, 0.0}, b};
        const auto acc = w.compute_accelerations();
        MALLOY_CHECK_VEC2_NEAR(acc[0], Vec2(0.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(acc[1], Vec2(0.0, 0.0), eps);

        // The energy diagnostics must be finite for the same configuration, and
        // must use the same guard so force and energy stay consistent.
        MALLOY_CHECK_NEAR(total_potential_energy(w.bodies(), 1.0, 0.0), 0.0, eps);
        MALLOY_CHECK_TRUE(is_finite(total_energy(w.bodies(), 1.0, 0.0)));

        // A full step must leave every body finite and report success.
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_TRUE(is_finite(w.bodies()[0].position));
        MALLOY_CHECK_TRUE(is_finite(w.bodies()[0].velocity));
        MALLOY_CHECK_TRUE(is_finite(w.bodies()[1].position));
        MALLOY_CHECK_TRUE(is_finite(w.bodies()[1].velocity));
    }

    // --- Issue #2: g == 0 is legal and means no gravity, so coincident bodies
    //     must give exactly zero acceleration rather than 0 * inf. ---
    {
        std::vector<Body2D> b = {Body2D{Vec2{0.0, 0.0}, Vec2{}, 1.0},
                                 Body2D{Vec2{0.0, 0.0}, Vec2{}, 1.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{0.0, 0.0}, b};
        const auto acc = w.compute_accelerations();
        MALLOY_CHECK_VEC2_NEAR(acc[0], Vec2(0.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(acc[1], Vec2(0.0, 0.0), eps);
        // specific_orbital_energy has no potential term when g == 0, so it must
        // not evaluate 0/0 for a coincident pair either.
        const Real e = malloy::nbody::specific_orbital_energy(b[0], b[1], 0.0);
        MALLOY_CHECK_NEAR(e, 0.0, eps);
    }

    // --- Issue #8: the integrator is semi-implicit (symplectic) Euler, not
    //     explicit Euler. One hand-computable step pins the order exactly:
    //     two unit masses 1 apart with G=1 and no softening give a = (+-1, 0),
    //     so with dt = 0.5 the velocities become (+-0.5, 0) and the positions
    //     move by velocity*dt AFTER that update. Explicit Euler would use the
    //     old zero velocities and leave both positions unchanged. ---
    {
        NBodyWorld w{SimulationSettings{0.5}, NBodySettings{1.0, 0.0}, two_unit_bodies()};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(0.5, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(-0.5, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].position, Vec2(0.25, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].position, Vec2(0.75, 0.0), eps);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{1});
    }

    // --- Issue #8: pin the acceleration MAGNITUDE, not just its symmetry.
    //     Force symmetry alone survives a dropped G, a wrong inverse power, and
    //     even a sign flip, so assert exact known values instead. ---
    {
        // m=2 at the origin, m=3 at (1,0), G=2, no softening.
        // a0 = G*m1/r^2 = 2*3 = 6 toward +x; a1 = G*m0/r^2 = 2*2 = 4 toward -x.
        std::vector<Body2D> b = {Body2D{Vec2{0.0, 0.0}, Vec2{}, 2.0},
                                 Body2D{Vec2{1.0, 0.0}, Vec2{}, 3.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{2.0, 0.0}, b};
        const auto acc = w.compute_accelerations();
        MALLOY_CHECK_VEC2_NEAR(acc[0], Vec2(6.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(acc[1], Vec2(-4.0, 0.0), eps);
    }

    // --- Issue #8: softening enters the denominator SQUARED (docs/04). With
    //     r = 1 and softening = 2, r2 = 1 + 4 = 5, so |a| = 5^-1.5. If softening
    //     were added unsquared, r2 would be 3 and |a| would be 3^-1.5. ---
    {
        std::vector<Body2D> b = {Body2D{Vec2{0.0, 0.0}, Vec2{}, 1.0},
                                 Body2D{Vec2{1.0, 0.0}, Vec2{}, 1.0}};
        NBodyWorld w{SimulationSettings{0.001}, NBodySettings{1.0, 2.0}, b};
        const auto acc = w.compute_accelerations();
        MALLOY_CHECK_VEC2_NEAR(acc[0], Vec2(0.0894427190999916, 0.0), 1e-15);
        MALLOY_CHECK_VEC2_NEAR(acc[1], Vec2(-0.0894427190999916, 0.0), 1e-15);
    }

    // --- Issue #8: exercise the second half of the angular-momentum formula.
    //     Every existing case has position.y * velocity.x == 0, so the
    //     "- y*vx" term never contributes and could be deleted undetected.
    //     Bodies are built as named locals because the check macros are
    //     function-like and would split a braced initializer on its commas. ---
    {
        // Single body at (0,1) moving along +x: L = m(x*vy - y*vx) = -2.
        const std::vector<Body2D> spin_y = {Body2D{Vec2{0.0, 1.0}, Vec2{2.0, 0.0}, 1.0}};
        MALLOY_CHECK_NEAR(total_angular_momentum(spin_y), -2.0, eps);

        // Both terms nonzero and of opposite sign: L = 1*(2*3 - 4*5) = -14.
        const std::vector<Body2D> spin_both = {Body2D{Vec2{2.0, 4.0}, Vec2{5.0, 3.0}, 1.0}};
        MALLOY_CHECK_NEAR(total_angular_momentum(spin_both), -14.0, eps);
    }

    // --- Issue #8: the is_finite clauses in validation were never exercised.
    //     inf >= 0 is true, so without them an infinite G or mass passes
    //     validation and silently turns the whole world into NaN. ---
    {
        const Body2D good{Vec2{1.0, 2.0}, Vec2{3.0, 4.0}, 1.0};
        const Body2D inf_mass{Vec2{}, Vec2{}, inf};
        const Body2D neg_inf_mass{Vec2{}, Vec2{}, -inf};
        const Body2D nan_mass{Vec2{}, Vec2{}, nan};
        const Body2D inf_position{Vec2{inf, 0.0}, Vec2{}, 1.0};
        const Body2D nan_velocity{Vec2{}, Vec2{0.0, nan}, 1.0};

        MALLOY_CHECK_TRUE(good.is_valid());
        MALLOY_CHECK_FALSE(inf_mass.is_valid());
        MALLOY_CHECK_FALSE(neg_inf_mass.is_valid());
        MALLOY_CHECK_FALSE(nan_mass.is_valid());
        MALLOY_CHECK_FALSE(inf_position.is_valid());
        MALLOY_CHECK_FALSE(nan_velocity.is_valid());

        const NBodySettings no_gravity{0.0, 0.0}; // G = 0 is legal
        const NBodySettings inf_g{inf, 0.0};
        const NBodySettings nan_g{nan, 0.0};
        const NBodySettings inf_softening{1.0, inf};
        const NBodySettings nan_softening{1.0, nan};

        MALLOY_CHECK_TRUE(no_gravity.is_valid());
        MALLOY_CHECK_FALSE(inf_g.is_valid());
        MALLOY_CHECK_FALSE(nan_g.is_valid());
        MALLOY_CHECK_FALSE(inf_softening.is_valid());
        MALLOY_CHECK_FALSE(nan_softening.is_valid());

        // And end to end, so the rule is enforced where it matters.
        NBodyWorld w{SimulationSettings{0.001}, inf_g, two_unit_bodies()};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
    }

    std::cout << "malloy_nbody_tests passed\n";
    return 0;
}
