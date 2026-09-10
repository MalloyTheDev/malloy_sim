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
        // Issue #19: this deviation is not physics, and 0.02 was 47x looser
        // than the thing it measures.
        //
        // Eliminating the velocity from the two update lines gives the
        // Stormer-Verlet position recurrence, so the STORED velocity is a
        // backward difference: it approximates the true velocity half a step
        // earlier than the time it is reported at. Supplying a tangential
        // velocity, as this world does, is therefore short of the discrete
        // circular orbit by a radial component of (dt/2)|a|, and that seeds an
        // epicycle of eccentricity dt*Omega/2.
        //
        // Here Omega is sqrt(G M / R^3) = 1.0000005, so the epicycle's
        // amplitude is R dt Omega / 2 = 5.000002e-04. The physical
        // eccentricity of these initial conditions is about 1e-6, so the wobble
        // is roughly 500 times the physics.
        //
        // The run does not reach that amplitude, and the shortfall is part of
        // the derivation rather than slack. The deviation grows as
        // R e sin(Omega t), and 1000 steps at dt = 1e-3 is Omega t = 1.0 rad,
        // about a sixth of an orbit. So the maximum over this run is
        // 5.000002e-04 * sin(1.0000005) = 4.207358e-04, against a measured
        // 4.210280e-04: agreement to 0.07 per cent, the remainder being the
        // O(dt^2) shape of the epicycle.
        //
        // Bounded from BOTH sides deliberately. Fixing the initial conditions
        // with a half-step kick would drop this into the 1e-6 range, and that
        // is a change to what an initial velocity MEANS across the whole
        // domain, so it should fail here and be updated on purpose rather than
        // pass unnoticed.
        const Real tolerance = 4.25e-4;
        const Real artefact_floor = 4.15e-4;
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
        MALLOY_CHECK_TRUE(max_deviation > artefact_floor);
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

        // These two were both asserted at 0.02, which is far looser than the
        // errors they are measuring and left them unable to fail.
        //
        // Energy is not conserved by this integrator, but it does not drift
        // secularly either: it oscillates about the shadow Hamiltonian with a
        // bounded amplitude. Measured here: 7.42e-08.
        MALLOY_CHECK_NEAR(total_energy(world.bodies(), 1.0, 0.001), energy0, 1e-6);

        // Angular momentum is EXACT for this integrator, to all orders in dt,
        // and not merely accurate. For pairwise central forces the (i,j) and
        // (j,i) torque contributions are x_i x x_j and x_j x x_i, which cancel
        // identically; softening does not break that because it changes only
        // the magnitude. So the only error possible is floating-point rounding,
        // and the measured value is 1.8e-15 on a quantity of order 2.3, about
        // 8 ulp. At the old 0.02 the integrator could have been replaced with
        // one that is not symplectic at all and this would still have passed.
        MALLOY_CHECK_NEAR(total_angular_momentum(world.bodies()), angular0, 1e-12);
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

    // --- Issue #3: a step that produces non-finite state must report failure
    //     and roll back, not return Ok. Overflow reaches here even when every
    //     input was finite and passed validate(). ---
    {
        // g * mass overflows to infinity: 1e10 * 1e300 exceeds the range of
        // double, so the acceleration is non-finite even though both bodies and
        // both settings are individually valid.
        const std::vector<Body2D> before = {Body2D{Vec2{0.0, 0.0}, Vec2{}, 1e300},
                                            Body2D{Vec2{1.0, 0.0}, Vec2{}, 1e300}};
        NBodyWorld w{SimulationSettings{1e10}, NBodySettings{1e10, 1.0}, before};
        MALLOY_CHECK_TRUE(w.validate() == StepStatus::Ok); // the inputs are fine

        const auto result = w.step();
        MALLOY_CHECK_TRUE(result.status == StepStatus::InvalidState);
        MALLOY_CHECK_FALSE(result.ok());

        // Rolled back: the documented "state is left unchanged on failure"
        // contract has to hold for a failure detected after the update too.
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
        for (std::size_t i = 0; i < before.size(); ++i)
        {
            MALLOY_CHECK_VEC2_NEAR(w.bodies()[i].position, before[i].position, 0.0);
            MALLOY_CHECK_VEC2_NEAR(w.bodies()[i].velocity, before[i].velocity, 0.0);
        }
    }

    // --- Issue #3: a normal step is unaffected by the post-check. ---
    {
        NBodyWorld w{SimulationSettings{0.5}, NBodySettings{1.0, 0.0}, two_unit_bodies()};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{1});
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].position, Vec2(0.25, 0.0), eps);
    }

    // --- Elapsed simulation time comes from malloy_time, so it is ticks * dt
    //     rather than an accumulated sum and cannot drift. dt = 0.1 is chosen
    //     because accumulating it ten times does NOT give exactly 1.0. ---
    {
        NBodyWorld w{SimulationSettings{0.1}, NBodySettings{1.0, 1e-6},
                     two_unit_bodies()};
        MALLOY_CHECK_NEAR(w.elapsed_time(), 0.0, 0.0); // exact before any step
        for (int i = 0; i < 10; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{10});
        MALLOY_CHECK_NEAR(w.elapsed_time(), 1.0, 1e-17);
    }
    {
        // An unusable dt leaves both at zero rather than reporting a bogus time.
        NBodyWorld w{SimulationSettings{0.0}, NBodySettings{1.0, 1e-6},
                     two_unit_bodies()};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
        MALLOY_CHECK_NEAR(w.elapsed_time(), 0.0, 0.0);
    }

    // --- The settings accessors are public API and were untested. Trivial, but
    //     a getter returning the wrong member is exactly the kind of thing only
    //     a test catches. ---
    {
        const SimulationSettings sim{0.0025};
        const NBodySettings nb{2.5, 0.125};
        NBodyWorld w{sim, nb, two_unit_bodies()};
        MALLOY_CHECK_NEAR(w.simulation_settings().dt, 0.0025, eps);
        MALLOY_CHECK_NEAR(w.nbody_settings().g, 2.5, eps);
        MALLOY_CHECK_NEAR(w.nbody_settings().softening, 0.125, eps);
    }

    // --- The FORCE side of the softening contract is pinned above. The
    //     ENERGY side was not: every existing call uses either softening 0 or a
    //     softening so much smaller than the separation that dropping it shifts
    //     the result by about 3e-7, under tolerances of 0.02.
    //
    //     diagnostics.hpp states the point of this function: it uses the same
    //     softening as the integrator so that total_energy is a quantity the
    //     dynamics conserve. Here softening is comparable to the separation, so
    //     dropping it is unmissable: r = sqrt(9 + 16) = 5 exactly, not 3. ---
    {
        const std::vector<Body2D> b = {Body2D{Vec2{1.0, -1.0}, Vec2{}, 2.0},
                                       Body2D{Vec2{4.0, -1.0}, Vec2{}, 3.0}};
        // U = -g m1 m2 / r = -(1.5)(2)(3)/5
        MALLOY_CHECK_NEAR(total_potential_energy(b, 1.5, 4.0), -1.8, eps);
        // Unsoftened the same pair would give -(1.5)(2)(3)/3 = -3.
        MALLOY_CHECK_NEAR(total_potential_energy(b, 1.5, 0.0), -3.0, eps);
    }

    // --- Issue #17: a softening large enough that its SQUARE overflows is
    //     refused, in both domains that have one.
    //
    //     It has to be refused rather than tolerated, because the failure is
    //     invisible in every diagnostic. A squared softening of infinity makes
    //     every separation infinite, so every force is exactly zero and every
    //     potential is exactly zero, and a simulation in which nothing happens
    //     conserves everything perfectly. Measured before the fix: two like
    //     charges that should fly apart sat still for 2000 steps while the
    //     energy column read 0.00000000e+00, when the correct value was 0.5. ---
    {
        NBodySettings ok;
        ok.softening = 1.0e154; // squares to 1e308, still finite
        MALLOY_CHECK_TRUE(ok.is_valid());

        NBodySettings over;
        over.softening = 1.4e154; // squares to just over DBL_MAX
        MALLOY_CHECK_FALSE(over.is_valid());

        NBodySettings absurd;
        absurd.softening = 1.0e200;
        MALLOY_CHECK_FALSE(absurd.is_valid());

        // Still refused for the older reasons.
        NBodySettings negative;
        negative.softening = -1.0;
        MALLOY_CHECK_FALSE(negative.is_valid());
        NBodySettings not_finite;
        not_finite.softening = inf;
        MALLOY_CHECK_FALSE(not_finite.is_valid());
    }

    // --- Issue #16: a position or velocity whose SQUARE overflows is refused.
    //     Finite is not enough. Everything that squares a vector reaches
    //     infinity above about 1.34e154, so accepting state up to 1.8e308 left
    //     a window in which a world validated clean while every energy it
    //     reported was inf. ---
    {
        const Real too_big = 1.4e154; // finite, and its square is not
        MALLOY_CHECK_TRUE((Body2D{Vec2{1.0, 2.0}, Vec2{}, 1.0}.is_valid()));
        MALLOY_CHECK_FALSE((Body2D{Vec2{too_big, 0.0}, Vec2{}, 1.0}.is_valid()));
        MALLOY_CHECK_FALSE((Body2D{Vec2{}, Vec2{0.0, too_big}, 1.0}.is_valid()));
    }

    std::cout << "malloy_nbody_tests passed\n";
    return 0;
}
