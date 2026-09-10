#include <malloy/springs/springs.hpp>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <test_check.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using malloy::math::Real;
using malloy::math::Vec2;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;
using malloy::springs::accumulate_spring_forces;
using malloy::springs::Spring;
using malloy::springs::SpringBody2D;
using malloy::springs::SpringNetwork;
using malloy::springs::SpringWorld;
using malloy::springs::total_elastic_energy;
using malloy::springs::total_kinetic_energy;
using malloy::springs::total_momentum;

namespace
{
Spring make_spring(std::size_t a, std::size_t b, Real rest, Real k, Real c)
{
    Spring s;
    s.a = a;
    s.b = b;
    s.rest_length = rest;
    s.stiffness = k;
    s.damping = c;
    return s;
}

SpringBody2D body_at(Vec2 position, Vec2 velocity, Real mass)
{
    SpringBody2D b;
    b.position = position;
    b.velocity = velocity;
    b.mass = mass;
    return b;
}
} // namespace

int main()
{
    const Real eps = 1e-12;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    // --- Spring validation, including the self-spring case. ---
    {
        MALLOY_CHECK_TRUE(make_spring(0, 1, 1.0, 10.0, 0.5).is_valid());
        MALLOY_CHECK_TRUE(make_spring(0, 1, 0.0, 0.0, 0.0).is_valid()); // all zero is legal
        // A body cannot be sprung to itself: the axis would be undefined forever.
        MALLOY_CHECK_FALSE(make_spring(2, 2, 1.0, 10.0, 0.0).is_valid());
        MALLOY_CHECK_FALSE(make_spring(0, 1, -1.0, 10.0, 0.0).is_valid());
        MALLOY_CHECK_FALSE(make_spring(0, 1, 1.0, -10.0, 0.0).is_valid());
        MALLOY_CHECK_FALSE(make_spring(0, 1, 1.0, 10.0, -0.5).is_valid());
        MALLOY_CHECK_FALSE(make_spring(0, 1, nan, 10.0, 0.0).is_valid());
        MALLOY_CHECK_FALSE(make_spring(0, 1, 1.0, inf, 0.0).is_valid());
    }

    // --- Body validation. ---
    {
        MALLOY_CHECK_TRUE(body_at(Vec2{1.0, 2.0}, Vec2{3.0, 4.0}, 2.5).is_valid());
        MALLOY_CHECK_FALSE(body_at(Vec2{}, Vec2{}, 0.0).is_valid());
        MALLOY_CHECK_FALSE(body_at(Vec2{}, Vec2{}, -1.0).is_valid());
        MALLOY_CHECK_FALSE(body_at(Vec2{}, Vec2{}, inf).is_valid());
        MALLOY_CHECK_FALSE(body_at(Vec2{nan, 0.0}, Vec2{}, 1.0).is_valid());
    }

    // --- Network add/remove is deterministic, and remove preserves the
    //     relative order of everything else. ---
    {
        SpringNetwork n;
        MALLOY_CHECK_EQ(n.add(make_spring(0, 1, 1.0, 1.0, 0.0)), std::size_t{0});
        MALLOY_CHECK_EQ(n.add(make_spring(1, 2, 2.0, 2.0, 0.0)), std::size_t{1});
        MALLOY_CHECK_EQ(n.add(make_spring(2, 3, 3.0, 3.0, 0.0)), std::size_t{2});
        MALLOY_CHECK_EQ(n.size(), std::size_t{3});

        MALLOY_CHECK_TRUE(n.remove(1));
        MALLOY_CHECK_EQ(n.size(), std::size_t{2});
        // The survivors keep their relative order: rest lengths 1 then 3.
        MALLOY_CHECK_NEAR(n.springs()[0].rest_length, 1.0, eps);
        MALLOY_CHECK_NEAR(n.springs()[1].rest_length, 3.0, eps);

        MALLOY_CHECK_FALSE(n.remove(99)); // out of range changes nothing
        MALLOY_CHECK_EQ(n.size(), std::size_t{2});

        // Endpoint ids are checked against the world's body count.
        MALLOY_CHECK_TRUE(n.is_valid_for(4));
        MALLOY_CHECK_FALSE(n.is_valid_for(3)); // body 3 no longer exists
    }

    // --- Hooke's law, with exact hand-computed values. Bodies are 3 apart on a
    //     spring with rest length 2, so the extension is 1 and k = 10 gives a
    //     force of 10 pulling them together. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{3.0, 0.0}, Vec2{}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 2.0, 10.0, 0.0));

        std::vector<Vec2> forces(2);
        accumulate_spring_forces(n, bodies, forces);
        // Stretched, so body 0 is pulled toward +x and body 1 toward -x.
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(10.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(-10.0, 0.0), eps);
        // Equal and opposite, exactly.
        MALLOY_CHECK_VEC2_NEAR(forces[0] + forces[1], Vec2(0.0, 0.0), 0.0);
    }

    // --- Compression pushes apart: the opposite direction from extension. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{1.0, 0.0}, Vec2{}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 3.0, 10.0, 0.0)); // rest 3, actual 1: compressed by 2

        std::vector<Vec2> forces(2);
        accumulate_spring_forces(n, bodies, forces);
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(-20.0, 0.0), eps); // pushed away from b
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(20.0, 0.0), eps);
    }

    // --- A spring at exactly its rest length with no relative motion produces
    //     no force at all. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{0.0, 0.0}, Vec2{1.0, 1.0}, 1.0),
            body_at(Vec2{2.0, 0.0}, Vec2{1.0, 1.0}, 1.0)}; // same velocity
        SpringNetwork n;
        n.add(make_spring(0, 1, 2.0, 100.0, 5.0));

        std::vector<Vec2> forces(2);
        accumulate_spring_forces(n, bodies, forces);
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(0.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(0.0, 0.0), eps);
    }

    // --- Zero stiffness produces no elastic force however stretched. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{50.0, 0.0}, Vec2{}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 0.0, 0.0));

        std::vector<Vec2> forces(2);
        accumulate_spring_forces(n, bodies, forces);
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(0.0, 0.0), eps);
    }

    // --- Damping uses the AXIAL relative velocity, not the full relative
    //     speed. The perpendicular component must contribute nothing. ---
    {
        // At rest length, so the elastic term is zero and only damping acts.
        // Relative velocity is (3, 4): axial component along +x is 3.
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{0.0, 0.0}, Vec2{0.0, 0.0}, 1.0),
            body_at(Vec2{2.0, 0.0}, Vec2{3.0, 4.0}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 2.0, 100.0, 2.0));

        std::vector<Vec2> forces(2);
        accumulate_spring_forces(n, bodies, forces);
        // c * v_axial = 2 * 3 = 6, along the axis only. Using |v_rel| = 5 would
        // give 10, and a perpendicular leak would give a nonzero y.
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(6.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(-6.0, 0.0), eps);
    }

    // --- Damping opposes relative axial motion, in both directions. ---
    {
        SpringNetwork n;
        n.add(make_spring(0, 1, 2.0, 0.0, 3.0)); // damper only

        // Separating: the force must pull them back together.
        const std::vector<SpringBody2D> apart = {
            body_at(Vec2{0.0, 0.0}, Vec2{-1.0, 0.0}, 1.0),
            body_at(Vec2{2.0, 0.0}, Vec2{1.0, 0.0}, 1.0)};
        std::vector<Vec2> f_apart(2);
        accumulate_spring_forces(n, apart, f_apart);
        MALLOY_CHECK_TRUE(f_apart[0].x > 0.0); // body 0 pulled toward +x
        MALLOY_CHECK_TRUE(f_apart[1].x < 0.0);

        // Approaching: the force must push them apart.
        const std::vector<SpringBody2D> closing = {
            body_at(Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, 1.0),
            body_at(Vec2{2.0, 0.0}, Vec2{-1.0, 0.0}, 1.0)};
        std::vector<Vec2> f_close(2);
        accumulate_spring_forces(n, closing, f_close);
        MALLOY_CHECK_TRUE(f_close[0].x < 0.0);
        MALLOY_CHECK_TRUE(f_close[1].x > 0.0);
    }

    // --- MANY-TO-ONE: several springs incident on one body must ACCUMULATE.
    //     An implementation that assigns instead of adding, or that stops after
    //     the first incident spring, gives a different answer here. ---
    {
        // Body 1 in the middle, pulled by two stretched springs from opposite
        // sides with DIFFERENT stiffnesses, so the two contributions cannot
        // cancel and neither can be mistaken for the total.
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{-3.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{4.0, 0.0}, Vec2{}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 2.0, 0.0)); // length 3, ext 2, k 2 -> 4
        n.add(make_spring(1, 2, 1.0, 5.0, 0.0)); // length 4, ext 3, k 5 -> 15

        std::vector<Vec2> forces(3);
        accumulate_spring_forces(n, bodies, forces);
        // Spring 0 pulls body 1 toward body 0, i.e. -x by 4.
        // Spring 1 pulls body 1 toward body 2, i.e. +x by 15.
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(11.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(4.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(forces[2], Vec2(-15.0, 0.0), eps);
        // Every spring being equal and opposite means the whole network sums to
        // zero, whatever its topology.
        MALLOY_CHECK_VEC2_NEAR(forces[0] + forces[1] + forces[2], Vec2(0.0, 0.0), 1e-14);
    }

    // --- accumulate_spring_forces ADDS into the buffer rather than assigning,
    //     which is what makes the many-to-one case work at all. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{3.0, 0.0}, Vec2{}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 2.0, 10.0, 0.0));

        std::vector<Vec2> forces = {Vec2{100.0, 7.0}, Vec2{-100.0, -7.0}};
        accumulate_spring_forces(n, bodies, forces);
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(110.0, 7.0), eps);
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(-110.0, -7.0), eps);
    }

    // --- Removing one spring changes only that spring's contribution. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{-3.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{4.0, 0.0}, Vec2{}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 2.0, 0.0));
        n.add(make_spring(1, 2, 1.0, 5.0, 0.0));

        std::vector<Vec2> both(3);
        accumulate_spring_forces(n, bodies, both);

        MALLOY_CHECK_TRUE(n.remove(0)); // drop the first spring
        std::vector<Vec2> one(3);
        accumulate_spring_forces(n, bodies, one);

        // Body 2 only ever touched the surviving spring, so it is unchanged.
        MALLOY_CHECK_VEC2_NEAR(one[2], both[2], 0.0);
        // Body 0 only touched the removed spring, so it now gets nothing.
        MALLOY_CHECK_VEC2_NEAR(one[0], Vec2(0.0, 0.0), eps);
        // Body 1 loses exactly the removed spring's contribution.
        MALLOY_CHECK_VEC2_NEAR(one[1], both[1] - Vec2(-4.0, 0.0), eps);
    }

    // --- Coincident endpoints have no axis, so they contribute nothing rather
    //     than producing NaN. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{1.0, 1.0}, Vec2{}, 1.0),
            body_at(Vec2{1.0, 1.0}, Vec2{2.0, 0.0}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 2.0, 100.0, 5.0));

        std::vector<Vec2> forces(2);
        accumulate_spring_forces(n, bodies, forces);
        MALLOY_CHECK_TRUE(malloy::math::is_finite(forces[0]));
        MALLOY_CHECK_TRUE(malloy::math::is_finite(forces[1]));
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(0.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(0.0, 0.0), eps);
    }

    // --- A mismatched force buffer is a caller error and changes nothing. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{3.0, 0.0}, Vec2{}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 2.0, 10.0, 0.0));

        std::vector<Vec2> wrong_size(5, Vec2{1.0, 1.0});
        accumulate_spring_forces(n, bodies, wrong_size);
        for (const Vec2& f : wrong_size)
        {
            MALLOY_CHECK_VEC2_NEAR(f, Vec2(1.0, 1.0), 0.0);
        }
    }

    // --- A spring referencing a body that does not exist makes the world
    //     invalid, rather than being silently skipped every step. ---
    {
        SpringNetwork n;
        n.add(make_spring(0, 5, 1.0, 10.0, 0.0));
        SpringWorld w{SimulationSettings{0.001}, n,
                      {body_at(Vec2{}, Vec2{}, 1.0), body_at(Vec2{1.0, 0.0}, Vec2{}, 1.0)}};
        MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidState);
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
    }
    {
        SpringNetwork n;
        n.add(make_spring(1, 1, 1.0, 10.0, 0.0)); // self-spring
        SpringWorld w{SimulationSettings{0.001}, n,
                      {body_at(Vec2{}, Vec2{}, 1.0), body_at(Vec2{1.0, 0.0}, Vec2{}, 1.0)}};
        MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidState);
    }

    // --- Failed validation leaves the state untouched. ---
    {
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 10.0, 0.0));
        const std::vector<SpringBody2D> before = {
            body_at(Vec2{0.0, 0.0}, Vec2{1.0, 2.0}, 1.0),
            body_at(Vec2{3.0, 0.0}, Vec2{}, 1.0)};
        SpringWorld w{SimulationSettings{0.0}, n, before}; // dt invalid
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].position, before[0].position, 0.0);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, before[0].velocity, 0.0);
    }

    // --- A stretched spring pulls its endpoints together: the end-to-end
    //     behaviour a force-only library could not test. ---
    {
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 20.0, 0.0));
        SpringWorld w{SimulationSettings{0.001}, n,
                      {body_at(Vec2{-2.0, 0.0}, Vec2{}, 1.0),
                       body_at(Vec2{2.0, 0.0}, Vec2{}, 1.0)}};
        const Real start = malloy::math::distance(w.bodies()[1].position,
                                                  w.bodies()[0].position);
        for (int i = 0; i < 100; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const Real after = malloy::math::distance(w.bodies()[1].position,
                                                  w.bodies()[0].position);
        MALLOY_CHECK_TRUE(after < start); // it contracted
    }

    // --- INVARIANT: a spring network cannot change total momentum, at any
    //     stiffness and any damping, because every contribution is equal and
    //     opposite. ---
    {
        for (const Real damping : {0.0, 1.5})
        {
            SpringNetwork n;
            n.add(make_spring(0, 1, 1.0, 30.0, damping));
            n.add(make_spring(1, 2, 2.0, 12.0, damping));
            n.add(make_spring(0, 2, 4.0, 7.0, damping));
            const std::vector<SpringBody2D> start = {
                body_at(Vec2{-1.5, 0.3}, Vec2{0.4, -0.2}, 1.0),
                body_at(Vec2{0.7, -0.4}, Vec2{-0.1, 0.5}, 2.5),
                body_at(Vec2{1.9, 1.1}, Vec2{0.25, 0.15}, 0.75)};
            const Vec2 p0 = total_momentum(start);

            SpringWorld w{SimulationSettings{0.0005}, n, start};
            for (int i = 0; i < 4000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_VEC2_NEAR(total_momentum(w.bodies()), p0, 1e-10);
        }
    }

    // --- Damping removes energy and an undamped network does not. Kinetic plus
    //     elastic is what a spring conserves, so both terms are needed. ---
    {
        SpringNetwork undamped;
        undamped.add(make_spring(0, 1, 1.0, 40.0, 0.0));
        SpringNetwork damped;
        damped.add(make_spring(0, 1, 1.0, 40.0, 2.0));

        const std::vector<SpringBody2D> start = {
            body_at(Vec2{-1.5, 0.0}, Vec2{}, 1.0),
            body_at(Vec2{1.5, 0.0}, Vec2{}, 1.0)};
        const Real e0 = total_kinetic_energy(start) + total_elastic_energy(undamped, start);

        SpringWorld a{SimulationSettings{0.0002}, undamped, start};
        SpringWorld b{SimulationSettings{0.0002}, damped, start};
        for (int i = 0; i < 5000; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        const Real ea = total_kinetic_energy(a.bodies()) +
                        total_elastic_energy(a.network(), a.bodies());
        const Real eb = total_kinetic_energy(b.bodies()) +
                        total_elastic_energy(b.network(), b.bodies());

        // Undamped: bounded by the integrator's own drift, not exactly equal.
        MALLOY_CHECK_NEAR(ea, e0, e0 * 0.01);
        // Damped: strictly and obviously less, not merely different.
        MALLOY_CHECK_TRUE(eb < e0 * 0.5);
    }

    // --- Elastic energy is zero at rest length and grows with the square of
    //     the extension. ---
    {
        SpringNetwork n;
        n.add(make_spring(0, 1, 2.0, 8.0, 0.0));
        const std::vector<SpringBody2D> at_rest = {
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0), body_at(Vec2{2.0, 0.0}, Vec2{}, 1.0)};
        const std::vector<SpringBody2D> stretched = {
            body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0), body_at(Vec2{5.0, 0.0}, Vec2{}, 1.0)};
        MALLOY_CHECK_NEAR(total_elastic_energy(n, at_rest), 0.0, eps);
        // (1/2) k e^2 = 0.5 * 8 * 9 = 36
        MALLOY_CHECK_NEAR(total_elastic_energy(n, stretched), 36.0, eps);
    }

    // --- With no springs at all, bodies coast unchanged: the guard that M13
    //     did not alter anything about free translational motion. ---
    {
        SpringWorld w{SimulationSettings{0.5}, SpringNetwork{},
                      {body_at(Vec2{1.0, 2.0}, Vec2{3.0, -4.0}, 2.0)}};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(3.0, -4.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].position, Vec2(2.5, 0.0), eps);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{1});
    }

    // --- The accumulator is cleared each step, so forces do not compound. ---
    {
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 10.0, 0.0));
        SpringWorld w{SimulationSettings{0.0001}, n,
                      {body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
                       body_at(Vec2{2.0, 0.0}, Vec2{}, 1.0)}};
        MALLOY_CHECK_TRUE(w.step().ok());
        const Vec2 first = w.forces()[0];
        MALLOY_CHECK_TRUE(w.step().ok());
        const Vec2 second = w.forces()[0];
        // dt is tiny, so the geometry barely moved: without clearing, the second
        // reading would be roughly twice the first.
        MALLOY_CHECK_NEAR(second.x, first.x, first.x * 0.01);
    }

    // --- Determinism and elapsed time. ---
    {
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 25.0, 0.4));
        n.add(make_spring(1, 2, 1.5, 9.0, 0.2));
        const std::vector<SpringBody2D> start = {
            body_at(Vec2{-1.0, 0.2}, Vec2{0.3, 0.0}, 1.0),
            body_at(Vec2{0.4, -0.3}, Vec2{0.0, 0.2}, 2.0),
            body_at(Vec2{1.7, 0.6}, Vec2{-0.2, 0.1}, 0.5)};
        SpringWorld a{SimulationSettings{0.001}, n, start};
        SpringWorld b{SimulationSettings{0.001}, n, start};
        for (int i = 0; i < 500; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        for (std::size_t i = 0; i < start.size(); ++i)
        {
            MALLOY_CHECK_VEC2_NEAR(a.bodies()[i].position, b.bodies()[i].position, 0.0);
            MALLOY_CHECK_VEC2_NEAR(a.bodies()[i].velocity, b.bodies()[i].velocity, 0.0);
        }
    }
    {
        SpringWorld w{SimulationSettings{0.1}, SpringNetwork{},
                      {body_at(Vec2{}, Vec2{}, 1.0)}};
        for (int i = 0; i < 10; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_NEAR(w.elapsed_time(), 1.0, 1e-17);
        MALLOY_CHECK_NEAR(w.simulation_settings().dt, 0.1, eps);
    }

    // --- Empty world. ---
    {
        SpringWorld w{SimulationSettings{0.01}, SpringNetwork{}, {}};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{1});
        MALLOY_CHECK_VEC2_NEAR(total_momentum({}), Vec2(0.0, 0.0), eps);
        MALLOY_CHECK_NEAR(total_kinetic_energy({}), 0.0, eps);
        MALLOY_CHECK_NEAR(total_elastic_energy(SpringNetwork{}, {}), 0.0, eps);
    }

    // --- Every force-value assertion above uses a HORIZONTAL spring, so the
    //     y component of the separation is never exercised: replacing
    //     dot(delta, delta) with delta.x * delta.x passes all of them. That is
    //     the asymmetric-configuration rule in docs/05, and this is the case
    //     that enforces it.
    //
    //     Off-axis and off-origin. delta = (3, 4), so length is exactly 5 and
    //     the axis is exactly (0.6, 0.8) with no rounding. ---
    {
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{1.0, 2.0}, Vec2{-1.0, 1.0}, 1.0),
            body_at(Vec2{4.0, 6.0}, Vec2{1.0, 2.0}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 3.0, 10.0, 3.0));

        // extension = 5 - 3 = 2, so the spring term is 10 * 2 = 20.
        // relative velocity = (2, 1), axial speed = dot((2,1),(0.6,0.8)) = 2,
        // so the damper term is 3 * 2 = 6. Total 26 along (0.6, 0.8).
        std::vector<Vec2> forces = {Vec2{}, Vec2{}};
        accumulate_spring_forces(n, bodies, forces);
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(15.6, 20.8), eps);
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(-15.6, -20.8), eps);
    }
    {
        // A purely VERTICAL spring, where delta.x is zero. Under the horizontal
        // -only mistake this hits the coincident-endpoint guard and produces
        // exactly nothing, which is the failure that would be silent.
        const std::vector<SpringBody2D> bodies = {
            body_at(Vec2{-2.0, 1.0}, Vec2{}, 1.0),
            body_at(Vec2{-2.0, 4.0}, Vec2{}, 1.0)};
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 5.0, 0.0));

        std::vector<Vec2> forces = {Vec2{}, Vec2{}};
        accumulate_spring_forces(n, bodies, forces);
        MALLOY_CHECK_VEC2_NEAR(forces[0], Vec2(0.0, 10.0), eps); // 5 * (3 - 1)
        MALLOY_CHECK_VEC2_NEAR(forces[1], Vec2(0.0, -10.0), eps);
    }

    // --- Issue #15, first half: a separation whose SQUARE overflows is refused
    //     rather than silently deleting the spring.
    //
    //     The force kernel computes dot(delta, delta), which reaches infinity
    //     once the endpoints are more than sqrt(DBL_MAX) apart, about 1.34e154,
    //     while both positions are still perfectly finite. The kernel's guard
    //     then cannot distinguish that from coincident endpoints and skips the
    //     spring, so the spring stops existing and the bodies coast apart
    //     forever with every step reporting Ok. Measured before the fix: two
    //     bodies 2e154 apart with stiffness 1 reported kinetic energy of
    //     exactly 0.00000000e+00 for 2000 steps. ---
    {
        SpringNetwork n;
        n.add(make_spring(0, 1, 1.0, 1.0, 0.0));

        // 2e150 apart: the square is 4e300, still representable, so this runs.
        SpringWorld near{SimulationSettings{0.001}, n,
                         {body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
                          body_at(Vec2{2.0e150, 0.0}, Vec2{}, 1.0)}};
        MALLOY_CHECK_TRUE(near.validate() == StepStatus::Ok);
        MALLOY_CHECK_TRUE(near.step().ok());
        // It really did feel the spring rather than being skipped.
        MALLOY_CHECK_TRUE(near.bodies()[0].velocity.x > 0.0);

        // 2e154 apart: the square overflows, and that is now reported.
        SpringWorld far{SimulationSettings{0.001}, n,
                        {body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
                         body_at(Vec2{2.0e154, 0.0}, Vec2{}, 1.0)}};
        MALLOY_CHECK_TRUE(far.validate() == StepStatus::InvalidState);
        MALLOY_CHECK_FALSE(far.step().ok());
        // And the state is untouched by the refusal.
        MALLOY_CHECK_NEAR(far.bodies()[1].position.x, 2.0e154, 0.0);
    }

    // --- Issue #15, second half: a timestep past a spring's stability limit is
    //     refused. Symplectic Euler on a spring pair diverges geometrically
    //     past it, and until the separation reaches the range above, nothing
    //     reports that.
    //
    //     For two unit masses the reduced mass is 0.5, so omega^2 = 2k and the
    //     undamped bound dt < 2/omega becomes k < 2/dt^2. At dt = 5e-4 that is
    //     exactly k = 8.0e6, and the boundary is checked from both sides. ---
    {
        const Real dt = 5.0e-4;
        const auto pair_world = [&](Real stiffness, Real damping) {
            SpringNetwork n;
            n.add(make_spring(0, 1, 1.0, stiffness, damping));
            return SpringWorld{SimulationSettings{dt}, n,
                               {body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
                                body_at(Vec2{1.0, 0.0}, Vec2{}, 1.0)}};
        };

        MALLOY_CHECK_TRUE(pair_world(7.9e6, 0.0).validate() == StepStatus::Ok);
        MALLOY_CHECK_TRUE(pair_world(8.1e6, 0.0).validate() ==
                          StepStatus::InvalidSettings);

        // Zero stiffness has no limit at all: there is no oscillation to be
        // unstable.
        MALLOY_CHECK_TRUE(pair_world(0.0, 0.0).validate() == StepStatus::Ok);

        // Damping TIGHTENS the bound rather than helping, which is the part
        // that would be got wrong by using the undamped formula. With
        // k = 2e6, omega is 2000 and the undamped limit is dt < 1e-3. Adding
        // c = 1000 gives gamma = 1000 and a limit of
        // 2 (sqrt(1e6 + 4e6) - 1000) / 4e6 = 6.180e-4, so a timestep of 8e-4
        // is inside the undamped bound and outside the real one.
        SpringNetwork damped;
        damped.add(make_spring(0, 1, 1.0, 2.0e6, 1000.0));
        SpringWorld w{SimulationSettings{8.0e-4}, damped,
                      {body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
                       body_at(Vec2{1.0, 0.0}, Vec2{}, 1.0)}};
        MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidSettings);

        SpringNetwork undamped;
        undamped.add(make_spring(0, 1, 1.0, 2.0e6, 0.0));
        SpringWorld u{SimulationSettings{8.0e-4}, undamped,
                      {body_at(Vec2{0.0, 0.0}, Vec2{}, 1.0),
                       body_at(Vec2{1.0, 0.0}, Vec2{}, 1.0)}};
        MALLOY_CHECK_TRUE(u.validate() == StepStatus::Ok);
    }

    // --- Issue #16: a position or velocity whose SQUARE overflows is refused.
    //     Finite is not enough. Everything that squares a vector reaches
    //     infinity above about 1.34e154, so accepting state up to 1.8e308 left
    //     a window in which a world validated clean while every energy it
    //     reported was inf. ---
    {
        const Real too_big = 1.4e154; // finite, and its square is not
        MALLOY_CHECK_TRUE(body_at(Vec2{1.0, 2.0}, Vec2{}, 1.0).is_valid());
        MALLOY_CHECK_FALSE(body_at(Vec2{too_big, 0.0}, Vec2{}, 1.0).is_valid());
        MALLOY_CHECK_FALSE(body_at(Vec2{}, Vec2{0.0, too_big}, 1.0).is_valid());
    }

    std::cout << "malloy_springs_tests passed\n";
    return 0;
}
