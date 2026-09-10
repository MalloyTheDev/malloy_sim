#include <malloy/rigid/rigid.hpp>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <test_check.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using malloy::collide::Aabb;
using malloy::collide::Circle;
using malloy::math::Real;
using malloy::math::Vec2;
using malloy::rigid::center_of_mass;
using malloy::rigid::mass_properties;
using malloy::rigid::MassProperties;
using malloy::rigid::RigidBody2D;
using malloy::rigid::RigidWorld;
using malloy::rigid::shift_inertia;
using malloy::rigid::to_local;
using malloy::rigid::to_world;
using malloy::rigid::total_angular_momentum;
using malloy::rigid::total_kinetic_energy;
using malloy::rigid::total_linear_momentum;
using malloy::rigid::velocity_at;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;

namespace
{
// The deliberately asymmetric body this milestone is tested against, per
// ADR 0007 and docs/05. Every one of these is true at once:
//   body origin is NOT the centre of mass
//   the angle is nonzero
//   inertia is not 1
//   mass is not 1
//   both velocity components are nonzero
//   the spin is nonzero
// A body missing any of these hides a whole class of mistake.
RigidBody2D awkward_body()
{
    RigidBody2D b;
    b.position = Vec2{3.0, -2.0};
    b.angle = 0.7;
    b.velocity = Vec2{1.5, -0.25};
    b.angular_velocity = 0.4;
    b.mass = 2.5;
    b.inertia = 3.75;
    b.local_center_of_mass = Vec2{0.6, -0.4};
    return b;
}
} // namespace

int main()
{
    const Real eps = 1e-12;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real pi = 3.14159265358979323846;

    // --- Validation ---
    {
        MALLOY_CHECK_TRUE(awkward_body().is_valid());

        RigidBody2D zero_mass = awkward_body();
        zero_mass.mass = 0.0;
        MALLOY_CHECK_FALSE(zero_mass.is_valid());

        RigidBody2D zero_inertia = awkward_body();
        zero_inertia.inertia = 0.0;
        MALLOY_CHECK_FALSE(zero_inertia.is_valid());

        // Infinity is now VALID and means immovable. M11 rejected it because
        // nothing needed statics until contact response did (ADR 0007); M14
        // needs them, so the contract changed deliberately rather than drifting.
        RigidBody2D inf_inertia = awkward_body();
        inf_inertia.inertia = inf;
        MALLOY_CHECK_TRUE(inf_inertia.is_valid());
        MALLOY_CHECK_TRUE(inf_inertia.is_static());

        // NaN is still invalid, because NaN > 0 is false.
        RigidBody2D nan_mass = awkward_body();
        nan_mass.mass = nan;
        MALLOY_CHECK_FALSE(nan_mass.is_valid());

        // A radius must be non-negative and finite.
        RigidBody2D bad_radius = awkward_body();
        bad_radius.radius = -1.0;
        MALLOY_CHECK_FALSE(bad_radius.is_valid());
        RigidBody2D inf_radius = awkward_body();
        inf_radius.radius = inf;
        MALLOY_CHECK_FALSE(inf_radius.is_valid());

        RigidBody2D nan_angle = awkward_body();
        nan_angle.angle = nan;
        MALLOY_CHECK_FALSE(nan_angle.is_valid());

        RigidBody2D nan_spin = awkward_body();
        nan_spin.angular_velocity = nan;
        MALLOY_CHECK_FALSE(nan_spin.is_valid());

        RigidBody2D bad_com = awkward_body();
        bad_com.local_center_of_mass = Vec2{inf, 0.0};
        MALLOY_CHECK_FALSE(bad_com.is_valid());
    }

    // --- Mass properties from geometry. Density is what turns area into mass,
    //     and it is not 1, so a dropped density factor shows. ---
    {
        // Disc R = 2, density 3: m = 3 * 4pi = 12pi, I = 3 * 8pi = 24pi.
        const MassProperties disc = mass_properties(Circle{Vec2{5.0, 1.0}, 2.0}, 3.0);
        MALLOY_CHECK_NEAR(disc.mass, 12.0 * pi, 1e-12);
        MALLOY_CHECK_NEAR(disc.inertia, 24.0 * pi, 1e-12);
        MALLOY_CHECK_VEC2_NEAR(disc.center_of_mass, Vec2(5.0, 1.0), eps);
        // Cross-check against the textbook disc relation I = m R^2 / 2.
        MALLOY_CHECK_NEAR(disc.inertia, disc.mass * 4.0 / 2.0, 1e-12);

        // Box 6 by 2, density 0.5: m = 6, I = 0.5 * 40 = 20.
        const MassProperties box =
            mass_properties(Aabb{Vec2{1.0, 5.0}, Vec2{7.0, 7.0}}, 0.5);
        MALLOY_CHECK_NEAR(box.mass, 6.0, eps);
        MALLOY_CHECK_NEAR(box.inertia, 20.0, eps);
        MALLOY_CHECK_VEC2_NEAR(box.center_of_mass, Vec2(4.0, 6.0), eps);
        // Textbook: I = m(w^2 + h^2)/12 = 6 * 40 / 12 = 20.
        MALLOY_CHECK_NEAR(box.inertia, box.mass * 40.0 / 12.0, eps);

        // A non-positive density or an invalid shape yields nothing usable.
        MALLOY_CHECK_NEAR(mass_properties(Circle{Vec2{}, 2.0}, 0.0).mass, 0.0, eps);
        MALLOY_CHECK_NEAR(mass_properties(Circle{Vec2{}, -1.0}, 3.0).mass, 0.0, eps);
        MALLOY_CHECK_NEAR(mass_properties(Circle{Vec2{}, 2.0}, nan).inertia, 0.0, eps);
    }

    // --- Parallel-axis theorem, explicit rather than emergent. ---
    {
        // I = I_com + m d^2 = 20 + 6 * 9 = 74.
        MALLOY_CHECK_NEAR(shift_inertia(20.0, 6.0, 3.0), 74.0, eps);
        // Zero distance leaves it alone.
        MALLOY_CHECK_NEAR(shift_inertia(20.0, 6.0, 0.0), 20.0, eps);
        // A bad distance returns the input rather than a plausible wrong answer.
        MALLOY_CHECK_NEAR(shift_inertia(20.0, 6.0, -3.0), 20.0, eps);
        MALLOY_CHECK_NEAR(shift_inertia(20.0, 6.0, nan), 20.0, eps);
    }

    // --- Pose conversions round-trip, on a body with a nonzero angle. ---
    {
        const RigidBody2D b = awkward_body();
        const Vec2 local{1.3, -0.9};
        MALLOY_CHECK_VEC2_NEAR(to_local(b, to_world(b, local)), local, 1e-12);

        // A quarter turn maps local +x onto world +y, which pins the rotation
        // direction: counter-clockwise positive.
        RigidBody2D q;
        q.angle = pi / 2.0;
        q.position = Vec2{0.0, 0.0};
        MALLOY_CHECK_VEC2_NEAR(to_world(q, Vec2{1.0, 0.0}), Vec2(0.0, 1.0), 1e-12);

        // The centre of mass is the origin plus the ROTATED local offset, so it
        // is not simply position + local_center_of_mass unless the angle is 0.
        RigidBody2D c;
        c.position = Vec2{10.0, 20.0};
        c.angle = pi / 2.0;
        c.local_center_of_mass = Vec2{2.0, 0.0};
        MALLOY_CHECK_VEC2_NEAR(center_of_mass(c), Vec2(10.0, 22.0), 1e-12);
    }

    // --- velocity_at measures its arm from the centre of mass, not the origin.
    //     With the two coinciding this test would pass either way, so the body
    //     used here deliberately offsets them. ---
    {
        RigidBody2D b;
        b.position = Vec2{0.0, 0.0};
        b.angle = 0.0;
        b.local_center_of_mass = Vec2{1.0, 0.0}; // COM at (1, 0)
        b.velocity = Vec2{0.0, 0.0};
        b.angular_velocity = 2.0;
        b.mass = 1.0;
        b.inertia = 1.0;

        // At the centre of mass the rotational term vanishes.
        MALLOY_CHECK_VEC2_NEAR(velocity_at(b, Vec2{1.0, 0.0}), Vec2(0.0, 0.0), eps);
        // One unit to the +x side of the COM: omega x r = (0, 2).
        MALLOY_CHECK_VEC2_NEAR(velocity_at(b, Vec2{2.0, 0.0}), Vec2(0.0, 2.0), eps);
        // At the body ORIGIN, one unit on the -x side: (0, -2). Measuring the
        // arm from the origin instead would wrongly give (0, 0) here.
        MALLOY_CHECK_VEC2_NEAR(velocity_at(b, Vec2{0.0, 0.0}), Vec2(0.0, -2.0), eps);
    }

    // --- Free motion: nothing accelerates a body, so the spin is unchanged and
    //     the centre of mass travels in a straight line at velocity * dt. ---
    {
        const RigidBody2D start = awkward_body();
        RigidWorld w{SimulationSettings{0.5}, {start}};
        const Vec2 com0 = center_of_mass(start);
        MALLOY_CHECK_TRUE(w.step().ok());

        const RigidBody2D& b = w.bodies()[0];
        MALLOY_CHECK_NEAR(b.angular_velocity, start.angular_velocity, eps);
        MALLOY_CHECK_VEC2_NEAR(b.velocity, start.velocity, eps);
        MALLOY_CHECK_NEAR(b.angle, start.angle + start.angular_velocity * 0.5, eps);
        MALLOY_CHECK_VEC2_NEAR(center_of_mass(b), com0 + start.velocity * 0.5, 1e-12);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{1});
    }

    // --- A spinning body whose origin is offset from its centre of mass must
    //     ORBIT that centre, not rotate about its own origin. With zero linear
    //     velocity the centre of mass must not move at all. ---
    {
        RigidBody2D b;
        b.position = Vec2{0.0, 0.0};
        b.angle = 0.0;
        b.local_center_of_mass = Vec2{1.0, 0.0};
        b.angular_velocity = 1.0;
        b.mass = 2.0;
        b.inertia = 1.0;

        RigidWorld w{SimulationSettings{pi / 2.0}, {b}};
        MALLOY_CHECK_TRUE(w.step().ok()); // a quarter turn

        // The centre of mass stayed at (1, 0)...
        MALLOY_CHECK_VEC2_NEAR(center_of_mass(w.bodies()[0]), Vec2(1.0, 0.0), 1e-12);
        // ...while the origin swung round it to (1, -1).
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].position, Vec2(1.0, -1.0), 1e-12);
    }

    // --- Angle is NOT wrapped by integration (ADR 0007): after many turns it
    //     keeps accumulating rather than being silently canonicalised. ---
    {
        RigidBody2D b = awkward_body();
        b.velocity = Vec2{};
        b.angular_velocity = 1.0;
        RigidWorld w{SimulationSettings{1.0}, {b}};
        for (int i = 0; i < 20; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_NEAR(w.bodies()[0].angle, b.angle + 20.0, 1e-12);
        MALLOY_CHECK_TRUE(w.bodies()[0].angle > 2.0 * pi); // not wrapped
    }

    // --- INVARIANT: with no impulses, linear and angular momentum and kinetic
    //     energy are all conserved exactly over a long run. ---
    {
        std::vector<RigidBody2D> bodies = {awkward_body()};
        RigidBody2D second = awkward_body();
        second.position = Vec2{-4.0, 6.0};
        second.angle = -1.1;
        second.velocity = Vec2{-0.7, 2.0};
        second.angular_velocity = -0.9;
        second.mass = 1.25;
        second.inertia = 0.8;
        second.local_center_of_mass = Vec2{-0.3, 0.5};
        bodies.push_back(second);

        const Vec2 p0 = total_linear_momentum(bodies);
        const Real l0 = total_angular_momentum(bodies);
        const Real k0 = total_kinetic_energy(bodies);

        RigidWorld w{SimulationSettings{0.001}, bodies};
        for (int i = 0; i < 5000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(w.bodies()), p0, 1e-12);
        MALLOY_CHECK_NEAR(total_angular_momentum(w.bodies()), l0, 1e-9);
        MALLOY_CHECK_NEAR(total_kinetic_energy(w.bodies()), k0, 1e-12);
    }

    // --- An impulse THROUGH the centre of mass changes linear velocity only.
    //     This is the cleanest check that the torque arm is measured from the
    //     centre of mass: measuring from the body origin would spin it up. ---
    {
        const RigidBody2D start = awkward_body();
        RigidWorld w{SimulationSettings{0.01}, {start}};
        const Vec2 impulse{2.0, -1.0};
        MALLOY_CHECK_TRUE(w.apply_impulse_at(0, impulse, center_of_mass(start)));

        const RigidBody2D& b = w.bodies()[0];
        MALLOY_CHECK_VEC2_NEAR(b.velocity, start.velocity + impulse / start.mass, 1e-12);
        MALLOY_CHECK_NEAR(b.angular_velocity, start.angular_velocity, 1e-12);
    }

    // --- An OFF-CENTRE impulse changes both, and the angular change is exactly
    //     (r x J) / I. Both components of J are nonzero and both components of
    //     r are nonzero, so neither half of the determinant can vanish. ---
    {
        const RigidBody2D start = awkward_body();
        RigidWorld w{SimulationSettings{0.01}, {start}};

        const Vec2 point = center_of_mass(start) + Vec2{0.8, -0.5};
        const Vec2 impulse{1.5, 2.25};
        MALLOY_CHECK_TRUE(w.apply_impulse_at(0, impulse, point));

        const Real expected_torque = 0.8 * 2.25 - (-0.5) * 1.5; // rx*Jy - ry*Jx
        MALLOY_CHECK_TRUE(std::abs(expected_torque) > 1e-6);    // genuinely nonzero
        const RigidBody2D& b = w.bodies()[0];
        MALLOY_CHECK_VEC2_NEAR(b.velocity, start.velocity + impulse / start.mass, 1e-12);
        MALLOY_CHECK_NEAR(b.angular_velocity,
                          start.angular_velocity + expected_torque / start.inertia,
                          1e-12);
    }

    // --- Both torque signs occur, and mirroring the arm mirrors the spin. ---
    {
        const RigidBody2D start = awkward_body();
        const Vec2 arm{0.8, -0.5};
        const Vec2 impulse{1.5, 2.25};

        RigidWorld plus{SimulationSettings{0.01}, {start}};
        MALLOY_CHECK_TRUE(
            plus.apply_impulse_at(0, impulse, center_of_mass(start) + arm));
        RigidWorld minus{SimulationSettings{0.01}, {start}};
        MALLOY_CHECK_TRUE(
            minus.apply_impulse_at(0, impulse, center_of_mass(start) - arm));

        const Real up = plus.bodies()[0].angular_velocity - start.angular_velocity;
        const Real down = minus.bodies()[0].angular_velocity - start.angular_velocity;
        MALLOY_CHECK_TRUE(up * down < 0.0); // opposite signs
        MALLOY_CHECK_NEAR(up, -down, 1e-12);
    }

    // --- An impulse changes total linear momentum by exactly J, and total
    //     angular momentum about the origin by exactly (p x J). ---
    {
        const RigidBody2D start = awkward_body();
        RigidWorld w{SimulationSettings{0.01}, {start}};
        const std::vector<RigidBody2D> before = {start};

        const Vec2 point = center_of_mass(start) + Vec2{0.8, -0.5};
        const Vec2 impulse{1.5, 2.25};
        MALLOY_CHECK_TRUE(w.apply_impulse_at(0, impulse, point));

        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(w.bodies()),
                               total_linear_momentum(before) + impulse, 1e-12);
        const Real expected = point.x * impulse.y - point.y * impulse.x;
        MALLOY_CHECK_NEAR(total_angular_momentum(w.bodies()) -
                              total_angular_momentum(before),
                          expected, 1e-9);
    }

    // --- apply_impulse_at refuses bad input rather than corrupting a body. ---
    {
        const RigidBody2D start = awkward_body();
        RigidWorld w{SimulationSettings{0.01}, {start}};
        MALLOY_CHECK_FALSE(w.apply_impulse_at(5, Vec2{1.0, 0.0}, Vec2{}));
        MALLOY_CHECK_FALSE(w.apply_impulse_at(0, Vec2{nan, 0.0}, Vec2{}));
        MALLOY_CHECK_FALSE(w.apply_impulse_at(0, Vec2{1.0, 0.0}, Vec2{inf, 0.0}));
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, start.velocity, 0.0);
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, start.angular_velocity, 0.0);
    }

    // --- Failed validation leaves the state untouched. ---
    {
        const RigidBody2D start = awkward_body();
        RigidWorld w{SimulationSettings{0.0}, {start}};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
        MALLOY_CHECK_NEAR(w.elapsed_time(), 0.0, 0.0);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].position, start.position, 0.0);
        MALLOY_CHECK_NEAR(w.bodies()[0].angle, start.angle, 0.0);
    }
    {
        RigidBody2D bad = awkward_body();
        bad.inertia = -1.0;
        RigidWorld w{SimulationSettings{0.01}, {bad}};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
    }

    // --- Elapsed time comes from malloy_time: ticks * dt, not accumulated. ---
    {
        RigidWorld w{SimulationSettings{0.1}, {awkward_body()}};
        for (int i = 0; i < 10; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_NEAR(w.elapsed_time(), 1.0, 1e-17);
    }

    // --- Determinism: identical worlds stepped identically stay identical. ---
    {
        const std::vector<RigidBody2D> start = {awkward_body()};
        RigidWorld a{SimulationSettings{0.003}, start};
        RigidWorld b{SimulationSettings{0.003}, start};
        for (int i = 0; i < 500; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        MALLOY_CHECK_VEC2_NEAR(a.bodies()[0].position, b.bodies()[0].position, 0.0);
        MALLOY_CHECK_NEAR(a.bodies()[0].angle, b.bodies()[0].angle, 0.0);
    }

    // --- Empty world: diagnostics are zero and a step is a no-op that still
    //     counts, matching the other domains. ---
    {
        const std::vector<RigidBody2D> none;
        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(none), Vec2(0.0, 0.0), eps);
        MALLOY_CHECK_NEAR(total_angular_momentum(none), 0.0, eps);
        MALLOY_CHECK_NEAR(total_kinetic_energy(none), 0.0, eps);

        RigidWorld w{SimulationSettings{0.01}, none};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{1});
    }

    // --- Settings accessor, previously untested public API. ---
    {
        RigidWorld w{SimulationSettings{0.0025}, {awkward_body()}};
        MALLOY_CHECK_NEAR(w.simulation_settings().dt, 0.0025, eps);
    }

    // --- M14: statics. Infinite mass and inertia give exactly zero inverses, so
    //     every impulse formula handles immovability without a special case. ---
    {
        RigidBody2D wall = awkward_body();
        wall.mass = inf;
        wall.inertia = inf;
        MALLOY_CHECK_TRUE(wall.is_valid());
        MALLOY_CHECK_TRUE(wall.is_static());
        MALLOY_CHECK_NEAR(wall.inverse_mass(), 0.0, 0.0);
        MALLOY_CHECK_NEAR(wall.inverse_inertia(), 0.0, 0.0);

        const RigidBody2D movable = awkward_body();
        MALLOY_CHECK_FALSE(movable.is_static());
        MALLOY_CHECK_NEAR(movable.inverse_mass(), 1.0 / 2.5, eps);
        MALLOY_CHECK_NEAR(movable.inverse_inertia(), 1.0 / 3.75, eps);
    }

    // --- An impulse cannot move an immovable body. ---
    {
        RigidBody2D wall = awkward_body();
        wall.mass = inf;
        wall.inertia = inf;
        wall.velocity = Vec2{};
        wall.angular_velocity = 0.0;
        RigidWorld w{SimulationSettings{0.01}, {wall}};
        MALLOY_CHECK_TRUE(w.apply_impulse_at(0, Vec2{1000.0, -500.0},
                                             center_of_mass(wall) + Vec2{1.0, 1.0}));
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(0.0, 0.0), 0.0);
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 0.0);
    }

    // --- Two discs approaching head on collide and separate. Their centres of
    //     mass are at their disc centres here, so the normal passes through
    //     both and neither should spin: the pure translational case. ---
    {
        RigidBody2D a;
        a.position = Vec2{-1.0, 0.0};
        a.velocity = Vec2{1.0, 0.0};
        a.mass = 1.0;
        a.inertia = 1.0;
        a.radius = 0.5;
        RigidBody2D b = a;
        b.position = Vec2{1.0, 0.0};
        b.velocity = Vec2{-1.0, 0.0};

        RigidWorld w{SimulationSettings{0.5}, {a, b}, 1.0};
        MALLOY_CHECK_TRUE(w.step().ok()); // brings them to +-0.5, exactly touching
        // Equal masses, elastic, head on: they exchange velocities.
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(-1.0, 0.0), 1e-12);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(1.0, 0.0), 1e-12);
        // Normal through both centres of mass means no torque at all.
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 1e-12);
        MALLOY_CHECK_NEAR(w.bodies()[1].angular_velocity, 0.0, 1e-12);
    }

    // --- THE M14 POINT: an off-centre contact SPINS the body. The disc is
    //     centred on the body origin and the centre of mass is offset from it,
    //     so the contact normal misses the centre of mass and generates torque.
    //     With the two coinciding this test would show no spin at all. ---
    {
        // Starting x is derived, not guessed: after one step at dt = 0.5 the
        // ball centre reaches 0.3, giving a centre gap of 0.7 against radii
        // summing to 1.0, so it overlaps by 0.3 and the contact arm is nonzero.
        RigidBody2D ball;
        ball.position = Vec2{-0.2, 0.0};
        ball.velocity = Vec2{1.0, 0.0};
        ball.angular_velocity = 0.0;
        ball.mass = 1.0;
        ball.inertia = 0.5;
        ball.radius = 0.5;
        ball.local_center_of_mass = Vec2{0.0, 0.4}; // offset, so the arm is nonzero

        RigidBody2D wall;
        wall.position = Vec2{1.0, 0.0};
        wall.mass = inf;
        wall.inertia = inf;
        wall.radius = 0.5;

        RigidWorld w{SimulationSettings{0.5}, {ball, wall}, 1.0};
        MALLOY_CHECK_TRUE(w.step().ok());

        // It bounced back off the immovable body...
        MALLOY_CHECK_TRUE(w.bodies()[0].velocity.x < 0.0);
        // ...and picked up spin from the off-centre contact, which is the whole
        // capability this milestone adds.
        MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angular_velocity) > 1e-6);
        // The wall did not budge.
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(0.0, 0.0), 0.0);
        MALLOY_CHECK_NEAR(w.bodies()[1].angular_velocity, 0.0, 0.0);
    }

    // --- Mirroring the offset mirrors the spin, so both torque signs occur. ---
    {
        RigidBody2D ball;
        ball.position = Vec2{-0.2, 0.0};
        ball.velocity = Vec2{1.0, 0.0};
        ball.mass = 1.0;
        ball.inertia = 0.5;
        ball.radius = 0.5;
        RigidBody2D wall;
        wall.position = Vec2{1.0, 0.0};
        wall.mass = inf;
        wall.inertia = inf;
        wall.radius = 0.5;

        RigidBody2D up = ball;
        up.local_center_of_mass = Vec2{0.0, 0.4};
        RigidBody2D down = ball;
        down.local_center_of_mass = Vec2{0.0, -0.4};

        RigidWorld a{SimulationSettings{0.5}, {up, wall}, 1.0};
        RigidWorld b{SimulationSettings{0.5}, {down, wall}, 1.0};
        MALLOY_CHECK_TRUE(a.step().ok());
        MALLOY_CHECK_TRUE(b.step().ok());
        const Real spin_a = a.bodies()[0].angular_velocity;
        const Real spin_b = b.bodies()[0].angular_velocity;
        MALLOY_CHECK_TRUE(spin_a * spin_b < 0.0); // opposite signs
        MALLOY_CHECK_NEAR(spin_a, -spin_b, 1e-12);
    }

    // --- A zero radius means the body does not take part in contacts. ---
    {
        RigidBody2D a;
        a.position = Vec2{-0.1, 0.0};
        a.velocity = Vec2{1.0, 0.0};
        a.mass = 1.0;
        a.inertia = 1.0;
        a.radius = 0.0; // no collision shape
        RigidBody2D b = a;
        b.position = Vec2{0.1, 0.0};
        b.velocity = Vec2{-1.0, 0.0};
        b.radius = 0.5;

        RigidWorld w{SimulationSettings{0.001}, {a, b}, 1.0};
        MALLOY_CHECK_TRUE(w.step().ok());
        // They pass straight through each other: velocities unchanged.
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(1.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(-1.0, 0.0), eps);
    }

    // --- INVARIANT: a contact between two free bodies conserves both linear
    //     and angular momentum, at any restitution, because the impulse is
    //     equal and opposite and acts at one shared point. ---
    {
        for (const Real restitution : {0.0, 0.4, 1.0})
        {
            RigidBody2D a;
            a.position = Vec2{-0.9, 0.15};
            a.velocity = Vec2{1.4, -0.3};
            a.angular_velocity = 0.5;
            a.mass = 1.25;
            a.inertia = 0.8;
            a.radius = 0.5;
            a.local_center_of_mass = Vec2{0.1, 0.2};

            RigidBody2D b;
            b.position = Vec2{0.9, -0.1};
            b.velocity = Vec2{-0.7, 0.2};
            b.angular_velocity = -0.3;
            b.mass = 2.75;
            b.inertia = 1.6;
            b.radius = 0.5;
            b.local_center_of_mass = Vec2{-0.15, 0.05};

            const std::vector<RigidBody2D> start = {a, b};
            const Vec2 p0 = total_linear_momentum(start);
            const Real l0 = total_angular_momentum(start);

            RigidWorld w{SimulationSettings{0.002}, start, restitution};
            for (int i = 0; i < 600; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            // Linear momentum is exact: the impulse is equal and opposite, and
            // the positional correction moves positions without touching any
            // velocity.
            MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(w.bodies()), p0, 1e-10);
            // Angular momentum is NOT exact, and the reason is the positional
            // correction rather than the impulse. Moving positions without
            // changing velocities shifts the orbital term m (r x v) by
            // c x (v_b - v_a), where c is the correction. Verified by disabling
            // the correction, which makes this hold to 1e-8. The impulse itself
            // is checked exactly in the zero-penetration case below.
            MALLOY_CHECK_NEAR(total_angular_momentum(w.bodies()), l0, 1e-2);
        }
    }

    // --- The IMPULSE conserves angular momentum exactly. Set up a contact at
    //     exactly zero penetration so the positional correction does nothing,
    //     which isolates the impulse from the correction artefact above. ---
    {
        RigidBody2D a;
        a.position = Vec2{-0.5, 0.0};
        a.velocity = Vec2{1.3, 0.4};
        a.angular_velocity = 0.6;
        a.mass = 1.25;
        a.inertia = 0.8;
        a.radius = 0.5;
        a.local_center_of_mass = Vec2{0.1, 0.2};

        RigidBody2D b;
        b.position = Vec2{0.5, 0.0};
        b.velocity = Vec2{-0.9, -0.2};
        b.angular_velocity = -0.35;
        b.mass = 2.75;
        b.inertia = 1.6;
        b.radius = 0.5;
        b.local_center_of_mass = Vec2{-0.15, 0.05};

        // Centres exactly 1.0 apart against radii summing to 1.0: touching, so
        // penetration is exactly zero and the correction is a no-op.
        const std::vector<RigidBody2D> start = {a, b};
        const Vec2 p0 = total_linear_momentum(start);
        const Real l0 = total_angular_momentum(start);

        // dt of zero would be invalid, so take the impulse on the first step
        // with a tiny dt. Integration still introduces about 2e-9 of
        // penetration, so the correction still fires and still perturbs the
        // angular momentum, just proportionally less.
        //
        // That proportionality is the evidence. The long run above drifts by
        // about 5e-4 with ordinary penetrations; here, with penetration nine
        // orders smaller, the drift is around 5e-9. The correction is the
        // cause, and the impulse is exact.
        RigidWorld w{SimulationSettings{1e-9}, start, 1.0};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(w.bodies()), p0, 1e-12);
        MALLOY_CHECK_NEAR(total_angular_momentum(w.bodies()), l0, 1e-7);
    }

    // --- Restitution 1 conserves kinetic energy across a contact, and below 1
    //     strictly removes some. Both terms of the kinetic energy matter,
    //     because an off-centre contact moves energy into spin. ---
    {
        RigidBody2D a;
        a.position = Vec2{-0.9, 0.15};
        a.velocity = Vec2{1.4, -0.3};
        a.mass = 1.25;
        a.inertia = 0.8;
        a.radius = 0.5;
        a.local_center_of_mass = Vec2{0.1, 0.2};
        RigidBody2D b;
        b.position = Vec2{0.9, -0.1};
        b.velocity = Vec2{-0.7, 0.2};
        b.mass = 2.75;
        b.inertia = 1.6;
        b.radius = 0.5;
        b.local_center_of_mass = Vec2{-0.15, 0.05};

        const std::vector<RigidBody2D> start = {a, b};
        const Real k0 = total_kinetic_energy(start);

        RigidWorld elastic{SimulationSettings{0.002}, start, 1.0};
        RigidWorld lossy{SimulationSettings{0.002}, start, 0.3};
        for (int i = 0; i < 600; ++i)
        {
            MALLOY_CHECK_TRUE(elastic.step().ok());
            MALLOY_CHECK_TRUE(lossy.step().ok());
        }
        MALLOY_CHECK_NEAR(total_kinetic_energy(elastic.bodies()), k0, 1e-9);
        MALLOY_CHECK_TRUE(total_kinetic_energy(lossy.bodies()) < k0);
    }

    // --- The relative velocity at a contact must include omega x r, not just
    //     the centre-of-mass velocities. A body with ZERO linear velocity that
    //     is spinning has a moving surface, and that surface can be closing on
    //     something even though its centre of mass is not. ---
    {
        RigidBody2D spinner;
        spinner.position = Vec2{-0.5, 0.0};
        spinner.velocity = Vec2{};          // centre of mass going nowhere
        spinner.angular_velocity = 3.0;     // but spinning hard
        spinner.mass = 1.0;
        spinner.inertia = 0.5;
        spinner.radius = 0.5;
        spinner.local_center_of_mass = Vec2{0.0, 0.4}; // arm at the contact

        RigidBody2D wall;
        wall.position = Vec2{0.5, 0.0};
        wall.mass = inf;
        wall.inertia = inf;
        wall.radius = 0.5;

        // Exactly touching, so the contact exists with zero penetration.
        RigidWorld w{SimulationSettings{1e-9}, {spinner, wall}, 1.0};
        MALLOY_CHECK_TRUE(w.step().ok());

        // Comparing only centre-of-mass velocities would give a closing speed of
        // exactly zero here, so no impulse at all and nothing would change. The
        // spin must therefore have been altered by the contact.
        MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angular_velocity - 3.0) > 1e-6);
    }

    // --- A contact pair already separating must not be impulsed again. Doing so
    //     conserves momentum and looks plausible, so only checking the
    //     velocities directly catches it. ---
    {
        RigidBody2D a;
        a.position = Vec2{-0.4, 0.0};
        a.velocity = Vec2{-1.0, 0.0};   // already moving apart
        a.mass = 1.0;
        a.inertia = 1.0;
        a.radius = 0.5;
        RigidBody2D b;
        b.position = Vec2{0.4, 0.0};
        b.velocity = Vec2{1.0, 0.0};
        b.mass = 1.0;
        b.inertia = 1.0;
        b.radius = 0.5;

        // Overlapping by 0.2, so a contact is found, but they are separating.
        RigidWorld w{SimulationSettings{1e-9}, {a, b}, 1.0};
        MALLOY_CHECK_TRUE(w.step().ok());
        // Only the positional correction should have acted.
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(-1.0, 0.0), 1e-9);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(1.0, 0.0), 1e-9);
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 1e-9);
    }

    // --- Diagnostics skip immovable bodies: infinite mass times zero velocity
    //     would be NaN, and an immovable body carries nothing to report. ---
    {
        RigidBody2D wall;
        wall.mass = inf;
        wall.inertia = inf;
        wall.radius = 1.0;
        RigidBody2D movable = awkward_body();
        const std::vector<RigidBody2D> mixed = {wall, movable};

        MALLOY_CHECK_TRUE(malloy::math::is_finite(total_linear_momentum(mixed)));
        MALLOY_CHECK_TRUE(std::isfinite(total_angular_momentum(mixed)));
        MALLOY_CHECK_TRUE(std::isfinite(total_kinetic_energy(mixed)));
        // The wall contributes exactly nothing.
        const std::vector<RigidBody2D> alone = {movable};
        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(mixed),
                               total_linear_momentum(alone), 0.0);
        MALLOY_CHECK_NEAR(total_kinetic_energy(mixed), total_kinetic_energy(alone), 0.0);
    }

    // --- Restitution is validated. ---
    {
        RigidWorld too_bouncy{SimulationSettings{0.01}, {awkward_body()}, 1.5};
        MALLOY_CHECK_TRUE(too_bouncy.validate() == StepStatus::InvalidSettings);
        RigidWorld negative{SimulationSettings{0.01}, {awkward_body()}, -0.1};
        MALLOY_CHECK_TRUE(negative.validate() == StepStatus::InvalidSettings);
        RigidWorld ok{SimulationSettings{0.01}, {awkward_body()}, 0.5};
        MALLOY_CHECK_TRUE(ok.validate() == StepStatus::Ok);
        MALLOY_CHECK_NEAR(ok.restitution(), 0.5, eps);
    }

    std::cout << "malloy_rigid_tests passed\n";
    return 0;
}
