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
using malloy::collide::Halfplane;
using malloy::rigid::RigidBody2D;
using malloy::rigid::RigidSettings;
using malloy::rigid::RigidWorld;
using malloy::rigid::shift_inertia;
using malloy::rigid::to_local;
using malloy::rigid::to_world;
using malloy::rigid::total_angular_momentum;
using malloy::rigid::total_potential_energy;
using malloy::rigid::total_energy;
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
        // NOT static: infinite inertia alone means it cannot spin, and it can
        // still be pushed. This assertion used to read is_static() == true,
        // which is the contract that changed and why it changed.
        MALLOY_CHECK_FALSE(inf_inertia.is_static());
        MALLOY_CHECK_TRUE(inf_inertia.has_infinite_inertia());
        MALLOY_CHECK_FALSE(inf_inertia.has_infinite_mass());

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

        RigidWorld w{SimulationSettings{0.5}, {a, b}, RigidSettings{1.0}};
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

        RigidWorld w{SimulationSettings{0.5}, {ball, wall}, RigidSettings{1.0}};
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

        RigidWorld a{SimulationSettings{0.5}, {up, wall}, RigidSettings{1.0}};
        RigidWorld b{SimulationSettings{0.5}, {down, wall}, RigidSettings{1.0}};
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

        RigidWorld w{SimulationSettings{0.001}, {a, b}, RigidSettings{1.0}};
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

            RigidWorld w{SimulationSettings{0.002}, start, RigidSettings{restitution}};
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
        // angular momentum, just far less.
        //
        // The measured figures are +5.02583e-04 for the long run above and
        // -1.20527e-09 here, a factor of 4.2e5.
        //
        // Do not read that ratio as the penetration ratio alone. The
        // perturbation is c x (v_b - v_a), so it scales with the RELATIVE
        // VELOCITY as well as with the correction, and these two setups differ
        // in both. Proportionality to penetration is established by holding a
        // setup fixed and varying dt, not by comparing these two numbers.
        //
        // The correction is the cause, and the impulse is exact.
        RigidWorld w{SimulationSettings{1e-9}, start, RigidSettings{1.0}};
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

        RigidWorld elastic{SimulationSettings{0.002}, start, RigidSettings{1.0}};
        RigidWorld lossy{SimulationSettings{0.002}, start, RigidSettings{0.3}};
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
        RigidWorld w{SimulationSettings{1e-9}, {spinner, wall}, RigidSettings{1.0}};
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
        RigidWorld w{SimulationSettings{1e-9}, {a, b}, RigidSettings{1.0}};
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
        RigidWorld too_bouncy{SimulationSettings{0.01}, {awkward_body()}, RigidSettings{1.5}};
        MALLOY_CHECK_TRUE(too_bouncy.validate() == StepStatus::InvalidSettings);
        RigidWorld negative{SimulationSettings{0.01}, {awkward_body()}, RigidSettings{-0.1}};
        MALLOY_CHECK_TRUE(negative.validate() == StepStatus::InvalidSettings);
        RigidWorld ok{SimulationSettings{0.01}, {awkward_body()}, RigidSettings{0.5}};
        MALLOY_CHECK_TRUE(ok.validate() == StepStatus::Ok);
        MALLOY_CHECK_NEAR(ok.settings().restitution, 0.5, eps);
    }

    // --- M15: settings validation. ---
    {
        MALLOY_CHECK_TRUE(RigidSettings{}.is_valid());
        MALLOY_CHECK_TRUE((RigidSettings{0.0, Vec2{0.0, -9.81}}.is_valid()));
        RigidSettings bouncy;
        bouncy.restitution = 1.5;
        MALLOY_CHECK_FALSE(bouncy.is_valid());
        RigidSettings bad_gravity;
        bad_gravity.gravity = Vec2{0.0, nan};
        MALLOY_CHECK_FALSE(bad_gravity.is_valid());
        RigidSettings inf_gravity;
        inf_gravity.gravity = Vec2{inf, 0.0};
        MALLOY_CHECK_FALSE(inf_gravity.is_valid());
    }

    // --- Zero gravity reproduces the pre-M15 behaviour EXACTLY, which is the
    //     guard for every scenario written before the setting existed. ---
    {
        const std::vector<RigidBody2D> start = {awkward_body()};
        RigidWorld unset{SimulationSettings{0.002}, start, RigidSettings{0.9}};
        RigidSettings explicit_zero;
        explicit_zero.restitution = 0.9;
        explicit_zero.gravity = Vec2{0.0, 0.0};
        RigidWorld zero{SimulationSettings{0.002}, start, explicit_zero};
        for (int i = 0; i < 300; ++i)
        {
            MALLOY_CHECK_TRUE(unset.step().ok());
            MALLOY_CHECK_TRUE(zero.step().ok());
        }
        MALLOY_CHECK_VEC2_NEAR(unset.bodies()[0].position, zero.bodies()[0].position, 0.0);
        MALLOY_CHECK_NEAR(unset.bodies()[0].angle, zero.bodies()[0].angle, 0.0);
        MALLOY_CHECK_NEAR(total_potential_energy(unset.bodies(), Vec2{}), 0.0, eps);
    }

    // --- Gravity is applied BEFORE the position update, which is what keeps
    //     this semi-implicit Euler. With g = (0,-10) and dt = 0.5 the velocity
    //     becomes -5 and the centre of mass moves by -2.5, not 0 (explicit
    //     Euler) and not -1.25 (the exact solution). ---
    {
        RigidBody2D b;
        b.mass = 2.0;
        b.inertia = 1.0;
        RigidSettings g;
        g.gravity = Vec2{0.0, -10.0};
        RigidWorld w{SimulationSettings{0.5}, {b}, g};
        const Vec2 com0 = center_of_mass(b);
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(0.0, -5.0), eps);
        MALLOY_CHECK_VEC2_NEAR(center_of_mass(w.bodies()[0]), com0 + Vec2(0.0, -2.5), eps);
        // Gravity acts through the centre of mass, so it alone cannot spin a
        // body however long it falls.
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 0.0);
    }

    // --- Gravity accelerates every body equally regardless of mass, and does
    //     NOT move a static one: it is an acceleration, so without an explicit
    //     skip an immovable wall would start falling. ---
    {
        RigidBody2D light;
        light.mass = 0.5;
        light.inertia = 1.0;
        RigidBody2D heavy;
        heavy.position = Vec2{50.0, 0.0};
        heavy.mass = 500.0;
        heavy.inertia = 1.0;
        RigidBody2D wall;
        wall.position = Vec2{-50.0, 0.0};
        wall.mass = inf;
        wall.inertia = inf;

        RigidSettings g;
        g.gravity = Vec2{3.0, -4.0};
        RigidWorld w{SimulationSettings{0.25}, {light, heavy, wall}, g};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(0.75, -1.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(0.75, -1.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[2].velocity, Vec2(0.0, 0.0), 0.0); // the wall
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[2].position, Vec2(-50.0, 0.0), 0.0);
    }

    // --- INVARIANT: in free fall the energy loss per step is exactly
    //     (1/2)(sum m)|g|^2 dt^2, the same figure derived for the particle
    //     domain in M12, because it is the same integrator on the same field. ---
    {
        RigidBody2D b;
        b.position = Vec2{0.0, 10.0};
        b.velocity = Vec2{1.0, 0.0};
        b.mass = 2.0;
        b.inertia = 1.0;
        RigidSettings g;
        g.gravity = Vec2{0.0, -3.0};
        const Real dt = 0.01;
        const int steps = 300;

        const std::vector<RigidBody2D> start = {b};
        const Real e0 = total_energy(start, g.gravity);
        RigidWorld w{SimulationSettings{dt}, start, g};
        for (int i = 0; i < steps; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const Real per_step = 0.5 * 2.0 * 9.0 * dt * dt;
        MALLOY_CHECK_NEAR(total_energy(w.bodies(), g.gravity),
                          e0 - per_step * static_cast<Real>(steps), 1e-10);
        MALLOY_CHECK_TRUE(per_step * static_cast<Real>(steps) > 1e-6);
    }

    // --- Potential energy: higher is more with gravity pointing down, measured
    //     at the CENTRE OF MASS rather than the body origin, and statics are
    //     excluded because infinite mass times a position is not a number. ---
    {
        const Vec2 g{0.0, -10.0};
        RigidBody2D high;
        high.position = Vec2{0.0, 3.0};
        high.mass = 2.0;
        high.inertia = 1.0;
        RigidBody2D low = high;
        low.position = Vec2{0.0, 1.0};
        MALLOY_CHECK_NEAR(total_potential_energy({high}, g), 60.0, eps);
        MALLOY_CHECK_NEAR(total_potential_energy({low}, g), 20.0, eps);

        // The offset centre of mass is what counts, not the origin.
        RigidBody2D offset = high;
        offset.local_center_of_mass = Vec2{0.0, 1.0}; // centre of mass at y = 4
        MALLOY_CHECK_NEAR(total_potential_energy({offset}, g), 80.0, eps);

        RigidBody2D wall;
        wall.mass = inf;
        wall.inertia = inf;
        MALLOY_CHECK_NEAR(total_potential_energy({wall}, g), 0.0, 0.0);
    }

    // --- A body dropped onto an immovable floor lands and stays there. ---
    {
        RigidBody2D ball;
        ball.position = Vec2{0.0, 6.0};
        ball.mass = 1.0;
        ball.inertia = 0.5;
        ball.radius = 0.5;
        RigidBody2D floor;
        floor.position = Vec2{0.0, -5.0};
        floor.mass = inf;
        floor.inertia = inf;
        floor.radius = 5.0;

        RigidSettings g;
        g.restitution = 0.4;
        g.gravity = Vec2{0.0, -9.81};
        RigidWorld w{SimulationSettings{0.001}, {ball, floor}, g};
        for (int i = 0; i < 20000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const RigidBody2D& b = w.bodies()[0];
        MALLOY_CHECK_TRUE(malloy::math::is_finite(b.position));
        // Resting on the floor disc: centres about 5.5 apart, and not sunk in.
        const Real gap = malloy::math::distance(b.position, Vec2{0.0, -5.0});
        MALLOY_CHECK_TRUE(gap > 5.0);
        MALLOY_CHECK_TRUE(gap < 5.7);
        MALLOY_CHECK_TRUE(std::abs(b.velocity.y) < 0.5);
    }

    // --- EMERGENT: a body whose centre of mass is offset from its disc centre,
    //     resting on a floor, ROCKS. Gravity acts through the centre of mass
    //     while the contact acts at the rim, so the two are not collinear and
    //     the pair generates torque. Nothing in M15 implements this: it falls
    //     out of M14's contact arm once there is a field to act against. ---
    {
        RigidBody2D wobbler;
        wobbler.position = Vec2{0.0, 1.0};
        wobbler.angle = 0.0;
        wobbler.mass = 1.0;
        wobbler.inertia = 0.05; // small, so the torque shows quickly
        wobbler.radius = 0.5;
        wobbler.local_center_of_mass = Vec2{0.35, 0.0}; // offset sideways

        RigidBody2D floor;
        floor.position = Vec2{0.0, -5.0};
        floor.mass = inf;
        floor.inertia = inf;
        floor.radius = 5.0;

        RigidSettings g;
        g.restitution = 0.1;
        g.gravity = Vec2{0.0, -9.81};
        RigidWorld w{SimulationSettings{0.0005}, {wobbler, floor}, g};
        for (int i = 0; i < 8000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        // It turned. With the centre of mass at the disc centre this would stay
        // at exactly zero however long it sat there.
        MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angle) > 1e-3);
        MALLOY_CHECK_TRUE(malloy::math::is_finite(w.bodies()[0].position));
    }

    // --- The claim scenarios/dropped_bodies.scn is built on: two bodies that
    //     differ only in mass fall identically, land on the same step, and
    //     bounce identically.
    //
    //     The FALL matches because gravity is an acceleration and never divides
    //     by a mass. The BOUNCE matches for a separate reason: the contact
    //     impulse divides by (1/m + (r x n)^2 / I), so the resulting change in
    //     velocity is -(1+e)(v.n) / (1 + (m/I)(r x n)^2). It depends on the
    //     RATIO m/I, not on m, and these two share it (0.5/0.08 = 5.0/0.80).
    //     Change one inertia and only the free fall would still agree.
    //
    //     The pair is translated rather than mirrored, so identical means
    //     identical here and not sign flipped. ---
    {
        const auto falling = [](Real mass, Real inertia, Real x) {
            RigidBody2D b;
            b.mass = mass;
            b.inertia = inertia;
            b.radius = 0.40;
            b.local_center_of_mass = Vec2{0.15, 0.0};
            b.position = Vec2{x, 3.0};
            return b;
        };
        const auto floor_at = [inf](Real x) {
            RigidBody2D f;
            f.mass = inf;
            f.inertia = inf;
            f.radius = 0.60;
            f.position = Vec2{x, -1.20};
            return f;
        };

        RigidSettings g;
        g.restitution = 0.45;
        g.gravity = Vec2{0.0, -9.81};
        RigidWorld w{SimulationSettings{0.001},
                     {falling(0.5, 0.08, -1.60), falling(5.0, 0.80, 1.40),
                      floor_at(-1.50), floor_at(1.50)},
                     g};
        for (int i = 0; i < 1200; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }

        const RigidBody2D& light = w.bodies()[0];
        const RigidBody2D& heavy = w.bodies()[1];
        MALLOY_CHECK_NEAR(light.position.y, heavy.position.y, 1e-12);
        MALLOY_CHECK_NEAR(light.position.x + 3.0, heavy.position.x, 1e-12);
        MALLOY_CHECK_NEAR(light.angle, heavy.angle, 1e-12);
        MALLOY_CHECK_NEAR(light.velocity.y, heavy.velocity.y, 1e-12);

        // And the run really did include a landing, so the agreement above is
        // about a contact and not merely about two bodies still in free fall.
        MALLOY_CHECK_TRUE(light.position.y < 0.0);
        MALLOY_CHECK_TRUE(std::abs(light.angle) > 1e-6);
        MALLOY_CHECK_TRUE(light.velocity.y > 0.0); // rebounding

        // The mass genuinely differs, so this is not two copies of one body.
        MALLOY_CHECK_NEAR(heavy.mass, 10.0 * light.mass, eps);
    }

    // --- Break the shared ratio and the bounce stops agreeing, which is what
    //     shows the agreement above is the ratio and not the mass. ---
    {
        const auto falling = [](Real mass, Real inertia, Real x) {
            RigidBody2D b;
            b.mass = mass;
            b.inertia = inertia;
            b.radius = 0.40;
            b.local_center_of_mass = Vec2{0.15, 0.0};
            b.position = Vec2{x, 3.0};
            return b;
        };
        const auto floor_at = [inf](Real x) {
            RigidBody2D f;
            f.mass = inf;
            f.inertia = inf;
            f.radius = 0.60;
            f.position = Vec2{x, -1.20};
            return f;
        };

        RigidSettings g;
        g.restitution = 0.45;
        g.gravity = Vec2{0.0, -9.81};
        RigidWorld w{SimulationSettings{0.001},
                     {falling(0.5, 0.08, -1.60), falling(5.0, 0.20, 1.40),
                      floor_at(-1.50), floor_at(1.50)},
                     g};
        for (int i = 0; i < 1200; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angle - w.bodies()[1].angle) > 1e-6);
    }

    // --- Infinity is per quantity, so a half-infinite body is immovable in
    //     exactly one sense and fully real in the other.
    //
    //     This was a defect, not a gap. is_static() used to be an OR, so a body
    //     with finite mass and infinite inertia called itself static: contacts
    //     moved it, while gravity and every diagnostic skipped it. The momentum
    //     handed to it vanished from the report, which falsified the exact
    //     linear-momentum conservation the header claims. ---
    {
        RigidBody2D slider; // can translate, cannot spin
        slider.mass = 2.0;
        slider.inertia = inf;
        slider.velocity = Vec2{3.0, -1.0};
        slider.position = Vec2{1.0, 4.0};
        MALLOY_CHECK_TRUE(slider.is_valid());
        MALLOY_CHECK_FALSE(slider.is_static());
        MALLOY_CHECK_NEAR(slider.inverse_mass(), 0.5, eps);
        MALLOY_CHECK_NEAR(slider.inverse_inertia(), 0.0, 0.0);

        RigidBody2D flywheel; // can spin, cannot translate
        flywheel.mass = inf;
        flywheel.inertia = 0.25;
        flywheel.angular_velocity = 4.0;
        MALLOY_CHECK_TRUE(flywheel.is_valid());
        MALLOY_CHECK_FALSE(flywheel.is_static());
        MALLOY_CHECK_NEAR(flywheel.inverse_mass(), 0.0, 0.0);
        MALLOY_CHECK_NEAR(flywheel.inverse_inertia(), 4.0, eps);

        RigidBody2D wall; // neither
        wall.mass = inf;
        wall.inertia = inf;
        MALLOY_CHECK_TRUE(wall.is_static());

        // The slider carries real linear momentum, which the OR version
        // discarded. 2 * (3, -1).
        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum({slider}), Vec2(6.0, -2.0), eps);
        // And real translational kinetic energy: 0.5 * 2 * (9 + 1) = 10.
        MALLOY_CHECK_NEAR(total_kinetic_energy({slider}), 10.0, eps);
        // Its spin term is dropped rather than becoming inf * 0 = NaN.
        MALLOY_CHECK_TRUE(malloy::math::is_finite(total_kinetic_energy({slider})));
        // Orbital angular momentum survives, spin does not: m * cross(r, v)
        // with r = (1, 4) and v = (3, -1) gives 2 * (1*-1 - 4*3) = -26.
        MALLOY_CHECK_NEAR(total_angular_momentum({slider}), -26.0, eps);
        // And it has gravitational potential energy, being of finite mass.
        MALLOY_CHECK_NEAR(total_potential_energy({slider}, Vec2{0.0, -10.0}), 80.0, eps);

        // The flywheel is the mirror: real spin energy, no translation terms.
        MALLOY_CHECK_NEAR(total_kinetic_energy({flywheel}), 2.0, eps); // 0.5*0.25*16
        MALLOY_CHECK_NEAR(total_angular_momentum({flywheel}), 1.0, eps); // 0.25*4
        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum({flywheel}), Vec2(0.0, 0.0), 0.0);
        MALLOY_CHECK_NEAR(total_potential_energy({flywheel}, Vec2{0.0, -10.0}), 0.0, 0.0);

        // The wall contributes nothing at all, and nothing is NaN.
        MALLOY_CHECK_NEAR(total_kinetic_energy({wall}), 0.0, 0.0);
        MALLOY_CHECK_NEAR(total_angular_momentum({wall}), 0.0, 0.0);
    }

    // --- Gravity asks about the MASS, not about is_static(). A body that
    //     merely cannot spin still falls; one of infinite mass does not. ---
    {
        RigidBody2D slider;
        slider.mass = 1.0;
        slider.inertia = inf;
        RigidBody2D flywheel;
        flywheel.position = Vec2{20.0, 0.0};
        flywheel.mass = inf;
        flywheel.inertia = 1.0;

        RigidSettings g;
        g.gravity = Vec2{0.0, -8.0};
        RigidWorld w{SimulationSettings{0.25}, {slider, flywheel}, g};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(0.0, -2.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(0.0, 0.0), 0.0);
        // The slider fell but did not turn: it cannot.
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 0.0);
    }

    // --- A contact against a half-infinite body conserves linear momentum
    //     EXACTLY, which is the claim the OR version falsified. The slider is
    //     pushed and its momentum change is visible in the total. ---
    {
        RigidBody2D ball;
        ball.mass = 1.0;
        ball.inertia = 0.5;
        ball.radius = 0.5;
        ball.position = Vec2{-0.6, 0.0};
        ball.velocity = Vec2{4.0, 0.0};

        RigidBody2D slider;
        slider.mass = 3.0;
        slider.inertia = inf; // cannot spin, but is pushed
        slider.radius = 0.5;
        slider.position = Vec2{0.3, 0.0};

        const std::vector<RigidBody2D> start = {ball, slider};
        const Vec2 before = total_linear_momentum(start);
        MALLOY_CHECK_VEC2_NEAR(before, Vec2(4.0, 0.0), eps);

        RigidWorld w{SimulationSettings{0.001}, start, RigidSettings{1.0}};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(w.bodies()), before, 1e-12);
        // The slider really did move, so this is not conservation by inaction.
        MALLOY_CHECK_TRUE(w.bodies()[1].velocity.x > 0.1);
        // And it did not start spinning, however off-centre the contact was.
        MALLOY_CHECK_NEAR(w.bodies()[1].angular_velocity, 0.0, 0.0);
    }

    // --- Two bodies of infinite mass and finite inertia are NOT static under
    //     the AND rule, so they now reach impulse code that the OR rule used to
    //     short-circuit. With the contact normal running through both centres
    //     of mass, every arm is parallel to it and the effective mass is
    //     exactly zero. The guard must return rather than divide by it. ---
    {
        RigidBody2D a;
        a.mass = inf;
        a.inertia = 1.0;
        a.radius = 0.5;
        a.position = Vec2{-0.4, 0.0};
        a.velocity = Vec2{1.0, 0.0}; // closing, so the separating exit is not taken
        RigidBody2D b = a;
        b.position = Vec2{0.4, 0.0};
        b.velocity = Vec2{-1.0, 0.0};
        MALLOY_CHECK_FALSE(a.is_static());

        RigidWorld w{SimulationSettings{0.001}, {a, b}, RigidSettings{1.0}};
        MALLOY_CHECK_TRUE(w.step().ok());
        // Nothing became NaN, and no impulse was applied, because there is no
        // finite effective mass for one to divide by.
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(1.0, 0.0), 0.0);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(-1.0, 0.0), 0.0);
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 0.0);
        MALLOY_CHECK_TRUE(malloy::math::is_finite(w.bodies()[0].position));
    }

    // --- An immovable body with a velocity is a KINEMATIC body: a moving
    //     platform. It is integrated like anything else, gravity does not
    //     accelerate it, no impulse can slow it, and it carries no momentum in
    //     the diagnostics, so it acts on the world without being part of it.
    //
    //     Nothing tested this in either direction, so adding a static skip to
    //     the pose loop used to pass the suite just as well. ---
    {
        RigidBody2D platform;
        platform.mass = inf;
        platform.inertia = inf;
        platform.radius = 0.5;
        platform.position = Vec2{0.0, 0.0};
        platform.velocity = Vec2{2.0, 0.0};

        RigidBody2D ball;
        ball.mass = 1.0;
        ball.inertia = 0.5;
        ball.radius = 0.5;
        ball.position = Vec2{0.9, 0.0};

        RigidSettings g;
        g.restitution = 1.0;
        g.gravity = Vec2{0.0, -9.81};
        RigidWorld w{SimulationSettings{0.01}, {platform, ball}, g};
        for (int i = 0; i < 50; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }

        // It moved exactly as its velocity says, undisturbed by gravity, by the
        // impact, or by anything else: 2.0 * 0.01 * 50.
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].position, Vec2(1.0, 0.0), 1e-12);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(2.0, 0.0), 0.0);
        // And it drove the ball ahead of it rather than passing through.
        MALLOY_CHECK_TRUE(w.bodies()[1].velocity.x > 0.5);
        MALLOY_CHECK_TRUE(w.bodies()[1].position.x > 0.9);
        // The platform contributes no momentum, being infinitely massive.
        MALLOY_CHECK_VEC2_NEAR(total_linear_momentum({w.bodies()[0]}), Vec2(0.0, 0.0),
                               0.0);
    }

    // --- M16: ground planes are validated with the settings, so a malformed
    //     one is refused before it can produce a NaN normal. ---
    {
        RigidSettings ok;
        ok.ground.push_back(Halfplane{Vec2{0.0, 1.0}, -2.0});
        ok.ground.push_back(Halfplane{Vec2{0.6, 0.8}, 3.5});
        MALLOY_CHECK_TRUE(ok.is_valid());

        RigidSettings bad;
        bad.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});
        bad.ground.push_back(Halfplane{Vec2{1.0, 1.0}, 0.0}); // not a unit normal
        MALLOY_CHECK_FALSE(bad.is_valid());

        RigidBody2D b;
        RigidWorld w{SimulationSettings{0.01}, {b}, bad};
        MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidSettings);
        MALLOY_CHECK_FALSE(w.step().ok());
    }

    // --- No ground reproduces the pre-M16 trajectory exactly, so every world
    //     written before planes existed is unaffected. ---
    {
        const std::vector<RigidBody2D> start = {awkward_body()};
        RigidSettings without;
        without.restitution = 0.8;
        RigidSettings empty_ground = without;
        empty_ground.ground.clear();

        RigidWorld a{SimulationSettings{0.002}, start, without};
        RigidWorld b{SimulationSettings{0.002}, start, empty_ground};
        for (int i = 0; i < 300; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        MALLOY_CHECK_VEC2_NEAR(a.bodies()[0].position, b.bodies()[0].position, 0.0);
        MALLOY_CHECK_NEAR(a.bodies()[0].angle, b.bodies()[0].angle, 0.0);
    }

    // --- THE POINT OF THE MILESTONE. A body sliding along a flat floor keeps
    //     its horizontal velocity EXACTLY and picks up no spin at all, because
    //     the contact normal is the plane's own and never turns.
    //
    //     A floor built from discs cannot do this. Its normal points at
    //     whichever disc centre is nearest, so it swings by up to 14 degrees
    //     across the floor in dropped_bodies.scn, which steers the body and
    //     spins it. The tolerance here is zero, not small. ---
    {
        RigidBody2D slider;
        slider.mass = 1.5;
        slider.inertia = 0.2;
        slider.radius = 0.4;
        slider.position = Vec2{-3.0, 0.4}; // resting on the floor at y = 0
        slider.velocity = Vec2{2.0, 0.0};

        RigidSettings g;
        g.restitution = 0.0; // no bouncing, so it stays in contact
        g.gravity = Vec2{0.0, -9.81};
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.001}, {slider}, g};
        for (int i = 0; i < 3000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const RigidBody2D& now = w.bodies()[0];
        MALLOY_CHECK_NEAR(now.velocity.x, 2.0, 0.0);       // exactly, not nearly
        MALLOY_CHECK_NEAR(now.angular_velocity, 0.0, 0.0); // no spin from a flat floor
        MALLOY_CHECK_NEAR(now.position.x, -3.0 + 2.0 * 3.0, 1e-12);
        // It stayed on the surface rather than sinking through or hopping.
        MALLOY_CHECK_TRUE(now.position.y > 0.39);
        MALLOY_CHECK_TRUE(now.position.y < 0.41);
    }

    // --- On a slanted plane a frictionless body accelerates exactly down the
    //     slope at g sin(theta), and this is derivable in closed form for one
    //     step. The plane is a 3-4-5 slope, normal (0.6, 0.8), so cos(theta) is
    //     0.8 and sin(theta) is 0.6, and the body starts EXACTLY touching it.
    //
    //     With g = (0, -10) and dt = 0.01, gravity gives v = (0, -0.1). A
    //     restitution-0 contact removes the normal component and leaves the
    //     tangential one, so
    //         v' = v - (v . n) n = (0, -0.1) - (-0.08)(0.6, 0.8)
    //            = (0.048, -0.036)
    //     whose magnitude is 0.06 = g dt sin(theta), pointing down the slope
    //     along (0.8, -0.6). ---
    {
        RigidBody2D block;
        block.mass = 1.0;
        block.inertia = 0.5;
        block.radius = 0.5;
        // 0.6 * 1.5 + 0.8 * -0.5 = 0.5, exactly one radius clear of the plane.
        block.position = Vec2{1.5, -0.5};

        RigidSettings g;
        g.restitution = 0.0;
        g.gravity = Vec2{0.0, -10.0};
        g.ground.push_back(Halfplane{Vec2{0.6, 0.8}, 0.0});

        RigidWorld w{SimulationSettings{0.01}, {block}, g};
        MALLOY_CHECK_TRUE(w.step().ok());
        const RigidBody2D& now = w.bodies()[0];
        MALLOY_CHECK_VEC2_NEAR(now.velocity, Vec2(0.048, -0.036), 1e-15);
        MALLOY_CHECK_NEAR(malloy::math::length(now.velocity), 0.06, 1e-15);

        // Frictionless and centred, so nothing turns it. The bound is rounding
        // level rather than exactly zero, and the reason is worth stating: the
        // arm is parallel to the impulse, so the true cross product is zero,
        // but cross() computes arm.x * impulse.y - arm.y * impulse.x, and on a
        // slanted normal those two products multiply 0.6 and 0.8 in opposite
        // orders and round differently. The residue is about 7e-18. On the
        // axis-aligned floor above, the same expression is exactly zero,
        // because one factor in each product is a literal zero.
        MALLOY_CHECK_NEAR(now.angular_velocity, 0.0, 1e-15);
    }

    // --- A perfectly elastic head-on bounce reverses the velocity exactly and
    //     produces no spin, which is the control case for the slope above. ---
    {
        RigidBody2D ball;
        ball.mass = 2.0;
        ball.inertia = 0.3;
        ball.radius = 0.5;
        ball.position = Vec2{4.0, 0.5}; // exactly touching
        ball.velocity = Vec2{0.0, -3.0};

        RigidSettings g;
        g.restitution = 1.0;
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.001}, {ball}, g};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(0.0, 3.0), 1e-15);
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 0.0);
    }

    // --- An off-centre contact against a plane DOES generate spin, so the zero
    //     above is a property of the geometry and not of planes being inert. ---
    {
        RigidBody2D wobbler;
        wobbler.mass = 1.0;
        wobbler.inertia = 0.05;
        wobbler.radius = 0.5;
        wobbler.local_center_of_mass = Vec2{0.3, 0.0}; // offset from the disc
        wobbler.position = Vec2{0.0, 0.45};
        wobbler.velocity = Vec2{0.0, -2.0};

        RigidSettings g;
        g.restitution = 0.6;
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.001}, {wobbler}, g};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angular_velocity) > 1e-6);

        // Mirroring the offset mirrors the spin, so both signs occur.
        RigidBody2D mirrored = wobbler;
        mirrored.local_center_of_mass = Vec2{-0.3, 0.0};
        RigidWorld m{SimulationSettings{0.001}, {mirrored}, g};
        MALLOY_CHECK_TRUE(m.step().ok());
        MALLOY_CHECK_NEAR(m.bodies()[0].angular_velocity,
                          -w.bodies()[0].angular_velocity, 1e-15);
    }

    // --- A body of zero radius takes part in no contacts, planes included, so
    //     it falls straight through the floor. ---
    {
        RigidBody2D ghost;
        ghost.mass = 1.0;
        ghost.inertia = 1.0;
        ghost.radius = 0.0;
        ghost.position = Vec2{0.0, 1.0};

        RigidSettings g;
        g.gravity = Vec2{0.0, -10.0};
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.01}, {ghost}, g};
        for (int i = 0; i < 200; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_TRUE(w.bodies()[0].position.y < -1.0);
    }

    // --- Two planes make a corner, and a body settles into it rather than
    //     being fought over by them. Also exercises more than one plane, since
    //     a loop that resolved only the first would pass every case above. ---
    {
        RigidBody2D ball;
        ball.mass = 1.0;
        ball.inertia = 0.4;
        ball.radius = 0.5;
        ball.position = Vec2{2.0, 3.0};

        RigidSettings g;
        g.restitution = 0.2;
        g.gravity = Vec2{0.0, -9.81};
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});   // floor at y = 0
        g.ground.push_back(Halfplane{Vec2{-1.0, 0.0}, -1.0}); // wall at x = 1

        RigidWorld w{SimulationSettings{0.001}, {ball}, g};
        for (int i = 0; i < 6000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const RigidBody2D& now = w.bodies()[0];
        MALLOY_CHECK_TRUE(malloy::math::is_finite(now.position));
        MALLOY_CHECK_TRUE(now.position.y > 0.49);  // resting on the floor
        MALLOY_CHECK_TRUE(now.position.y < 0.55);
        MALLOY_CHECK_TRUE(now.position.x < 0.51);  // and left of the wall
        MALLOY_CHECK_TRUE(std::abs(now.velocity.y) < 0.2);
    }

    // --- M17: friction validation. Negative and non-finite are refused; a
    //     coefficient above 1 is NOT, because rubber on rubber really does
    //     exceed it and clamping would silently change the caller's model. ---
    {
        RigidSettings ok;
        ok.friction = 1.7;
        MALLOY_CHECK_TRUE(ok.is_valid());
        MALLOY_CHECK_TRUE(RigidSettings{}.is_valid()); // default 0, frictionless

        RigidSettings negative;
        negative.friction = -0.1;
        MALLOY_CHECK_FALSE(negative.is_valid());
        RigidSettings not_a_number;
        not_a_number.friction = nan;
        MALLOY_CHECK_FALSE(not_a_number.is_valid());
        RigidSettings unbounded;
        unbounded.friction = inf;
        MALLOY_CHECK_FALSE(unbounded.is_valid());
    }

    // --- Zero friction reproduces the pre-M17 trajectory exactly, so every
    //     world written before this existed is unaffected. ---
    {
        RigidBody2D slider;
        slider.mass = 1.5;
        slider.inertia = 0.2;
        slider.radius = 0.4;
        slider.position = Vec2{-3.0, 0.4};
        slider.velocity = Vec2{2.0, 0.0};

        RigidSettings without;
        without.restitution = 0.0;
        without.gravity = Vec2{0.0, -9.81};
        without.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});
        RigidSettings explicit_zero = without;
        explicit_zero.friction = 0.0;

        RigidWorld a{SimulationSettings{0.001}, {slider}, without};
        RigidWorld b{SimulationSettings{0.001}, {slider}, explicit_zero};
        for (int i = 0; i < 2000; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        MALLOY_CHECK_VEC2_NEAR(a.bodies()[0].velocity, b.bodies()[0].velocity, 0.0);
        // And it really did keep sliding, since nothing was there to stop it.
        MALLOY_CHECK_NEAR(a.bodies()[0].velocity.x, 2.0, 0.0);
    }

    // --- THE MILESTONE'S INVARIANT: a uniform disc launched with no spin onto
    //     a flat floor settles into rolling at exactly two thirds of its launch
    //     speed, and that figure depends on NEITHER the friction coefficient
    //     NOR gravity. Both only decide how long it takes to get there.
    //
    //     The reason it is exact rather than approximate: a tangential impulse
    //     at the contact point changes v and omega together in a fixed ratio,
    //     so `I*omega - m*R*v` is unchanged by it however large or small the
    //     impulse is, clamped or not. Starting from omega = 0 that constant is
    //     -m*R*v0, and rolling means omega = -v/R, which gives
    //     -(3/2)*m*R*v_roll, so v_roll = (2/3)*v0. The discretisation cannot
    //     move it, which is why this is not a tolerance-tuned test.
    //
    //     I = m*R^2/2 for a uniform disc is not asserted here by hand: it comes
    //     out of collide::second_moment_of_area and rigid::mass_properties,
    //     which have their own tests against closed-form values. ---
    {
        struct Run
        {
            Real friction;
            Real gravity;
        };
        const Run runs[] = {{0.30, 9.81}, {0.80, 9.81}, {0.15, 4.00}};
        const Real v0 = 3.0;
        const Real radius = 0.5;
        const Real mass = 2.0;

        for (const Run& run : runs)
        {
            RigidBody2D disc;
            disc.mass = mass;
            disc.inertia = mass * radius * radius / 2.0; // uniform disc
            disc.radius = radius;
            disc.position = Vec2{-8.0, radius}; // exactly touching the floor
            disc.velocity = Vec2{v0, 0.0};
            disc.angular_velocity = 0.0;

            RigidSettings g;
            g.restitution = 0.0;
            g.friction = run.friction;
            g.gravity = Vec2{0.0, -run.gravity};
            g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

            RigidWorld w{SimulationSettings{0.0005}, {disc}, g};
            for (int i = 0; i < 6000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }

            const RigidBody2D& now = w.bodies()[0];
            MALLOY_CHECK_NEAR(now.velocity.x, 2.0 * v0 / 3.0, 1e-5);

            // Rolling without slipping: the contact point is not moving. The
            // bound is not zero, and the residual is derivable rather than
            // tuned. In steady contact the body sinks by g*dt^2/(1 + e) per
            // step before the correction, which is g*dt^2 here because e is 0,
            // and the contact POINT is reported midway through that overlap.
            // The lever arm is therefore R - depth/2 rather than R, so rolling
            // is established about a point just inside the surface and the
            // contact creeps BACKWARD by v * (depth/2) / R.
            //
            // For the first run that is
            // -2.0 * (9.81 * 0.0005^2 / 2) / 0.5 = -4.905e-6, and the measured
            // residual is -4.905e-6. The sign matters and was once recorded
            // backwards: the assertion below is a magnitude bound, so it
            // cannot catch a sign error on its own.
            //
            // The second run has the SAME residual, not a smaller one, because
            // the penetration is set by gravity, dt and restitution and never
            // by the friction coefficient. Only the third run, at g = 4.0,
            // is smaller, at -2.0e-6.
            // Derived from this run's own gravity rather than hard coded, and
            // SIGNED. A magnitude bound cannot tell R - depth/2 from
            // R + depth/2, which is exactly how the sign came to be recorded
            // backwards in the first place: both give the same 4.905e-6.
            const Real depth = run.gravity * 0.0005 * 0.0005;
            const Real expected_creep = -now.velocity.x * (depth / 2.0) / radius;
            const Vec2 contact{now.position.x, 0.0};
            MALLOY_CHECK_NEAR(velocity_at(now, contact).x, expected_creep, 1e-9);
            // Which for this geometry means omega = -v/R, to the same bound
            // divided by R.
            MALLOY_CHECK_NEAR(now.angular_velocity, -now.velocity.x / radius, 2e-5);
        }
    }

    // --- And exactly one third of the launch kinetic energy is gone, which is
    //     the same statement in energy terms: (1/2)m v0^2 becomes
    //     (1/2)m v^2 + (1/2)I omega^2 = (1/3) m v0^2. ---
    {
        const Real v0 = 3.0;
        const Real radius = 0.5;
        const Real mass = 2.0;

        RigidBody2D disc;
        disc.mass = mass;
        disc.inertia = mass * radius * radius / 2.0;
        disc.radius = radius;
        disc.position = Vec2{-8.0, radius};
        disc.velocity = Vec2{v0, 0.0};

        RigidSettings g;
        g.restitution = 0.0;
        g.friction = 0.4;
        g.gravity = Vec2{0.0, -9.81};
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.0005}, {disc}, g};
        for (int i = 0; i < 6000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const RigidBody2D& now = w.bodies()[0];
        const Real rolling_energy =
            0.5 * mass * now.velocity.x * now.velocity.x +
            0.5 * disc.inertia * now.angular_velocity * now.angular_velocity;
        MALLOY_CHECK_NEAR(rolling_energy, mass * v0 * v0 / 3.0, 1e-3);
    }

    // --- Backspin reversal. A disc launched forward while spinning the wrong
    //     way walks forward, stops, and comes BACK. Nothing in the project
    //     could do this before: it needs a tangential impulse to convert spin
    //     into translation. ---
    {
        RigidBody2D disc;
        disc.mass = 1.0;
        disc.inertia = 0.125; // uniform disc, R = 0.5
        disc.radius = 0.5;
        disc.position = Vec2{0.0, 0.5};
        disc.velocity = Vec2{2.0, 0.0};
        disc.angular_velocity = 14.0; // backspin: same sign as forward motion

        RigidSettings g;
        g.restitution = 0.0;
        g.friction = 0.5;
        g.gravity = Vec2{0.0, -9.81};
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.0005}, {disc}, g};
        Real furthest = disc.position.x;
        for (int i = 0; i < 6000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
            furthest = std::max(furthest, w.bodies()[0].position.x);
        }
        // It went out, and then came back past where it turned around.
        MALLOY_CHECK_TRUE(furthest > 0.1);
        MALLOY_CHECK_TRUE(w.bodies()[0].velocity.x < 0.0);
        MALLOY_CHECK_TRUE(w.bodies()[0].position.x < furthest - 0.05);
    }

    // --- Spin-down: a body dropped with spin and no translation loses the
    //     spin to the floor, and gains translation from it. ---
    {
        RigidBody2D top;
        top.mass = 1.0;
        top.inertia = 0.125;
        top.radius = 0.5;
        top.position = Vec2{0.0, 0.5};
        top.angular_velocity = -20.0;

        RigidSettings g;
        g.restitution = 0.0;
        g.friction = 0.4;
        g.gravity = Vec2{0.0, -9.81};
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.0005}, {top}, g};
        for (int i = 0; i < 4000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angular_velocity) < 20.0);
        MALLOY_CHECK_TRUE(w.bodies()[0].velocity.x > 0.1); // spin became motion
    }

    // --- Static holding on a slope. A body that CANNOT roll, because its
    //     inertia is infinite, stays put when the friction coefficient reaches
    //     the tangent of the slope angle and slides when it does not. The
    //     3-4-5 plane has tan(theta) = 0.6/0.8 = 0.75.
    //
    //     Infinite inertia meaning "can translate but cannot spin" is the
    //     per-quantity rule; under the old OR-based is_static() this body would
    //     have called itself immovable and the test would be vacuous. ---
    {
        // An explicit return type, and step failures reported through a flag:
        // MALLOY_CHECK_TRUE expands to `return 1`, which cannot live in a
        // lambda that returns a Vec2.
        bool every_step_ok = true;
        const auto slide_test = [&](Real friction) -> Vec2 {
            RigidBody2D block;
            block.mass = 1.0;
            block.inertia = inf; // cannot roll, so this is pure sliding
            block.radius = 0.5;
            block.position = Vec2{1.5, -0.5}; // exactly touching the 3-4-5 plane

            RigidSettings g;
            g.restitution = 0.0;
            g.friction = friction;
            g.gravity = Vec2{0.0, -9.81};
            g.ground.push_back(Halfplane{Vec2{0.6, 0.8}, 0.0});

            RigidWorld w{SimulationSettings{0.001}, {block}, g};
            for (int i = 0; i < 2000; ++i)
            {
                if (!w.step().ok())
                {
                    every_step_ok = false;
                    break;
                }
            }
            return w.bodies()[0].velocity;
        };

        // Above the friction angle: it holds. Two seconds of gravity and it has
        // gone essentially nowhere.
        const Vec2 held = slide_test(0.9);
        MALLOY_CHECK_TRUE(every_step_ok);
        MALLOY_CHECK_TRUE(malloy::math::length(held) < 0.02);

        // Below it: it accelerates at exactly g(sin - mu cos). With mu = 0.3
        // that is 9.81 * (0.6 - 0.3 * 0.8) = 3.5316, over 2 seconds giving
        // 7.0632 down the slope along (0.8, -0.6).
        const Vec2 sliding = slide_test(0.3);
        MALLOY_CHECK_TRUE(every_step_ok);
        const Real expected = 9.81 * (0.6 - 0.3 * 0.8) * 2.0;
        MALLOY_CHECK_NEAR(malloy::math::length(sliding), expected, 1e-2);
        MALLOY_CHECK_NEAR(sliding.x, expected * 0.8, 1e-2);
        MALLOY_CHECK_NEAR(sliding.y, -expected * 0.6, 1e-2);
    }

    // --- Friction can only ever REMOVE kinetic energy. An unclamped or
    //     sign-flipped tangential impulse pumps it in, and this is the
    //     assertion that catches that at every step rather than at the end. ---
    {
        for (const Real friction : {0.1, 0.6, 2.0})
        {
            RigidBody2D a;
            a.mass = 1.0;
            a.inertia = 0.125;
            a.radius = 0.5;
            a.position = Vec2{-1.2, 0.5};
            a.velocity = Vec2{4.0, 0.0};
            a.angular_velocity = -3.0;

            RigidBody2D b;
            b.mass = 2.0;
            b.inertia = 0.4;
            b.radius = 0.5;
            b.position = Vec2{1.2, 0.5};
            b.velocity = Vec2{-1.0, 0.0};
            b.angular_velocity = 5.0;

            RigidSettings g;
            g.restitution = 0.9;
            g.friction = friction;
            g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

            RigidWorld w{SimulationSettings{0.0005}, {a, b}, g};
            Real previous = total_kinetic_energy(w.bodies());
            for (int i = 0; i < 4000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                const Real now = total_kinetic_energy(w.bodies());
                MALLOY_CHECK_TRUE(now <= previous + 1e-12);
                previous = now;
            }
        }
    }

    // --- Friction is gated on the NORMAL impulse, so a body touching nothing
    //     cannot be slowed by it. Without that gating a falling body would be
    //     dragged sideways by a friction term with no contact to justify it. ---
    {
        RigidBody2D falling;
        falling.mass = 1.0;
        falling.inertia = 0.125;
        falling.radius = 0.5;
        falling.position = Vec2{0.0, 40.0}; // far above the floor
        falling.velocity = Vec2{5.0, 0.0};

        RigidSettings g;
        g.friction = 2.0; // enormous, and still irrelevant in mid-air
        g.gravity = Vec2{0.0, -9.81};
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.001}, {falling}, g};
        for (int i = 0; i < 500; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_NEAR(w.bodies()[0].velocity.x, 5.0, 0.0); // exactly, still
        MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 0.0);
    }

    // --- The ordering in the step contract is observable, and this is what
    //     makes it so. A body whose centre of mass is offset from its disc
    //     drops STRAIGHT DOWN onto the floor with no spin, so the tangential
    //     relative velocity before the normal impulse is exactly zero. The
    //     normal impulse then acts at a point that is not below the centre of
    //     mass, which spins the body, and THAT gives the contact point a
    //     sideways velocity for friction to resist.
    //
    //     So friction computed from the pre-impulse velocity does nothing here
    //     and the body picks up no horizontal motion at all, while friction
    //     computed from what remains does. The mutation that reads the
    //     pre-impulse velocity escaped the entire suite until this existed. ---
    {
        RigidBody2D lopsided;
        lopsided.mass = 1.0;
        lopsided.inertia = 0.05;
        lopsided.radius = 0.5;
        lopsided.local_center_of_mass = Vec2{0.3, 0.0}; // offset from the disc
        lopsided.position = Vec2{0.0, 0.5};             // exactly touching
        lopsided.velocity = Vec2{0.0, -2.0};            // straight down, no spin
        lopsided.angular_velocity = 0.0;

        RigidSettings g;
        g.restitution = 0.0;
        g.friction = 0.5;
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.001}, {lopsided}, g};
        MALLOY_CHECK_TRUE(w.step().ok());

        // It arrived with no horizontal velocity and no spin whatsoever, and
        // both are now nonzero: the spin came from the off-centre normal
        // impulse, and the horizontal motion came from friction resisting the
        // slide that spin produced.
        MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angular_velocity) > 1e-6);
        MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].velocity.x) > 1e-6);
    }

    // --- A body genuinely at rest on the ground has EXACTLY zero tangential
    //     relative velocity, so there is no direction for friction to act
    //     along. Dividing by that zero to build a tangent produces NaN and
    //     poisons the body, which validate() would then reject.
    //
    //     Every other test here slides by some tiny residual and so never
    //     reaches the exact zero. This one does. ---
    {
        RigidBody2D resting;
        resting.mass = 1.0;
        resting.inertia = 0.125;
        resting.radius = 0.5;
        resting.position = Vec2{0.0, 0.5}; // exactly touching, motionless
        // velocity and angular_velocity are both exactly zero

        RigidSettings g;
        g.restitution = 0.0;
        g.friction = 0.7;
        g.gravity = Vec2{0.0, -9.81};
        g.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});

        RigidWorld w{SimulationSettings{0.001}, {resting}, g};
        for (int i = 0; i < 500; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const RigidBody2D& now = w.bodies()[0];
        MALLOY_CHECK_TRUE(malloy::math::is_finite(now.position));
        MALLOY_CHECK_TRUE(malloy::math::is_finite(now.velocity));
        // And it did not wander sideways, which is what an invented tangent
        // direction would have made it do.
        MALLOY_CHECK_NEAR(now.velocity.x, 0.0, 0.0);
        MALLOY_CHECK_NEAR(now.position.x, 0.0, 0.0);
        MALLOY_CHECK_NEAR(now.angular_velocity, 0.0, 0.0);
    }

    std::cout << "malloy_rigid_tests passed\n";
    return 0;
}
