#include <malloy/rigid/rigid.hpp>

#include <malloy/collide/collide.hpp>
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

    // --- Issue #16: a position or velocity whose SQUARE overflows is refused.
    //     Finite is not enough. Everything that squares a vector reaches
    //     infinity above about 1.34e154, so accepting state up to 1.8e308 left
    //     a window in which a world validated clean while every energy it
    //     reported was inf. ---
    {
        const Real too_big = 1.4e154; // finite, and its square is not
        RigidBody2D ok;
        ok.position = Vec2{1.0, 2.0};
        MALLOY_CHECK_TRUE(ok.is_valid());

        RigidBody2D far = ok;
        far.position = Vec2{too_big, 0.0};
        MALLOY_CHECK_FALSE(far.is_valid());

        RigidBody2D fast = ok;
        fast.velocity = Vec2{0.0, too_big};
        MALLOY_CHECK_FALSE(fast.is_valid());

        // The local centre of mass is squared too, through center_of_mass.
        RigidBody2D offset = ok;
        offset.local_center_of_mass = Vec2{too_big, 0.0};
        MALLOY_CHECK_FALSE(offset.is_valid());
    }

    // --- Issue #18: the angle accumulates by repeated addition, and this pins
    //     how far that drifts.
    //
    //     malloy_time computes elapsed time as tick_count * dt precisely so it
    //     "cannot drift the way repeated floating-point addition would"
    //     (fixed_step.hpp). The angle is the same pattern and does NOT get the
    //     same treatment, and it cannot: angular velocity is changed by
    //     contacts, so there is no constant increment to multiply.
    //
    //     What is left is to know the size of it. Summing N terms of magnitude
    //     theta accumulates at most u*theta*N/2 with u = 2^-53, which for the
    //     100000 steps below is 3.11e-09 on a final angle of 560 radians. The
    //     comparison is against the exactly computed N*omega*dt.
    //
    //     ADR 0007 takes an explicit position on this area, so the behaviour is
    //     documented rather than quietly changed. ---
    {
        RigidBody2D spinner;
        spinner.mass = 1.0;
        spinner.inertia = 1.0;
        spinner.radius = 0.0; // no contacts, so omega really is constant
        spinner.angular_velocity = 1.4;

        const Real dt = 0.004;
        const int steps = 100000;
        RigidWorld w{SimulationSettings{dt}, {spinner}, RigidSettings{}};
        for (int i = 0; i < steps; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }

        const Real exact = static_cast<Real>(steps) * spinner.angular_velocity * dt;
        const Real bound =
            (1.0 / 9007199254740992.0) * exact * static_cast<Real>(steps) / 2.0;
        MALLOY_CHECK_NEAR(w.bodies()[0].angle, exact, bound);

        // The angle is NOT wrapped, which is the documented contract, so the
        // comparison above is against 560 radians rather than a canonicalized
        // remainder.
        MALLOY_CHECK_TRUE(w.bodies()[0].angle > 500.0);
    }

    // --- M20: rotation in three dimensions.
    //
    //     This is the one part of 3D that cannot be reached from the 2D code by
    //     replacing types. There the inertia is a scalar and omega lies along a
    //     fixed axis, so `omega x (I omega)` is identically zero and a free body
    //     spins forever about one axis. Here that term is the whole subject.
    //
    //     The principal moments are (1, 2, 3) wherever three distinct ones are
    //     wanted: asymmetric on purpose (docs/05), so no two coefficients in
    //     Euler's equations coincide and the intermediate axis is y. ---
    {
        using malloy::math::Quat;
        using malloy::math::Vec3;
        using malloy::rigid::spin_angular_momentum;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::rotational_energy;

        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();
        const Vec3 asymmetric{1.0, 2.0, 3.0};

        // --- Validation, including the malformed input a scenario file can
        //     produce. ---
        {
            RigidBody3D body;
            MALLOY_CHECK_TRUE(body.is_valid()); // the defaults are a usable body

            for (const Real bad_mass : {0.0, -1.0, inf3, nan3})
            {
                RigidBody3D b;
                b.mass = bad_mass;
                MALLOY_CHECK_FALSE(b.is_valid());
            }

            // Every principal moment must be positive and finite, and each is
            // checked separately: a body with two good moments and one bad one
            // is not a body.
            for (const Real bad : {0.0, -1.0, inf3, nan3})
            {
                RigidBody3D bx;
                bx.inertia = Vec3{bad, 1.0, 1.0};
                MALLOY_CHECK_FALSE(bx.is_valid());
                RigidBody3D by;
                by.inertia = Vec3{1.0, bad, 1.0};
                MALLOY_CHECK_FALSE(by.is_valid());
                RigidBody3D bz;
                bz.inertia = Vec3{1.0, 1.0, bad};
                MALLOY_CHECK_FALSE(bz.is_valid());
            }

            // The orientation must be a UNIT quaternion, not merely a finite
            // one. Every rotation formula assumes it, and the all-zero
            // quaternion is the shape a half-written scenario file produces.
            RigidBody3D unnormalized;
            unnormalized.orientation = Quat{2.0, Vec3{}};
            MALLOY_CHECK_FALSE(unnormalized.is_valid());
            RigidBody3D zero_quat;
            zero_quat.orientation = Quat{0.0, Vec3{}};
            MALLOY_CHECK_FALSE(zero_quat.is_valid());
            RigidBody3D nan_quat;
            nan_quat.orientation = Quat{nan3, Vec3{}};
            MALLOY_CHECK_FALSE(nan_quat.is_valid());

            // A quaternion built from an axis and an angle is always unit, so
            // the honest way to write a tilted body stays valid.
            RigidBody3D tilted;
            tilted.orientation =
                malloy::math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 2.2);
            MALLOY_CHECK_TRUE(tilted.is_valid());

            // Squarable, not merely finite (docs/04): 1.4e154 is a perfectly
            // good double whose square is not, and a body carrying one reports
            // infinite energy while looking sound.
            RigidBody3D huge_position;
            huge_position.position = Vec3{0.0, 1.4e154, 0.0};
            MALLOY_CHECK_TRUE(malloy::math::is_finite(huge_position.position));
            MALLOY_CHECK_FALSE(huge_position.is_valid());
            RigidBody3D huge_velocity;
            huge_velocity.velocity = Vec3{1.4e154, 0.0, 0.0};
            MALLOY_CHECK_FALSE(huge_velocity.is_valid());
            RigidBody3D huge_spin;
            huge_spin.angular_velocity = Vec3{0.0, 0.0, 1.4e154};
            MALLOY_CHECK_FALSE(huge_spin.is_valid());
        }

        // --- A world with bad settings or a bad body reports it and changes
        //     nothing. ---
        {
            RigidBody3D start;
            start.inertia = asymmetric;
            start.angular_velocity = Vec3{0.7, -1.3, 0.4};
            start.velocity = Vec3{0.2, 0.3, -0.1};

            Rigid3DWorld bad_dt{SimulationSettings{0.0}, {start}};
            MALLOY_CHECK_TRUE(bad_dt.validate() == StepStatus::InvalidSettings);
            MALLOY_CHECK_TRUE(bad_dt.step().status == StepStatus::InvalidSettings);
            MALLOY_CHECK_EQ(bad_dt.tick_count(), std::uint64_t{0});
            MALLOY_CHECK_NEAR(bad_dt.elapsed_time(), 0.0, 0.0);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                bad_dt.bodies()[0].angular_velocity, start.angular_velocity, 0.0));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                bad_dt.bodies()[0].position, start.position, 0.0));

            RigidBody3D broken = start;
            broken.inertia = Vec3{1.0, -2.0, 3.0};
            Rigid3DWorld bad_body{SimulationSettings{0.001}, {broken}};
            MALLOY_CHECK_TRUE(bad_body.validate() == StepStatus::InvalidState);
            MALLOY_CHECK_TRUE(bad_body.step().status == StepStatus::InvalidState);
            MALLOY_CHECK_EQ(bad_body.tick_count(), std::uint64_t{0});

            Rigid3DWorld fine{SimulationSettings{0.001}, {start}};
            MALLOY_CHECK_TRUE(fine.validate() == StepStatus::Ok);
            MALLOY_CHECK_TRUE(fine.step().ok());
            MALLOY_CHECK_EQ(fine.tick_count(), std::uint64_t{1});

            // An empty world is valid and steps: no bodies is not an error.
            Rigid3DWorld empty{SimulationSettings{0.001}, {}};
            MALLOY_CHECK_TRUE(empty.step().ok());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::rigid::total_angular_momentum3d(empty.bodies()), Vec3{}, 0.0));
        }

        // --- A sphere cannot precess.
        //
        //     Every right-hand side in Euler's equations is a DIFFERENCE of
        //     principal moments, so when all three are equal the angular
        //     velocity is bit-for-bit constant no matter which way the body is
        //     spinning. Written with a generic omega, along no principal axis,
        //     so the only thing making the rate vanish is the difference. ---
        {
            RigidBody3D ball;
            ball.inertia = Vec3{2.5, 2.5, 2.5};
            ball.angular_velocity = Vec3{0.7, -1.3, 0.4};

            Rigid3DWorld world{SimulationSettings{0.001}, {ball}};
            const Vec3 l0 = spin_angular_momentum(world.bodies()[0]);
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                world.bodies()[0].angular_velocity, ball.angular_velocity, 0.0));

            // The body spins about omega and its angular momentum is parallel
            // to omega, so the world-frame vector is fixed as well, to rounding.
            // Rounding only: the quaternion's axis stays parallel to omega
            // to within a few ulp per step, so the bound is the same
            // N half-ulps of the magnitude that the 2D angle accumulates.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                spin_angular_momentum(world.bodies()[0]), l0,
                (1.0 / 9007199254740992.0) * malloy::math::length(l0) *
                    2000.0 / 2.0));
            MALLOY_CHECK_TRUE(malloy::math::is_unit(world.bodies()[0].orientation));
        }

        // --- Spin about a single principal axis: exactly constant, and the
        //     orientation follows the same discrete polygon the 2D angle does.
        //
        //     With omega along a principal axis, `L x omega` is zero, so the
        //     drift law below predicts no drift at all and the state is
        //     bit-identical after thousands of steps. ---
        {
            const Real spin = 1.3;
            const Real dt = 0.001;
            const int steps = 5000;
            const Vec3 axes[3] = {Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0},
                                  Vec3{0.0, 0.0, 1.0}};
            for (int a = 0; a < 3; ++a)
            {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = axes[a] * spin;

                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                const Vec3 l0 = spin_angular_momentum(world.bodies()[0]);
                const Real t0 = rotational_energy(world.bodies()[0]);
                for (int i = 0; i < steps; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                }
                const RigidBody3D& out = world.bodies()[0];
                MALLOY_CHECK_TRUE(
                    malloy::math::approx_equal(out.angular_velocity, body.angular_velocity, 0.0));
                MALLOY_CHECK_TRUE(
                    malloy::math::approx_equal(spin_angular_momentum(out), l0, 0.0));
                MALLOY_CHECK_NEAR(rotational_energy(out), t0, 0.0);
                MALLOY_CHECK_TRUE(malloy::math::is_unit(out.orientation));

                // The quaternion update is a rotation by omega/2 in the
                // (w, axis) plane stepped with explicit Euler, so each step
                // turns the HALF angle by exactly atan(omega dt / 2) and the
                // renormalization removes the magnitude growth. The body has
                // therefore turned by twice that, N times over.
                const Real half = static_cast<Real>(steps) * std::atan(spin * dt * 0.5);
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    out.orientation, malloy::math::from_axis_angle(axes[a], 2.0 * half),
                    1e-14));

                // It really did turn: 6.5 radians, not a body sitting still.
                MALLOY_CHECK_TRUE(2.0 * half > 6.0);
            }
        }

        // --- A symmetric top precesses at a rate that is exactly predictable.
        //
        //     With Ix == Iy the third equation has a zero coefficient, so wz is
        //     bit-constant, and the first two are a plane rotation at
        //     Omega = wz (Iz - Ix) / Ix stepped with explicit Euler. That is the
        //     same discrete polygon as the cyclotron in M18: per step the
        //     magnitude grows by exactly sqrt(1 + (Omega dt)^2) and the angle
        //     advances by exactly atan(Omega dt), never by Omega dt. ---
        {
            const Real dt = 0.001;
            const int steps = 2000;
            const Real across = 0.4;
            const Real along = 1.5;

            RigidBody3D top;
            top.inertia = Vec3{1.0, 1.0, 2.0};
            top.angular_velocity = Vec3{across, 0.0, along};

            Rigid3DWorld world{SimulationSettings{dt}, {top}};
            for (int i = 0; i < steps; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }

            const Real omega = along * (2.0 - 1.0) / 1.0;
            const Real growth =
                std::pow(1.0 + omega * omega * dt * dt, static_cast<Real>(steps) / 2.0);
            const Real turned = static_cast<Real>(steps) * std::atan(omega * dt);
            const Vec3& w = world.bodies()[0].angular_velocity;

            // A few ulp per step over 2000 steps of a component near 0.4.
            const Real bound = 1e-12;
            MALLOY_CHECK_NEAR(w.x, across * growth * std::cos(turned), bound);
            MALLOY_CHECK_NEAR(w.y, across * growth * std::sin(turned), bound);
            MALLOY_CHECK_NEAR(w.z, along, 0.0); // (Ix - Iy) is exactly zero

            // The precession is real, not a rounding artefact: the transverse
            // part has swung nearly 3 radians round.
            MALLOY_CHECK_TRUE(turned > 2.9);
            MALLOY_CHECK_TRUE(w.x < 0.0); // past a quarter turn, into the far half
        }

        // --- What semi-implicit Euler does to the conserved quantities, derived
        //     rather than tolerated.
        //
        //     The update is omega' = omega + dt I^-1 u with u = L x omega, and
        //     both L . u and omega . u vanish identically, so the first-order
        //     terms cancel and exactly the second-order ones survive:
        //
        //         |L'|^2 = |L|^2 + dt^2 |u|^2
        //         T'     = T     + dt^2 (I^-1 u) . u / 2
        //
        //     Neither is an inequality or a bound. They are equalities, they are
        //     checked every step, and they explain the exact cases above: when
        //     omega lies along a principal axis L is parallel to omega, u is
        //     zero, and nothing drifts at all. ---
        {
            const Real dt = 0.01;
            RigidBody3D body;
            body.inertia = asymmetric;
            body.angular_velocity = Vec3{0.9, 1.4, -0.6};

            Rigid3DWorld world{SimulationSettings{dt}, {body}};
            Real worst_momentum = 0.0;
            Real worst_energy = 0.0;
            Real worst_frame = 0.0;
            for (int i = 0; i < 200; ++i)
            {
                const Vec3 w = world.bodies()[0].angular_velocity;
                const Vec3 l{asymmetric.x * w.x, asymmetric.y * w.y, asymmetric.z * w.z};
                const Vec3 u = malloy::math::cross(l, w);
                const Vec3 iu{u.x / asymmetric.x, u.y / asymmetric.y,
                              u.z / asymmetric.z};
                const Real momentum_before = malloy::math::length_squared(l);
                const Real energy_before = rotational_energy(world.bodies()[0]);

                MALLOY_CHECK_TRUE(world.step().ok());

                const RigidBody3D& out = world.bodies()[0];
                const Vec3 w1 = out.angular_velocity;
                const Vec3 l1{asymmetric.x * w1.x, asymmetric.y * w1.y,
                              asymmetric.z * w1.z};
                const Real momentum_after = malloy::math::length_squared(l1);

                worst_momentum = std::fmax(
                    worst_momentum,
                    std::abs((momentum_after - momentum_before) -
                             dt * dt * malloy::math::length_squared(u)));
                worst_energy = std::fmax(
                    worst_energy,
                    std::abs((rotational_energy(out) - energy_before) -
                             0.5 * dt * dt * malloy::math::dot(iu, u)));

                // The world-frame vector has the same length as the body-frame
                // one, which is the whole reason the orientation must stay a
                // unit quaternion: a scaled one would not preserve length.
                worst_frame =
                    std::fmax(worst_frame, std::abs(malloy::math::length(
                                                        spin_angular_momentum(out)) -
                                                    std::sqrt(momentum_after)));
            }

            // |L|^2 is about 12 here, so one ulp of it is 1.8e-15, and each
            // comparison costs a few: the difference of two such numbers plus
            // the rounding in the predicted term. The drift being measured is
            // 3.5e-4, so this pins it to eleven significant figures.
            MALLOY_CHECK_TRUE(worst_momentum < 1e-13);
            MALLOY_CHECK_TRUE(worst_energy < 1e-13);
            MALLOY_CHECK_TRUE(worst_frame < 1e-13);

            // Both quantities GROW. They do not merely fail to be constant, and
            // the sign is not a matter of luck: |u|^2 and (I^-1 u) . u are both
            // sums of squares over positive moments.
            MALLOY_CHECK_TRUE(rotational_energy(world.bodies()[0]) >
                              rotational_energy(body));
        }

        // --- The intermediate-axis theorem.
        //
        //     The headline invariant of M20, and the thing 2D cannot express at
        //     all. With moments (1, 2, 3) a body spun about x or z returns to
        //     where it started; spun about y, the intermediate axis, it flips
        //     over. The same equations, the same code, the same perturbation:
        //     only the axis differs.
        //
        //     Linearizing about a spin W along one axis gives a second-order
        //     equation for the two transverse components. With the moments
        //     SORTED as I1 < I2 < I3, the three rates are
        //
        //         about I2 (intermediate):  +W^2 (I2-I1)(I3-I2) / (I1 I3)
        //         about I1 (smallest):      -W^2 (I2-I1)(I3-I1) / (I2 I3)
        //         about I3 (largest):       -W^2 (I3-I1)(I3-I2) / (I1 I2)
        //
        //     Every bracket is a larger moment minus a smaller one and so is
        //     positive, which puts the whole theorem in the leading sign:
        //     exponential about the intermediate axis, oscillatory about the
        //     other two.
        //
        //     For (1, 2, 3) those are +W^2/3, -W^2/3 and -W^2, so the rate is
        //     W/sqrt(3) about the intermediate and smallest axes and W about
        //     the largest. The two stable cases below therefore check two
        //     DIFFERENT numbers rather than the same one twice. ---
        {
            const Real dt = 0.001;
            const int steps = 20000;
            const Real spin = 2.0;
            const Real nudge = 1.0e-3;
            const Real rate = spin / std::sqrt(3.0);

            // Spun about x, the SMALLEST moment. The transverse pair conserves
            // wy^2 + 3 wz^2, so it traces an ellipse and comes back.
            {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = Vec3{spin, nudge, 0.0};

                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                Real worst = 0.0;
                Real lowest_spin = spin;
                for (int i = 0; i < steps; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const Vec3& w = world.bodies()[0].angular_velocity;
                    worst = std::fmax(worst, std::sqrt(w.y * w.y + 3.0 * w.z * w.z));
                    lowest_spin = std::fmin(lowest_spin, w.x);
                }

                // Explicit Euler on a rotation grows the amplitude by exactly
                // sqrt(1 + (mu dt)^2) per step, the same polygon factor as the
                // cyclotron, so this is an equality and not a fudge. What is
                // left over is the O(nudge^2) the linearization drops.
                const Real predicted =
                    nudge * std::pow(1.0 + rate * rate * dt * dt,
                                     static_cast<Real>(steps) / 2.0);
                MALLOY_CHECK_TRUE(worst <= predicted * (1.0 + 1e-6));
                MALLOY_CHECK_TRUE(worst >= predicted * 0.99);

                // The spin itself never wavers: it is disturbed only at second
                // order in the nudge.
                MALLOY_CHECK_TRUE(std::abs(lowest_spin - spin) < 1e-5);
            }

            // Spun about z, the LARGEST moment. Here the transverse pair
            // conserves wx^2 + wy^2 exactly, a circle rather than an ellipse.
            {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = Vec3{0.0, nudge, spin};

                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                Real worst = 0.0;
                Real lowest_spin = spin;
                for (int i = 0; i < steps; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const Vec3& w = world.bodies()[0].angular_velocity;
                    worst = std::fmax(worst, std::sqrt(w.x * w.x + w.y * w.y));
                    lowest_spin = std::fmin(lowest_spin, w.z);
                }

                // mu is W rather than W/sqrt(3) on this axis, which is a
                // different number and therefore a real check on the
                // coefficients rather than a repeat of the case above.
                const Real predicted =
                    nudge * std::pow(1.0 + spin * spin * dt * dt,
                                     static_cast<Real>(steps) / 2.0);
                MALLOY_CHECK_TRUE(worst <= predicted * (1.0 + 1e-6));
                MALLOY_CHECK_TRUE(worst >= predicted * 0.99);
                MALLOY_CHECK_TRUE(std::abs(lowest_spin - spin) < 1e-5);
            }

            // Spun about y, the INTERMEDIATE moment. The body flips.
            {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = Vec3{nudge, spin, 0.0};

                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                Real lowest = spin;
                Real fastest = 0.0;
                for (int i = 0; i < steps; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const Vec3& w = world.bodies()[0].angular_velocity;
                    lowest = std::fmin(lowest, w.y);
                    fastest = std::fmax(fastest, malloy::math::length(w));
                }

                // It does not wobble: it reverses. A nudge of a thousandth
                // turns into a full change of sign of the spin itself, which is
                // four orders of magnitude larger than anything the two stable
                // axes did with the identical nudge.
                MALLOY_CHECK_TRUE(lowest < -0.9 * spin);

                // And it reverses along the path the conserved quantities
                // allow. Starting on the separatrix fixes |L|^2 = (I2 W)^2 and
                // 2T = I2 W^2; where the body passes through wy = 0 those two
                // give wz^2 = I2 (I2 - I1) W^2 / (I3 (I3 - I1)) = W^2/3 and
                // wx^2 = W^2, so the fastest the body ever turns is
                //
                //     |omega| = 2 W / sqrt(3)
                //
                // The gap from that is the accumulated per-step drift measured
                // above, which over 20000 steps is a few parts in ten thousand.
                const Real peak = 2.0 * spin / std::sqrt(3.0);
                MALLOY_CHECK_TRUE(fastest > peak);
                MALLOY_CHECK_TRUE(fastest < peak * 1.001);
            }

            // The exponential growth rate itself, against the exact discrete
            // linearization rather than its continuum limit. The transverse map
            // has eigenvalues 1 +/- rate*dt, so starting from (nudge, 0) the
            // growing component is exactly
            //
            //     (nudge/2) [ (1 + rate dt)^n + (1 - rate dt)^n ]
            //
            // A small nudge keeps the dropped O(nudge^2) terms near 1e-17, so
            // this is asserted to fifteen significant figures.
            {
                const Real tiny = 1.0e-6;
                const int short_run = 2000;
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = Vec3{tiny, spin, 0.0};

                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                for (int i = 0; i < short_run; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                }
                const Real predicted =
                    0.5 * tiny *
                    (std::pow(1.0 + rate * dt, short_run) +
                     std::pow(1.0 - rate * dt, short_run));
                MALLOY_CHECK_NEAR(world.bodies()[0].angular_velocity.x, predicted, 1e-15);

                // Five times its starting size in two seconds, and still only
                // five parts per million of the spin: this is the early phase,
                // where the linearization is what governs.
                MALLOY_CHECK_TRUE(predicted > 5.0 * tiny);
            }
        }

        // --- The orientation is built from the angular velocity AFTER Euler's
        //     equations have been applied, not before.
        //
        //     That is what makes this semi-implicit, the same "velocities
        //     first" ordering every other world here uses, and it is otherwise
        //     invisible: swapping the two changes the result only at O(dt), the
        //     same order as the scheme's own error, so no convergence test can
        //     see it.
        //
        //     It is visible as a DIRECTION. Starting from the identity, the
        //     quaternion increment is (dt/2)(0, omega), so after one step the
        //     orientation's vector part points exactly along whichever omega
        //     was used. Here the two differ by rate * dt, which is large enough
        //     to separate them by eleven orders of magnitude. ---
        {
            const Real dt = 0.01;
            RigidBody3D body;
            body.inertia = asymmetric;
            body.angular_velocity = Vec3{0.9, 1.4, -0.6};

            Rigid3DWorld world{SimulationSettings{dt}, {body}};
            MALLOY_CHECK_TRUE(world.step().ok());

            const RigidBody3D& out = world.bodies()[0];
            MALLOY_CHECK_TRUE(malloy::math::is_unit(out.orientation));

            // Parallel to the updated omega, to rounding.
            MALLOY_CHECK_TRUE(
                malloy::math::length(
                    malloy::math::cross(out.orientation.v, out.angular_velocity)) < 1e-16);

            // And measurably NOT parallel to the one the step started with, so
            // this distinguishes the two rather than holding either way.
            MALLOY_CHECK_TRUE(
                malloy::math::length(
                    malloy::math::cross(out.orientation.v, body.angular_velocity)) > 1e-5);

            // The step really did change omega, which is what gives the two
            // directions something to disagree about.
            MALLOY_CHECK_FALSE(malloy::math::approx_equal(
                out.angular_velocity, body.angular_velocity, 1e-4));
        }

        // --- Angular momentum is spin PLUS orbital, as it is in two dimensions.
        //
        //     `spin_angular_momentum` is R (I omega) and nothing else.
        //     `total_angular_momentum3d` adds m (r x v) about the world origin,
        //     because dropping it makes a whole class of motion look conserved
        //     when it is not: a body flying past the origin carries angular
        //     momentum about it whether or not it is spinning.
        //
        //     Written with a velocity NOT parallel to the position, so the
        //     orbital term is nonzero in all three components. A configuration
        //     with r parallel to v would have an orbital term of zero and would
        //     pass whether or not it was summed. ---
        {
            RigidBody3D flyer;
            flyer.mass = 2.0;
            flyer.inertia = asymmetric;
            flyer.position = Vec3{1.0, 2.0, -1.0};
            flyer.velocity = Vec3{0.3, -0.1, 0.5};
            flyer.angular_velocity = Vec3{0.4, 1.1, -0.9};

            // 2 * ((1,2,-1) x (0.3,-0.1,0.5)) = 2 * (0.9,-0.8,-0.7), by hand.
            const Vec3 orbital_start{1.8, -1.6, -1.4};
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                flyer.mass * malloy::math::cross(flyer.position, flyer.velocity),
                orbital_start, 1e-15));

            // Large enough that omitting it could not hide in a tolerance: the
            // orbital term here is comparable to the spin term, not a
            // correction to it.
            MALLOY_CHECK_TRUE(malloy::math::length(orbital_start) > 2.0);

            Rigid3DWorld world{SimulationSettings{0.001}, {flyer}};
            const Vec3 total_start =
                malloy::rigid::total_angular_momentum3d(world.bodies());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                total_start, spin_angular_momentum(flyer) + orbital_start, 1e-15));

            // The two are genuinely different vectors, so a total that returned
            // only the spin would fail here rather than pass by coincidence.
            MALLOY_CHECK_FALSE(malloy::math::approx_equal(
                total_start, spin_angular_momentum(flyer), 1e-6));

            const int steps = 4000;
            Real worst_orbital = 0.0;
            for (int i = 0; i < steps; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                const RigidBody3D& out = world.bodies()[0];
                const Vec3 orbital =
                    out.mass * malloy::math::cross(out.position, out.velocity);

                // The total is exactly the two parts added, every step.
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    malloy::rigid::total_angular_momentum3d(world.bodies()),
                    spin_angular_momentum(out) + orbital, 0.0));

                worst_orbital = std::fmax(
                    worst_orbital, malloy::math::length(orbital - orbital_start));
            }

            // The orbital term is exactly constant in exact arithmetic, since
            // d(r x p)/dt = v x mv = 0 for a free body. What is left is the
            // rounding in r, which accumulates over N steps of a quantity of
            // this size: the same N half-ulp bound the 2D angle uses.
            MALLOY_CHECK_TRUE(worst_orbital <
                              (1.0 / 9007199254740992.0) *
                                  malloy::math::length(orbital_start) *
                                  static_cast<Real>(steps) / 2.0);

            // The body really did travel, so the constancy above is a
            // cancellation and not a body sitting at the origin.
            MALLOY_CHECK_TRUE(malloy::math::distance(world.bodies()[0].position,
                                                     flyer.position) > 1.0);
        }

        // --- The zero-drift condition is an EIGENSPACE condition, not a
        //     principal-axis one.
        //
        //     The drift vanishes exactly when u = L x omega is zero, that is
        //     when I omega is parallel to omega. For three distinct moments the
        //     only such directions are the three principal axes. For an
        //     AXISYMMETRIC body two moments are equal, that eigenspace is a
        //     plane, and every direction in it qualifies: not just the two
        //     principal axes lying in it.
        //
        //     The angular velocity below is in the degenerate plane but along
        //     neither axis of it, so it is a case the principal-axis statement
        //     does not cover and the eigenspace one does. ---
        {
            RigidBody3D top;
            top.inertia = Vec3{1.5, 1.5, 4.0}; // axisymmetric about z
            top.angular_velocity = Vec3{0.6, -0.8, 0.0};

            // In the eigenspace: I omega is a multiple of omega, so u is zero.
            const Vec3& w0 = top.angular_velocity;
            const Vec3 iw{top.inertia.x * w0.x, top.inertia.y * w0.y,
                          top.inertia.z * w0.z};
            // Zero to rounding: the eigenspace relation is exact, but the
            // hand-chosen components are not bit-exact in binary.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::math::cross(iw, w0), Vec3{}, 1e-15));

            // And along no principal axis: both in-plane components are
            // nonzero, which is what makes this a different case.
            MALLOY_CHECK_TRUE(std::abs(w0.x) > 0.1);
            MALLOY_CHECK_TRUE(std::abs(w0.y) > 0.1);

            Rigid3DWorld world{SimulationSettings{0.001}, {top}};
            const Real energy_start = rotational_energy(world.bodies()[0]);
            for (int i = 0; i < 3000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                world.bodies()[0].angular_velocity, Vec3(0.6, -0.8, 0.0), 0.0));
            MALLOY_CHECK_NEAR(rotational_energy(world.bodies()[0]), energy_start, 0.0);

            // Tilt the same body out of that plane and it drifts, so the
            // exactness above is the eigenspace and not the inertia.
            RigidBody3D tilted = top;
            tilted.angular_velocity = Vec3{0.6, -0.8, 0.5};
            Rigid3DWorld drifting{SimulationSettings{0.001}, {tilted}};
            const Real tilted_start = rotational_energy(drifting.bodies()[0]);
            for (int i = 0; i < 3000; ++i)
            {
                MALLOY_CHECK_TRUE(drifting.step().ok());
            }
            MALLOY_CHECK_TRUE(rotational_energy(drifting.bodies()[0]) > tilted_start);
        }

        // --- The world frame is not the body frame, and that is the whole
        //     reason `angular_momentum` rotates.
        //
        //     Under torque-free rotation the WORLD-frame angular momentum is
        //     fixed while the BODY-frame vector I omega tumbles with the body.
        //     Measuring both at two step sizes separates the two claims from
        //     each other: the world-frame deviation is integration error and
        //     falls with dt, the body-frame swing is the physics and does not.
        //
        //     The orientation is integrated with a first-order scheme, so the
        //     former is O(dt) and a tenfold smaller step buys a tenfold smaller
        //     deviation. That ratio is asserted, because "small" alone would
        //     also be satisfied by an implementation that never rotated
        //     anything at all. ---
        {
            Real deviation[2] = {0.0, 0.0};
            Real tumble[2] = {0.0, 0.0};
            const Real steps[2] = {1e-3, 1e-4};
            for (int k = 0; k < 2; ++k)
            {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = Vec3{0.9, 1.4, -0.6};

                Rigid3DWorld world{SimulationSettings{steps[k]}, {body}};
                const Vec3 world_start = spin_angular_momentum(world.bodies()[0]);
                const Vec3 body_start{asymmetric.x * body.angular_velocity.x,
                                      asymmetric.y * body.angular_velocity.y,
                                      asymmetric.z * body.angular_velocity.z};

                const int count = static_cast<int>(2.0 / steps[k]);
                for (int i = 0; i < count; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const RigidBody3D& out = world.bodies()[0];
                    const Vec3 in_body{asymmetric.x * out.angular_velocity.x,
                                       asymmetric.y * out.angular_velocity.y,
                                       asymmetric.z * out.angular_velocity.z};
                    deviation[k] = std::fmax(
                        deviation[k],
                        malloy::math::length(spin_angular_momentum(out) - world_start));
                    tumble[k] =
                        std::fmax(tumble[k], malloy::math::length(in_body - body_start));
                }
            }

            // First order: ten times the step, ten times the error.
            const Real ratio = deviation[0] / deviation[1];
            MALLOY_CHECK_TRUE(ratio > 9.5);
            MALLOY_CHECK_TRUE(ratio < 10.5);

            // The body really is tumbling, by more than |L| itself, and by the
            // same amount at both step sizes because it is not an error term.
            MALLOY_CHECK_TRUE(tumble[0] > 5.0);
            MALLOY_CHECK_TRUE(tumble[1] > 5.0);
            MALLOY_CHECK_NEAR(tumble[0], tumble[1], 0.01);

            // Three orders of magnitude between the two, which is what makes
            // returning the unrotated vector from `angular_momentum` a failure
            // rather than a rounding difference.
            MALLOY_CHECK_TRUE(tumble[0] > 1000.0 * deviation[0]);
        }

        // --- Nothing acts on the translation, and the diagnostics are sums.
        //
        //     A tumbling body still flies straight, which is the part of this
        //     that is NOT new: the rotation is uncoupled from the centre of
        //     mass because no contact or force couples them yet. ---
        {
            const Real dt = 0.001;
            const int steps = 4000;

            RigidBody3D first;
            first.mass = 1.5;
            first.inertia = asymmetric;
            first.velocity = Vec3{0.3, -0.7, 0.2};
            first.angular_velocity = Vec3{0.4, 1.1, -0.9};

            RigidBody3D second;
            second.mass = 2.5;
            second.inertia = Vec3{2.0, 0.5, 4.0};
            second.position = Vec3{1.0, 2.0, -1.0};
            second.velocity = Vec3{-0.1, 0.4, 0.6};
            second.angular_velocity = Vec3{-1.2, 0.3, 0.8};
            second.orientation = malloy::math::from_axis_angle(Vec3{1.0, 1.0, 1.0}, 0.9);

            const std::vector<RigidBody3D> start = {first, second};
            const Vec3 momentum = malloy::rigid::total_linear_momentum3d(start);

            Rigid3DWorld world{SimulationSettings{dt}, start};
            for (int i = 0; i < steps; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }

            // No force acts, so the velocities and hence the linear momentum
            // are bit-identical, not merely close.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::rigid::total_linear_momentum3d(world.bodies()), momentum, 0.0));

            // The positions are a sum of N identical increments, so the error is
            // N/2 half-ulps of the final value.
            for (std::size_t i = 0; i < start.size(); ++i)
            {
                const Vec3 exact =
                    start[i].position + start[i].velocity * (dt * static_cast<Real>(steps));
                const Real bound = (1.0 / 9007199254740992.0) *
                                   malloy::math::length(exact) *
                                   static_cast<Real>(steps) / 2.0;
                MALLOY_CHECK_TRUE(
                    malloy::math::approx_equal(world.bodies()[i].position, exact, bound));
                MALLOY_CHECK_TRUE(malloy::math::is_unit(world.bodies()[i].orientation));
            }

            // Both bodies really are tumbling, so the run above is not a pair of
            // trivially exact principal-axis spins in disguise.
            for (const RigidBody3D& body : world.bodies())
            {
                MALLOY_CHECK_TRUE(rotational_energy(body) > 0.0);
            }
            MALLOY_CHECK_FALSE(malloy::math::approx_equal(
                world.bodies()[0].angular_velocity, first.angular_velocity, 1e-6));

            // The energy is the translational part plus the rotational part, to
            // rounding, and the translational part alone is exactly conserved.
            Real expected = 0.0;
            for (const RigidBody3D& body : world.bodies())
            {
                expected += 0.5 * body.mass * malloy::math::dot(body.velocity, body.velocity);
                expected += rotational_energy(body);
            }
            MALLOY_CHECK_NEAR(malloy::rigid::total_kinetic_energy3d(world.bodies()),
                              expected, 1e-15);
        }
    }

    // --- M21: a constant applied torque, in the world frame.
    //
    //     M20 was torque-free: nothing acted on the rotation, so a body
    //     spinning about a principal axis stayed there forever. M21 adds a
    //     torque, and with it the one thing 2D cannot do with a torque, which
    //     is move the spin axis rather than only speed the spin up.
    //
    //     The torque is a WORLD-frame setting. That is what makes the headline
    //     invariant exact: dL/dt = torque holds in the world frame, so the
    //     world-frame angular momentum grows along a straight line. ---
    {
        using malloy::math::Quat;
        using malloy::math::Vec3;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::rotational_energy;
        using malloy::rigid::spin_angular_momentum;
        using malloy::rigid::total_angular_momentum3d;

        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();
        const Vec3 asymmetric{1.0, 2.0, 3.0};

        // --- The torque setting is validated, and squarable rather than merely
        //     finite: it enters |u|^2 in the drift law, so one whose square
        //     overflows is refused up front, the bound softening carries. ---
        {
            MALLOY_CHECK_TRUE(Rigid3DSettings{}.is_valid());          // zero
            MALLOY_CHECK_TRUE((Rigid3DSettings{Vec3{1.0, -2.0, 3.0}}.is_valid()));
            MALLOY_CHECK_FALSE((Rigid3DSettings{Vec3{inf3, 0.0, 0.0}}.is_valid()));
            MALLOY_CHECK_FALSE((Rigid3DSettings{Vec3{0.0, nan3, 0.0}}.is_valid()));
            MALLOY_CHECK_FALSE((Rigid3DSettings{Vec3{0.0, 0.0, 1.4e154}}.is_valid()));

            RigidBody3D body;
            body.inertia = asymmetric;
            Rigid3DWorld bad{SimulationSettings{0.001}, {body},
                             Rigid3DSettings{Vec3{inf3, 0.0, 0.0}}};
            MALLOY_CHECK_TRUE(bad.validate() == StepStatus::InvalidSettings);
            MALLOY_CHECK_TRUE(bad.step().status == StepStatus::InvalidSettings);
            MALLOY_CHECK_EQ(bad.tick_count(), std::uint64_t{0});
        }

        // --- Zero torque is the M20 world, bit for bit. A defaulted settings
        //     and an explicit zero must both reproduce the torque-free run
        //     exactly, or the torque term is doing something when it should do
        //     nothing. ---
        {
            RigidBody3D body;
            body.inertia = asymmetric;
            body.angular_velocity = Vec3{0.4, 1.1, -0.9};
            body.orientation = malloy::math::from_axis_angle(Vec3{1.0, 1.0, 1.0}, 0.9);

            Rigid3DWorld defaulted{SimulationSettings{0.001}, {body}};
            Rigid3DWorld zeroed{SimulationSettings{0.001}, {body},
                                Rigid3DSettings{Vec3{}}};
            for (int i = 0; i < 3000; ++i)
            {
                MALLOY_CHECK_TRUE(defaulted.step().ok());
                MALLOY_CHECK_TRUE(zeroed.step().ok());
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    defaulted.bodies()[0].angular_velocity,
                    zeroed.bodies()[0].angular_velocity, 0.0));
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    defaulted.bodies()[0].orientation, zeroed.bodies()[0].orientation, 0.0));
            }
        }

        // --- Spin-up about a principal axis, EXACTLY.
        //
        //     A body at rest, unrotated, with the torque along a principal
        //     axis. The torque is world-frame, but the body only ever rotates
        //     about that same axis, and a rotation about an axis fixes it, so
        //     the body-frame torque stays equal to the world one bit for bit.
        //     The gyroscopic term is zero because omega stays along the axis.
        //     What is left is I wx' = T, integrated exactly by forward Euler
        //     because the right side is constant:
        //
        //       wx(n) = n dt T / Ix,   wy = wz = 0,   L_world = (T t, 0, 0).
        //
        //     The transverse components are not just small, they are zero, and
        //     the axis is not just nearly fixed, it is fixed. ---
        {
            const Real T = 0.5;
            const Real dt = 0.001;
            const int steps = 5000;
            RigidBody3D body;
            body.inertia = asymmetric;
            Rigid3DWorld world{SimulationSettings{dt}, {body},
                               Rigid3DSettings{Vec3{T, 0.0, 0.0}}};

            for (int n = 1; n <= steps; ++n)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                const RigidBody3D& out = world.bodies()[0];

                // Transverse components exactly zero: the axis never moves.
                MALLOY_CHECK_NEAR(out.angular_velocity.y, 0.0, 0.0);
                MALLOY_CHECK_NEAR(out.angular_velocity.z, 0.0, 0.0);

                // The spin grows linearly. The bound is the N-half-ulp
                // accumulation of a repeated constant addition, the same one
                // the 2D angle test uses.
                const Real expected = static_cast<Real>(n) * dt * T / asymmetric.x;
                const Real bound =
                    (1.0 / 9007199254740992.0) * expected * static_cast<Real>(n) / 2.0;
                MALLOY_CHECK_NEAR(out.angular_velocity.x, expected, bound);

                // The world-frame angular momentum is (T t, 0, 0) exactly: the
                // transverse parts are zero, and the x part is Ix wx unrotated.
                const Vec3 lw = total_angular_momentum3d(world.bodies());
                MALLOY_CHECK_NEAR(lw.x, T * static_cast<Real>(n) * dt,
                                  asymmetric.x * bound);
                MALLOY_CHECK_NEAR(lw.y, 0.0, 0.0);
                MALLOY_CHECK_NEAR(lw.z, 0.0, 0.0);
            }

            // It really spun up: from rest to 2.5 rad/s, not a body sitting
            // still under a torque that quietly did nothing.
            MALLOY_CHECK_TRUE(world.bodies()[0].angular_velocity.x > 2.0);
        }

        // --- The forced per-step laws, derived and asserted every step.
        //
        //     With L = I omega and u = L x omega + T_body, the forward-Euler
        //     step omega' = omega + dt I^-1 u gives, exactly,
        //
        //       |L'|^2 = |L|^2 + 2 dt (L . T_body) + dt^2 |u|^2
        //       T'     = T     +     dt (omega . T_body) + dt^2 (I^-1 u).u / 2
        //
        //     The cross term drops out of both first-order parts because it is
        //     orthogonal to L and to omega; only the torque survives there.
        //     With T_body = 0 these are the M20 laws. A tumbling body, a tilted
        //     start and an off-axis torque, so no term is accidentally zero. ---
        {
            const Vec3 torque{0.7, -0.4, 0.9};
            const Real dt = 0.01;
            RigidBody3D body;
            body.inertia = asymmetric;
            body.angular_velocity = Vec3{0.9, 1.4, -0.6};
            body.orientation = malloy::math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.8);

            Rigid3DWorld world{SimulationSettings{dt}, {body}, Rigid3DSettings{torque}};
            Real worst_momentum = 0.0;
            Real worst_energy = 0.0;
            for (int i = 0; i < 200; ++i)
            {
                const RigidBody3D& s = world.bodies()[0];
                const Vec3 tb =
                    malloy::math::rotate(malloy::math::conjugate(s.orientation), torque);
                const Vec3 w = s.angular_velocity;
                const Vec3 l{asymmetric.x * w.x, asymmetric.y * w.y, asymmetric.z * w.z};
                const Vec3 u = malloy::math::cross(l, w) + tb;
                const Vec3 iu{u.x / asymmetric.x, u.y / asymmetric.y, u.z / asymmetric.z};
                const Real momentum_before = malloy::math::length_squared(l);
                const Real energy_before = rotational_energy(s);

                MALLOY_CHECK_TRUE(world.step().ok());

                const RigidBody3D& s1 = world.bodies()[0];
                const Vec3 w1 = s1.angular_velocity;
                const Vec3 l1{asymmetric.x * w1.x, asymmetric.y * w1.y,
                              asymmetric.z * w1.z};

                const Real pred_momentum =
                    momentum_before + 2.0 * dt * malloy::math::dot(l, tb) +
                    dt * dt * malloy::math::length_squared(u);
                const Real pred_energy = energy_before + dt * malloy::math::dot(w, tb) +
                                         0.5 * dt * dt * malloy::math::dot(iu, u);
                worst_momentum = std::fmax(
                    worst_momentum, std::abs(malloy::math::length_squared(l1) - pred_momentum));
                worst_energy =
                    std::fmax(worst_energy, std::abs(rotational_energy(s1) - pred_energy));
            }
            MALLOY_CHECK_TRUE(worst_momentum < 1e-13);
            MALLOY_CHECK_TRUE(worst_energy < 1e-13);

            // The torque did net work here, so the energy is not merely
            // drifting at second order as it did in M20.
            MALLOY_CHECK_TRUE(rotational_energy(world.bodies()[0]) !=
                              rotational_energy(body));
        }

        // --- The world-frame linear law, for a body that is NOT spinning up
        //     cleanly: a fast symmetric top with a torque across its spin.
        //
        //     dL/dt = torque holds in the world frame no matter what the body
        //     does, so L_world(t) = L_world(0) + torque t is the continuum law
        //     and the discrete scheme tracks it to FIRST order. Measuring the
        //     deviation at two step sizes separates the law (exact) from the
        //     integration error (falls with dt), the same way M20 separated the
        //     world frame from the body frame. ---
        {
            const Real T = 0.3;
            const Real spin = 8.0;
            Real deviation[2] = {0.0, 0.0};
            Real tilt[2] = {0.0, 0.0};
            const Real steps[2] = {1e-3, 1e-4};
            for (int k = 0; k < 2; ++k)
            {
                RigidBody3D top;
                top.inertia = Vec3{1.0, 1.0, 2.0};
                top.angular_velocity = Vec3{0.0, 0.0, spin};
                Rigid3DWorld world{SimulationSettings{steps[k]}, {top},
                                   Rigid3DSettings{Vec3{T, 0.0, 0.0}}};
                const Vec3 start = total_angular_momentum3d(world.bodies());

                const int count = static_cast<int>(4.0 / steps[k]);
                for (int n = 1; n <= count; ++n)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const Vec3 lw = total_angular_momentum3d(world.bodies());
                    const Vec3 exact{start.x + T * static_cast<Real>(n) * steps[k],
                                     start.y, start.z};
                    deviation[k] = std::fmax(deviation[k], malloy::math::length(lw - exact));
                    tilt[k] = std::fmax(tilt[k], std::abs(lw.x));
                }
            }

            // First order: ten times the step, ten times the deviation.
            const Real ratio = deviation[0] / deviation[1];
            MALLOY_CHECK_TRUE(ratio > 9.5);
            MALLOY_CHECK_TRUE(ratio < 10.5);

            // The angular momentum really did tilt into x, by more than a
            // radian's worth, so the law is being tracked over real motion and
            // not a body sitting still.
            MALLOY_CHECK_TRUE(tilt[0] > 1.0);
        }

        // --- The gyroscopic response: a torque across the spin MOVES the spin
        //     axis, where in 2D a torque can only change the rate.
        //
        //     Two runs of the same fast top. One torque is along the spin axis
        //     and only speeds it up; the transverse spin stays put. The other
        //     is across the spin and swings the axis over. Same body, same
        //     torque magnitude, the only difference is direction, and the
        //     responses differ by orders of magnitude. ---
        {
            const Real magnitude = 0.5;
            const Real spin = 10.0;
            const Real dt = 0.001;
            const int steps = 3000;

            const auto transverse_swing = [&](const Vec3& torque) {
                RigidBody3D top;
                top.inertia = Vec3{1.0, 1.0, 2.0};
                top.angular_velocity = Vec3{0.0, 0.0, spin};
                Rigid3DWorld world{SimulationSettings{dt}, {top},
                                   Rigid3DSettings{torque}};
                Real worst = 0.0;
                for (int i = 0; i < steps; ++i)
                {
                    world.step();
                    const Vec3 lw = spin_angular_momentum(world.bodies()[0]);
                    worst = std::fmax(worst, std::sqrt(lw.x * lw.x + lw.y * lw.y));
                }
                return worst;
            };

            // Along z, the spin axis: the momentum stays along z, so the
            // transverse part barely moves.
            const Real along = transverse_swing(Vec3{0.0, 0.0, magnitude});
            // Across the spin, along x: the axis swings over.
            const Real across = transverse_swing(Vec3{magnitude, 0.0, 0.0});

            MALLOY_CHECK_TRUE(across > 1.0);        // a full radian of momentum
            MALLOY_CHECK_TRUE(along < 0.05);        // essentially unmoved
            MALLOY_CHECK_TRUE(across > 50.0 * along); // and the two are not close
        }
    }

    // --- M22: gravity and restitution contacts against a ground plane.
    //
    //     The 3D echo of M10's colliding particles, not M14's rigid contacts: a
    //     contact on a CENTRED sphere passes through the centre of mass, so it
    //     imparts no spin, and the response is a pure normal impulse. The sharp
    //     invariants are therefore translational, and the milestone's own claim,
    //     that the contact is rotationally inert, is itself a test below. ---
    {
        using malloy::collide::Plane3;
        using malloy::collide::Sphere;
        using malloy::math::Vec3;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::total_kinetic_energy3d;
        using malloy::rigid::total_potential_energy3d;

        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();

        // --- Validation: the radius, and the new settings. ---
        {
            RigidBody3D body;
            MALLOY_CHECK_TRUE(body.is_valid()); // default radius 0 is fine
            for (const Real bad : {-1.0, inf3, nan3})
            {
                RigidBody3D b;
                b.radius = bad;
                MALLOY_CHECK_FALSE(b.is_valid());
            }
            RigidBody3D collides;
            collides.radius = 0.5;
            MALLOY_CHECK_TRUE(collides.is_valid());

            MALLOY_CHECK_TRUE(Rigid3DSettings{}.is_valid());
            Rigid3DSettings good;
            good.restitution = 0.5;
            good.gravity = Vec3{0.0, 0.0, -9.81};
            good.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
            MALLOY_CHECK_TRUE(good.is_valid());

            for (const Real bad : {-0.1, 1.5, nan3, inf3})
            {
                Rigid3DSettings s;
                s.restitution = bad;
                MALLOY_CHECK_FALSE(s.is_valid());
            }
            Rigid3DSettings bad_gravity;
            bad_gravity.gravity = Vec3{0.0, 0.0, 1.0e200}; // square overflows
            MALLOY_CHECK_FALSE(bad_gravity.is_valid());
            Rigid3DSettings bad_plane;
            bad_plane.ground = {Plane3{Vec3{0.0, 0.0, 2.0}, 0.0}}; // not unit
            MALLOY_CHECK_FALSE(bad_plane.is_valid());

            RigidBody3D b;
            b.radius = 0.5;
            Rigid3DWorld world{SimulationSettings{0.001}, {b}, bad_gravity};
            MALLOY_CHECK_TRUE(world.validate() == StepStatus::InvalidSettings);
            MALLOY_CHECK_TRUE(world.step().status == StepStatus::InvalidSettings);
        }

        // --- Velocity restitution, EXACT, on a tilted plane with no gravity so
        //     the impulse is the only thing that changes velocity.
        //
        //     The normal component reverses to -e times itself; the tangential
        //     component is untouched (no friction); and the spin is untouched
        //     (a normal contact on a centred sphere has no moment arm). Tilted
        //     so the normal has all three components, and the body carries a
        //     tangential velocity and a spin so their invariance is a real
        //     check and not a zero. ---
        for (const Real e : {1.0, 0.5, 0.0})
        {
            const Vec3 n = malloy::math::normalize(Vec3{1.0, 2.0, 2.0});
            Rigid3DSettings st;
            st.restitution = e;
            st.ground = {Plane3{n, 0.0}};
            RigidBody3D b;
            b.mass = 2.0;
            b.radius = 1.0;
            b.inertia = Vec3{2.0, 2.0, 2.0}; // a sphere: isotropic, so no tumble
            b.position = n * 0.98; // penetrating by 0.02
            const Vec3 tangential{0.7, 0.0, -0.35}; // perpendicular to n
            b.velocity = n * (-3.0) + tangential;
            b.angular_velocity = Vec3{0.9, -1.1, 0.4};

            Rigid3DWorld world{SimulationSettings{0.001}, {b}, st};
            const Real vn_before = malloy::math::dot(b.velocity, n);
            MALLOY_CHECK_TRUE(world.step().ok());
            const RigidBody3D& out = world.bodies()[0];

            // Normal component reversed and scaled by e, to rounding.
            const Real vn_after = malloy::math::dot(out.velocity, n);
            MALLOY_CHECK_NEAR(vn_after, -e * vn_before, 1e-14);

            // Tangential component untouched: no friction.
            const Vec3 tangential_after =
                out.velocity - n * malloy::math::dot(out.velocity, n);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(tangential_after, tangential, 1e-14));

            // Spin untouched by the contact. The inertia is isotropic (a real
            // sphere), so with no torque there is no tumble either, and the
            // angular velocity is bit-identical to its start. The tumbling case
            // is covered by the rotationally-inert test below.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                out.angular_velocity, b.angular_velocity, 0.0));
        }

        // --- Energy across a single straight-down bounce: KE_after = e^2 KE.
        //
        //     A sphere moving straight into a floor, no gravity, so the only
        //     velocity is normal and all the kinetic energy is in it. The bounce
        //     scales the speed by e, so the energy scales by e^2, exactly. ---
        for (const Real e : {1.0, 0.6})
        {
            Rigid3DSettings st;
            st.restitution = e;
            st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
            RigidBody3D b;
            b.mass = 1.5;
            b.radius = 0.5;
            b.position = Vec3{0.0, 0.0, 0.49}; // penetrating
            b.velocity = Vec3{0.0, 0.0, -4.0};

            Rigid3DWorld world{SimulationSettings{0.001}, {b}, st};
            const Real ke_before = total_kinetic_energy3d(world.bodies());
            MALLOY_CHECK_TRUE(world.step().ok());
            const Real ke_after = total_kinetic_energy3d(world.bodies());
            MALLOY_CHECK_NEAR(ke_after, e * e * ke_before, 1e-13);
        }

        // --- Free-flight energy shed = (1/2)(sum m)|g|^2 dt^2 per step, the
        //     same derived constant as M12 and M15, now in 3D. No ground, so it
        //     is pure free flight. An equality, asserted every step. ---
        {
            const Vec3 g{0.4, -0.5, -9.81}; // gravity NOT axis-aligned
            const Real dt = 0.001;
            Rigid3DSettings st;
            st.gravity = g;
            RigidBody3D b;
            b.mass = 1.5;
            b.position = Vec3{2.0, -1.0, 100.0};
            b.velocity = Vec3{0.3, 0.7, 0.0};

            Rigid3DWorld world{SimulationSettings{dt}, {b}, st};
            const Real predicted =
                0.5 * b.mass * malloy::math::length_squared(g) * dt * dt;
            Real worst = 0.0;
            for (int i = 0; i < 500; ++i)
            {
                const Real before = total_kinetic_energy3d(world.bodies()) +
                                    total_potential_energy3d(world.bodies(), g);
                MALLOY_CHECK_TRUE(world.step().ok());
                const Real after = total_kinetic_energy3d(world.bodies()) +
                                   total_potential_energy3d(world.bodies(), g);
                worst = std::fmax(worst, std::abs((before - after) - predicted));
            }
            MALLOY_CHECK_TRUE(worst < 1e-11);
            // The shed is real and one-signed: energy falls, it does not wander.
            MALLOY_CHECK_TRUE(predicted > 0.0);
        }

        // --- The contact is rotationally INERT, which is M22's scope claim
        //     stated as a test. A spinning sphere under a torque, dropped onto a
        //     floor, has exactly the same orientation and angular velocity as
        //     the same body with no floor at all: the bounce changes where it is
        //     and how fast it moves, never how it spins. Bit for bit. ---
        {
            Rigid3DSettings grounded;
            grounded.restitution = 0.8;
            grounded.gravity = Vec3{0.0, 0.0, -9.81};
            grounded.torque = Vec3{0.3, -0.2, 0.1};
            grounded.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};

            Rigid3DSettings floating = grounded;
            floating.ground = {}; // identical but for the floor

            RigidBody3D b;
            b.radius = 0.5;
            b.inertia = Vec3{1.0, 2.0, 3.0};
            b.position = Vec3{0.0, 0.0, 1.2};
            b.velocity = Vec3{0.1, 0.0, -2.0};
            b.angular_velocity = Vec3{1.0, 0.5, -0.3};

            Rigid3DWorld with_floor{SimulationSettings{0.001}, {b}, grounded};
            Rigid3DWorld without{SimulationSettings{0.001}, {b}, floating};
            int bounces = 0;
            for (int i = 0; i < 6000; ++i)
            {
                const Real vz0 = with_floor.bodies()[0].velocity.z;
                MALLOY_CHECK_TRUE(with_floor.step().ok());
                MALLOY_CHECK_TRUE(without.step().ok());
                if (vz0 < 0.0 && with_floor.bodies()[0].velocity.z > 0.0) ++bounces;

                // The rotational state is identical, floor or no floor.
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    with_floor.bodies()[0].angular_velocity,
                    without.bodies()[0].angular_velocity, 0.0));
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    with_floor.bodies()[0].orientation,
                    without.bodies()[0].orientation, 0.0));
            }
            // The floor really did act: the grounded body bounced at least once,
            // and it stayed above the floor while the floating one fell far
            // below it.
            MALLOY_CHECK_TRUE(bounces >= 1);
            MALLOY_CHECK_TRUE(with_floor.bodies()[0].position.z > -0.6);
            MALLOY_CHECK_TRUE(without.bodies()[0].position.z < -5.0);
        }

        // --- A sphere at rest on the floor stays there: it neither sinks
        //     through nor is flung off. The steady penetration is a scheme
        //     characterization, g dt^2 / (1 + e) in magnitude, so it is pinned
        //     as a bound rather than a golden number. ---
        {
            const Real g = 9.81, dt = 0.001, R = 0.3, e = 0.5;
            Rigid3DSettings st;
            st.restitution = e;
            st.gravity = Vec3{0.0, 0.0, -g};
            st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
            RigidBody3D b;
            b.radius = R;
            b.position = Vec3{0.0, 0.0, R}; // resting exactly on the floor
            Rigid3DWorld world{SimulationSettings{dt}, {b}, st};

            Real lowest = 0.0, highest = 0.0; // bottom starts at 0
            for (int i = 0; i < 20000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                const Real bottom = world.bodies()[0].position.z - R;
                lowest = std::fmin(lowest, bottom);
                highest = std::fmax(highest, bottom);
            }
            // Never sinks more than a few times the characteristic penetration,
            // and never climbs meaningfully: it rests.
            // The end-of-step position is pinned to the surface: the
            // positional correction removes the penetration each step, so the
            // bottom stays within a hair of the floor and never tunnels far
            // below it or is flung off. The band is generous (R is 0.3), so
            // this is a stability characterization, not a golden number.
            MALLOY_CHECK_TRUE(lowest > -1e-3);
            MALLOY_CHECK_TRUE(highest < 1e-3);
            // It did rest, not drift away: the final bottom is essentially zero.
            MALLOY_CHECK_TRUE(std::abs(world.bodies()[0].position.z - R) < 1e-3);
        }

        // --- Two planes make a corner, resolved in a fixed order. A sphere
        //     pushed into the corner of a floor and a wall ends up outside both,
        //     which a single-plane pass would not guarantee. ---
        {
            Rigid3DSettings st;
            st.restitution = 0.0; // land and stay
            st.gravity = Vec3{0.0, 0.0, -9.81};
            st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0},   // floor z >= 0
                         Plane3{Vec3{1.0, 0.0, 0.0}, 0.0}};  // wall  x >= 0
            RigidBody3D b;
            b.radius = 0.5;
            b.position = Vec3{0.2, 0.0, 1.0};
            b.velocity = Vec3{-1.0, 0.0, 0.0}; // drifting into the wall
            Rigid3DWorld world{SimulationSettings{0.001}, {b}, st};
            for (int i = 0; i < 3000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            const RigidBody3D& out = world.bodies()[0];
            // Outside both planes, to within the steady penetration.
            MALLOY_CHECK_TRUE(out.position.z > 0.5 - 1e-3);
            MALLOY_CHECK_TRUE(out.position.x > 0.5 - 1e-3);
        }

        // --- A zero-radius body does NOT collide: radius 0 means no shape, so
        //     it falls straight through a plane as if the plane were not there.
        //     Checked against the identical body with no ground, which must
        //     move bit for bit the same. This pins the "radius 0 = no collision"
        //     contract, which every pre-M22 body relies on. ---
        {
            Rigid3DSettings grounded;
            grounded.gravity = Vec3{0.0, 0.0, -9.81};
            grounded.restitution = 1.0;
            grounded.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
            Rigid3DSettings floating = grounded;
            floating.ground = {};

            RigidBody3D b;
            b.radius = 0.0; // no collision shape
            b.position = Vec3{0.0, 0.0, 1.0};
            b.velocity = Vec3{0.0, 0.0, -1.0};

            Rigid3DWorld with_floor{SimulationSettings{0.001}, {b}, grounded};
            Rigid3DWorld without{SimulationSettings{0.001}, {b}, floating};
            for (int i = 0; i < 3000; ++i)
            {
                MALLOY_CHECK_TRUE(with_floor.step().ok());
                MALLOY_CHECK_TRUE(without.step().ok());
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    with_floor.bodies()[0].position,
                    without.bodies()[0].position, 0.0));
            }
            // It really did fall through: far below the floor, not resting on it.
            MALLOY_CHECK_TRUE(with_floor.bodies()[0].position.z < -3.0);
        }
    }

    // --- M23: Coulomb friction for 3D contacts, the first contact here that
    //     imparts spin.
    //
    //     A normal impulse on a centred sphere has no lever arm (r x n = 0), so
    //     M22 was rotationally trivial. Friction is tangential, r x t is not
    //     zero, and the rotational effective-mass term finally does work. The
    //     headline is the same shape as the 2D rolling ratio: a sliding sphere
    //     rolls without slipping at v_roll = 5/7 v0 for a solid sphere,
    //     independent of the friction coefficient and of gravity. ---
    {
        using malloy::collide::Plane3;
        using malloy::math::Vec3;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;

        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();
        const Real R = 0.5;
        const Real m = 2.0;
        const Real inertia = 0.4 * m * R * R; // (2/5) m R^2, a solid sphere

        // --- Validation: friction is non-negative and finite. ---
        {
            MALLOY_CHECK_TRUE(Rigid3DSettings{}.is_valid()); // default 0
            Rigid3DSettings s;
            s.friction = 0.7;
            MALLOY_CHECK_TRUE(s.is_valid());
            s.friction = 2.5; // above 1 is physically real, not capped
            MALLOY_CHECK_TRUE(s.is_valid());
            for (const Real bad : {-0.1, nan3, inf3})
            {
                Rigid3DSettings b;
                b.friction = bad;
                MALLOY_CHECK_FALSE(b.is_valid());
            }
        }

        // --- Rolling without slipping, EXACT, in one step. A large coefficient
        //     arrests the slide fully, so a single impulse takes the sphere to
        //     the rolling state. The slide is DIAGONAL, so the tangent has two
        //     components and the axis it spins about is not a coordinate axis:
        //     a real check on the r x t machinery, not an axis-aligned accident.
        //
        //       v_roll = 5/7 v0,  and the contact point comes to rest. ---
        {
            const Real speed = 3.0;
            const Vec3 dir = malloy::math::normalize(Vec3{1.0, 2.0, 0.0});
            Rigid3DSettings st;
            st.restitution = 0.0;
            st.friction = 100.0; // effectively unclamped
            st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};

            RigidBody3D b;
            b.mass = m;
            b.radius = R;
            b.inertia = Vec3{inertia, inertia, inertia};
            b.position = Vec3{0.0, 0.0, R - 0.01}; // penetrating, so a contact fires
            b.velocity = dir * speed + Vec3{0.0, 0.0, -1.0}; // sliding + into floor

            Rigid3DWorld world{SimulationSettings{0.001}, {b}, st};
            MALLOY_CHECK_TRUE(world.step().ok());
            const RigidBody3D& out = world.bodies()[0];

            // Horizontal speed is exactly 5/7 of the slide speed, still along
            // the slide direction.
            const Vec3 horizontal{out.velocity.x, out.velocity.y, 0.0};
            MALLOY_CHECK_NEAR(malloy::math::length(horizontal), 5.0 * speed / 7.0, 1e-13);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::math::normalize(horizontal), dir, 1e-13));

            // The contact point is at rest: rolling without slipping. arm points
            // from the centre of mass to the sphere's lowest point.
            const Vec3 arm{0.0, 0.0, -R};
            const Vec3 omega_world = malloy::math::rotate(out.orientation, out.angular_velocity);
            const Vec3 contact_velocity = out.velocity + malloy::math::cross(omega_world, arm);
            const Vec3 contact_tangential{contact_velocity.x, contact_velocity.y, 0.0};
            MALLOY_CHECK_TRUE(malloy::math::length(contact_tangential) < 1e-13);

            // It really did spin up: friction turned a non-spinning slide into a
            // spin, which a normal contact (M22) could never do.
            MALLOY_CHECK_TRUE(malloy::math::length(out.angular_velocity) > 1.0);
            // Spin magnitude is v_roll / R.
            MALLOY_CHECK_NEAR(malloy::math::length(out.angular_velocity),
                              (5.0 * speed / 7.0) / R, 1e-13);
        }

        // --- The rolling ratio is INDEPENDENT of the friction coefficient.
        //     Two different coefficients, both large enough to reach rolling
        //     under gravity, settle at the SAME 5/7 v0. That independence is the
        //     content of the theorem, the 3D echo of the 2D rolling ratio. ---
        {
            const Real v0 = 4.0;
            const auto rolled_speed = [&](Real mu) -> Real {
                Rigid3DSettings st;
                st.restitution = 0.0;
                st.friction = mu;
                st.gravity = Vec3{0.0, 0.0, -9.81};
                st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
                RigidBody3D b;
                b.mass = m;
                b.radius = R;
                b.inertia = Vec3{inertia, inertia, inertia};
                b.position = Vec3{0.0, 0.0, R};
                b.velocity = Vec3{v0, 0.0, 0.0};
                Rigid3DWorld world{SimulationSettings{0.0005}, {b}, st};
                for (int i = 0; i < 40000; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                }
                return world.bodies()[0].velocity.x;
            };

            const Real slow = rolled_speed(0.2);
            const Real fast = rolled_speed(0.6);
            MALLOY_CHECK_NEAR(slow, 5.0 * v0 / 7.0, 1e-4);
            MALLOY_CHECK_NEAR(fast, 5.0 * v0 / 7.0, 1e-4);
            // The same to a much tighter bound than the distance from v0: the
            // coefficient sets how FAST it rolls, not the speed it rolls at.
            MALLOY_CHECK_NEAR(slow, fast, 1e-6);
            MALLOY_CHECK_TRUE(std::abs(slow - v0) > 1.0);
        }

        // --- Coulomb's cone: a small coefficient cannot arrest the slide in one
        //     step. The tangential impulse is clamped to friction * normal
        //     impulse, so the horizontal speed drops by exactly mu * (1+e) |vz|,
        //     not by the 2/7 v0 a full arrest would give. ---
        {
            const Real v0 = 3.0, vz = 2.0, mu = 0.1;
            Rigid3DSettings st;
            st.restitution = 0.0;
            st.friction = mu;
            st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
            RigidBody3D b;
            b.mass = m;
            b.radius = R;
            b.inertia = Vec3{inertia, inertia, inertia};
            b.position = Vec3{0.0, 0.0, R - 0.01};
            b.velocity = Vec3{v0, 0.0, -vz};
            Rigid3DWorld world{SimulationSettings{0.001}, {b}, st};
            MALLOY_CHECK_TRUE(world.step().ok());

            const Real drop = v0 - world.bodies()[0].velocity.x;
            // Normal impulse magnitude is (1 + e) |vz| m; clamp is mu times it;
            // the velocity drop is that over m.
            const Real clamped_drop = mu * (1.0 + 0.0) * vz; // = mu (1+e) |vz|
            MALLOY_CHECK_NEAR(drop, clamped_drop, 1e-14);
            // A full arrest would have dropped it by 2/7 v0, much more.
            MALLOY_CHECK_TRUE(2.0 * v0 / 7.0 > 3.0 * clamped_drop);
        }

        // --- Friction cannot act in mid-air: no normal impulse means no
        //     tangential impulse. A moving, spinning body high above the floor
        //     keeps its velocity and spin bit for bit. ---
        {
            Rigid3DSettings st;
            st.friction = 0.5;
            st.gravity = Vec3{0.0, 0.0, -9.81};
            st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
            RigidBody3D b;
            b.mass = m;
            b.radius = R;
            b.inertia = Vec3{inertia, inertia, inertia};
            b.position = Vec3{0.0, 0.0, 5.0}; // well clear of the floor
            const Vec3 v_start{3.0, -1.0, 0.0};
            const Vec3 w_start{0.0, 1.0, 0.0};
            b.velocity = v_start;
            b.angular_velocity = w_start;
            Rigid3DWorld world{SimulationSettings{0.001}, {b}, st};
            for (int i = 0; i < 100; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                // The horizontal velocity is untouched: only gravity acts, on vz.
                MALLOY_CHECK_NEAR(world.bodies()[0].velocity.x, v_start.x, 0.0);
                MALLOY_CHECK_NEAR(world.bodies()[0].velocity.y, v_start.y, 0.0);
                // Isotropic inertia and no torque: spin is bit-constant.
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    world.bodies()[0].angular_velocity, w_start, 0.0));
            }
        }

        // --- Friction off (the default) is M22 exactly. A sliding sphere on a
        //     frictionless floor keeps its horizontal velocity and never spins
        //     up: bit-identical to the same world with friction left at zero. ---
        {
            Rigid3DSettings rough;
            rough.restitution = 0.5;
            rough.gravity = Vec3{0.0, 0.0, -9.81};
            rough.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
            rough.friction = 0.0; // explicit zero
            Rigid3DSettings smooth = rough; // friction defaulted to 0 as well

            RigidBody3D b;
            b.mass = m;
            b.radius = R;
            b.inertia = Vec3{inertia, inertia, inertia};
            b.position = Vec3{0.0, 0.0, R};
            b.velocity = Vec3{3.0, 0.0, 0.0};

            Rigid3DWorld a{SimulationSettings{0.001}, {b}, rough};
            Rigid3DWorld c{SimulationSettings{0.001}, {b}, smooth};
            for (int i = 0; i < 4000; ++i)
            {
                MALLOY_CHECK_TRUE(a.step().ok());
                MALLOY_CHECK_TRUE(c.step().ok());
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    a.bodies()[0].velocity, c.bodies()[0].velocity, 0.0));
            }
            // Frictionless: the sphere slides forever at its start speed and
            // never spins up.
            MALLOY_CHECK_NEAR(a.bodies()[0].velocity.x, 3.0, 0.0);
            MALLOY_CHECK_TRUE(malloy::math::length(a.bodies()[0].angular_velocity) < 1e-15);
        }

        // --- The inverse inertia is rotated into the world frame, which only
        //     matters for a NON-isotropic body at a non-identity orientation. A
        //     body with distinct principal moments is tilted 90 degrees about x,
        //     so its third principal axis (moment I3) now lies along the world y
        //     axis the sphere rolls about. The rolling speed therefore uses I3,
        //     not the body-y moment I2:
        //
        //       v_roll = v0 / (1 + I3 / (m R^2)),
        //
        //     which is a different number from the I2 value. Getting it right
        //     requires the R I^-1 R^T rotation; using the body-frame inertia
        //     directly would give the I2 answer. A sphere's isotropic inertia
        //     hides this, so it needs a body whose moments differ. ---
        {
            const Real v0 = 2.2;
            const Vec3 moments{0.1, 0.2, 0.3}; // distinct on purpose
            const Real Ithree = 0.3;
            const Real roll_mass = 1.0;
            const Real roll_radius = 0.5;

            Rigid3DSettings st;
            st.restitution = 0.0;
            st.friction = 0.6;
            st.gravity = Vec3{0.0, 0.0, -9.81};
            st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};

            RigidBody3D b;
            b.mass = roll_mass;
            b.radius = roll_radius;
            b.inertia = moments;
            // 90 degrees about x sends body z onto world -y, so the world y roll
            // axis sees the I3 moment.
            b.orientation = malloy::math::from_axis_angle(Vec3{1.0, 0.0, 0.0},
                                                          1.57079632679489661923);
            b.position = Vec3{0.0, 0.0, roll_radius};
            b.velocity = Vec3{v0, 0.0, 0.0};

            Rigid3DWorld world{SimulationSettings{0.0005}, {b}, st};
            for (int i = 0; i < 60000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            const Real mR2 = roll_mass * roll_radius * roll_radius;
            const Real expected = v0 / (1.0 + Ithree / mR2);
            MALLOY_CHECK_NEAR(world.bodies()[0].velocity.x, expected, 1e-4);
            // And NOT the body-y value, which is what ignoring the rotation gives.
            const Real wrong = v0 / (1.0 + 0.2 / mR2);
            MALLOY_CHECK_TRUE(std::abs(expected - wrong) > 0.2); // they are far apart
            MALLOY_CHECK_TRUE(std::abs(world.bodies()[0].velocity.x - wrong) > 0.15);
        }

        // --- Friction with NOTHING to resist: a sphere dropped straight down,
        //     no horizontal velocity and no spin, has zero sliding at the
        //     contact, so friction must do nothing rather than divide by that
        //     zero. It bounces normally and never spins up, every step valid. ---
        {
            Rigid3DSettings st;
            st.restitution = 0.7;
            st.friction = 0.9; // high, but there is no sliding for it to act on
            st.gravity = Vec3{0.0, 0.0, -9.81};
            st.ground = {Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}};
            RigidBody3D b;
            b.mass = m;
            b.radius = R;
            b.inertia = Vec3{inertia, inertia, inertia};
            b.position = Vec3{0.0, 0.0, 2.0};
            b.velocity = Vec3{0.0, 0.0, 0.0}; // straight down under gravity
            Rigid3DWorld world{SimulationSettings{0.001}, {b}, st};
            for (int i = 0; i < 5000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok()); // no NaN from 0/0
                // No horizontal motion is ever created, and no spin.
                MALLOY_CHECK_NEAR(world.bodies()[0].velocity.x, 0.0, 0.0);
                MALLOY_CHECK_NEAR(world.bodies()[0].velocity.y, 0.0, 0.0);
                MALLOY_CHECK_TRUE(
                    malloy::math::length(world.bodies()[0].angular_velocity) < 1e-15);
            }
        }
    }

    // --- M24: sphere against sphere, two movable bodies.
    //
    //     M22 and M23 bounced a sphere off an immovable plane, which carries
    //     momentum away silently. With two movable bodies the sharp invariant
    //     is that TOTAL momentum is conserved: the impulse is equal and
    //     opposite, so what one body gains the other loses. The contact goes
    //     through the same impulse core as the ground, with the plane replaced
    //     by a second real body. ---
    {
        using malloy::collide::Sphere;
        using malloy::math::Vec3;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::total_kinetic_energy3d;
        using malloy::rigid::total_linear_momentum3d;

        // Two unit-radius spheres overlapping slightly along x (centres 1.9
        // apart), approaching head-on. A failed step leaves velocities
        // unchanged, so the exact post-collision checks below would catch it.
        const auto head_on = [](Real ma, Real mb, Real va, Real vb,
                                Real e) -> Rigid3DWorld {
            Rigid3DSettings st;
            st.restitution = e;
            RigidBody3D a;
            a.mass = ma; a.radius = 1.0; a.position = Vec3{0.0, 0.0, 0.0};
            a.velocity = Vec3{va, 0.0, 0.0};
            RigidBody3D b;
            b.mass = mb; b.radius = 1.0; b.position = Vec3{1.9, 0.0, 0.0};
            b.velocity = Vec3{vb, 0.0, 0.0};
            Rigid3DWorld w{SimulationSettings{0.001}, {a, b}, st};
            w.step();
            return w;
        };

        // Equal masses, elastic, head-on: the classic exchange of velocities.
        {
            const auto w = head_on(1.0, 1.0, 2.0, -2.0, 1.0);
            MALLOY_CHECK_NEAR(w.bodies()[0].velocity.x, -2.0, 1e-14);
            MALLOY_CHECK_NEAR(w.bodies()[1].velocity.x, 2.0, 1e-14);
        }
        // Equal masses, one at rest: the mover stops, the target leaves with it.
        {
            const auto w = head_on(1.0, 1.0, 3.0, 0.0, 1.0);
            MALLOY_CHECK_NEAR(w.bodies()[0].velocity.x, 0.0, 1e-14);
            MALLOY_CHECK_NEAR(w.bodies()[1].velocity.x, 3.0, 1e-14);
        }
        // Unequal masses, elastic: the exact 1D elastic result.
        {
            const Real ma = 3.0, mb = 1.0, va = 2.0, vb = -2.0;
            const auto w = head_on(ma, mb, va, vb, 1.0);
            const Real sum = ma + mb;
            MALLOY_CHECK_NEAR(w.bodies()[0].velocity.x,
                              ((ma - mb) * va + 2.0 * mb * vb) / sum, 1e-14); // 0
            MALLOY_CHECK_NEAR(w.bodies()[1].velocity.x,
                              ((mb - ma) * vb + 2.0 * ma * va) / sum, 1e-14); // 4
        }
        // Inelastic (e = 0), equal masses head-on at equal speed: both stop.
        {
            const auto w = head_on(1.0, 1.0, 2.0, -2.0, 0.0);
            MALLOY_CHECK_NEAR(w.bodies()[0].velocity.x, 0.0, 1e-14);
            MALLOY_CHECK_NEAR(w.bodies()[1].velocity.x, 0.0, 1e-14);
        }

        // --- Total momentum is conserved (both e), and energy never grows and
        //     here strictly falls: friction dissipates while the spheres slide
        //     across each other, and for e < 1 the bounce removes energy too.
        //     A general collision: unequal masses and radii, off-axis
        //     velocities, spin, and friction, so no term is trivially zero. ---
        for (const Real e : {1.0, 0.5})
        {
            Rigid3DSettings st;
            st.restitution = e;
            st.friction = 0.4;
            RigidBody3D a;
            a.mass = 1.5; a.radius = 0.8; a.inertia = Vec3{0.2, 0.3, 0.25};
            a.position = Vec3{0.0, 0.0, 0.0}; a.velocity = Vec3{2.0, 0.4, -0.3};
            a.angular_velocity = Vec3{0.5, -0.2, 0.1};
            RigidBody3D b;
            b.mass = 2.5; b.radius = 1.0; b.inertia = Vec3{0.4, 0.35, 0.5};
            b.position = Vec3{1.7, 0.2, 0.1}; b.velocity = Vec3{-1.0, 0.1, 0.2};
            b.angular_velocity = Vec3{-0.3, 0.4, -0.1};

            const std::vector<RigidBody3D> start = {a, b};
            const Vec3 momentum = total_linear_momentum3d(start);
            const Real energy = total_kinetic_energy3d(start);

            Rigid3DWorld world{SimulationSettings{0.001}, start, st};
            Real worst_momentum = 0.0;
            for (int i = 0; i < 400; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                worst_momentum = std::fmax(
                    worst_momentum,
                    malloy::math::length(total_linear_momentum3d(world.bodies()) - momentum));
            }
            // Equal and opposite impulse, no external force: total momentum is
            // conserved to rounding across every contact.
            MALLOY_CHECK_TRUE(worst_momentum < 1e-13);
            const Real final_energy = total_kinetic_energy3d(world.bodies());
            MALLOY_CHECK_TRUE(final_energy < energy);        // strictly lost
            MALLOY_CHECK_TRUE(final_energy > energy - 10.0); // but not absurdly
        }

        // --- Frictionless and elastic: total kinetic energy is conserved
        //     across the collision, to the scheme's contact drift (the overlap
        //     resolves over a few steps, each applying an impulse, ~1e-7 here,
        //     the same class as M14's positional-correction artefact). The exact
        //     single-impulse conservation is the head-on swap above. ---
        {
            Rigid3DSettings st;
            st.restitution = 1.0; st.friction = 0.0;
            RigidBody3D a;
            a.mass = 1.5; a.radius = 0.8; a.inertia = Vec3{0.2, 0.3, 0.25};
            a.position = Vec3{0.0, 0.0, 0.0}; a.velocity = Vec3{2.0, 0.4, -0.3};
            a.angular_velocity = Vec3{0.5, -0.2, 0.1};
            RigidBody3D b;
            b.mass = 2.5; b.radius = 1.0; b.inertia = Vec3{0.4, 0.35, 0.5};
            b.position = Vec3{1.7, 0.2, 0.1}; b.velocity = Vec3{-1.0, 0.1, 0.2};
            b.angular_velocity = Vec3{-0.3, 0.4, -0.1};

            const std::vector<RigidBody3D> start = {a, b};
            const Real energy = total_kinetic_energy3d(start);
            const Vec3 momentum = total_linear_momentum3d(start);

            Rigid3DWorld world{SimulationSettings{0.001}, start, st};
            for (int i = 0; i < 400; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_NEAR(total_kinetic_energy3d(world.bodies()), energy, 1e-6);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                total_linear_momentum3d(world.bodies()), momentum, 1e-13));
            MALLOY_CHECK_FALSE(malloy::math::approx_equal(
                world.bodies()[0].velocity, a.velocity, 1e-6));
        }

        // --- Restitution reverses the relative normal velocity by -e, off-axis
        //     so the normal is not a coordinate axis. ---
        {
            const Vec3 n = malloy::math::normalize(Vec3{1.0, 1.0, 0.0});
            Rigid3DSettings st;
            st.restitution = 0.6;
            RigidBody3D a;
            a.mass = 2.0; a.radius = 1.0; a.position = Vec3{0.0, 0.0, 0.0};
            a.velocity = n * 2.0;
            RigidBody3D b;
            b.mass = 1.0; b.radius = 1.0; b.position = n * 1.9; b.velocity = n * -1.0;
            const Real approach = malloy::math::dot(b.velocity - a.velocity, n);
            Rigid3DWorld world{SimulationSettings{0.001}, {a, b}, st};
            MALLOY_CHECK_TRUE(world.step().ok());
            const Real separate = malloy::math::dot(
                world.bodies()[1].velocity - world.bodies()[0].velocity, n);
            MALLOY_CHECK_NEAR(separate, -0.6 * approach, 1e-13);
        }

        // --- Coincident centres: fallback normal +x, no NaN, pushed apart. ---
        {
            Rigid3DSettings st;
            st.restitution = 1.0;
            RigidBody3D a;
            a.mass = 1.0; a.radius = 1.0; a.position = Vec3{0.0, 0.0, 0.0};
            RigidBody3D b;
            b.mass = 1.0; b.radius = 1.0; b.position = Vec3{0.0, 0.0, 0.0};
            Rigid3DWorld world{SimulationSettings{0.001}, {a, b}, st};
            MALLOY_CHECK_TRUE(world.step().ok());
            MALLOY_CHECK_TRUE(world.bodies()[1].position.x - world.bodies()[0].position.x >
                              1.9);
            MALLOY_CHECK_NEAR(world.bodies()[0].position.y, 0.0, 0.0);
            MALLOY_CHECK_NEAR(world.bodies()[0].position.z, 0.0, 0.0);
        }

        // --- Friction between two spheres spins both from rest. ---
        {
            Rigid3DSettings st;
            st.restitution = 0.5; st.friction = 0.5;
            RigidBody3D a;
            a.mass = 1.0; a.radius = 1.0; a.inertia = Vec3{0.4, 0.4, 0.4};
            a.position = Vec3{0.0, 0.0, 0.0}; a.velocity = Vec3{3.0, 0.0, 0.0};
            RigidBody3D b;
            b.mass = 1.0; b.radius = 1.0; b.inertia = Vec3{0.4, 0.4, 0.4};
            b.position = Vec3{1.6, 0.9, 0.0}; b.velocity = Vec3{0.0, 0.0, 0.0};
            Rigid3DWorld world{SimulationSettings{0.001}, {a, b}, st};
            MALLOY_CHECK_TRUE(world.step().ok());
            MALLOY_CHECK_TRUE(malloy::math::length(world.bodies()[0].angular_velocity) > 1e-3);
            MALLOY_CHECK_TRUE(malloy::math::length(world.bodies()[1].angular_velocity) > 1e-3);
        }

        // --- Radius zero does not collide, on EITHER side. First the striker
        //     has no shape, then the target has none: both pass through. ---
        for (int who = 0; who < 2; ++who)
        {
            Rigid3DSettings st;
            st.restitution = 1.0;
            RigidBody3D a;
            a.mass = 1.0; a.radius = (who == 0) ? 0.0 : 1.0;
            a.position = Vec3{-2.0, 0.0, 0.0}; a.velocity = Vec3{1.0, 0.0, 0.0};
            RigidBody3D b;
            b.mass = 1.0; b.radius = (who == 0) ? 1.0 : 0.0;
            b.position = Vec3{0.0, 0.0, 0.0};
            Rigid3DWorld world{SimulationSettings{0.001}, {a, b}, st};
            for (int i = 0; i < 4000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_NEAR(world.bodies()[0].velocity.x, 1.0, 0.0); // drove through
            MALLOY_CHECK_TRUE(world.bodies()[0].position.x > 1.0);
        }

        // --- The positional correction is split by inverse mass: the heavier
        //     body moves less. Two overlapping spheres at rest, mass 1 and 3,
        //     no restitution to keep it about the correction alone. The light
        //     one is pushed out three times as far. ---
        {
            Rigid3DSettings st;
            st.restitution = 0.0;
            RigidBody3D light;
            light.mass = 1.0; light.radius = 1.0; light.position = Vec3{0.0, 0.0, 0.0};
            RigidBody3D heavy;
            heavy.mass = 3.0; heavy.radius = 1.0; heavy.position = Vec3{1.5, 0.0, 0.0};
            Rigid3DWorld world{SimulationSettings{0.001}, {light, heavy}, st};
            MALLOY_CHECK_TRUE(world.step().ok());
            const Real light_move = 0.0 - world.bodies()[0].position.x;
            const Real heavy_move = world.bodies()[1].position.x - 1.5;
            MALLOY_CHECK_TRUE(light_move > 0.0);
            MALLOY_CHECK_TRUE(heavy_move > 0.0);
            MALLOY_CHECK_NEAR(light_move, 3.0 * heavy_move, 1e-12);
        }

        // --- A glancing collision between UNEQUAL spheres with friction spins
        //     both, and the induced spins pin the contact arms. Each body's arm
        //     is its radius along the normal toward the contact: a's forward,
        //     b's back. Getting b's sign wrong flips its spin; using the wrong
        //     radius for a changes its spin. The two values below are from a
        //     real run and fix both. Unequal radii on purpose. ---
        {
            Rigid3DSettings st;
            st.restitution = 0.5; st.friction = 0.5;
            RigidBody3D a;
            a.mass = 1.0; a.radius = 0.6; a.inertia = Vec3{0.2, 0.3, 0.25};
            a.position = Vec3{0.0, 0.0, 0.0}; a.velocity = Vec3{3.0, 0.0, 0.0};
            RigidBody3D b;
            b.mass = 2.0; b.radius = 1.0; b.inertia = Vec3{0.5, 0.6, 0.7};
            b.position = Vec3{1.4, 0.55, 0.0}; // centres 1.504 apart, radii sum 1.6
            Rigid3DWorld world{SimulationSettings{0.001}, {a, b}, st};
            for (int i = 0; i < 300; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            // The collision is in the xy plane, so both spin about z. Signs and
            // magnitudes are a real-run pin of the two-body arm machinery.
            MALLOY_CHECK_NEAR(world.bodies()[0].angular_velocity.z, 0.60376549, 1e-6);
            MALLOY_CHECK_NEAR(world.bodies()[1].angular_velocity.z, 0.35938422, 1e-6);
            MALLOY_CHECK_TRUE(world.bodies()[1].angular_velocity.z > 0.1); // definite sign
        }
    }

    // --- M25: 3D mass properties. A body's mass, centre of mass and principal
    //     moments of inertia COMPUTED from its geometry (solid spheres) and
    //     density, rather than typed in. The 3D sibling of the 2D
    //     mass_properties, and the piece that will let a generated shape carry
    //     real physics. Validated against closed forms, then run through M20's
    //     tumbling to check the whole chain. ---
    {
        using malloy::math::Vec3;
        using malloy::rigid::mass_properties_3d;
        using malloy::rigid::MassProperties3D;
        using malloy::rigid::rigid_body_from;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::SolidSphere;

        // --- A single solid sphere: the closed forms. Off the origin, so the
        //     centre of mass is the sphere's centre and not zero by accident. ---
        {
            const Real R = 0.5, density = 1000.0;
            const MassProperties3D mp =
                mass_properties_3d({SolidSphere{Vec3{2.0, -1.0, 3.0}, R, density}});
            const Real mass = density * (4.0 / 3.0) * pi * R * R * R;
            MALLOY_CHECK_NEAR(mp.mass, mass, 1e-9);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(mp.center_of_mass,
                                                         Vec3{2.0, -1.0, 3.0}, 1e-12));
            // Isotropic: all three moments are (2/5) m R^2.
            const Real moment = 0.4 * mass * R * R;
            MALLOY_CHECK_NEAR(mp.inertia.x, moment, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.y, moment, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.z, moment, 1e-9);
        }

        // --- A dumbbell: two equal spheres on the x axis at +/- d. The moment
        //     about the axis through both centres is just the two spheres' own,
        //     with no parallel-axis shift; the two transverse moments each gain
        //     m d^2 per sphere. So
        //       I_axial = 2 (2/5) m R^2,
        //       I_transverse = 2 (2/5) m R^2 + 2 m d^2,
        //     and the transverse pair is equal (an axisymmetric body). ---
        {
            const Real R = 0.5, density = 3.0, d = 2.0;
            const MassProperties3D mp = mass_properties_3d(
                {SolidSphere{Vec3{d, 0.0, 0.0}, R, density},
                 SolidSphere{Vec3{-d, 0.0, 0.0}, R, density}});
            const Real m = density * (4.0 / 3.0) * pi * R * R * R;
            const Real own = 0.4 * m * R * R;
            MALLOY_CHECK_NEAR(mp.mass, 2.0 * m, 1e-9);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(mp.center_of_mass, Vec3{}, 1e-12));
            // Ascending: the axial moment is smallest.
            MALLOY_CHECK_NEAR(mp.inertia.x, 2.0 * own, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.y, 2.0 * own + 2.0 * m * d * d, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.z, 2.0 * own + 2.0 * m * d * d, 1e-9);
        }

        // --- The orientation is the principal-axis frame. A dumbbell along a
        //     DIAGONAL direction has its smallest-moment axis along that
        //     diagonal, so the returned orientation rotates the body x axis onto
        //     it. This is where the diagonalization's rotation is exercised, not
        //     just its eigenvalues. ---
        {
            const Vec3 axis = malloy::math::normalize(Vec3{1.0, 1.0, 0.0});
            const Real R = 0.4, density = 2.0, d = 1.5;
            const MassProperties3D mp = mass_properties_3d(
                {SolidSphere{axis * d, R, density}, SolidSphere{axis * -d, R, density}});
            // The body x axis (the smallest principal moment) maps onto the
            // dumbbell axis, up to sign.
            const Vec3 mapped = malloy::math::rotate(mp.orientation, Vec3{1.0, 0.0, 0.0});
            const Real alignment = std::abs(malloy::math::dot(mapped, axis));
            MALLOY_CHECK_NEAR(alignment, 1.0, 1e-12);
            // And the moments are the dumbbell's, whatever direction it points.
            const Real m = density * (4.0 / 3.0) * pi * R * R * R;
            const Real own = 0.4 * m * R * R;
            MALLOY_CHECK_NEAR(mp.inertia.x, 2.0 * own, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.y, 2.0 * own + 2.0 * m * d * d, 1e-9);
        }

        // --- Different densities: a dense small core and a light large shell,
        //     centres apart. The centre of mass sits nearer the dense one, and
        //     the total mass is the sum. ---
        {
            const MassProperties3D mp = mass_properties_3d(
                {SolidSphere{Vec3{0.0, 0.0, 0.0}, 0.4, 8000.0},   // dense core
                 SolidSphere{Vec3{2.0, 0.0, 0.0}, 0.5, 1000.0}}); // light shell
            const Real dense = 8000.0 * (4.0 / 3.0) * pi * 0.4 * 0.4 * 0.4;
            const Real light = 1000.0 * (4.0 / 3.0) * pi * 0.5 * 0.5 * 0.5;
            MALLOY_CHECK_NEAR(mp.mass, dense + light, 1e-6);
            // COM on x at (dense*0 + light*2)/(dense+light).
            MALLOY_CHECK_NEAR(mp.center_of_mass.x, light * 2.0 / (dense + light), 1e-9);
            MALLOY_CHECK_TRUE(mp.center_of_mass.x < 1.0); // nearer the dense core
        }

        // --- Validation: an empty body or an invalid part yields nothing
        //     usable, rather than a meaningless number. ---
        {
            MALLOY_CHECK_NEAR(mass_properties_3d({}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                                  {SolidSphere{Vec3{}, -1.0, 1.0}}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                                  {SolidSphere{Vec3{}, 1.0, 0.0}}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                                  {SolidSphere{Vec3{0.0, nan, 0.0}, 1.0, 1.0}}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                                  {SolidSphere{Vec3{}, inf, 1.0}}).mass, 0.0, 0.0);
            // A valid part alongside an invalid one rejects the WHOLE body, not
            // just the bad part: one bad radius or density spoils it. (A single
            // bad part is also caught by the positive-mass guard, so it takes a
            // good part beside it to pin the per-part check.)
            MALLOY_CHECK_NEAR(mass_properties_3d(
                {SolidSphere{Vec3{0.0, 0.0, 0.0}, 2.0, 3.0},   // a heavy good part
                 SolidSphere{Vec3{4.0, 0.0, 0.0}, -0.5, 1.0}}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                {SolidSphere{Vec3{0.0, 0.0, 0.0}, 1.0, 1.0},
                 SolidSphere{Vec3{2.0, 0.0, 0.0}, 1.0, 0.0}}).mass, 0.0, 0.0);
        }

        // --- rigid_body_from copies the principal frame, not just the moments.
        //     A diagonal-dumbbell body has a non-identity orientation, and the
        //     built body must carry it, or it would spin about the wrong axes. ---
        {
            const Vec3 axis = malloy::math::normalize(Vec3{1.0, 1.0, 0.0});
            const MassProperties3D mp = mass_properties_3d(
                {SolidSphere{axis * 1.5, 0.4, 2.0}, SolidSphere{axis * -1.5, 0.4, 2.0}});
            MALLOY_CHECK_FALSE(malloy::math::approx_equal(mp.orientation,
                                                          malloy::math::Quat{}, 1e-6));
            const RigidBody3D body = rigid_body_from(mp);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(body.orientation,
                                                         mp.orientation, 0.0));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(body.inertia, mp.inertia, 0.0));
            MALLOY_CHECK_NEAR(body.mass, mp.mass, 0.0);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(body.position,
                                                         mp.center_of_mass, 0.0));
            MALLOY_CHECK_TRUE(body.is_valid());
        }

        // --- The whole chain, through M20's dynamics: a compound body with
        //     three DISTINCT principal moments, built from spheres, tumbles
        //     about its intermediate principal axis exactly as the intermediate
        //     axis theorem demands. This is the strongest test: the derived
        //     moments, their ordering, the orientation, and rigid_body_from all
        //     have to be right for the flip to land on the middle axis. ---
        {
            using malloy::rigid::Rigid3DWorld;
            const Real R = 0.4, density = 2.0;
            // A cross: a long pair on x, a shorter pair on y, so the three
            // moments are distinct and the principal axes are the coordinate
            // axes (the tensor is already diagonal, orientation the identity).
            const MassProperties3D mp = mass_properties_3d(
                {SolidSphere{Vec3{1.5, 0.0, 0.0}, R, density},
                 SolidSphere{Vec3{-1.5, 0.0, 0.0}, R, density},
                 SolidSphere{Vec3{0.0, 0.8, 0.0}, R, density},
                 SolidSphere{Vec3{0.0, -0.8, 0.0}, R, density}});
            // Distinct, ascending. y is the intermediate axis.
            MALLOY_CHECK_TRUE(mp.inertia.x < mp.inertia.y - 0.1);
            MALLOY_CHECK_TRUE(mp.inertia.y < mp.inertia.z - 0.1);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(mp.orientation,
                                                         malloy::math::Quat{}, 1e-12));

            const Real spin = 5.0, nudge = 1e-3;
            // Spun about y (intermediate): it flips.
            {
                RigidBody3D body = rigid_body_from(mp);
                body.angular_velocity = Vec3{nudge, spin, 0.0};
                Rigid3DWorld world{SimulationSettings{0.0005}, {body}};
                Real lowest = spin;
                for (int i = 0; i < 40000; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    lowest = std::fmin(lowest, world.bodies()[0].angular_velocity.y);
                }
                MALLOY_CHECK_TRUE(lowest < -0.9 * spin); // reversed: a real flip
            }
            // Spun about x (smallest): it does NOT flip.
            {
                RigidBody3D body = rigid_body_from(mp);
                body.angular_velocity = Vec3{spin, nudge, 0.0};
                Rigid3DWorld world{SimulationSettings{0.0005}, {body}};
                Real lowest = spin;
                for (int i = 0; i < 40000; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    lowest = std::fmin(lowest, world.bodies()[0].angular_velocity.x);
                }
                MALLOY_CHECK_TRUE(lowest > 0.9 * spin); // held: stable axis
            }
        }
    }

    std::cout << "malloy_rigid_tests passed\n";
    return 0;
}
