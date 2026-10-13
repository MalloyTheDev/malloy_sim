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
#include <utility>
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

// --- References for the 3D free rotation (issue #88), written from the
//     physics and the documented step, never from the world's own code. ---

// The angle, in radians, of the rotation taking orientation b to a: the
// distance between two orientations. Blind to the sign of a quaternion, since
// q and -q are the same rotation.
Real rotation_between(const malloy::math::Quat& a, const malloy::math::Quat& b)
{
    const malloy::math::Quat d = a * malloy::math::conjugate(b);
    return 2.0 * std::atan2(malloy::math::length(d.v), std::abs(d.w));
}

// Jacobi's elliptic functions of (u | m), 0 <= m < 1, by the arithmetic-
// geometric mean and its descending recurrence (Abramowitz and Stegun 16.4):
// a_0 = 1, b_0 = sqrt(1 - m), c_0 = sqrt(m), each step the means and the half
// difference, until c vanishes; then phi_N = 2^N a_N u and
// phi_(n-1) = (phi_n + asin(c_n sin(phi_n) / a_n)) / 2 give sn = sin(phi_0) and
// cn = cos(phi_0). `quarter` is K(m) = pi / (2 a_N), the quarter period.
struct Jacobi
{
    Real sn;
    Real cn;
    Real dn;
    Real quarter;
};

Jacobi jacobi(Real u, Real m)
{
    Real a[16];
    Real c[16];
    a[0] = 1.0;
    c[0] = std::sqrt(m);
    Real b = std::sqrt(1.0 - m);
    int n = 0;
    while (n < 15 && c[n] > 1e-17 * a[n]) // quadratic: five or six rounds
    {
        a[n + 1] = 0.5 * (a[n] + b);
        c[n + 1] = 0.5 * (a[n] - b);
        b = std::sqrt(a[n] * b);
        ++n;
    }
    Real phi = std::ldexp(a[n] * u, n);
    for (int k = n; k > 0; --k)
    {
        phi = 0.5 * (phi + std::asin(c[k] * std::sin(phi) / a[k]));
    }
    const Real sn = std::sin(phi);
    return Jacobi{sn, std::cos(phi), std::sqrt(1.0 - m * sn * sn),
                  1.5707963267948966 / a[n]};
}

// The closed-form torque-free motion of an asymmetric top (Landau and
// Lifshitz, Mechanics, section 37): moments I1 < I2 < I3 along x, y, z, a
// start (w1, 0, w3) with w1, w3 > 0 and |L|^2 > 2 T I2, so the body tumbles
// about its z axis. Then, with lambda and m below,
//
//     w(t) = (w1 cn(lambda t), B sn(lambda t), w3 dn(lambda t)),
//     B = w1 sqrt(I1 (I3 - I1) / (I2 (I3 - I2))),
//     lambda = w3 sqrt((I3 - I1) (I3 - I2) / (I1 I2)),
//     m = (I2 - I1) I1 w1^2 / ((I3 - I2) I3 w3^2),
//
// the book's amplitudes with |L|^2 and 2T substituted from the start, so no
// difference of nearly equal numbers is formed. The body frame, which the
// world stores angular velocity in. The reduced period is 4 K(m) / lambda.
malloy::math::Vec3 free_top(const malloy::math::Vec3& I, const malloy::math::Vec3& w0, Real t,
                            Real* period = nullptr)
{
    const Real lambda = w0.z * std::sqrt((I.z - I.x) * (I.z - I.y) / (I.x * I.y));
    const Real m = (I.y - I.x) * I.x * w0.x * w0.x / ((I.z - I.y) * I.z * w0.z * w0.z);
    const Real across = w0.x * std::sqrt(I.x * (I.z - I.x) / (I.y * (I.z - I.y)));
    const Jacobi f = jacobi(lambda * t, m);
    if (period != nullptr)
    {
        *period = 4.0 * f.quarter / lambda;
    }
    return malloy::math::Vec3{w0.x * f.cn, across * f.sn, w0.z * f.dn};
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
        // M40 detects contacts at the START of a step, so the first step only
        // brings them to +-0.5, exactly touching, still approaching (before
        // M40 the same step also resolved the hit); the second resolves it.
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_EQ(w.bodies()[0].position.x, -0.5);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(1.0, 0.0), 0.0);
        MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(-1.0, 0.0), 0.0);
        MALLOY_CHECK_TRUE(w.step().ok());
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
        // M40 finds the contact at the start of the NEXT step, at that pose,
        // so the bounce is resolved in the second step.
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
        for (int i = 0; i < 2; ++i) // to contact, then the bounce (M40 order)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
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

            // Rolling without slipping: the contact point is not moving.
            // Before M40 it crept backward by v (depth/2) / R (-4.905e-6 in the
            // first run): the disc sank g dt^2 each step before the correction,
            // and rolling was established about the contact point reported
            // midway through that overlap, R - depth/2 below the centre. M40's
            // rest phase holds a touching point at a normal velocity of 0 (the
            // gravity is cancelled within the step), so the disc never sinks:
            // it stays exactly at height R, the contact point is exactly R
            // below the centre, and the creep is 0. One point, so its normal
            // and friction rows are each solved exactly by one sweep, and what
            // is left is rounding of a 2 m/s velocity, far inside 1e-9.
            MALLOY_CHECK_NEAR(now.position.y, radius, 1e-12);
            const Vec2 contact{now.position.x, 0.0};
            MALLOY_CHECK_NEAR(velocity_at(now, contact).x, 0.0, 1e-9);
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

    // --- M36: box against box in the 2D rigid domain, the separating-axis
    //     contact fed through the SAME impulse core as the discs (ADR 0008),
    //     reduced to one contact at the centroid of the clipped manifold: a
    //     bounce, not a stack. Total linear momentum stays exact, and an aligned
    //     head-on hit makes no spin (the review-pass block below pins that). ---
    {
        using malloy::math::Vec2;
        const Vec2 half{0.5, 0.5};
        const auto box_body = [&](Real mass, const Vec2& position, const Vec2& velocity) {
            RigidBody2D b;
            b.mass = mass;
            b.inertia = mass * (half.x * half.x + half.y * half.y) / 3.0; // solid box
            b.half_extents = half;
            b.position = position;
            b.velocity = velocity;
            return b;
        };

        // Head-on along x, overlapping 0.1, unequal masses so the conserved
        // momentum is a nonzero value. Total linear momentum holds at every
        // restitution; the heavier body is slowed, so a real impulse fired.
        for (const Real restitution : {0.0, 0.5, 1.0})
        {
            const std::vector<RigidBody2D> start{
                box_body(2.0, Vec2{0.0, 0.0}, Vec2{1.0, 0.0}),
                box_body(1.0, Vec2{0.9, 0.0}, Vec2{-1.0, 0.0})};
            const Vec2 p0 = total_linear_momentum(start); // (1, 0)
            RigidWorld w{SimulationSettings{0.002}, start, RigidSettings{restitution}};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(w.bodies()), p0, 1e-10);
            MALLOY_CHECK_TRUE(w.bodies()[0].velocity.x < 1.0 - 1e-6);
        }

        // Frictionless restitution 1 conserves total kinetic energy (2D rotation
        // is driftless, so this is tight); restitution below 1 strictly loses
        // some.
        {
            const std::vector<RigidBody2D> start{
                box_body(1.0, Vec2{0.0, 0.0}, Vec2{1.5, 0.0}),
                box_body(1.0, Vec2{0.9, 0.0}, Vec2{-1.5, 0.0})};
            const Real k0 = total_kinetic_energy(start);
            RigidWorld elastic{SimulationSettings{0.002}, start, RigidSettings{1.0}};
            RigidWorld lossy{SimulationSettings{0.002}, start, RigidSettings{0.4}};
            for (int i = 0; i < 30; ++i)
            {
                MALLOY_CHECK_TRUE(elastic.step().ok());
                MALLOY_CHECK_TRUE(lossy.step().ok());
            }
            MALLOY_CHECK_NEAR(total_kinetic_energy(elastic.bodies()), k0, 1e-9);
            MALLOY_CHECK_TRUE(total_kinetic_energy(lossy.bodies()) < k0 - 0.1);
        }

        // A box against a DISC is deferred: the mix is skipped, so an overlapping
        // box and disc pass through each other. This guards that scope boundary
        // and will change when box/disc lands.
        {
            RigidBody2D box = box_body(1.0, Vec2{0.0, 0.0}, Vec2{1.0, 0.0});
            RigidBody2D disc;
            disc.mass = 1.0;
            disc.inertia = 1.0;
            disc.radius = 0.5; // a disc, not a box
            disc.position = Vec2{0.6, 0.0};
            disc.velocity = Vec2{-1.0, 0.0};
            RigidWorld w{SimulationSettings{0.002}, {box, disc}, RigidSettings{1.0}};
            MALLOY_CHECK_TRUE(w.step().ok());
            MALLOY_CHECK_NEAR(w.bodies()[0].velocity.x, 1.0, 1e-12);
            MALLOY_CHECK_NEAR(w.bodies()[1].velocity.x, -1.0, 1e-12);
        }

        // Two boxes clear of each other never touch: no impulse, no change.
        {
            RigidWorld w{SimulationSettings{0.002},
                         {box_body(1.0, Vec2{0.0, 0.0}, Vec2{}),
                          box_body(1.0, Vec2{3.0, 0.0}, Vec2{})},
                         RigidSettings{1.0}};
            for (int i = 0; i < 100; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_NEAR(malloy::math::length(w.bodies()[0].velocity), 0.0, 1e-12);
            MALLOY_CHECK_NEAR(malloy::math::length(w.bodies()[1].velocity), 0.0, 1e-12);
        }
    }

    // --- M37: a 2D box on a ground plane (Halfplane), the 2D sibling of M30. A
    //     box body's contact with the floor is the manifold of its penetrating
    //     corners reduced to a single centroid contact, so a flat drop stays flat
    //     (no spurious torque) and a corner landing tips. ---
    {
        using malloy::math::Vec2;
        const Vec2 half{0.5, 0.5};
        const auto box_body = [&](const Vec2& position, Real angle) {
            RigidBody2D b;
            b.mass = 1.0;
            b.inertia = 1.0 * (half.x * half.x + half.y * half.y) / 3.0; // solid box
            b.half_extents = half;
            b.position = position;
            b.angle = angle;
            return b;
        };
        const auto floor_settings = [](Real restitution, Real friction) {
            RigidSettings s;
            s.gravity = Vec2{0.0, -4.0};
            s.restitution = restitution;
            s.friction = friction;
            s.ground.push_back(malloy::collide::Halfplane{Vec2{0.0, 1.0}, 0.0});
            return s;
        };

        // A flat, symmetric drop induces NO spin and NO horizontal drift: the
        // centroid contact lies below the centre of mass, so its normal impulse
        // makes no torque. The square settles with its centre at half a side
        // height. This is the headline invariant.
        {
            RigidBody2D body = box_body(Vec2{0.0, 2.0}, 0.0);
            RigidWorld w{SimulationSettings{0.002}, {body}, floor_settings(0.5, 0.0)};
            Real worst_spin = 0.0;
            Real worst_drift = 0.0;
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                const RigidBody2D& b = w.bodies()[0];
                worst_spin = std::fmax(worst_spin, std::abs(b.angular_velocity));
                worst_drift = std::fmax(worst_drift, std::abs(b.velocity.x));
            }
            MALLOY_CHECK_TRUE(worst_spin < 1e-12);  // symmetric: no torque
            MALLOY_CHECK_TRUE(worst_drift < 1e-12); // symmetric: no drift
            const RigidBody2D& b = w.bodies()[0];
            MALLOY_CHECK_NEAR(b.position.y, 0.5, 1e-2); // resting at half height
            MALLOY_CHECK_TRUE(malloy::math::length(b.velocity) < 0.1);
        }

        // The box collider is what stops it: a box body (radius 0) rests on the
        // floor, while a body with neither a box nor a disc falls through.
        {
            RigidBody2D solid = box_body(Vec2{0.0, 1.5}, 0.0);
            RigidBody2D ghost = box_body(Vec2{0.0, 1.5}, 0.0);
            ghost.half_extents = Vec2{}; // no box, and radius 0: does not collide
            RigidWorld w{SimulationSettings{0.002}, {solid, ghost}, floor_settings(0.0, 0.0)};
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_TRUE(w.bodies()[0].position.y > 0.4);  // caught on the floor
            MALLOY_CHECK_TRUE(w.bodies()[1].position.y < -3.0); // fell through
        }

        // A corner-first landing DOES spin the box up: the contact is off the
        // centre of mass, so it has a lever arm. It tips and comes to rest above
        // the floor.
        {
            RigidBody2D body = box_body(Vec2{0.0, 2.0}, 0.4); // tilted
            RigidWorld w{SimulationSettings{0.002}, {body}, floor_settings(0.2, 0.5)};
            Real peak_spin = 0.0;
            for (int i = 0; i < 3000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                if (i > 400) // after it has reached the floor
                {
                    peak_spin = std::fmax(peak_spin, std::abs(w.bodies()[0].angular_velocity));
                }
            }
            MALLOY_CHECK_TRUE(peak_spin > 0.1);              // the corner torqued it
            MALLOY_CHECK_TRUE(w.bodies()[0].position.y > 0.0); // never sank through
        }

        // A box placed exactly at rest stays put: no drift, no spin, no sinking.
        {
            RigidBody2D body = box_body(Vec2{0.0, 0.5}, 0.0);
            RigidWorld w{SimulationSettings{0.002}, {body}, floor_settings(0.5, 0.0)};
            for (int i = 0; i < 1500; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            const RigidBody2D& b = w.bodies()[0];
            MALLOY_CHECK_NEAR(b.position.y, 0.5, 1e-2);
            MALLOY_CHECK_TRUE(malloy::math::length(b.velocity) < 0.05);
            MALLOY_CHECK_TRUE(std::abs(b.angular_velocity) < 1e-12);
        }

        // Validation: a non-positive half-extent is an invalid body, so the world
        // refuses to step.
        {
            RigidBody2D bad = box_body(Vec2{}, 0.0);
            bad.half_extents = Vec2{-1.0, 0.5};
            RigidWorld w{SimulationSettings{0.01}, {bad}, RigidSettings{}};
            MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
        }
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
        //     so the only thing making the rate vanish is the difference.
        //
        //     Since issue #88 the free rotation is the regrouped splitting of
        //     the step contract: for a sphere both axis coefficients are 0, so
        //     only the |L|^2 rotation acts, a turn of the body about omega by
        //     |omega| dt that leaves the body-frame momentum alone. The exact
        //     motion is therefore q(t) = q0 * (rotation by |omega| t about
        //     omega), which holds to rounding: at most 48 eps of angle per step
        //     (one turn, one product, one renormalization, and the reference's
        //     own rounding), a sine-and-cosine error of a few eps of each
        //     step's angle, and 16 eps for the comparison. Started tilted, so a
        //     turn composed on the wrong side (r * q for q * r) shows. ---
        {
            RigidBody3D ball;
            ball.inertia = Vec3{2.5, 2.5, 2.5};
            ball.angular_velocity = Vec3{0.7, -1.3, 0.4};
            ball.orientation = malloy::math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.8);

            const Real dt = 0.001;
            const int steps = 2000;
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            const Real rate = malloy::math::length(ball.angular_velocity);
            Rigid3DWorld world{SimulationSettings{dt}, {ball}};
            const Vec3 l0 = spin_angular_momentum(world.bodies()[0]);
            Real worst_turn = 0.0;
            Real worst_norm = 0.0;
            for (int n = 1; n <= steps; ++n)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                const Real t = static_cast<Real>(n) * dt;
                const malloy::math::Quat exact =
                    ball.orientation * malloy::math::from_axis_angle(ball.angular_velocity, rate * t);
                worst_turn = std::fmax(worst_turn,
                                       rotation_between(world.bodies()[0].orientation, exact) /
                                           ((48.0 * n + 8.0 * rate * t + 16.0) * machine_eps));
                worst_norm = std::fmax(worst_norm,
                                       std::abs(malloy::math::norm(world.bodies()[0].orientation) - 1.0));
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                world.bodies()[0].angular_velocity, ball.angular_velocity, 0.0));
            MALLOY_CHECK_TRUE(worst_turn <= 1.0);

            // A unit quaternion at every step, to the 4 eps renormalization
            // leaves: the turns are unit, and step (3) removes their rounding.
            MALLOY_CHECK_TRUE(worst_norm <= 4.0 * machine_eps);

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

            // It really turned: 3 radians, not a body sitting still.
            MALLOY_CHECK_TRUE(rate * dt * steps > 3.0);
        }

        // --- Spin about a single principal axis: exactly constant, and the
        //     orientation turns by exactly the spin times the time.
        //
        //     With omega along a principal axis, `L x omega` is zero, so the
        //     body-frame momentum is a fixed point of every piece of the free
        //     rotation: each axis turn either has angle 0 (skipped) or turns
        //     about the spin axis itself, whose two other components are 0, and
        //     the |L|^2 turn leaves the momentum alone. The angular velocity is
        //     therefore bit-identical after thousands of steps, on all three
        //     axes, the intermediate one included. ---
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

                // Every turn of every step is about the spin axis, and the
                // angles of one step add to (c_a I_a + I_a / I_i) spin dt =
                // spin dt whichever axis it is, so the exact motion is a turn
                // by spin t. Forward Euler turned the half angle by
                // atan(spin dt / 2) per step instead, 9.2e-7 rad short here.
                // The bound is test 1's: 48 eps per step, the reference's
                // sine and cosine, and the comparison.
                const Real machine_eps = std::numeric_limits<Real>::epsilon();
                const Real turned = static_cast<Real>(steps) * dt * spin;
                MALLOY_CHECK_TRUE(rotation_between(out.orientation,
                                                   malloy::math::from_axis_angle(axes[a], turned)) <=
                                  (48.0 * steps + 8.0 * turned + 16.0) * machine_eps);

                // It really did turn: 6.5 radians, not a body sitting still.
                MALLOY_CHECK_TRUE(turned > 6.0);
            }
        }

        // --- A symmetric top precesses at a rate that is exactly predictable.
        //
        //     With Ix == Iy the third equation has a zero coefficient, so wz is
        //     constant, and the first two are a plane rotation at
        //     Omega = wz (Iz - Ix) / Ix. Forward Euler stepped that rotation as
        //     a polygon, growing the magnitude by sqrt(1 + (Omega dt)^2) per
        //     step (issue #88). The regrouped flow is EXACT for an axisymmetric
        //     body: its axis coefficient for x is 0, and the one for z turns
        //     (wx, wy) by exactly Omega dt per step while R2 leaves wz's bits
        //     alone. So the transverse part is (across cos, across sin) of
        //     Omega t, to rounding: per step a few eps of angle from the turn
        //     and of magnitude from the rotation's rounding, within the scope's
        //     (64 N + 8 |w| t + 16) eps of angle, times the radius. ---
        {
            const Real dt = 0.001;
            const int steps = 2000;
            const Real across = 0.4;
            const Real along = 1.5;
            const Real machine_eps = std::numeric_limits<Real>::epsilon();

            RigidBody3D top;
            top.inertia = Vec3{1.0, 1.0, 2.0};
            top.angular_velocity = Vec3{across, 0.0, along};

            Rigid3DWorld world{SimulationSettings{dt}, {top}};
            for (int i = 0; i < steps; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }

            const Real omega = along * (2.0 - 1.0) / 1.0;
            const Real t = static_cast<Real>(steps) * dt;
            const Real turned = omega * t;
            const Vec3& w = world.bodies()[0].angular_velocity;

            const Real bound =
                (64.0 * steps + 8.0 * malloy::math::length(top.angular_velocity) * t + 16.0) * machine_eps * across;
            MALLOY_CHECK_NEAR(w.x, across * std::cos(turned), bound);
            MALLOY_CHECK_NEAR(w.y, across * std::sin(turned), bound);
            MALLOY_CHECK_NEAR(w.z, along, 0.0); // (Ix - Iy) is exactly zero

            // The precession is real, not a rounding artefact: the transverse
            // part has swung 3 radians round.
            MALLOY_CHECK_TRUE(turned > 2.9);
            MALLOY_CHECK_TRUE(w.x < 0.0); // past a quarter turn, into the far half
        }

        // --- Oblate and prolate tops whose symmetry axis is NOT z, started
        //     tilted: the whole motion in closed form (scope test 2).
        //
        //     For moments I_s about the symmetry axis e_s and I_t across it,
        //     H = |L|^2 / (2 I_t) + (1/I_s - 1/I_t) L_s^2 / 2, and the two terms
        //     commute (|L| is a Casimir and L_s is constant), so the exact flow
        //     is the product of both exact flows:
        //
        //       q(t) = q0 * turn(L_b0, |L| t / I_t) * turn(e_s, nu t),
        //       nu = (1/I_s - 1/I_t) I_s w_s,
        //
        //     with the body-frame momentum turned by -nu t about e_s: w_s
        //     constant, and the cross-axis part rotating at
        //     Omega = (I_s - I_t) / I_t w_s (the textbook rate, signed). The
        //     oblate top (2, 3.5, 2) has e_s = y as its LARGEST moment and the
        //     prolate one (0.75, 2.5, 2.5) has e_s = x as its SMALLEST, so they
        //     run through the two different axis turns of the step. The start
        //     orientation is tilted, so a turn composed on the wrong side or
        //     about the wrong axis shows. Tolerance: (64 N + 8 |w| t + 16) eps of
        //     angle, and that times |w| on the angular velocity. ---
        {
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            const Real dt = 0.002;
            const int steps = 2500;
            struct Top
            {
                Vec3 inertia;
                Vec3 axis; // the symmetry axis e_s
                Vec3 spin;
            };
            const Top tops[2] = {{Vec3{2.0, 3.5, 2.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.6, 1.1, -0.45}},
                                 {Vec3{0.75, 2.5, 2.5}, Vec3{1.0, 0.0, 0.0}, Vec3{1.3, -0.5, 0.35}}};
            for (const Top& top : tops)
            {
                RigidBody3D body;
                body.inertia = top.inertia;
                body.angular_velocity = top.spin;
                body.orientation = malloy::math::from_axis_angle(Vec3{0.4, 1.0, -0.7}, 1.1);

                const Real i_s = malloy::math::dot(top.inertia, top.axis);
                const Real i_t = top.axis.x > 0.0 ? top.inertia.y : top.inertia.x;
                const Vec3 l0{top.inertia.x * top.spin.x, top.inertia.y * top.spin.y,
                              top.inertia.z * top.spin.z};
                const Real w_s = malloy::math::dot(top.spin, top.axis);
                const Real nu = (1.0 / i_s - 1.0 / i_t) * i_s * w_s;
                const Real size = malloy::math::length(top.spin);

                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                Real worst_turn = 0.0;
                Real worst_spin = 0.0;
                for (int n = 1; n <= steps; ++n)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const RigidBody3D& out = world.bodies()[0];
                    const Real t = static_cast<Real>(n) * dt;
                    const Real tol = (64.0 * n + 8.0 * size * t + 16.0) * machine_eps;
                    const malloy::math::Quat exact =
                        body.orientation *
                        malloy::math::from_axis_angle(l0, malloy::math::length(l0) * t / i_t) *
                        malloy::math::from_axis_angle(top.axis, nu * t);
                    worst_turn = std::fmax(worst_turn, rotation_between(out.orientation, exact) / tol);
                    const Vec3 spin = malloy::math::rotate(
                        malloy::math::from_axis_angle(top.axis, -nu * t), top.spin);
                    worst_spin = std::fmax(worst_spin,
                                           malloy::math::length(out.angular_velocity - spin) / (tol * size));
                    // The symmetry-axis component is never touched (R2).
                    MALLOY_CHECK_NEAR(malloy::math::dot(out.angular_velocity, top.axis), w_s, 0.0);
                }
                MALLOY_CHECK_TRUE(worst_turn <= 1.0);
                MALLOY_CHECK_TRUE(worst_spin <= 1.0);

                // The cross-axis spin really did go round, by more than a
                // radian, and the body turned by several.
                MALLOY_CHECK_TRUE(std::abs(nu) * dt * steps > 1.0);
                MALLOY_CHECK_TRUE(malloy::math::length(l0) * dt * steps / i_t > 3.0);
            }
        }

        // --- What the free rotation does to the conserved quantities, derived
        //     rather than tolerated.
        //
        //     Before issue #88 the update was forward Euler,
        //     omega' = omega + dt I^-1 (L x omega), which grew |L|^2 by exactly
        //     dt^2 |L x omega|^2 every step and the energy with it: a tumbling
        //     body sped up without limit. Now every piece of the step is an
        //     exact flow of one term of the energy: each turns the body-frame
        //     momentum L_b by a rotation and the body by the matching turn, so
        //     |L_b| is unchanged and the WORLD momentum R L_b is unchanged, for
        //     any inertia. Those are checked every step to rounding, a few eps
        //     of |L| per step: (64 n + 16) eps |L| after n steps.
        //
        //     The energy is NOT a sub-flow invariant for three distinct
        //     moments; its bounded second-order error is measured against an
        //     independent reference in the asymmetric-top block below. ---
        {
            const Real dt = 0.01;
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            RigidBody3D body;
            body.inertia = asymmetric;
            body.angular_velocity = Vec3{0.9, 1.4, -0.6};
            body.orientation = malloy::math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.8);

            Rigid3DWorld world{SimulationSettings{dt}, {body}};
            const Vec3 world_start = spin_angular_momentum(world.bodies()[0]);
            const Real size = malloy::math::length(world_start);
            Real worst_momentum = 0.0;
            Real worst_world = 0.0;
            Real worst_frame = 0.0;
            Real worst_norm = 0.0;
            for (int n = 1; n <= 200; ++n)
            {
                MALLOY_CHECK_TRUE(world.step().ok());

                const RigidBody3D& out = world.bodies()[0];
                const Vec3 w1 = out.angular_velocity;
                const Vec3 l1{asymmetric.x * w1.x, asymmetric.y * w1.y,
                              asymmetric.z * w1.z};
                const Real bound = (64.0 * n + 16.0) * machine_eps * size;
                worst_momentum =
                    std::fmax(worst_momentum, std::abs(malloy::math::length(l1) - size) / bound);
                worst_world = std::fmax(
                    worst_world,
                    malloy::math::length(spin_angular_momentum(out) - world_start) / bound);

                // The world-frame vector has the same length as the body-frame
                // one, which is the whole reason the orientation must stay a
                // unit quaternion: a scaled one would not preserve length.
                worst_frame =
                    std::fmax(worst_frame, std::abs(malloy::math::length(
                                                        spin_angular_momentum(out)) -
                                                    malloy::math::length(l1)));
                worst_norm =
                    std::fmax(worst_norm, std::abs(malloy::math::norm(out.orientation) - 1.0));
            }

            MALLOY_CHECK_TRUE(worst_momentum <= 1.0);
            MALLOY_CHECK_TRUE(worst_world <= 1.0);
            MALLOY_CHECK_TRUE(worst_frame < 1e-13);
            MALLOY_CHECK_TRUE(worst_norm <= 4.0 * machine_eps); // scope test 7

            // The body really tumbled: the body-frame spin moved by a large
            // fraction of itself, so the conservation above is not a frozen
            // state.
            MALLOY_CHECK_FALSE(malloy::math::approx_equal(
                world.bodies()[0].angular_velocity, body.angular_velocity, 0.3));
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
        //     DIFFERENT numbers rather than the same one twice.
        //
        //     Forward Euler stepped those motions as polygons, growing every
        //     transverse amplitude by sqrt(1 + (mu dt)^2) per step (issue #88).
        //     Below, each claim is checked against what the regrouped flow does
        //     EXACTLY: its linearization about a spin is a product of its own
        //     pieces, plane rotations where a turn is about the spin axis and
        //     shears where a small component turns the large one, and each
        //     prediction is written from those pieces in closed form. ---
        {
            const Real dt = 0.001;
            const int steps = 20000;
            const Real spin = 2.0;
            const Real nudge = 1.0e-3;
            const Real rate = spin / std::sqrt(3.0);
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            const Real c_j = 1.0 / asymmetric.x - 1.0 / asymmetric.y; //  1/2
            const Real c_k = 1.0 / asymmetric.z - 1.0 / asymmetric.y; // -1/6

            // About a STABLE axis the continuum transverse amplitude is exactly
            // constant, not merely in the linearization: with |L| and the
            // energy both conserved, wy^2 + 3 wz^2 = (|L|^2 - 2T) / 2 about x and
            // wx^2 + wy^2 = (6T - |L|^2) / 2 about z. The step keeps |L|; what it
            // does to the amplitude is its linear map M of the two transverse
            // momenta per step. M is symmetric (equal diagonal, determinant 1),
            // and such a map keeps the quadratic form -m21 a^2 + m12 b^2, so the
            // amplitude runs round one fixed ellipse between two values, each
            // reached on an axis of it. Both are asserted to 1e-10 of the nudge:
            // the rounding of 20000 steps of a few eps each (2e-11), the dropped
            // O(nudge^2) terms (1e-13) and the sampling of the extremes (1e-13).
            // The departure from the continuum is 3e-8 about x and 4e-7 about z,
            // so this pins the step's own map, and forward Euler's growth, a
            // factor 1.01 to 1.04 by the end, is far outside it.

            // Spun about x, the SMALLEST moment. With (y, z) the transverse
            // momenta, A_j(dt/2) is the plane rotation by psi = (dt/2) c_j Lx,
            // R = [[c, s], [-s, c]], and A_k(dt) turns about z by dt c_k z, which
            // to first order adds -Lx dt c_k z to y, the shear [[1, tau], [0, 1]]
            // with tau = -dt c_k Lx. So M = R S R, m12 = sin 2psi + tau c^2 and
            // m21 = -sin 2psi + tau s^2. From (y, z) = (2 nudge, 0) the amplitude
            // is the nudge where z = 0, and nudge sqrt(4 (-m21) / (3 m12)) where
            // y = 0.
            {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = Vec3{spin, nudge, 0.0};

                const Real psi = 0.5 * dt * c_j * asymmetric.x * spin;
                const Real tau = -dt * c_k * asymmetric.x * spin;
                const Real c = std::cos(psi);
                const Real s = std::sin(psi);
                const Real m12 = 2.0 * s * c + tau * c * c;
                const Real m21 = -2.0 * s * c + tau * s * s;
                const Real other = nudge * std::sqrt(4.0 * -m21 / (3.0 * m12));

                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                Real worst = 0.0;
                Real least = spin;
                Real lowest_spin = spin;
                for (int i = 0; i < steps; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const Vec3& w = world.bodies()[0].angular_velocity;
                    const Real amplitude = std::sqrt(w.y * w.y + 3.0 * w.z * w.z);
                    worst = std::fmax(worst, amplitude);
                    least = std::fmin(least, amplitude);
                    lowest_spin = std::fmin(lowest_spin, w.x);
                }

                MALLOY_CHECK_NEAR(worst, std::fmax(nudge, other), 1e-10 * nudge);
                MALLOY_CHECK_NEAR(least, std::fmin(nudge, other), 1e-10 * nudge);
                // The ellipse is really not a circle at this step: the two
                // extremes differ by 3e-8 of the nudge, 300 times the tolerance.
                MALLOY_CHECK_TRUE(std::abs(other - nudge) > 1e-8 * nudge);

                // The spin itself never wavers: it is disturbed only at second
                // order in the nudge.
                MALLOY_CHECK_TRUE(std::abs(lowest_spin - spin) < 1e-5);
            }

            // Spun about z, the LARGEST moment. With (a, b) = (Lx, Ly),
            // A_j(dt/2) turns about x by (dt/2) c_j a, which to first order adds
            // Lz (dt/2) c_j a to b, the shear [[1, 0], [sigma, 1]] with
            // sigma = (dt/2) c_j Lz, and A_k(dt) is the plane rotation of (a, b)
            // by phi = dt c_k Lz. So M = S R S, m12 = sin phi and
            // m21 = 2 sigma cos phi - sin phi + sigma^2 sin phi. From
            // (a, b) = (0, 2 nudge) the amplitude is the nudge where a = 0 and
            // 2 nudge sqrt(m12 / -m21) where b = 0. The coefficients differ from
            // the case above (mu is W here, W/sqrt(3) there), so this is a
            // second check, not the same one twice.
            {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = Vec3{0.0, nudge, spin};

                const Real sigma = 0.5 * dt * c_j * asymmetric.z * spin;
                const Real phi = dt * c_k * asymmetric.z * spin;
                const Real m12 = std::sin(phi);
                const Real m21 = 2.0 * sigma * std::cos(phi) - std::sin(phi) +
                                 sigma * sigma * std::sin(phi);
                const Real other = 2.0 * nudge * std::sqrt(m12 / -m21);

                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                Real worst = 0.0;
                Real least = spin;
                Real lowest_spin = spin;
                for (int i = 0; i < steps; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const Vec3& w = world.bodies()[0].angular_velocity;
                    const Real amplitude = std::sqrt(w.x * w.x + w.y * w.y);
                    worst = std::fmax(worst, amplitude);
                    least = std::fmin(least, amplitude);
                    lowest_spin = std::fmin(lowest_spin, w.z);
                }

                MALLOY_CHECK_NEAR(worst, std::fmax(nudge, other), 1e-10 * nudge);
                MALLOY_CHECK_NEAR(least, std::fmin(nudge, other), 1e-10 * nudge);
                MALLOY_CHECK_TRUE(std::abs(other - nudge) > 1e-7 * nudge);
                MALLOY_CHECK_TRUE(std::abs(lowest_spin - spin) < 1e-5);
            }

            // Spun about y, the INTERMEDIATE moment. The body flips.
            {
                const auto flip = [&](Real h, Real& lowest, Real& fastest, Real& drift,
                                      Real& start_energy, Real& start_momentum) {
                    RigidBody3D body;
                    body.inertia = asymmetric;
                    body.angular_velocity = Vec3{nudge, spin, 0.0};
                    Rigid3DWorld world{SimulationSettings{h}, {body}};
                    start_energy = rotational_energy(world.bodies()[0]);
                    start_momentum =
                        malloy::math::length_squared(spin_angular_momentum(world.bodies()[0]));
                    lowest = spin;
                    fastest = 0.0;
                    drift = 0.0;
                    const int count = static_cast<int>(std::lround(20.0 / h));
                    for (int i = 0; i < count; ++i)
                    {
                        if (!world.step().ok())
                        {
                            return false;
                        }
                        const RigidBody3D& out = world.bodies()[0];
                        lowest = std::fmin(lowest, out.angular_velocity.y);
                        fastest = std::fmax(fastest, malloy::math::length(out.angular_velocity));
                        drift = std::fmax(drift, std::abs(rotational_energy(out) - start_energy));
                    }
                    return true;
                };
                Real lowest = 0.0, fastest = 0.0, drift = 0.0, energy = 0.0, momentum = 0.0;
                MALLOY_CHECK_TRUE(flip(dt, lowest, fastest, drift, energy, momentum));

                // It does not wobble: it reverses. A nudge of a thousandth
                // turns into a full change of sign of the spin itself, which is
                // four orders of magnitude larger than anything the two stable
                // axes did with the identical nudge.
                MALLOY_CHECK_TRUE(lowest < -0.9 * spin);

                // And it reverses along the path the conserved quantities
                // allow. With moments (1, 2, 3), for ANY state,
                //
                //     |omega|^2 = (8 T - |L|^2) / 3 - wy^2 / 3,
                //
                // so where the body passes through wy = 0 it turns as fast as
                // its energy and |L| allow. On the separatrix that is
                // |omega|^2 = 4 W^2 / 3; the start's own value, peak below, is
                // (8 T0 - |L0|^2) / 3 = 4 W^2 / 3 + nudge^2, the nudge's
                // 1e-6 above it. The step keeps |L| to
                // rounding, and the energy moves by the splitting's
                // second-order error, measured in the same run as `drift`; so
                // the fastest sample lies between
                //
                //     peak^2 - (8/3) drift - dt^2 W^4 / 36   and
                //     peak^2 + (8/3) drift,
                //
                // the low side's last term the sample nearest the crossing
                // (wy moves there at most W^2 / sqrt(3) per unit time, half a
                // step away), plus the rounding of |L|^2, (64 N + 16) eps of it
                // twice. Forward Euler, whose energy grew, ran past the peak by
                // parts in ten thousand.
                const Real peak_sq = (8.0 * energy - momentum) / 3.0;
                const Real rounding = 2.0 * (64.0 * steps + 16.0) * machine_eps * momentum / 3.0;
                MALLOY_CHECK_TRUE(fastest * fastest <= peak_sq + 8.0 / 3.0 * drift + rounding);
                MALLOY_CHECK_TRUE(fastest * fastest >=
                                  peak_sq - 8.0 / 3.0 * drift -
                                      dt * dt * spin * spin * spin * spin / 36.0 - rounding);
                MALLOY_CHECK_NEAR(peak_sq, 4.0 * spin * spin / 3.0 + nudge * nudge, 1e-14);

                // The energy's error is the second-order one, not a drift:
                // halving the step quarters it (the scope's order criterion,
                // p in [1.8, 2.2]), and it is far below a part in a million.
                Real lowest2 = 0.0, fastest2 = 0.0, drift2 = 0.0, energy2 = 0.0, momentum2 = 0.0;
                MALLOY_CHECK_TRUE(flip(0.5 * dt, lowest2, fastest2, drift2, energy2, momentum2));
                const Real order = std::log2(drift / drift2);
                MALLOY_CHECK_TRUE(order >= 1.8 && order <= 2.2);
                MALLOY_CHECK_TRUE(drift < 1e-6 * energy);
                MALLOY_CHECK_TRUE(lowest2 < -0.9 * spin);
            }

            // The exponential growth rate itself, against the exact discrete
            // linearization rather than its continuum limit. With (a, c) =
            // (Lx, Lz) small and Ly = I2 W, A_j(dt/2) adds -Ly (dt/2) c_j a to c
            // and A_k(dt) adds Ly dt c_k c to a: three shears,
            // M = S1 S2 S1 with S1 = [[1, 0], [-s1, 1]], S2 = [[1, s2], [0, 1]],
            // s1 = (dt/2) c_j Ly and s2 = dt c_k Ly. Its diagonal is
            // 1 - s1 s2 = 1 + (W dt)^2 / 6 = cosh(theta), so its eigenvalues are
            // e^(+/-theta) with theta = 2 asinh(rate dt / 2), and from (tiny, 0)
            // the growing component is exactly
            //
            //     tiny cosh(n theta),
            //
            // where the continuum has tiny cosh(rate t): theta is rate dt to
            // second order. A small nudge keeps the dropped O(nudge^2) terms
            // near 1e-17, so this is asserted to fifteen significant figures,
            // as the forward-Euler law it replaces was.
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
                const Real theta = 2.0 * std::asinh(0.5 * rate * dt);
                const Real predicted = tiny * std::cosh(short_run * theta);
                MALLOY_CHECK_NEAR(world.bodies()[0].angular_velocity.x, predicted, 1e-15);

                // The continuum value differs by 6e-13, 600 times the
                // tolerance, so the check sees the step and not just the physics.
                MALLOY_CHECK_TRUE(std::abs(tiny * std::cosh(rate * short_run * dt) - predicted) >
                                  5e-13);

                // Five times its starting size in two seconds, and still only
                // five parts per million of the spin: this is the early phase,
                // where the linearization is what governs.
                MALLOY_CHECK_TRUE(predicted > 5.0 * tiny);
            }
        }

        // --- The ORDER of the step under a torque, against an analytic law
        //     (scope test 5): one full kick, then the free rotation.
        //
        //     A sphere of moment I under a world torque tau parallel to its
        //     spin turns about that fixed axis by the continuum angle
        //
        //       theta(t) = (|L0| t + |tau| t^2 / 2) / I.
        //
        //     A step kicks the spin by tau dt / I and then turns the body by dt
        //     times the KICKED spin, so after n steps the angle is the sum of
        //     dt w_m over m = 1..n, w_m = w0 + m dt |tau| / I, which is
        //
        //       theta_n = theta(t_n) + |tau| dt t_n / (2 I),
        //
        //     ahead of the continuum by a first-order offset. Half a kick on
        //     each side of the turn would put it exactly on theta(t_n), and the
        //     turn before the kick the same offset behind, so the law separates
        //     the three orders.
        //
        //     Dyadic data make every kick exact: I = 2, tau = 0.25 e3,
        //     w0 = 0.25 e3, dt = 1/64, so tau dt / I = 2^-9 and the spin is bit
        //     for bit 0.25 + n / 512. The orientation is checked to the sphere
        //     test's angle tolerance. The start is turned about e3 (not the
        //     identity), which keeps the body-frame torque bit-equal to the
        //     world one; the torque's frame is tested on a tilted body in M21.
        //     Halving dt halves the offset (first order under forcing), and
        //     with no torque it is zero (the free flow is exact for a sphere). ---
        {
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            const Real moment = 2.0;
            const Vec3 axis{0.0, 0.0, 1.0};
            const malloy::math::Quat start = malloy::math::from_axis_angle(axis, 0.7);
            const auto run = [&](Real dt, Real torque, int steps, Real& offset) {
                RigidBody3D ball;
                ball.inertia = Vec3{moment, moment, moment};
                ball.angular_velocity = Vec3{0.0, 0.0, 0.25};
                ball.orientation = start;
                malloy::rigid::Rigid3DSettings settings;
                settings.torque = Vec3{0.0, 0.0, torque};
                Rigid3DWorld world{SimulationSettings{dt}, {ball}, settings};
                Real worst = 0.0;
                for (int n = 1; n <= steps; ++n)
                {
                    if (!world.step().ok())
                    {
                        return 2.0;
                    }
                    const RigidBody3D& out = world.bodies()[0];
                    const Real t = static_cast<Real>(n) * dt;
                    const Real spin = 0.25 + static_cast<Real>(n) * dt * torque / moment;
                    if (!(out.angular_velocity.x == 0.0 && out.angular_velocity.y == 0.0 &&
                          out.angular_velocity.z == spin))
                    {
                        return 2.0; // the kicks are exact: the spin is bit for bit
                    }
                    const Real continuum = (moment * 0.25 * t + torque * t * t / 2.0) / moment;
                    const Real kd = continuum + torque * dt * t / (2.0 * moment);
                    worst = std::fmax(worst, rotation_between(out.orientation,
                                                              start * malloy::math::from_axis_angle(axis, kd)) /
                                                 ((48.0 * n + 8.0 * spin * t + 16.0) * machine_eps));
                    offset = rotation_between(out.orientation,
                                              start * malloy::math::from_axis_angle(axis, continuum));
                }
                return worst;
            };
            Real offset = 0.0;
            Real halved = 0.0;
            Real none = 0.0;
            MALLOY_CHECK_TRUE(run(1.0 / 64.0, 0.25, 640, offset) <= 1.0);
            MALLOY_CHECK_TRUE(run(1.0 / 128.0, 0.25, 1280, halved) <= 1.0);
            MALLOY_CHECK_TRUE(run(1.0 / 64.0, 0.0, 640, none) <= 1.0);

            // At t = 10: an offset of 0.25 (1/64) 10 / 4 rad, then half that,
            // then none, each to the angle tolerance.
            const Real tol = (48.0 * 1280.0 + 8.0 * 1.5 * 10.0 + 16.0) * machine_eps;
            MALLOY_CHECK_NEAR(offset, 0.25 / 64.0 * 10.0 / 4.0, tol);
            MALLOY_CHECK_NEAR(halved, 0.5 * offset, tol);
            MALLOY_CHECK_NEAR(none, 0.0, tol);
        }

        // --- A free ASYMMETRIC top against its closed form (scope test 6).
        //
        //     Spheres and symmetric tops are exact under the free rotation; for
        //     three distinct moments it is second order. That is measured
        //     against an independent reference, the elliptic-function solution
        //     (`free_top` above), not against the method itself: a
        //     self-convergence test passes an error made the same way at every
        //     step size.
        //
        //       1. Order: the body-frame spin's error at t = 10 s at
        //          dt = 1/30, 1/60 and 1/120 falls by four per halving,
        //          p in [1.8, 2.2] on both pairs.
        //       2. Energy envelope: e(h), the largest |T_n - T_0| / T_0 over
        //          12 s (more than two reduced periods), is second order,
        //          p in [1.8, 2.2], and C = max e(h) / h^2 is its constant.
        //       3. No secular drift over 48 reduced periods at dt = 1/30 and
        //          1/60. The late maximum stays within 2 C h^2: a modified
        //          energy exists for small h, so the error is periodic and one
        //          period already samples its maximum. And with s the least-
        //          squares slope of the signed error against time,
        //          D = |s| T_long / (C h^2) <= 0.25: a bounded periodic error
        //          of amplitude C h^2 gives D <= 12 P / T_long = 0.25 at
        //          T_long = 48 P, and a drift of a quarter of the envelope over
        //          the run fails.
        //       4. |L| and the world momentum to rounding throughout, and a
        //          unit quaternion at every step.
        //
        //     Moments (1, 2, 3) started at (0.6, 0, 1.2): |L|^2 = 13.32 is above
        //     2 T I2 = 9.36, so it tumbles about z; m = 1/12, lambda = 1.2, and
        //     the reduced period P = 4 K(m) / lambda is 5.35 s. ---
        {
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            const Vec3 w0{0.6, 0.0, 1.2};
            Real period = 0.0;
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(free_top(asymmetric, w0, 0.0, &period), w0,
                                                         1e-15));
            MALLOY_CHECK_TRUE(period > 5.3 && period < 5.4);

            // 1. Order against the closed form.
            Real error[3] = {0.0, 0.0, 0.0};
            const Real rates[3] = {30.0, 60.0, 120.0};
            for (int k = 0; k < 3; ++k)
            {
                const Real dt = 1.0 / rates[k];
                const int steps = static_cast<int>(10.0 * rates[k]);
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = w0;
                body.orientation = malloy::math::from_axis_angle(Vec3{0.3, -1.0, 0.6}, 2.0);
                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                for (int i = 0; i < steps; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                }
                error[k] = malloy::math::length(world.bodies()[0].angular_velocity -
                                                free_top(asymmetric, w0, steps * dt));
            }
            const Real p1 = std::log2(error[0] / error[1]);
            const Real p2 = std::log2(error[1] / error[2]);
            MALLOY_CHECK_TRUE(p1 >= 1.8 && p1 <= 2.2);
            MALLOY_CHECK_TRUE(p2 >= 1.8 && p2 <= 2.2);
            MALLOY_CHECK_TRUE(error[0] > 1e-6); // a real error, not rounding

            // 2. and 3. The energy, short and long, with |L|, L and |q|.
            const auto energy_run = [&](Real dt, Real span, Real& peak, Real& slope, Real& worst_l,
                                        Real& worst_q) {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = w0;
                body.orientation = malloy::math::from_axis_angle(Vec3{0.3, -1.0, 0.6}, 2.0);
                Rigid3DWorld world{SimulationSettings{dt}, {body}};
                const Real t0 = rotational_energy(world.bodies()[0]);
                const Vec3 l0 = spin_angular_momentum(world.bodies()[0]);
                const int steps = static_cast<int>(std::ceil(span / dt));
                Real st = 0.0, se = 0.0, stt = 0.0, ste = 0.0;
                peak = 0.0;
                worst_l = 0.0;
                worst_q = 0.0;
                for (int n = 1; n <= steps; ++n)
                {
                    if (!world.step().ok())
                    {
                        return false;
                    }
                    const RigidBody3D& out = world.bodies()[0];
                    const Real t = n * dt;
                    const Real e = (rotational_energy(out) - t0) / t0;
                    peak = std::fmax(peak, std::abs(e));
                    st += t;
                    se += e;
                    stt += t * t;
                    ste += t * e;
                    worst_l = std::fmax(worst_l, malloy::math::length(spin_angular_momentum(out) - l0) /
                                                     ((64.0 * n + 16.0) * machine_eps * malloy::math::length(l0)));
                    worst_q = std::fmax(worst_q, std::abs(malloy::math::norm(out.orientation) - 1.0));
                }
                const Real count = static_cast<Real>(steps);
                slope = (ste - st * se / count) / (stt - st * st / count);
                return true;
            };
            Real e30 = 0.0, e60 = 0.0, unused = 0.0, worst_l = 0.0, worst_q = 0.0;
            MALLOY_CHECK_TRUE(energy_run(1.0 / 30.0, 12.0, e30, unused, worst_l, worst_q));
            MALLOY_CHECK_TRUE(energy_run(1.0 / 60.0, 12.0, e60, unused, worst_l, worst_q));
            const Real p_energy = std::log2(e30 / e60);
            MALLOY_CHECK_TRUE(p_energy >= 1.8 && p_energy <= 2.2);
            const Real envelope = std::fmax(e30 * 30.0 * 30.0, e60 * 60.0 * 60.0); // C

            for (const Real rate : {30.0, 60.0})
            {
                const Real h = 1.0 / rate;
                const Real span = 48.0 * period;
                Real late = 0.0, slope = 0.0;
                MALLOY_CHECK_TRUE(energy_run(h, span, late, slope, worst_l, worst_q));
                MALLOY_CHECK_TRUE(late <= 2.0 * envelope * h * h);
                MALLOY_CHECK_TRUE(std::abs(slope) * span / (envelope * h * h) <= 0.25);
                MALLOY_CHECK_TRUE(worst_l <= 1.0);
                MALLOY_CHECK_TRUE(worst_q <= 4.0 * machine_eps);
            }
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
        //     The spin is constant exactly when u = L x omega is zero, that is
        //     when I omega is parallel to omega, and there the step keeps it bit
        //     for bit (R1 to R3: every axis turn has angle 0 or turns two zero
        //     components). For three distinct moments the only such directions
        //     are the three principal axes. For an AXISYMMETRIC body two moments
        //     are equal, that eigenspace is a plane, and every direction in it
        //     qualifies: not just the two principal axes lying in it.
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

            // Tilt the same body out of that plane and it precesses, so the
            // exactness above is the eigenspace and not the inertia. The flow
            // is exact for an axisymmetric body: wz is kept bit for bit, and
            // (wx, wy) turn about z by Omega t, Omega = (Iz - Ix) / Ix wz, here
            // 2.5 rad in 3 s (forward Euler grew the energy instead). The
            // energy is kept to rounding, (64 N + 16) eps of itself.
            RigidBody3D tilted = top;
            tilted.angular_velocity = Vec3{0.6, -0.8, 0.5};
            Rigid3DWorld drifting{SimulationSettings{0.001}, {tilted}};
            const Real tilted_start = rotational_energy(drifting.bodies()[0]);
            for (int i = 0; i < 3000; ++i)
            {
                MALLOY_CHECK_TRUE(drifting.step().ok());
            }
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            const Vec3& w = drifting.bodies()[0].angular_velocity;
            const Real omega = (4.0 - 1.5) / 1.5 * 0.5;
            const Real turned = omega * 3.0;
            const Vec3 exact{0.6 * std::cos(turned) + 0.8 * std::sin(turned),
                             0.6 * std::sin(turned) - 0.8 * std::cos(turned), 0.5};
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(w, exact, (64.0 * 3000 + 16.0) * machine_eps));
            MALLOY_CHECK_NEAR(w.z, 0.5, 0.0);
            MALLOY_CHECK_NEAR(rotational_energy(drifting.bodies()[0]), tilted_start,
                              (64.0 * 3000 + 16.0) * machine_eps * tilted_start);
            MALLOY_CHECK_TRUE(turned > 2.0); // it really went round
        }

        // --- The world frame is not the body frame, and that is the whole
        //     reason `angular_momentum` rotates.
        //
        //     Under torque-free rotation the WORLD-frame angular momentum is
        //     fixed while the BODY-frame vector I omega tumbles with the body.
        //     Measuring both at two step sizes separates the two claims from
        //     each other: the world-frame deviation is the integrator's, the
        //     body-frame swing is the physics and does not depend on dt.
        //
        //     Since issue #88 every piece of the free rotation keeps the world
        //     momentum exactly, so the deviation is rounding alone at BOTH step
        //     sizes, (64 N + 16) eps |L| after N steps, where forward Euler's
        //     was a first-order integration error. A bound this tight
        //     is not satisfied by an implementation that never rotated
        //     anything: its world vector would swing with the body's. ---
        {
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            Real deviation[2] = {0.0, 0.0};
            Real tumble[2] = {0.0, 0.0};
            Real allowed[2] = {0.0, 0.0};
            const Real steps[2] = {1e-3, 1e-4};
            for (int k = 0; k < 2; ++k)
            {
                RigidBody3D body;
                body.inertia = asymmetric;
                body.angular_velocity = Vec3{0.9, 1.4, -0.6};
                body.orientation = malloy::math::from_axis_angle(Vec3{-0.5, 1.0, 2.0}, 1.3);

                Rigid3DWorld world{SimulationSettings{steps[k]}, {body}};
                const Vec3 world_start = spin_angular_momentum(world.bodies()[0]);
                const Vec3 body_start{asymmetric.x * body.angular_velocity.x,
                                      asymmetric.y * body.angular_velocity.y,
                                      asymmetric.z * body.angular_velocity.z};

                const int count = static_cast<int>(2.0 / steps[k]);
                allowed[k] = (64.0 * count + 16.0) * machine_eps * malloy::math::length(world_start);
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

            // Rounding only, at both step sizes.
            MALLOY_CHECK_TRUE(deviation[0] <= allowed[0]);
            MALLOY_CHECK_TRUE(deviation[1] <= allowed[1]);

            // The body really is tumbling, by more than |L| itself, and by the
            // same amount at both step sizes because it is not an error term.
            MALLOY_CHECK_TRUE(tumble[0] > 5.0);
            MALLOY_CHECK_TRUE(tumble[1] > 5.0);
            MALLOY_CHECK_NEAR(tumble[0], tumble[1], 0.01);

            // Ten orders of magnitude between the two, which is what makes
            // returning the unrotated vector from `angular_momentum` a failure
            // rather than a rounding difference.
            MALLOY_CHECK_TRUE(tumble[0] > 1e10 * deviation[0]);
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
        //     finite: it enters the momentum law as dt^2 |torque|^2, so one
        //     whose square overflows is refused up front, the bound softening
        //     carries. ---
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
        //     The gyroscopic term is zero because omega stays along the axis,
        //     and every turn of the free rotation is about it. What is left is
        //     I wx' = T, integrated exactly by the step's torque kick because
        //     the right side is constant:
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
        //     A step kicks the world momentum by exactly torque dt, and the free
        //     rotation that follows keeps the world momentum and |L| (issue
        //     #88), so for ANY inertia
        //
        //       L_world(n+1) = L_world(n) + torque dt,
        //       |L(n+1)|^2   = |L(n)|^2 + 2 dt (L . T_body) + dt^2 |torque|^2,
        //
        //     with L . T_body = L_world . torque. Forward Euler's last term was
        //     dt^2 |L x omega + T_body|^2, the gyroscopic part a growth of its
        //     own. The kick changes the energy by exactly
        //
        //       dt (omega . T_body) + dt^2 (T_body . I^-1 T_body) / 2,
        //
        //     and the free rotation keeps the energy of a sphere or an
        //     axisymmetric body, so for those that is the whole step; for three
        //     distinct moments the rotation adds its own second-order error,
        //     measured against the closed form in M20. A tumbling body, a tilted
        //     start and an off-axis torque, so no term is accidentally zero. The
        //     per-step laws are checked to a few ulp of |L|^2 (about 12, one ulp
        //     1.8e-15). The world law accumulated, L_world(t_n) =
        //     L_world(0) + torque t_n, is checked within (64 n + 16) eps of the
        //     largest |L| (scope test 4): from a tilted start, so a torque
        //     applied in the body frame fails it. ---
        {
            const Vec3 torque{0.7, -0.4, 0.9};
            const Real dt = 0.01;
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            for (const Vec3& inertia : {asymmetric, Vec3{1.5, 2.5, 1.5}})
            {
                const bool axisymmetric = inertia.x == inertia.z;
                RigidBody3D body;
                body.inertia = inertia;
                body.angular_velocity = Vec3{0.9, 1.4, -0.6};
                body.orientation = malloy::math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.8);

                Rigid3DWorld world{SimulationSettings{dt}, {body}, Rigid3DSettings{torque}};
                const Vec3 world_start = spin_angular_momentum(world.bodies()[0]);
                Real largest = malloy::math::length(world_start);
                Real worst_momentum = 0.0;
                Real worst_energy = 0.0;
                Real worst_world = 0.0;
                for (int n = 1; n <= 200; ++n)
                {
                    const RigidBody3D& s = world.bodies()[0];
                    const Vec3 tb =
                        malloy::math::rotate(malloy::math::conjugate(s.orientation), torque);
                    const Vec3 w = s.angular_velocity;
                    const Vec3 l{inertia.x * w.x, inertia.y * w.y, inertia.z * w.z};
                    const Vec3 itb{tb.x / inertia.x, tb.y / inertia.y, tb.z / inertia.z};
                    const Real momentum_before = malloy::math::length_squared(l);
                    const Real energy_before = rotational_energy(s);

                    MALLOY_CHECK_TRUE(world.step().ok());

                    const RigidBody3D& s1 = world.bodies()[0];
                    const Vec3 w1 = s1.angular_velocity;
                    const Vec3 l1{inertia.x * w1.x, inertia.y * w1.y, inertia.z * w1.z};

                    const Real pred_momentum = momentum_before + 2.0 * dt * malloy::math::dot(l, tb) +
                                               dt * dt * malloy::math::length_squared(torque);
                    const Real pred_energy = energy_before + dt * malloy::math::dot(w, tb) +
                                             0.5 * dt * dt * malloy::math::dot(tb, itb);
                    worst_momentum = std::fmax(
                        worst_momentum, std::abs(malloy::math::length_squared(l1) - pred_momentum));
                    worst_energy =
                        std::fmax(worst_energy, std::abs(rotational_energy(s1) - pred_energy));

                    const Vec3 lw = spin_angular_momentum(s1);
                    largest = std::fmax(largest, malloy::math::length(lw));
                    worst_world = std::fmax(
                        worst_world,
                        malloy::math::length(lw - (world_start + torque * (static_cast<Real>(n) * dt))) /
                            ((64.0 * n + 16.0) * machine_eps * largest));
                }
                MALLOY_CHECK_TRUE(worst_momentum < 1e-13);
                MALLOY_CHECK_TRUE(worst_world <= 1.0);
                if (axisymmetric)
                {
                    MALLOY_CHECK_TRUE(worst_energy < 1e-13);
                }

                // The torque did net work here, so the energy did not merely
                // wander at second order.
                MALLOY_CHECK_TRUE(rotational_energy(world.bodies()[0]) !=
                                  rotational_energy(body));
            }
        }

        // --- The world-frame linear law, for a body that is NOT spinning up
        //     cleanly: a fast symmetric top with a torque across its spin.
        //
        //     dL/dt = torque holds in the world frame no matter what the body
        //     does, so L_world(t) = L_world(0) + torque t is the continuum law.
        //     The step reproduces it exactly at every step boundary (the kick
        //     adds torque dt and the rotation keeps the world momentum), so the
        //     deviation is rounding at BOTH step sizes, (64 N + 16) eps of the
        //     largest |L|, where forward Euler's was a first-order error. ---
        {
            const Real T = 0.3;
            const Real spin = 8.0;
            const Real machine_eps = std::numeric_limits<Real>::epsilon();
            Real deviation[2] = {0.0, 0.0};
            Real allowed[2] = {0.0, 0.0};
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
                Real largest = malloy::math::length(start);
                for (int n = 1; n <= count; ++n)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    const Vec3 lw = total_angular_momentum3d(world.bodies());
                    const Vec3 exact{start.x + T * static_cast<Real>(n) * steps[k],
                                     start.y, start.z};
                    deviation[k] = std::fmax(deviation[k], malloy::math::length(lw - exact));
                    tilt[k] = std::fmax(tilt[k], std::abs(lw.x));
                    largest = std::fmax(largest, malloy::math::length(lw));
                }
                allowed[k] = (64.0 * count + 16.0) * machine_eps * largest;
            }

            MALLOY_CHECK_TRUE(deviation[0] <= allowed[0]);
            MALLOY_CHECK_TRUE(deviation[1] <= allowed[1]);

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

        // --- Coincident centres: fallback normal +x, no NaN, pushed apart.
        //     M40's projection removes the depth (2, both radii) a fifth at a
        //     time, at most 0.2 of the radius per push and 3 pushes per step,
        //     so the first step separates them by exactly 3 x 0.2 = 0.6, split
        //     equally (equal masses), and later steps take the rest: capped
        //     while the depth exceeds 1, then shrinking by 0.8 per push, until
        //     only the allowed 1e-6 of the radius remains. Positions only: at
        //     rest, they stay at rest (before M40 one step pushed them more
        //     than 1.9 apart). ---
        {
            Rigid3DSettings st;
            st.restitution = 1.0;
            RigidBody3D a;
            a.mass = 1.0; a.radius = 1.0; a.position = Vec3{0.0, 0.0, 0.0};
            RigidBody3D b;
            b.mass = 1.0; b.radius = 1.0; b.position = Vec3{0.0, 0.0, 0.0};
            Rigid3DWorld world{SimulationSettings{0.001}, {a, b}, st};
            MALLOY_CHECK_TRUE(world.step().ok());
            MALLOY_CHECK_NEAR(world.bodies()[0].position.x, -0.3, 1e-15);
            MALLOY_CHECK_NEAR(world.bodies()[1].position.x, 0.3, 1e-15);
            for (int i = 1; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            // The depth left after step 3 (0.41) shrinks by 0.8^3 a step, so
            // by step 60 only the slop remains, to rounding.
            MALLOY_CHECK_NEAR(world.bodies()[1].position.x - world.bodies()[0].position.x,
                              2.0 - 1e-6, 1e-12);
            MALLOY_CHECK_NEAR(world.bodies()[0].position.y, 0.0, 0.0);
            MALLOY_CHECK_NEAR(world.bodies()[0].position.z, 0.0, 0.0);
            MALLOY_CHECK_NEAR(malloy::math::length(world.bodies()[0].velocity), 0.0, 0.0);
            MALLOY_CHECK_NEAR(malloy::math::length(world.bodies()[1].velocity), 0.0, 0.0);
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
        //     both, and the induced spins pin the contact arms. Both arms reach
        //     the one contact point, midway between the surfaces (issue #35: they
        //     used to be each radius along the normal, two different points when
        //     the spheres overlap). Getting b's sign wrong flips its spin; a
        //     wrong arm for a changes its spin. The two values below are from a
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
            // magnitudes are a real-run pin of the two-body arm machinery,
            // cross-checked against an independent closed-form one-point
            // impact at the START pose (K P = dv for the 3x3 effective-mass
            // matrix, restitution 0.5 on the normal; it sticks, |P_t| / P_n =
            // 0.098 < 0.5), which agrees to ten digits. Before M40 they were
            // 0.60457214 and 0.37275515: the contact was resolved at the pose
            // one step of motion later, which moved the arms.
            MALLOY_CHECK_NEAR(world.bodies()[0].angular_velocity.z, 0.6034797849, 1e-6);
            MALLOY_CHECK_NEAR(world.bodies()[1].angular_velocity.z, 0.3716858050, 1e-6);
            MALLOY_CHECK_TRUE(world.bodies()[1].angular_velocity.z > 0.1); // definite sign
        }

        // --- Issue #35: the contact impulse conserves total angular momentum,
        //     however deep the overlap. Equal and opposite impulses at ONE point
        //     form no couple; at two points (each radius along the normal) they
        //     did, and friction's tangential part changed L by about the depth
        //     times the impulse: 0.0255 at depth 0.1. One step at a tiny dt, so
        //     the motion is negligible, with L taken at the ORIGINAL positions
        //     so the positional correction (a documented, separate artefact) is
        //     left out and only the impulse is measured. ---
        {
            Rigid3DSettings st;
            st.restitution = 0.5;
            st.friction = 0.5;
            for (const Real depth : {0.0, 0.1, 0.2})
            {
                RigidBody3D a;
                a.mass = 1.0; a.radius = 0.6; a.inertia = Vec3{0.2, 0.3, 0.25};
                a.velocity = Vec3{3.0, 0.0, 0.0};
                RigidBody3D b;
                b.mass = 2.0; b.radius = 1.0; b.inertia = Vec3{0.5, 0.6, 0.7};
                b.position = malloy::math::normalize(Vec3{1.0, 0.4, 0.0}) * (1.6 - depth);
                Rigid3DWorld world{SimulationSettings{1e-12}, {a, b}, st};
                const Vec3 before = malloy::rigid::total_angular_momentum3d(world.bodies());
                MALLOY_CHECK_TRUE(world.step().ok());
                std::vector<RigidBody3D> after = world.bodies();
                after[0].position = a.position;
                after[1].position = b.position;
                const Vec3 change = malloy::rigid::total_angular_momentum3d(after) - before;
                MALLOY_CHECK_NEAR(malloy::math::length(change), 0.0, 1e-9);
                // and the friction really acted: a picks up spin
                MALLOY_CHECK_TRUE(malloy::math::length(after[0].angular_velocity) > 0.5);
            }
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
            MALLOY_CHECK_NEAR(mass_properties_3d(std::vector<SolidSphere>{}).mass,
                              0.0, 0.0);
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

    // --- M26: box mass properties and combine. A solid box is the second mass
    //     primitive (three distinct moments and a real orientation, where a
    //     sphere had one moment and no orientation), and `combine` assembles a
    //     compound body from primitives of any kind. Both are validated against
    //     closed forms; combine is then pinned against the M25 list path, which
    //     never calls it, and the box is run through M20's dynamics. ---
    {
        using malloy::math::approx_equal;
        using malloy::math::diagonal3;
        using malloy::math::from_axis_angle;
        using malloy::math::Mat3;
        using malloy::math::Quat;
        using malloy::math::to_mat3;
        using malloy::math::transpose;
        using malloy::math::Vec3;
        using malloy::rigid::combine;
        using malloy::rigid::mass_properties_3d;
        using malloy::rigid::MassProperties3D;
        using malloy::rigid::rigid_body_from;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::SolidBox;
        using malloy::rigid::SolidSphere;

        // The full inertia tensor in the lab frame, reconstructed from a
        // result's principal moments and frame: R diag(I) R^T. Two results
        // describe the SAME body exactly when their masses, centres and these
        // tensors agree, no matter how the eigenvalues happened to be ordered,
        // so this is what commutativity and consistency should compare.
        const auto tensor_of = [](const MassProperties3D& mp) {
            const Mat3 rotation = to_mat3(mp.orientation);
            return rotation * diagonal3(mp.inertia) * transpose(rotation);
        };

        // --- A single axis-aligned box: the closed forms. Half-extents made
        //     strictly decreasing (hx > hy > hz), so the box-frame tensor is
        //     ALREADY ascending on the diagonal; the principal frame is then the
        //     identity and each moment lands on a known axis. Off the origin, so
        //     the centre of mass is the box centre and not zero by accident. ---
        {
            const Real hx = 3.0, hy = 2.0, hz = 1.0, density = 1.0;
            const MassProperties3D mp = mass_properties_3d(
                SolidBox{Vec3{2.0, -1.0, 0.5}, Vec3{hx, hy, hz}, Quat{}, density});
            const Real mass = density * 8.0 * hx * hy * hz; // 48
            MALLOY_CHECK_NEAR(mp.mass, mass, 1e-9);
            MALLOY_CHECK_TRUE(
                approx_equal(mp.center_of_mass, Vec3{2.0, -1.0, 0.5}, 1e-12));
            // I_i = (1/3) m (h_j^2 + h_k^2): the two OTHER half-extents, so the
            // axis with the smallest extent carries the LARGEST moment.
            const Real third = 1.0 / 3.0;
            const Real ix = third * mass * (hy * hy + hz * hz); // 80
            const Real iy = third * mass * (hx * hx + hz * hz); // 160
            const Real iz = third * mass * (hx * hx + hy * hy); // 208
            MALLOY_CHECK_NEAR(mp.inertia.x, ix, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.y, iy, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.z, iz, 1e-9);
            MALLOY_CHECK_TRUE(approx_equal(mp.orientation, Quat{}, 1e-12));
        }

        // --- The same box, turned by a non-trivial rotation. A rotation cannot
        //     change the principal MOMENTS, but the inertia tensor in the lab
        //     frame is now genuinely off-diagonal and the returned frame has to
        //     carry it. The invariant: the reconstructed tensor equals
        //     R_tilt diag(own) R_tilt^T built independently from the closed-form
        //     moments. A builder that ignored the orientation would return a
        //     diagonal tensor and fail this. ---
        {
            const Real hx = 3.0, hy = 2.0, hz = 1.0, density = 1.0;
            const Real mass = density * 8.0 * hx * hy * hz;
            const Real third = 1.0 / 3.0;
            // Already ascending (80 < 160 < 208), so this IS the sorted result.
            const Vec3 own{third * mass * (hy * hy + hz * hz),
                           third * mass * (hx * hx + hz * hz),
                           third * mass * (hx * hx + hy * hy)};
            const Quat tilt = from_axis_angle(Vec3{1.0, 2.0, -0.5}, 0.7);
            const MassProperties3D mp = mass_properties_3d(
                SolidBox{Vec3{}, Vec3{hx, hy, hz}, tilt, density});
            MALLOY_CHECK_NEAR(mp.mass, mass, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.x, own.x, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.y, own.y, 1e-9);
            MALLOY_CHECK_NEAR(mp.inertia.z, own.z, 1e-9);
            const Mat3 rotation = to_mat3(tilt);
            const Mat3 expected =
                rotation * diagonal3(own) * transpose(rotation);
            MALLOY_CHECK_TRUE(approx_equal(tensor_of(mp), expected, 1e-9));
            // The configuration is non-trivial: the lab tensor really is not the
            // box-frame diagonal, so the check above exercises the rotation.
            MALLOY_CHECK_TRUE(!approx_equal(expected, diagonal3(own), 1e-6));
        }

        // --- combine against the M25 list path: two spheres combined one at a
        //     time must equal the same two built together in one list. The list
        //     path never calls combine, so this pins combine's tensor
        //     reconstruction and its parallel-axis shift against an independent
        //     route. The spheres are apart and unequal, so a dropped shift or a
        //     mis-weighted centre shows. ---
        {
            const SolidSphere a{Vec3{0.0, 0.0, 0.0}, 1.0, 1.0};
            const SolidSphere b{Vec3{3.0, 0.0, 0.0}, 0.7, 4.0};
            const MassProperties3D together =
                mass_properties_3d(std::vector<SolidSphere>{a, b});
            const MassProperties3D piecewise =
                combine(mass_properties_3d(std::vector<SolidSphere>{a}),
                        mass_properties_3d(std::vector<SolidSphere>{b}));
            MALLOY_CHECK_NEAR(piecewise.mass, together.mass, 1e-9);
            MALLOY_CHECK_TRUE(approx_equal(piecewise.center_of_mass,
                                           together.center_of_mass, 1e-12));
            MALLOY_CHECK_TRUE(
                approx_equal(tensor_of(piecewise), tensor_of(together), 1e-9));
        }

        // --- combine is commutative up to eigenvalue ordering, and a zero-mass
        //     body is its identity (so a fold over parts can start from
        //     nothing). A box and a sphere, unequal masses and centres apart, so
        //     the mass sum and the mass-weighted centre are both pinned. ---
        {
            const MassProperties3D a = mass_properties_3d(
                SolidBox{Vec3{-1.0, 0.0, 0.0}, Vec3{1.0, 2.0, 0.5}, Quat{}, 3.0});
            const MassProperties3D b =
                mass_properties_3d(std::vector<SolidSphere>{
                    SolidSphere{Vec3{2.0, 1.0, 0.0}, 0.8, 2.0}});
            const MassProperties3D ab = combine(a, b);
            const MassProperties3D ba = combine(b, a);
            // Closed forms for the compound's mass and centre of mass.
            const Real ma = 3.0 * 8.0 * 1.0 * 2.0 * 0.5;               // 24
            const Real mb = 2.0 * (4.0 / 3.0) * pi * 0.8 * 0.8 * 0.8;  // sphere
            MALLOY_CHECK_NEAR(ab.mass, ma + mb, 1e-9);
            MALLOY_CHECK_NEAR(ab.center_of_mass.x,
                              (ma * (-1.0) + mb * 2.0) / (ma + mb), 1e-9);
            MALLOY_CHECK_NEAR(ab.center_of_mass.y, (mb * 1.0) / (ma + mb), 1e-9);
            // Commutative: the two orders describe the same body.
            MALLOY_CHECK_NEAR(ab.mass, ba.mass, 1e-12);
            MALLOY_CHECK_TRUE(
                approx_equal(ab.center_of_mass, ba.center_of_mass, 1e-12));
            MALLOY_CHECK_TRUE(approx_equal(tensor_of(ab), tensor_of(ba), 1e-9));
            MALLOY_CHECK_TRUE(rigid_body_from(ab).is_valid());
            // Identity: a zero-mass operand on either side returns the other
            // unchanged, exactly (not diagonalized again).
            const MassProperties3D zero{};
            const MassProperties3D left = combine(zero, a);
            const MassProperties3D right = combine(a, zero);
            MALLOY_CHECK_NEAR(left.mass, a.mass, 0.0);
            MALLOY_CHECK_TRUE(approx_equal(left.center_of_mass, a.center_of_mass, 0.0));
            MALLOY_CHECK_TRUE(approx_equal(left.inertia, a.inertia, 0.0));
            MALLOY_CHECK_TRUE(approx_equal(left.orientation, a.orientation, 0.0));
            MALLOY_CHECK_NEAR(right.mass, a.mass, 0.0);
            MALLOY_CHECK_TRUE(approx_equal(right.center_of_mass, a.center_of_mass, 0.0));
            MALLOY_CHECK_TRUE(approx_equal(right.inertia, a.inertia, 0.0));
        }

        // --- combine is associative up to eigenvalue ordering: three parts of
        //     different kinds fold to the same body either way. A tilted box in
        //     the middle, so the reconstruction is not a diagonal special
        //     case. ---
        {
            const MassProperties3D a =
                mass_properties_3d(std::vector<SolidSphere>{
                    SolidSphere{Vec3{0.0, 0.0, 0.0}, 1.0, 1.0}});
            const MassProperties3D b = mass_properties_3d(
                SolidBox{Vec3{2.0, 0.0, 0.0}, Vec3{0.5, 1.0, 1.5},
                         from_axis_angle(Vec3{0.0, 0.0, 1.0}, 0.5), 2.0});
            const MassProperties3D c =
                mass_properties_3d(std::vector<SolidSphere>{
                    SolidSphere{Vec3{0.0, 3.0, 1.0}, 0.6, 5.0}});
            const MassProperties3D left = combine(combine(a, b), c);
            const MassProperties3D right = combine(a, combine(b, c));
            MALLOY_CHECK_NEAR(left.mass, right.mass, 1e-9);
            MALLOY_CHECK_TRUE(
                approx_equal(left.center_of_mass, right.center_of_mass, 1e-12));
            MALLOY_CHECK_TRUE(approx_equal(tensor_of(left), tensor_of(right), 1e-9));
        }

        // --- Validation: a box with a negative or non-finite half-extent or
        //     density, a non-unit orientation, or a non-squarable centre yields
        //     nothing usable (zero mass), rather than a meaningless number. A
        //     NEGATIVE extent on each axis in turn pins all three positivity
        //     checks: a dropped one lets a negative width through and the mass
        //     goes negative rather than staying zero. ---
        {
            const Vec3 ok{1.0, 1.0, 1.0};
            MALLOY_CHECK_NEAR(mass_properties_3d(
                SolidBox{Vec3{}, Vec3{-1.0, 1.0, 1.0}, Quat{}, 1.0}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                SolidBox{Vec3{}, Vec3{1.0, -1.0, 1.0}, Quat{}, 1.0}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                SolidBox{Vec3{}, Vec3{1.0, 1.0, -1.0}, Quat{}, 1.0}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                SolidBox{Vec3{}, Vec3{1.0, inf, 1.0}, Quat{}, 1.0}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                SolidBox{Vec3{}, ok, Quat{}, -1.0}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                SolidBox{Vec3{}, ok, Quat{}, inf}).mass, 0.0, 0.0);
            MALLOY_CHECK_NEAR(mass_properties_3d(
                SolidBox{Vec3{0.0, nan, 0.0}, ok, Quat{}, 1.0}).mass, 0.0, 0.0);
            // A non-unit orientation is malformed input: every rotation formula
            // assumes a unit quaternion, so it is rejected, not rescaled.
            MALLOY_CHECK_NEAR(mass_properties_3d(
                SolidBox{Vec3{}, ok, Quat{2.0, Vec3{0.0, 0.0, 0.0}}, 1.0}).mass,
                0.0, 0.0);
        }

        // --- The whole chain through M20's dynamics. A solid box already has
        //     three distinct principal moments, so it is a natural
        //     intermediate-axis object. Built axis-aligned, its body y is the
        //     intermediate axis (moment 160, between 80 and 208): spun there it
        //     flips, spun about the largest axis it holds. The derived moments,
        //     their ordering and rigid_body_from all have to be right for the
        //     flip to land on the middle axis. ---
        {
            // Issue #21: box_inertia gives the moments in the box's OWN axis
            // order, unsorted, for a body whose half_extents keep that order.
            // Every permutation of (1, 2, 3) at mass 48 (16 times the sums of
            // the other two squares), and the sorted principal moments of the
            // same box are a permutation of them.
            using malloy::rigid::box_inertia;
            const Vec3 ascending = box_inertia(Vec3{1.0, 2.0, 3.0}, 48.0);
            MALLOY_CHECK_TRUE(approx_equal(ascending, Vec3{208.0, 160.0, 80.0}, 1e-9));
            MALLOY_CHECK_TRUE(approx_equal(box_inertia(Vec3{2.0, 3.0, 1.0}, 48.0),
                                           Vec3{160.0, 80.0, 208.0}, 1e-9));
            MALLOY_CHECK_TRUE(approx_equal(box_inertia(Vec3{3.0, 1.0, 2.0}, 48.0),
                                           Vec3{80.0, 208.0, 160.0}, 1e-9));
            const MassProperties3D sorted = mass_properties_3d(
                SolidBox{Vec3{}, Vec3{1.0, 2.0, 3.0}, Quat{}, 48.0 / 48.0});
            MALLOY_CHECK_TRUE(approx_equal(sorted.inertia, Vec3{80.0, 160.0, 208.0}, 1e-9));
            // Invalid input gives zero moments, which no body accepts.
            const Real bad_inf = std::numeric_limits<Real>::infinity();
            MALLOY_CHECK_TRUE(approx_equal(box_inertia(Vec3{1.0, 0.0, 1.0}, 1.0), Vec3{}, 0.0));
            MALLOY_CHECK_TRUE(approx_equal(box_inertia(Vec3{1.0, -1.0, 1.0}, 1.0), Vec3{}, 0.0));
            MALLOY_CHECK_TRUE(approx_equal(box_inertia(Vec3{1.0, bad_inf, 1.0}, 1.0), Vec3{}, 0.0));
            MALLOY_CHECK_TRUE(approx_equal(box_inertia(Vec3{1.0, 1.0, 1.0}, 0.0), Vec3{}, 0.0));
            MALLOY_CHECK_TRUE(approx_equal(box_inertia(Vec3{1.0, 1.0, 1.0}, -2.0), Vec3{}, 0.0));
        }
        {
            using malloy::rigid::Rigid3DWorld;
            const MassProperties3D mp = mass_properties_3d(
                SolidBox{Vec3{}, Vec3{3.0, 2.0, 1.0}, Quat{}, 1.0});
            MALLOY_CHECK_TRUE(mp.inertia.x < mp.inertia.y - 1.0);
            MALLOY_CHECK_TRUE(mp.inertia.y < mp.inertia.z - 1.0);
            MALLOY_CHECK_TRUE(approx_equal(mp.orientation, Quat{}, 1e-12));

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
            // Spun about z (largest): it does NOT flip.
            {
                RigidBody3D body = rigid_body_from(mp);
                body.angular_velocity = Vec3{nudge, 0.0, spin};
                Rigid3DWorld world{SimulationSettings{0.0005}, {body}};
                Real lowest = spin;
                for (int i = 0; i < 40000; ++i)
                {
                    MALLOY_CHECK_TRUE(world.step().ok());
                    lowest = std::fmin(lowest, world.bodies()[0].angular_velocity.z);
                }
                MALLOY_CHECK_TRUE(lowest > 0.9 * spin); // held: stable axis
            }
        }
    }

    // --- M27: mesh mass properties. The general primitive: a solid uniform
    //     body bounded by a closed triangle mesh, with mass, centre of mass and
    //     inertia computed as signed-tetrahedron volume integrals. The sphere
    //     and box were special cases; this handles any closed shape, and an
    //     axis-aligned box mesh must reproduce the M26 box exactly. Validated
    //     against closed forms, then run through M20's dynamics. ---
    {
        using malloy::math::approx_equal;
        using malloy::math::diagonal3;
        using malloy::math::from_axis_angle;
        using malloy::math::Mat3;
        using malloy::math::Quat;
        using malloy::math::to_mat3;
        using malloy::math::transpose;
        using malloy::math::Vec3;
        using malloy::rigid::combine;
        using malloy::rigid::mass_properties_3d;
        using malloy::rigid::MassProperties3D;
        using malloy::rigid::MeshTriangle;
        using malloy::rigid::rigid_body_from;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::SolidBox;
        using malloy::rigid::SolidMesh;
        using malloy::rigid::SolidSphere;

        // The full inertia tensor in the lab frame, reconstructed from a
        // result's principal moments and frame (as in M26): the ordering-proof
        // way to compare two mass-property sets.
        const auto tensor_of = [](const MassProperties3D& mp) {
            const Mat3 rotation = to_mat3(mp.orientation);
            return rotation * diagonal3(mp.inertia) * transpose(rotation);
        };

        // Orient a convex mesh's triangles outward: flip any whose normal points
        // toward the mesh centroid. Lets the builders below list faces without
        // fussing over winding, and the builder proper still rejects a mesh that
        // is genuinely inside out.
        const auto orient_outward = [](SolidMesh& mesh) {
            Vec3 center{};
            for (const Vec3& v : mesh.vertices)
            {
                center += v;
            }
            center = center / static_cast<Real>(mesh.vertices.size());
            for (MeshTriangle& t : mesh.triangles)
            {
                const Vec3& a = mesh.vertices[t.v0];
                const Vec3& b = mesh.vertices[t.v1];
                const Vec3& c = mesh.vertices[t.v2];
                const Vec3 normal = malloy::math::cross(b - a, c - a);
                const Vec3 outward = (a + b + c) * (1.0 / 3.0) - center;
                if (malloy::math::dot(normal, outward) < 0.0)
                {
                    const std::size_t tmp = t.v1;
                    t.v1 = t.v2;
                    t.v2 = tmp;
                }
            }
        };

        // A box as a mesh: eight corners, twelve triangles, wound outward.
        const auto box_mesh = [&](const Vec3& center, const Vec3& he, const Quat& q,
                                  Real density) {
            SolidMesh mesh;
            mesh.density = density;
            for (int i = 0; i < 8; ++i)
            {
                const Real sx = (i & 1) ? 1.0 : -1.0;
                const Real sy = (i & 2) ? 1.0 : -1.0;
                const Real sz = (i & 4) ? 1.0 : -1.0;
                mesh.vertices.push_back(
                    center +
                    malloy::math::rotate(q, Vec3{sx * he.x, sy * he.y, sz * he.z}));
            }
            const std::size_t face[6][4] = {{0, 2, 6, 4}, {1, 3, 7, 5},
                                            {0, 1, 5, 4}, {2, 3, 7, 6},
                                            {0, 1, 3, 2}, {4, 5, 7, 6}};
            for (const auto& f : face)
            {
                mesh.triangles.push_back(MeshTriangle{f[0], f[1], f[2]});
                mesh.triangles.push_back(MeshTriangle{f[0], f[2], f[3]});
            }
            orient_outward(mesh);
            return mesh;
        };

        const auto tet_mesh = [&](const Vec3& p0, const Vec3& p1, const Vec3& p2,
                                  const Vec3& p3, Real density) {
            SolidMesh mesh;
            mesh.density = density;
            mesh.vertices = {p0, p1, p2, p3};
            mesh.triangles = {MeshTriangle{0, 1, 2}, MeshTriangle{0, 1, 3},
                              MeshTriangle{0, 2, 3}, MeshTriangle{1, 2, 3}};
            orient_outward(mesh);
            return mesh;
        };

        // --- An axis-aligned box mesh reproduces the M26 box exactly: same
        //     mass, centre, principal moments and full tensor. This is the
        //     anchor, tying the general integral to a known closed form. ---
        {
            const Vec3 center{2.0, -1.0, 0.5};
            const Vec3 he{3.0, 2.0, 1.0};
            const MassProperties3D mesh = mass_properties_3d(box_mesh(center, he, Quat{}, 1.0));
            const MassProperties3D box =
                mass_properties_3d(SolidBox{center, he, Quat{}, 1.0});
            MALLOY_CHECK_NEAR(mesh.mass, box.mass, 1e-9);
            MALLOY_CHECK_TRUE(approx_equal(mesh.center_of_mass, box.center_of_mass, 1e-9));
            MALLOY_CHECK_TRUE(approx_equal(mesh.inertia, box.inertia, 1e-7));
            MALLOY_CHECK_TRUE(approx_equal(tensor_of(mesh), tensor_of(box), 1e-7));
        }

        // --- A cube away from the origin: the centre of mass comes out right
        //     and the parallel-axis shift leaves the correct central, isotropic
        //     moment. (13 sizes out, which was too close to show issue #45; the
        //     block below goes much further.) ---
        {
            const MassProperties3D mp =
                mass_properties_3d(box_mesh(Vec3{10.0, -5.0, 7.0}, Vec3{1.0, 1.0, 1.0},
                                            Quat{}, 2.0));
            MALLOY_CHECK_NEAR(mp.mass, 16.0, 1e-9); // 2 * 8
            MALLOY_CHECK_TRUE(approx_equal(mp.center_of_mass, Vec3{10.0, -5.0, 7.0}, 1e-8));
            const Real iso = (2.0 / 3.0) * 16.0; // (1/3) m (h^2 + h^2), h = 1
            MALLOY_CHECK_NEAR(mp.inertia.x, iso, 1e-6);
            MALLOY_CHECK_NEAR(mp.inertia.y, iso, 1e-6);
            MALLOY_CHECK_NEAR(mp.inertia.z, iso, 1e-6);
        }

        // --- Issue #45: a small mesh FAR from the origin. Spanned from the
        //     world origin, the signed tetrahedra cancel terms of size d^3 and
        //     d^5 down to an answer of size s^3 and s^5, so the relative error
        //     grew like eps (d / s)^3. At the two positions below, the old sums
        //     had the moments off by 4.7e-3 at 10 m, and off by -7.6e7 at 1 km,
        //     negative, so the body was refused (in the issue's own
        //     reproduction, another 1 km position was ACCEPTED with moments off
        //     by 3.7e6 and the centre of mass 0.15 m away). Spanned from the
        //     vertex mean it must match the same cube at the origin.
        //
        //     The tolerances are what the vertices themselves allow. A vertex 10
        //     m out is only known to about eps * 10 m = 2e-15 m, a relative 2e-13
        //     of the 1 cm edge; the volume is third order in the vertices and the
        //     moments fifth, so they can be off by about 1e-12, and 1e-11 leaves
        //     ten times that (measured: 1e-13). 1 km out it scales by 100: 1e-9
        //     (measured: 6e-12). The centre of mass is pinned to the rounding of
        //     its own coordinates: 1e-13 m at 10 m, 1e-11 m at 1 km (it was off
        //     by 7e-9 m and 0.15 m). A 1 cm cube of density 1000, turned about a
        //     generic axis, in a generic direction. ---
        {
            const Quat turn = from_axis_angle(Vec3{0.3, -0.7, 0.5}, 0.9);
            const Vec3 half{0.005, 0.005, 0.005};
            const MassProperties3D home =
                mass_properties_3d(box_mesh(Vec3{}, half, turn, 1000.0));
            MALLOY_CHECK_NEAR(home.mass, 1e-3, 1e-15); // 1000 * (0.01)^3
            const Vec3 direction = malloy::math::normalize(Vec3{0.37, -0.81, 0.45});
            struct Far
            {
                Real distance;
                Real tolerance;          // relative, on mass and moments
                Real position_tolerance; // absolute, on the centre of mass
            };
            for (const Far& f : {Far{10.0, 1e-11, 1e-13}, Far{1000.0, 1e-9, 1e-11}})
            {
                const Vec3 center = direction * f.distance;
                const MassProperties3D far =
                    mass_properties_3d(box_mesh(center, half, turn, 1000.0));
                MALLOY_CHECK_NEAR(far.mass / home.mass, 1.0, f.tolerance);
                MALLOY_CHECK_TRUE(
                    approx_equal(far.center_of_mass, center, f.position_tolerance));
                MALLOY_CHECK_NEAR(far.inertia.x / home.inertia.x, 1.0, f.tolerance);
                MALLOY_CHECK_NEAR(far.inertia.y / home.inertia.y, 1.0, f.tolerance);
                MALLOY_CHECK_NEAR(far.inertia.z / home.inertia.z, 1.0, f.tolerance);
                MALLOY_CHECK_TRUE(approx_equal(tensor_of(far), tensor_of(home),
                                               f.tolerance * home.inertia.x));
                MALLOY_CHECK_TRUE(rigid_body_from(far).is_valid());
            }
        }

        // --- A tetrahedron: the volume is det/6 and the centre of mass is the
        //     mean of the four vertices. One vertex sits at the origin; the
        //     tetrahedra are spanned from the vertex mean, an interior point, so
        //     all four faces contribute and the signed sum still adds up to the
        //     one solid. ---
        {
            const MassProperties3D mp = mass_properties_3d(tet_mesh(
                Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0}, Vec3{0.0, 3.0, 0.0},
                Vec3{0.0, 0.0, 4.0}, 3.0));
            MALLOY_CHECK_NEAR(mp.mass, 12.0, 1e-9); // 3 * (2*3*4/6) = 3 * 4
            MALLOY_CHECK_TRUE(
                approx_equal(mp.center_of_mass, Vec3{0.5, 0.75, 1.0}, 1e-12));
            MALLOY_CHECK_TRUE(rigid_body_from(mp).is_valid());
        }

        // --- A tilted box mesh: the principal moments are unchanged (a rotation
        //     cannot change them) and the full tensor matches the M26 box turned
        //     the same way. Ties the mesh path to M26's tilted case. ---
        {
            const Vec3 he{3.0, 2.0, 1.0};
            const Quat tilt = from_axis_angle(Vec3{1.0, 2.0, -0.5}, 0.7);
            const MassProperties3D mesh = mass_properties_3d(box_mesh(Vec3{}, he, tilt, 1.0));
            const MassProperties3D axis =
                mass_properties_3d(SolidBox{Vec3{}, he, Quat{}, 1.0});
            MALLOY_CHECK_TRUE(approx_equal(mesh.inertia, axis.inertia, 1e-7));
            const MassProperties3D box =
                mass_properties_3d(SolidBox{Vec3{}, he, tilt, 1.0});
            MALLOY_CHECK_TRUE(approx_equal(tensor_of(mesh), tensor_of(box), 1e-7));
        }

        // --- Winding carries the sign: a mesh wound inward describes negative
        //     volume and is rejected, rather than yielding a body inside out. ---
        {
            SolidMesh reversed = box_mesh(Vec3{}, Vec3{1.0, 1.0, 1.0}, Quat{}, 1.0);
            for (MeshTriangle& t : reversed.triangles)
            {
                const std::size_t tmp = t.v1;
                t.v1 = t.v2;
                t.v2 = tmp;
            }
            MALLOY_CHECK_NEAR(mass_properties_3d(reversed).mass, 0.0, 0.0);
        }

        // --- A mesh result composes with combine like any other primitive: a
        //     box mesh and a sphere, unequal masses, centres apart. The box mesh
        //     mass also cross-checks the box closed form (3 * 8 * 1 * 2 * 0.5). ---
        {
            const MassProperties3D box =
                mass_properties_3d(box_mesh(Vec3{-1.0, 0.0, 0.0}, Vec3{1.0, 2.0, 0.5},
                                            Quat{}, 3.0));
            const MassProperties3D sphere = mass_properties_3d(
                std::vector<SolidSphere>{SolidSphere{Vec3{2.0, 1.0, 0.0}, 0.8, 2.0}});
            const Real mbox = 3.0 * 8.0 * 1.0 * 2.0 * 0.5;             // 24
            const Real msph = 2.0 * (4.0 / 3.0) * pi * 0.8 * 0.8 * 0.8;
            MALLOY_CHECK_NEAR(box.mass, mbox, 1e-9);
            const MassProperties3D both = combine(box, sphere);
            MALLOY_CHECK_NEAR(both.mass, mbox + msph, 1e-9);
            MALLOY_CHECK_NEAR(both.center_of_mass.x,
                              (mbox * (-1.0) + msph * 2.0) / (mbox + msph), 1e-9);
            MALLOY_CHECK_NEAR(both.center_of_mass.y, (msph * 1.0) / (mbox + msph), 1e-9);
        }

        // --- Validation: a malformed or degenerate mesh yields nothing usable. ---
        {
            const SolidMesh good = box_mesh(Vec3{}, Vec3{1.0, 1.0, 1.0}, Quat{}, 1.0);
            MALLOY_CHECK_TRUE(mass_properties_3d(good).mass > 0.0);
            // Fewer than four vertices, or four faces: not a closed solid.
            SolidMesh few_vertices = good;
            few_vertices.vertices.resize(3);
            MALLOY_CHECK_NEAR(mass_properties_3d(few_vertices).mass, 0.0, 0.0);
            SolidMesh few_faces = good;
            few_faces.triangles.resize(3);
            MALLOY_CHECK_NEAR(mass_properties_3d(few_faces).mass, 0.0, 0.0);
            // A triangle index past the vertex list.
            SolidMesh bad_index = good;
            bad_index.triangles[0].v0 = good.vertices.size();
            MALLOY_CHECK_NEAR(mass_properties_3d(bad_index).mass, 0.0, 0.0);
            // A non-finite vertex, a non-positive and a non-finite density.
            SolidMesh nan_vertex = good;
            nan_vertex.vertices[0] = Vec3{0.0, nan, 0.0};
            MALLOY_CHECK_NEAR(mass_properties_3d(nan_vertex).mass, 0.0, 0.0);
            SolidMesh negative_density = good;
            negative_density.density = -1.0;
            MALLOY_CHECK_NEAR(mass_properties_3d(negative_density).mass, 0.0, 0.0);
            SolidMesh infinite_density = good;
            infinite_density.density = inf;
            MALLOY_CHECK_NEAR(mass_properties_3d(infinite_density).mass, 0.0, 0.0);
            // A flat (coplanar) mesh encloses no volume.
            SolidMesh flat = good;
            for (Vec3& v : flat.vertices)
            {
                v.z = 0.0;
            }
            MALLOY_CHECK_NEAR(mass_properties_3d(flat).mass, 0.0, 0.0);

            // Issue #34: a mesh must be CLOSED. Without its top face (the last
            // two triangles) a unit cube used to weigh 0.67 at the origin, 0.83
            // half a unit down and 1.67 three units down: the signed volume
            // stayed positive, so nothing refused it. Now it is refused
            // wherever it sits, while the closed cube weighs 8 everywhere.
            for (const Real z : {0.0, -0.5, -3.0})
            {
                const SolidMesh closed =
                    box_mesh(Vec3{0.0, 0.0, z}, Vec3{1.0, 1.0, 1.0}, Quat{}, 1.0);
                MALLOY_CHECK_NEAR(mass_properties_3d(closed).mass, 8.0, 1e-12);
                SolidMesh open = closed;
                open.triangles.resize(open.triangles.size() - 2);
                MALLOY_CHECK_NEAR(mass_properties_3d(open).mass, 0.0, 0.0);
            }
            // Wound inconsistently: one face flipped traverses an edge the same
            // way as its neighbour does.
            SolidMesh flipped = good;
            std::swap(flipped.triangles[0].v1, flipped.triangles[0].v2);
            MALLOY_CHECK_NEAR(mass_properties_3d(flipped).mass, 0.0, 0.0);
            // A face listed twice, and a triangle that repeats a vertex.
            SolidMesh doubled = good;
            doubled.triangles.push_back(doubled.triangles[0]);
            MALLOY_CHECK_NEAR(mass_properties_3d(doubled).mass, 0.0, 0.0);
            SolidMesh repeated = good;
            repeated.triangles[0].v1 = repeated.triangles[0].v0;
            MALLOY_CHECK_NEAR(mass_properties_3d(repeated).mass, 0.0, 0.0);
            // A sliver between two corners no edge joins (0 and 7, opposite
            // corners of the cube) pairs its own edges, 0->7 with 7->0, so only
            // the repeated-vertex rule refuses it.
            SolidMesh sliver = good;
            sliver.triangles.push_back(MeshTriangle{0, 7, 0});
            MALLOY_CHECK_NEAR(mass_properties_3d(sliver).mass, 0.0, 0.0);
        }

        // --- The whole chain through M20's dynamics. A box mesh has three
        //     distinct principal moments, so like the M26 box it is an
        //     intermediate-axis object: spun about its middle (body y) axis it
        //     flips. The derived moments and their ordering all have to be right
        //     for the flip to land there. ---
        {
            using malloy::rigid::Rigid3DWorld;
            const MassProperties3D mp =
                mass_properties_3d(box_mesh(Vec3{}, Vec3{3.0, 2.0, 1.0}, Quat{}, 1.0));
            MALLOY_CHECK_TRUE(mp.inertia.x < mp.inertia.y - 1.0);
            MALLOY_CHECK_TRUE(mp.inertia.y < mp.inertia.z - 1.0);
            const Real spin = 5.0, nudge = 1e-3;
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
    }

    // --- M28: a constant applied FORCE, the translational half of a wrench.
    //     M21 put a torque on a body; M28 puts a force. It enters as the
    //     acceleration F/m before the position update, so it scales with mass
    //     (a force, not gravity's acceleration) and acts through the centre of
    //     mass, making no torque. Together with the torque setting it expresses
    //     any constant wrench, and the two halves do not interfere. ---
    {
        using malloy::math::approx_equal;
        using malloy::math::Vec3;
        using malloy::rigid::mass_properties_3d;
        using malloy::rigid::rigid_body_from;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::SolidSphere;
        using malloy::rigid::total_force_potential3d;
        using malloy::rigid::total_kinetic_energy3d;
        using malloy::rigid::total_linear_momentum3d;

        // --- Validation: the force must be SQUARABLE, like the torque. ---
        {
            Rigid3DSettings good;
            good.force = Vec3{1.0, -2.0, 3.0};
            MALLOY_CHECK_TRUE(good.is_valid());
            Rigid3DSettings infinite;
            infinite.force = Vec3{inf, 0.0, 0.0};
            MALLOY_CHECK_FALSE(infinite.is_valid());
            Rigid3DSettings not_a_number;
            not_a_number.force = Vec3{0.0, nan, 0.0};
            MALLOY_CHECK_FALSE(not_a_number.is_valid());
            Rigid3DSettings unsquarable; // finite, but its square overflows
            unsquarable.force = Vec3{0.0, 0.0, 1.4e154};
            MALLOY_CHECK_FALSE(unsquarable.is_valid());
            // A world with a bad force refuses to step and leaves state alone.
            RigidBody3D body;
            Rigid3DWorld world{SimulationSettings{0.01}, {body}, infinite};
            MALLOY_CHECK_TRUE(world.step().status == StepStatus::InvalidSettings);
        }

        // --- Newton's second law, discretely exact in momentum. A body from
        //     rest under a constant force: the velocity reaches (F/m) t and the
        //     momentum reaches F t, because each step adds the impulse F dt
        //     exactly (m times the acceleration F/m). ---
        {
            const Real dt = 0.001;
            const Vec3 force{2.0, -3.0, 1.5};
            RigidBody3D body = rigid_body_from(mass_properties_3d(
                std::vector<SolidSphere>{SolidSphere{Vec3{}, 0.5, 1000.0}}));
            const Real m = body.mass;
            Rigid3DSettings settings;
            settings.force = force;
            Rigid3DWorld world{SimulationSettings{dt}, {body}, settings};
            const int steps = 5000;
            for (int i = 0; i < steps; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            const Real t = static_cast<Real>(steps) * dt;
            MALLOY_CHECK_TRUE(
                approx_equal(world.bodies()[0].velocity, force * (t / m), 1e-9));
            MALLOY_CHECK_TRUE(
                approx_equal(total_linear_momentum3d(world.bodies()), force * t, 1e-6));
        }

        // --- The force scales with mass, which is what makes it a force and not
        //     gravity. Two spheres of the same size but different density feel
        //     the same force; the lighter accelerates more, in inverse
        //     proportion to mass, and both gain the same momentum (equal
        //     impulse). This reads the mass computed in M25 to M27. ---
        {
            const Real dt = 0.001;
            const Vec3 force{5.0, 0.0, 0.0};
            RigidBody3D light = rigid_body_from(mass_properties_3d(
                std::vector<SolidSphere>{SolidSphere{Vec3{}, 0.5, 1000.0}}));
            RigidBody3D heavy = rigid_body_from(mass_properties_3d(
                std::vector<SolidSphere>{SolidSphere{Vec3{}, 0.5, 4000.0}}));
            Rigid3DSettings settings;
            settings.force = force;
            Rigid3DWorld world{SimulationSettings{dt}, {light, heavy}, settings};
            for (int i = 0; i < 1000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            const RigidBody3D& l = world.bodies()[0];
            const RigidBody3D& h = world.bodies()[1];
            MALLOY_CHECK_TRUE(l.velocity.x > h.velocity.x); // lighter is faster
            // Equal impulse: m v is the same for both (each is F.x t).
            MALLOY_CHECK_NEAR(l.mass * l.velocity.x, h.mass * h.velocity.x, 1e-9);
            // The speed ratio is the inverse mass ratio (4x denser is 4x slower).
            MALLOY_CHECK_NEAR(l.velocity.x / h.velocity.x, h.mass / l.mass, 1e-9);
        }

        // --- A centre-of-mass force makes NO torque: its lever arm is zero. A
        //     tumbling body run with a force has bit-for-bit the same angular
        //     velocity and orientation as the same body run free, while its
        //     linear velocity differs. The linear and rotational halves are
        //     orthogonal. ---
        {
            const Real dt = 0.005;
            RigidBody3D body;
            body.inertia = Vec3{2.0, 3.0, 5.0}; // asymmetric, so it tumbles
            body.angular_velocity = Vec3{0.5, 4.0, 0.2};
            body.orientation = malloy::math::from_axis_angle(Vec3{1.0, 1.0, 1.0}, 0.5);
            body.mass = 3.0;
            Rigid3DSettings pushed;
            pushed.force = Vec3{3.0, -1.0, 2.0};
            Rigid3DWorld forced{SimulationSettings{dt}, {body}, pushed};
            Rigid3DWorld free{SimulationSettings{dt}, {body}, Rigid3DSettings{}};
            for (int i = 0; i < 500; ++i)
            {
                MALLOY_CHECK_TRUE(forced.step().ok());
                MALLOY_CHECK_TRUE(free.step().ok());
                MALLOY_CHECK_TRUE(approx_equal(forced.bodies()[0].angular_velocity,
                                               free.bodies()[0].angular_velocity, 0.0));
                MALLOY_CHECK_TRUE(approx_equal(forced.bodies()[0].orientation,
                                               free.bodies()[0].orientation, 0.0));
            }
            // The force did act, so the linear velocities are not equal.
            MALLOY_CHECK_FALSE(approx_equal(forced.bodies()[0].velocity,
                                            free.bodies()[0].velocity, 1e-9));
        }

        // --- The energy law. Kinetic plus the force potential decreases by
        //     exactly (1/2) |F|^2 / m dt^2 per step, semi-implicit Euler's
        //     second-order shed for a constant acceleration a = F/m (the same
        //     shape as the gravity drift, carrying the 1/m the acceleration
        //     does not). A non-spinning sphere, so all the energy is
        //     translational. ---
        {
            const Real dt = 0.002;
            const Vec3 force{4.0, -2.0, 1.0};
            RigidBody3D body = rigid_body_from(mass_properties_3d(
                std::vector<SolidSphere>{SolidSphere{Vec3{}, 0.4, 2000.0}}));
            const Real m = body.mass;
            Rigid3DSettings settings;
            settings.force = force;
            Rigid3DWorld world{SimulationSettings{dt}, {body}, settings};
            const Real per_step =
                0.5 * malloy::math::length_squared(force) / m * dt * dt;
            Real worst = 0.0;
            for (int i = 0; i < 300; ++i)
            {
                const Real before = total_kinetic_energy3d(world.bodies()) +
                                    total_force_potential3d(world.bodies(), force);
                MALLOY_CHECK_TRUE(world.step().ok());
                const Real after = total_kinetic_energy3d(world.bodies()) +
                                   total_force_potential3d(world.bodies(), force);
                worst = std::fmax(worst, std::abs((after - before) + per_step));
            }
            MALLOY_CHECK_TRUE(worst < 1e-12);
        }

        // --- The complete wrench: force AND torque, at once. The body's angular
        //     state matches a torque-only run (the force adds no torque) and its
        //     linear velocity matches a force-only run (the torque adds no
        //     force), both bit for bit. The halves do not interfere, which is
        //     what lets an off-centre push be written as force plus couple. ---
        {
            const Real dt = 0.002;
            RigidBody3D body;
            body.inertia = Vec3{2.0, 3.0, 5.0};
            body.angular_velocity = Vec3{0.4, 0.1, -0.3};
            body.orientation = malloy::math::from_axis_angle(Vec3{0.0, 1.0, 0.0}, 0.3);
            body.mass = 5.0;
            const Vec3 force{2.0, 1.0, -1.0};
            const Vec3 torque{0.3, -0.5, 0.2};
            Rigid3DSettings wrench_settings;
            wrench_settings.force = force;
            wrench_settings.torque = torque;
            Rigid3DSettings force_only;
            force_only.force = force;
            Rigid3DSettings torque_only;
            torque_only.torque = torque;
            Rigid3DWorld wrench{SimulationSettings{dt}, {body}, wrench_settings};
            Rigid3DWorld linear{SimulationSettings{dt}, {body}, force_only};
            Rigid3DWorld angular{SimulationSettings{dt}, {body}, torque_only};
            for (int i = 0; i < 300; ++i)
            {
                MALLOY_CHECK_TRUE(wrench.step().ok());
                MALLOY_CHECK_TRUE(linear.step().ok());
                MALLOY_CHECK_TRUE(angular.step().ok());
                MALLOY_CHECK_TRUE(approx_equal(wrench.bodies()[0].angular_velocity,
                                               angular.bodies()[0].angular_velocity, 0.0));
                MALLOY_CHECK_TRUE(approx_equal(wrench.bodies()[0].velocity,
                                               linear.bodies()[0].velocity, 0.0));
            }
        }

        // --- The force potential: the sum of -(F . r), no mass factor (the
        //     force is uniform, not mass-scaled like gravity's m g), and zero
        //     for a zero force. ---
        {
            RigidBody3D a;
            a.position = Vec3{1.0, 2.0, 3.0};
            a.mass = 2.0;
            RigidBody3D b;
            b.position = Vec3{-1.0, 0.0, 4.0};
            b.mass = 5.0;
            const Vec3 force{10.0, 0.0, -2.0};
            const Real expected =
                -(malloy::math::dot(force, a.position) + malloy::math::dot(force, b.position));
            MALLOY_CHECK_NEAR(total_force_potential3d({a, b}, force), expected, 1e-12);
            MALLOY_CHECK_NEAR(total_force_potential3d({a, b}, Vec3{}), 0.0, 0.0);
        }
    }

    // --- M30: an oriented box against a ground plane. The first non-sphere 3D
    //     collision. A box gets a box collider (half_extents), and its contact
    //     with a plane is the manifold of its penetrating corners, reduced to a
    //     single centroid contact. A flat drop must stay flat (no spurious
    //     torque); a corner landing must tip (a real lever arm). ---
    {
        using malloy::math::from_axis_angle;
        using malloy::math::Quat;
        using malloy::math::Vec3;
        using malloy::rigid::mass_properties_3d;
        using malloy::rigid::rigid_body_from;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::SolidBox;

        // A box body whose inertia is computed from its geometry (M26), as the
        // scenario loader does, then given the box collider and a pose.
        const auto box_body = [&](const Vec3& half, Real mass, const Vec3& position,
                                  const Quat& orientation) {
            const Real density = mass / (8.0 * half.x * half.y * half.z);
            RigidBody3D body =
                rigid_body_from(mass_properties_3d(SolidBox{Vec3{}, half, Quat{}, density}));
            body.half_extents = half;
            body.position = position;
            body.orientation = orientation;
            return body;
        };
        const auto floor_settings = [](Real restitution, Real friction) {
            Rigid3DSettings s;
            s.gravity = Vec3{0.0, 0.0, -4.0};
            s.restitution = restitution;
            s.friction = friction;
            s.ground.push_back(malloy::collide::Plane3{Vec3{0.0, 0.0, 1.0}, 0.0});
            return s;
        };
        const Vec3 half{0.5, 0.5, 0.5};

        // --- A flat, symmetric drop induces NO spin and NO horizontal drift:
        //     the centroid contact lies below the centre of mass, so its normal
        //     impulse makes no torque. The cube settles resting with its centre
        //     at half a side height. This is the headline invariant. ---
        {
            RigidBody3D body = box_body(half, 1.0, Vec3{0.0, 0.0, 2.0}, Quat{});
            Rigid3DWorld world{SimulationSettings{0.002}, {body}, floor_settings(0.5, 0.0)};
            Real worst_spin = 0.0;
            Real worst_horizontal = 0.0;
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                const RigidBody3D& b = world.bodies()[0];
                worst_spin = std::fmax(worst_spin, malloy::math::length(b.angular_velocity));
                worst_horizontal = std::fmax(
                    worst_horizontal, std::abs(b.velocity.x) + std::abs(b.velocity.y));
            }
            MALLOY_CHECK_TRUE(worst_spin < 1e-12);       // symmetric: no torque
            MALLOY_CHECK_TRUE(worst_horizontal < 1e-12); // symmetric: no drift
            const RigidBody3D& b = world.bodies()[0];
            MALLOY_CHECK_NEAR(b.position.z, 0.5, 1e-2);      // resting at half height
            MALLOY_CHECK_TRUE(malloy::math::length(b.velocity) < 0.1);
        }

        // --- The box collider is what stops it: a box body (radius 0) rests on
        //     the plane, while a body with neither a box nor a sphere collider
        //     falls straight through. ---
        {
            RigidBody3D solid = box_body(half, 1.0, Vec3{0.0, 0.0, 1.5}, Quat{});
            RigidBody3D ghost = box_body(half, 1.0, Vec3{0.0, 0.0, 1.5}, Quat{});
            ghost.half_extents = Vec3{}; // no box, and radius is 0: does not collide
            Rigid3DWorld world{SimulationSettings{0.002},
                               {solid, ghost}, floor_settings(0.0, 0.0)};
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(world.bodies()[0].position.z > 0.4);   // caught on the floor
            MALLOY_CHECK_TRUE(world.bodies()[1].position.z < -3.0);  // fell through
        }

        // --- A corner-first landing DOES spin the box up: the contact is now
        //     off the centre of mass, so it has a lever arm. It tips and never
        //     sinks through the floor. ---
        {
            RigidBody3D body =
                box_body(half, 1.0, Vec3{0.0, 0.0, 2.0}, from_axis_angle(Vec3{1.0, 1.0, 0.0}, 0.6));
            Rigid3DWorld world{SimulationSettings{0.002}, {body}, floor_settings(0.2, 0.5)};
            Real peak_spin = 0.0;
            Real lowest_center = 10.0;
            bool landed = false;
            for (int i = 0; i < 3000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                const RigidBody3D& b = world.bodies()[0];
                if (i > 400) // after it has reached the floor
                {
                    peak_spin = std::fmax(peak_spin, malloy::math::length(b.angular_velocity));
                    lowest_center = std::fmin(lowest_center, b.position.z);
                }
            }
            landed = world.bodies()[0].position.z > 0.3;
            MALLOY_CHECK_TRUE(peak_spin > 0.1);      // the corner torqued it
            MALLOY_CHECK_TRUE(lowest_center > 0.0);  // never sank through
            MALLOY_CHECK_TRUE(landed);               // came to rest above the floor
        }

        // --- A box placed exactly at rest stays put: no drift, no spin, no
        //     sinking, over a long run. ---
        {
            RigidBody3D body = box_body(half, 1.0, Vec3{0.0, 0.0, 0.5}, Quat{});
            Rigid3DWorld world{SimulationSettings{0.002}, {body}, floor_settings(0.5, 0.0)};
            for (int i = 0; i < 1500; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            const RigidBody3D& b = world.bodies()[0];
            MALLOY_CHECK_NEAR(b.position.z, 0.5, 1e-2);
            MALLOY_CHECK_TRUE(malloy::math::length(b.velocity) < 0.05);
            MALLOY_CHECK_TRUE(malloy::math::length(b.angular_velocity) < 1e-12);
        }

        // --- Validation: a non-positive or non-finite half-extent is an invalid
        //     body, so the world refuses to step. ---
        {
            RigidBody3D bad = box_body(half, 1.0, Vec3{}, Quat{});
            bad.half_extents = Vec3{-1.0, 0.5, 0.5};
            Rigid3DWorld world{SimulationSettings{0.01}, {bad}, Rigid3DSettings{}};
            MALLOY_CHECK_TRUE(world.step().status == StepStatus::InvalidState);
        }
    }

    // --- M33: box against box, the first 3D body-vs-body collision that is not
    //     two spheres. It runs through the SAME impulse core (ADR 0008), with
    //     the separating-axis contact standing in for the sphere one. The
    //     clipped manifold is reduced to one contact at its centroid, so this is
    //     a bounce, not a stack (the iterative solver a stack needs is deferred).
    //     What stays EXACT is total linear momentum (the impulse is equal and
    //     opposite) and, at frictionless restitution 1, total kinetic energy; an
    //     aligned head-on hit makes no spin (the review-pass block pins that). ---
    {
        using malloy::math::from_axis_angle;
        using malloy::math::normalize;
        using malloy::math::Quat;
        using malloy::math::Vec3;
        using malloy::rigid::mass_properties_3d;
        using malloy::rigid::rigid_body_from;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::SolidBox;
        using malloy::rigid::total_kinetic_energy3d;
        using malloy::rigid::total_linear_momentum3d;

        const Vec3 half{0.5, 0.5, 0.5};
        // A unit cube whose inertia comes from geometry (M26), as the loader
        // builds one, then given the box collider, a pose and a velocity.
        const auto cube = [&](Real mass, const Vec3& position, const Vec3& velocity,
                              const Quat& orientation = Quat{}) {
            const Real density = mass / (8.0 * half.x * half.y * half.z);
            RigidBody3D body =
                rigid_body_from(mass_properties_3d(SolidBox{Vec3{}, half, Quat{}, density}));
            body.half_extents = half;
            body.position = position;
            body.velocity = velocity;
            body.orientation = orientation;
            return body;
        };

        // Head-on along x, overlapping 0.1, unequal masses so the conserved
        // momentum is a nonzero value rather than a trivial zero. Total LINEAR
        // momentum is conserved exactly at every restitution, because nothing but
        // the equal-and-opposite contact impulse acts (no gravity, no ground). A
        // skipped or one-sided resolve would break it, and the "it collided"
        // check rules out a resolve that quietly did nothing.
        for (const Real restitution : {0.0, 0.5, 1.0})
        {
            Rigid3DSettings st;
            st.restitution = restitution; // friction defaults to 0
            const std::vector<RigidBody3D> start{
                cube(2.0, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}),
                cube(1.0, Vec3{0.9, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0})};
            const Vec3 p0 = total_linear_momentum3d(start); // (1, 0, 0)
            Rigid3DWorld world{SimulationSettings{0.002}, start, st};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                total_linear_momentum3d(world.bodies()), p0, 1e-10));
            // The heavier body was slowed by the hit, so a real impulse fired.
            MALLOY_CHECK_TRUE(world.bodies()[0].velocity.x < 1.0 - 1e-6);
        }

        // Frictionless restitution 1 conserves TOTAL kinetic energy (translation
        // plus any rotation); restitution below 1 strictly removes some. Kept
        // short so the torque-free integrator's slow energy drift stays far under
        // the tolerance.
        {
            Rigid3DSettings elastic_st;
            elastic_st.restitution = 1.0;
            Rigid3DSettings lossy_st;
            lossy_st.restitution = 0.4;
            const std::vector<RigidBody3D> start{
                cube(1.0, Vec3{0.0, 0.0, 0.0}, Vec3{1.5, 0.0, 0.0}),
                cube(1.0, Vec3{0.9, 0.0, 0.0}, Vec3{-1.5, 0.0, 0.0})};
            const Real k0 = total_kinetic_energy3d(start);
            Rigid3DWorld elastic{SimulationSettings{0.002}, start, elastic_st};
            Rigid3DWorld lossy{SimulationSettings{0.002}, start, lossy_st};
            for (int i = 0; i < 30; ++i)
            {
                MALLOY_CHECK_TRUE(elastic.step().ok());
                MALLOY_CHECK_TRUE(lossy.step().ok());
            }
            MALLOY_CHECK_NEAR(total_kinetic_energy3d(elastic.bodies()), k0, 1e-6);
            MALLOY_CHECK_TRUE(total_kinetic_energy3d(lossy.bodies()) < k0 - 0.1);
        }

        // An edge-edge approach also resolves through the same core: a is
        // axis-aligned, b is tilted 45 degrees about (1,1,0) and driven back down
        // that diagonal into a. Total momentum is still conserved exactly, and
        // both bodies are pushed apart, so the edge branch of the contact point
        // is exercised in the solver and not just in the geometry test.
        {
            Rigid3DSettings st;
            st.restitution = 0.5;
            const Vec3 diag = normalize(Vec3{1.0, 1.0, 0.0});
            const std::vector<RigidBody3D> start{
                cube(1.0, Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 0.0}),
                cube(1.0, diag * 1.25, diag * -1.0,
                     from_axis_angle(Vec3{1.0, 1.0, 0.0}, 0.78539816339744830961))};
            const Vec3 p0 = total_linear_momentum3d(start);
            Rigid3DWorld world{SimulationSettings{0.002}, start, st};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                total_linear_momentum3d(world.bodies()), p0, 1e-10));
            // b was moving toward a along -diag; after the bounce it no longer is.
            MALLOY_CHECK_TRUE(malloy::math::dot(world.bodies()[1].velocity, diag) > -1.0 + 1e-6);
        }

        // (A box against a SPHERE was deferred at M33; M34 resolves it, and the
        // M34 block below covers that mix.)

        // Two boxes clear of each other never touch: no impulse, no change over a
        // long run (a guard against a false-positive overlap).
        {
            Rigid3DSettings st;
            st.restitution = 1.0;
            Rigid3DWorld world{SimulationSettings{0.002},
                               {cube(1.0, Vec3{0.0, 0.0, 0.0}, Vec3{}),
                                cube(1.0, Vec3{3.0, 0.0, 0.0}, Vec3{})},
                               st};
            for (int i = 0; i < 100; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_NEAR(malloy::math::length(world.bodies()[0].velocity), 0.0, 1e-12);
            MALLOY_CHECK_NEAR(malloy::math::length(world.bodies()[1].velocity), 0.0, 1e-12);
        }
    }

    // --- M34: box against sphere, the last movable pair in 3D, resolved through
    //     the same impulse core. A head-on hit along the line of centres lands
    //     the contact on that line, so unlike the box-box corner it imparts no
    //     spin: total linear momentum stays exact and, at frictionless
    //     restitution 1, so does total kinetic energy. The mix works in either
    //     argument order (the rigid branch flips the normal to run a to b). ---
    {
        using malloy::math::Quat;
        using malloy::math::Vec3;
        using malloy::rigid::mass_properties_3d;
        using malloy::rigid::rigid_body_from;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::SolidBox;
        using malloy::rigid::total_kinetic_energy3d;
        using malloy::rigid::total_linear_momentum3d;

        const Vec3 half{0.5, 0.5, 0.5};
        const auto box_body = [&](Real mass, const Vec3& pos, const Vec3& vel) {
            const Real density = mass / (8.0 * half.x * half.y * half.z);
            RigidBody3D b =
                rigid_body_from(mass_properties_3d(SolidBox{Vec3{}, half, Quat{}, density}));
            b.half_extents = half;
            b.position = pos;
            b.velocity = vel;
            return b;
        };
        const auto sphere_body = [&](Real mass, Real radius, const Vec3& pos,
                                     const Vec3& vel) {
            RigidBody3D b;
            b.mass = mass;
            const Real i = 0.4 * mass * radius * radius; // solid sphere
            b.inertia = Vec3{i, i, i};
            b.radius = radius;
            b.position = pos;
            b.velocity = vel;
            return b;
        };

        // Head-on: a box (mass 2, moving +x) into a sphere (mass 1, moving -x),
        // overlapping 0.1 along x. Total momentum is +1 and stays there at every
        // restitution; the box is slowed, so a real impulse fired.
        for (const Real restitution : {0.0, 0.5, 1.0})
        {
            Rigid3DSettings st;
            st.restitution = restitution;
            const std::vector<RigidBody3D> start{
                box_body(2.0, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}),
                sphere_body(1.0, 0.5, Vec3{0.9, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0})};
            const Vec3 p0 = total_linear_momentum3d(start); // (1, 0, 0)
            Rigid3DWorld world{SimulationSettings{0.002}, start, st};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                total_linear_momentum3d(world.bodies()), p0, 1e-10));
            MALLOY_CHECK_TRUE(world.bodies()[0].velocity.x < 1.0 - 1e-6);
        }

        // A head-on hit runs through both centres, so no spin is imparted:
        // frictionless restitution 1 conserves total kinetic energy exactly, and
        // restitution below 1 strictly loses some.
        {
            Rigid3DSettings elastic_st;
            elastic_st.restitution = 1.0;
            Rigid3DSettings lossy_st;
            lossy_st.restitution = 0.4;
            const std::vector<RigidBody3D> start{
                box_body(1.0, Vec3{0.0, 0.0, 0.0}, Vec3{1.5, 0.0, 0.0}),
                sphere_body(1.0, 0.5, Vec3{0.9, 0.0, 0.0}, Vec3{-1.5, 0.0, 0.0})};
            const Real k0 = total_kinetic_energy3d(start);
            Rigid3DWorld elastic{SimulationSettings{0.002}, start, elastic_st};
            Rigid3DWorld lossy{SimulationSettings{0.002}, start, lossy_st};
            for (int i = 0; i < 30; ++i)
            {
                MALLOY_CHECK_TRUE(elastic.step().ok());
                MALLOY_CHECK_TRUE(lossy.step().ok());
            }
            MALLOY_CHECK_NEAR(total_kinetic_energy3d(elastic.bodies()), k0, 1e-9);
            MALLOY_CHECK_TRUE(total_kinetic_energy3d(lossy.bodies()) < k0 - 0.1);
        }

        // The mix resolves in either argument order: here the SPHERE is body a
        // and the box is body b, which exercises the normal flip in the rigid
        // branch. The head-on outcome is exact: sphere mass 2 at +1 into box
        // mass 1 at -1 with restitution 0.5 leaves the sphere at rest and the box
        // at +1. Pinning those (not merely "it slowed") is what catches a normal
        // that was not flipped to run a toward b, which would tunnel the bodies
        // through and settle them at the wrong velocities.
        {
            Rigid3DSettings st;
            st.restitution = 0.5;
            const std::vector<RigidBody3D> start{
                sphere_body(2.0, 0.5, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}),
                box_body(1.0, Vec3{0.9, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0})};
            const Vec3 p0 = total_linear_momentum3d(start);
            Rigid3DWorld world{SimulationSettings{0.002}, start, st};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                total_linear_momentum3d(world.bodies()), p0, 1e-10));
            MALLOY_CHECK_TRUE(
                malloy::math::approx_equal(world.bodies()[0].velocity, Vec3{0.0, 0.0, 0.0}, 1e-9));
            MALLOY_CHECK_TRUE(
                malloy::math::approx_equal(world.bodies()[1].velocity, Vec3{1.0, 0.0, 0.0}, 1e-9));
        }

        // A box and a sphere clear of each other never touch: no change.
        {
            Rigid3DSettings st;
            st.restitution = 1.0;
            Rigid3DWorld world{SimulationSettings{0.002},
                               {box_body(1.0, Vec3{0.0, 0.0, 0.0}, Vec3{}),
                                sphere_body(1.0, 0.5, Vec3{3.0, 0.0, 0.0}, Vec3{})},
                               st};
            for (int i = 0; i < 100; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_NEAR(malloy::math::length(world.bodies()[0].velocity), 0.0, 1e-12);
            MALLOY_CHECK_NEAR(malloy::math::length(world.bodies()[1].velocity), 0.0, 1e-12);
        }
    }

    // --- Review pass after M37: box-box contact at the manifold CENTROID. Two
    //     aligned boxes meeting face-on must give exactly the one-dimensional
    //     answer, because the face contact's net impulse acts along the line of
    //     centres. M33 and M36 pushed at a single face corner instead, which
    //     invented a lever arm: this setup ended at kinetic energy 0.320 (3D) and
    //     0.425 (2D) instead of 0.5, with energy leaking into spin that should
    //     not exist. Masses 2 and 1 closing at +-1 with restitution 0.5 must
    //     leave the heavy box at rest and the light one moving at +1. ---
    {
        using malloy::math::Quat;
        using malloy::math::Vec2;
        using malloy::math::Vec3;
        using malloy::rigid::mass_properties_3d;
        using malloy::rigid::rigid_body_from;
        using malloy::rigid::RigidBody3D;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::SolidBox;

        const auto square = [](Real mass, const Vec2& position, const Vec2& velocity) {
            RigidBody2D b;
            b.mass = mass;
            b.inertia = mass * (0.25 + 0.25) / 3.0; // solid unit square
            b.half_extents = Vec2{0.5, 0.5};
            b.position = position;
            b.velocity = velocity;
            return b;
        };
        const auto cube = [](Real mass, const Vec3& position, const Vec3& velocity) {
            const Vec3 half{0.5, 0.5, 0.5};
            RigidBody3D b = rigid_body_from(
                mass_properties_3d(SolidBox{Vec3{}, half, Quat{}, mass / 1.0}));
            b.half_extents = half;
            b.position = position;
            b.velocity = velocity;
            return b;
        };

        // 2D, aligned: the exact 1D outcome, with no spin at all.
        {
            RigidWorld w{SimulationSettings{0.002},
                         {square(2.0, Vec2{0.0, 0.0}, Vec2{1.0, 0.0}),
                          square(1.0, Vec2{0.9, 0.0}, Vec2{-1.0, 0.0})},
                         RigidSettings{0.5}};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_VEC2_NEAR(w.bodies()[0].velocity, Vec2(0.0, 0.0), 1e-12);
            MALLOY_CHECK_VEC2_NEAR(w.bodies()[1].velocity, Vec2(1.0, 0.0), 1e-12);
            MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 1e-12);
            MALLOY_CHECK_NEAR(w.bodies()[1].angular_velocity, 0.0, 1e-12);
            MALLOY_CHECK_NEAR(total_kinetic_energy(w.bodies()), 0.5, 1e-12);
        }

        // 2D, offset by 0.3 across the hit: the overlap is lopsided, so the
        // centroid sits off the line of centres and the boxes MUST spin, while
        // linear momentum stays exact. (A contact pinned to the line of centres
        // would wrongly keep them still.)
        {
            const std::vector<RigidBody2D> start{
                square(2.0, Vec2{0.0, 0.0}, Vec2{1.0, 0.0}),
                square(1.0, Vec2{0.9, 0.3}, Vec2{-1.0, 0.0})};
            const Vec2 p0 = total_linear_momentum(start);
            RigidWorld w{SimulationSettings{0.002}, start, RigidSettings{0.5}};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_VEC2_NEAR(total_linear_momentum(w.bodies()), p0, 1e-10);
            MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angular_velocity) > 1e-3);
        }

        // 3D, aligned: the same exact 1D outcome, with no spin on any axis.
        {
            Rigid3DSettings st;
            st.restitution = 0.5;
            Rigid3DWorld w{SimulationSettings{0.002},
                           {cube(2.0, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}),
                            cube(1.0, Vec3{0.9, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0})},
                           st};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_TRUE(
                malloy::math::approx_equal(w.bodies()[0].velocity, Vec3{0.0, 0.0, 0.0}, 1e-12));
            MALLOY_CHECK_TRUE(
                malloy::math::approx_equal(w.bodies()[1].velocity, Vec3{1.0, 0.0, 0.0}, 1e-12));
            MALLOY_CHECK_NEAR(malloy::math::length(w.bodies()[0].angular_velocity), 0.0, 1e-12);
            MALLOY_CHECK_NEAR(malloy::math::length(w.bodies()[1].angular_velocity), 0.0, 1e-12);
            MALLOY_CHECK_NEAR(malloy::rigid::total_kinetic_energy3d(w.bodies()), 0.5, 1e-12);
        }

        // 3D, offset across the hit: the boxes must spin, momentum stays exact.
        {
            Rigid3DSettings st;
            st.restitution = 0.5;
            const std::vector<RigidBody3D> start{
                cube(2.0, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}),
                cube(1.0, Vec3{0.9, 0.3, 0.0}, Vec3{-1.0, 0.0, 0.0})};
            const Vec3 p0 = malloy::rigid::total_linear_momentum3d(start);
            Rigid3DWorld w{SimulationSettings{0.002}, start, st};
            for (int i = 0; i < 60; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::rigid::total_linear_momentum3d(w.bodies()), p0, 1e-10));
            MALLOY_CHECK_TRUE(malloy::math::length(w.bodies()[0].angular_velocity) > 1e-3);
        }
    }

    // --- Issues #22 and #23: the head-on and flat-landing cases above, turned to
    //     GENERIC orientations. Axis-aligned, every projection is exact; turned,
    //     a bare touch rounded clipped points off the face (a one-corner or
    //     lopsided manifold), clipping added a vertex partway along an edge, and
    //     a twisted face-on-face landing let an in-plane edge axis win. Each put
    //     torque into a hit that has none. ---
    {
        using malloy::math::from_axis_angle;
        using malloy::math::normalize;
        using malloy::math::Quat;
        using malloy::math::rotate;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;
        using malloy::math::Vec3;

        // 2D: box_collision2d's head-on hit (masses 2 and 1, restitution 0.5),
        // both squares turned and moving along their own x axis. The heavy one
        // stops, the light one leaves at 1, no spin, kinetic energy 0.5.
        for (int k = 0; k < 200; ++k)
        {
            const Real angle = 0.1 + 0.007 * k;
            const Vec2 u{std::cos(angle), std::sin(angle)};
            RigidBody2D a;
            a.mass = 2.0;
            a.inertia = 2.0 * 2.0 / 12.0;
            a.half_extents = Vec2{0.5, 0.5};
            a.angle = angle;
            a.position = u * -0.6;
            a.velocity = u;
            RigidBody2D b = a;
            b.mass = 1.0;
            b.inertia = 2.0 / 12.0;
            b.position = u * 0.6;
            b.velocity = u * -1.0;
            RigidWorld w{SimulationSettings{0.001}, {a, b}, RigidSettings{0.5}};
            for (int i = 0; i < 300; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_NEAR(w.bodies()[0].angular_velocity, 0.0, 1e-9);
            MALLOY_CHECK_NEAR(w.bodies()[1].angular_velocity, 0.0, 1e-9);
            MALLOY_CHECK_NEAR(malloy::rigid::total_kinetic_energy(w.bodies()), 0.5, 1e-9);
        }

        // 3D: box_collision's head-on hit, turned about a generic axis. 0.114,
        // 0.702 and 1.479 are angles that spun before the fix.
        for (const Real angle : {0.114, 0.3, 0.702, 1.0, 1.479, 2.2})
        {
            const Quat r = from_axis_angle(Vec3{0.3, -0.7, 0.2}, angle);
            const Vec3 u = rotate(r, Vec3{1.0, 0.0, 0.0});
            const Real moment = 2.0 * (0.25 + 0.25) / 3.0;
            RigidBody3D a;
            a.mass = 2.0;
            a.inertia = Vec3{moment, moment, moment};
            a.half_extents = Vec3{0.5, 0.5, 0.5};
            a.orientation = r;
            a.position = u * -1.2;
            a.velocity = u;
            RigidBody3D b = a;
            b.mass = 1.0;
            b.inertia = a.inertia * 0.5;
            b.position = u * 1.2;
            b.velocity = u * -1.0;
            Rigid3DSettings st;
            st.restitution = 0.5;
            Rigid3DWorld w{SimulationSettings{0.002}, {a, b}, st};
            for (int i = 0; i < 400; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_NEAR(malloy::math::length(w.bodies()[0].angular_velocity), 0.0, 1e-9);
            MALLOY_CHECK_NEAR(malloy::math::length(w.bodies()[1].angular_velocity), 0.0, 1e-9);
            MALLOY_CHECK_NEAR(malloy::rigid::total_kinetic_energy3d(w.bodies()), 0.5, 1e-9);
        }

        // 3D: a small box dropped flat onto a big one, turned about the shared
        // normal by every whole degree from 1 to 89, both boxes in a generic
        // orientation. Before the fix some twists resolved at one in-plane
        // edge point up to 0.9 away and left b spinning at nearly 2 rad/s.
        //
        // M40 resolves all four corners of b's face (issue #73), at the
        // default restitution 1. They approach a's face at the same speed, and
        // the per-point Newton law sends each away at that speed, so the
        // relative normal velocity is the same +1 at four points of a plane:
        // the RELATIVE angular velocity has no in-plane part, b turning with
        // a. (Before M40 the one contact under b's centre of mass left b with
        // no spin at all while a spun: the four corners did not leave
        // together.) Every point reversing at e = 1 conserves kinetic energy
        // exactly, sum lambda (u- + u+) / 2 = 0, and the impulses are equal
        // and opposite, so momentum is kept. Checked at the step that resolves
        // the landing, which ends with each body's free rotation over dt:
        //
        //   - Each body's normal is one of its principal axes (z), so a normal
        //     impulse gives it no spin about it: w = (wx, wy, 0) in each body
        //     frame. Then L x w lies along the normal, so the free rotation
        //     moves the in-plane spin only at SECOND order,
        //     (dt^2 / 2) d2w/dt2, with, in the body frame and u = L x w,
        //
        //       d2w/dt2 = w x I^-1 u - I^-1 (w x u) + I^-1 (L x I^-1 u),
        //
        //     the derivative of dw/dt = I^-1 u along the motion (the world
        //     inertia turns with the body). Forward Euler, a first-order map,
        //     had no such term, and this check bounded the raw in-plane part by
        //     the impact tolerance; the exact flow moves it by up to 2.5e-10
        //     here (issue #88). The term is subtracted for both bodies, and
        //     what is left is held to the same bound as before: corners at
        //     least 0.6 apart, each within 1e-12 of the velocity scale (1) of
        //     its target, give 2e-12 / 0.6 per axis, with the third-order
        //     remainder, (dt |w|)^3 |w| ~ 4e-13 at |w| = 0.14, inside it.
        //   - b is axisymmetric about its x axis (Iy = Iz), so its rotation is
        //     the exact precession: wx kept, (wy, wz) turned about x by
        //     Omega dt, Omega = (Iy - Ix) / Iy wx. From wz = 0 that leaves
        //     wz = -wy tan(Omega dt) after the step (forward Euler's first-order
        //     dt (Ix - Iy) / Iz wx wy differs from it by 1e-15).
        //   - The energy to the impact tolerance times the impulse (about
        //     1.8e-12): b's rotation keeps its energy exactly, and a's changes
        //     it by the splitting's third-order error, (dt |w|)^3 T ~ 4e-14.
        {
            const Quat r = from_axis_angle(Vec3{0.3, -0.7, 0.2}, 0.9);
            const Vec3 up = rotate(r, Vec3{0.0, 0.0, 1.0});
            const Real dt = 0.001;
            // (dt^2 / 2) d2w/dt2 of the torque-free motion, in the world frame.
            const auto precession = [dt](const RigidBody3D& body) {
                const Vec3& w = body.angular_velocity;
                const Vec3& m = body.inertia;
                const auto inverse = [&m](const Vec3& v) { return Vec3{v.x / m.x, v.y / m.y, v.z / m.z}; };
                const Vec3 l{m.x * w.x, m.y * w.y, m.z * w.z};
                const Vec3 u = malloy::math::cross(l, w);
                const Vec3 rate = inverse(u);
                const Vec3 second = malloy::math::cross(w, rate) - inverse(malloy::math::cross(w, u)) +
                                    inverse(malloy::math::cross(l, rate));
                return rotate(body.orientation, second) * (0.5 * dt * dt);
            };
            for (int degrees = 1; degrees < 90; ++degrees)
            {
                const Real twist = degrees * 3.14159265358979323846 / 180.0;
                RigidBody3D a;
                a.mass = 10.0;
                a.inertia = Vec3{1.0, 1.5, 2.0};
                a.half_extents = Vec3{1.0, 0.7, 0.5};
                a.orientation = r;
                RigidBody3D b;
                b.mass = 1.0;
                b.inertia = Vec3{0.06, 0.083, 0.083};
                b.half_extents = Vec3{0.4, 0.3, 0.3};
                b.orientation = normalize(r * from_axis_angle(Vec3{0.0, 0.0, 1.0}, twist));
                b.position = rotate(r, Vec3{0.1, 0.05, 0.0}) + up * 0.8053;
                b.velocity = up * -1.0;
                Rigid3DWorld w{SimulationSettings{dt}, {a, b}};
                const Real energy = malloy::rigid::total_kinetic_energy3d(w.bodies());
                const Vec3 momentum = malloy::rigid::total_linear_momentum3d(w.bodies());
                int steps = 0;
                while (steps < 40 && malloy::math::length(w.bodies()[0].velocity) == 0.0)
                {
                    MALLOY_CHECK_TRUE(w.step().ok());
                    ++steps;
                }
                MALLOY_CHECK_TRUE(steps < 40); // it landed
                const RigidBody3D& big = w.bodies()[0];
                const RigidBody3D& small = w.bodies()[1];
                MALLOY_CHECK_TRUE(malloy::math::dot(small.velocity - big.velocity, up) > 0.9);
                const Vec3 relative = rotate(small.orientation, small.angular_velocity) -
                                      rotate(big.orientation, big.angular_velocity) -
                                      (precession(small) - precession(big));
                MALLOY_CHECK_NEAR(malloy::math::length(relative - up * malloy::math::dot(relative, up)),
                                  0.0, 1e-11);
                const Vec3& ws = small.angular_velocity;
                const Vec3& is = small.inertia;
                MALLOY_CHECK_NEAR(ws.z, -ws.y * std::tan(dt * (is.y - is.x) / is.y * ws.x), 1e-12);
                MALLOY_CHECK_NEAR(malloy::rigid::total_kinetic_energy3d(w.bodies()), energy, 1e-10);
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                    malloy::rigid::total_linear_momentum3d(w.bodies()), momentum, 1e-14));
            }
        }
    }

    // --- Issue #41 (and #68): a box contact manifold was reduced to ONE
    //     contact at its PENETRATION-WEIGHTED centroid, with the deepest
    //     penetration. The plain mean of a tilted face's corners sits toward
    //     the raised side of the centre of mass, so the impulse tipped the box
    //     further, and with the translation-only correction lifting it every
    //     step a box on the floor rocked forever and gained energy even at
    //     restitution 0. M40 (issue #73) replaces the reduction: every
    //     manifold point gets its own impulse, and the weighted centroid
    //     survives only as the prepass's first guess. These tests were
    //     re-derived for that: tilted landings now stop every loaded point,
    //     and the issue's energy bound below still holds. ---
    {
        using malloy::math::conjugate;
        using malloy::math::cross;
        using malloy::math::dot;
        using malloy::math::from_axis_angle;
        using malloy::math::normalize;
        using malloy::math::Quat;
        using malloy::math::rotate;
        using malloy::math::Vec3;
        using malloy::rigid::box_inertia;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;

        // The world-frame inverse inertia applied to x, R diag(1/I) R^T x,
        // written out here rather than borrowed from the library under test.
        const auto inverse_inertia_world = [](const RigidBody3D& b, const Vec3& x) {
            const Vec3 local = rotate(conjugate(b.orientation), x);
            return rotate(b.orientation,
                          Vec3{local.x / b.inertia.x, local.y / b.inertia.y,
                               local.z / b.inertia.z});
        };
        // --- 2D, a tilted landing on the ground. No gravity, restitution 0 and
        //     no friction, so the step is one rest-phase solve that drives each
        //     loaded corner's normal velocity to 0. The box starts with both
        //     bottom corners in, about 0.005 and 0.025 deep (the pose the single
        //     weighted contact was pinned at before M40), moving down at 1.
        //     Tilted clockwise, so the deeper corner is the second in the
        //     corner order. Two loaded corners at different x fix both
        //     unknowns: v_y + omega x_k = 0 for k = 1, 2 gives v_y = omega = 0.
        //     Both ARE loaded: stopping the box without turning it takes
        //     lambda_k = m x_other / (x_other - x_k) at each, positive because
        //     the centre of mass lies between the corners. So the box stops
        //     dead, where the single contact left it falling at -1 + j / m and
        //     spinning at 0.36 rad/s. The bound: each corner ends within
        //     r = 1e-8 v_scale of 0, and v_scale (the speed after the centroid
        //     prepass, which only removes energy here, so |v| <= 1 and
        //     |omega| <= sqrt(m / I)) is at most 1 + sqrt(m / I) h, h = 0.5.
        //     With the corners x_2 - x_1 apart on either side of the centre of
        //     mass, |v_y| <= r and |omega| <= 2 r / (x_2 - x_1). The projection
        //     then rights it, positions only. ---
        {
            const Real angle = -0.02;
            const Vec2 half{0.5, 0.25};
            const Real mass = 1.5;
            const Real inertia = mass * (half.x * half.x + half.y * half.y) / 3.0;
            RigidBody2D box;
            box.mass = mass;
            box.inertia = inertia;
            box.half_extents = half;
            box.angle = angle;
            box.position = Vec2{0.3, 0.235};
            box.velocity = Vec2{0.0, -1.0};
            RigidSettings settings;
            settings.restitution = 0.0;
            settings.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});
            RigidWorld w{SimulationSettings{0.01}, {box}, settings};

            const auto depth = [&half](const RigidBody2D& b, Real sx) {
                const Vec2 u0{std::cos(b.angle), std::sin(b.angle)};
                const Vec2 u1{-std::sin(b.angle), std::cos(b.angle)};
                return -(b.position + u0 * (sx * half.x) + u1 * (-half.y)).y;
            };
            MALLOY_CHECK_TRUE(depth(box, -1.0) > 0.0); // both bottom corners are in
            MALLOY_CHECK_TRUE(depth(box, 1.0) - depth(box, -1.0) > 0.015); // genuinely tilted
            const Real deepest = depth(box, 1.0);

            MALLOY_CHECK_TRUE(w.step().ok());
            const RigidBody2D& b = w.bodies()[0];
            const Real r = 1e-8 * (1.0 + std::sqrt(mass / inertia) * half.x);
            MALLOY_CHECK_NEAR(b.velocity.x, 0.0, 0.0); // frictionless: normal only
            MALLOY_CHECK_NEAR(b.velocity.y, 0.0, r);
            MALLOY_CHECK_NEAR(b.angular_velocity, 0.0, 2.0 * r / (2.0 * half.x * std::cos(angle)));
            // The projection rights it: less tilt, and the deep corner lifted.
            MALLOY_CHECK_TRUE(std::abs(b.angle) < std::abs(angle));
            MALLOY_CHECK_TRUE(std::fmax(depth(b, -1.0), depth(b, 1.0)) < deepest);
        }

        // --- 3D, the same tilted landing on a plane: a box with three distinct
        //     half-widths, tilted 0.03 rad about a horizontal axis that is not a
        //     box axis, started with all four bottom corners in at four
        //     different depths (about 0.006 to 0.044). Four loaded corners of a
        //     plane fix v_z and both horizontal spin components, so the box
        //     stops dead (its centre of mass is over its face, so the corner
        //     impulses that do it are all positive), where the one contact at
        //     the weighted centroid left it spinning. Frictionless on a level
        //     plane, so the horizontal velocity is untouched exactly. The
        //     bound: each corner within r = 1e-8 v_scale of 0, v_scale <=
        //     1 + sqrt(m / I_min) h as in 2D; corners 2 x 0.5 by 2 x 0.4 apart
        //     bound the horizontal spin by r sqrt(1/0.5^2 + 1/0.4^2) < 4 r,
        //     and v_z by r plus that spin times the bottom face's offset from
        //     under the centre of mass (0.25 sin 0.03): within 2 r. ---
        {
            const Vec3 half{0.5, 0.4, 0.25};
            RigidBody3D box;
            box.mass = 2.0;
            box.inertia = box_inertia(half, box.mass);
            box.half_extents = half;
            box.orientation = from_axis_angle(normalize(Vec3{1.0, 2.0, 0.0}), 0.03);
            box.position = Vec3{0.2, -0.1, 0.225};
            box.velocity = Vec3{0.0, 0.0, -1.0};
            Rigid3DSettings settings;
            settings.restitution = 0.0;
            settings.ground.push_back(malloy::collide::Plane3{Vec3{0.0, 0.0, 1.0}, 0.0});
            Rigid3DWorld w{SimulationSettings{0.01}, {box}, settings};
            const auto deepest = [&half](const RigidBody3D& b, int& touching, Real& shallowest) {
                Real most = 0.0;
                touching = 0;
                shallowest = 1.0;
                for (int i = 0; i < 8; ++i)
                {
                    const Vec3 local{(i & 1) ? half.x : -half.x, (i & 2) ? half.y : -half.y,
                                     (i & 4) ? half.z : -half.z};
                    const Real z = (b.position + rotate(b.orientation, local)).z;
                    if (z <= 0.0)
                    {
                        ++touching;
                        most = std::fmax(most, -z);
                        shallowest = std::fmin(shallowest, -z);
                    }
                }
                return most;
            };
            int touching = 0;
            Real shallowest = 0.0;
            const Real before = deepest(box, touching, shallowest);
            MALLOY_CHECK_EQ(touching, 4);
            MALLOY_CHECK_TRUE(before - shallowest > 0.03); // four distinct depths

            MALLOY_CHECK_TRUE(w.step().ok());
            const RigidBody3D& b = w.bodies()[0];
            const Real r = 1e-8 * (1.0 + std::sqrt(box.mass / box.inertia.x) * half.x);
            MALLOY_CHECK_NEAR(b.velocity.x, 0.0, 0.0);
            MALLOY_CHECK_NEAR(b.velocity.y, 0.0, 0.0);
            MALLOY_CHECK_NEAR(b.velocity.z, 0.0, 2.0 * r);
            MALLOY_CHECK_NEAR(malloy::math::length(b.angular_velocity), 0.0, 4.0 * r);
            MALLOY_CHECK_TRUE(deepest(b, touching, shallowest) < before); // projected out
        }

        // --- 2D, box against box: a small tilted box landing on a wide one,
        //     both movable, no gravity, restitution 0, no friction, started
        //     where the collide layer gives two points (the small box's bottom
        //     corners, about 0.010 and 0.050 deep). Each loaded point ends at a
        //     relative normal velocity of 0, so the boxes leave the step
        //     turning together (the one weighted contact left the two points
        //     moving apart and together at once). The bound: within
        //     r = 1e-8 v_scale, v_scale the larger of each body's |v| +
        //     |omega| h, at most what the kinetic energy 0.5 allows it. The
        //     impulses are equal and opposite, so momentum is kept to
        //     rounding; the projection lifts the deeper corner. ---
        {
            const Real dt = 0.01;
            const Real angle = 0.05;
            RigidBody2D lower;
            lower.mass = 3.0;
            lower.half_extents = Vec2{1.0, 0.25};
            lower.inertia = 3.0 * (1.0 + 0.0625) / 3.0;
            RigidBody2D upper;
            upper.mass = 1.0;
            upper.half_extents = Vec2{0.4, 0.2};
            upper.inertia = 1.0 * (0.16 + 0.04) / 3.0;
            upper.angle = angle;
            upper.position = Vec2{0.3, 0.419742};
            upper.velocity = Vec2{0.0, -1.0};
            RigidWorld w{SimulationSettings{dt}, {lower, upper}, RigidSettings{0.0}};

            const std::vector<malloy::collide::Contact> hits = malloy::collide::contacts(
                malloy::collide::Obb2{lower.position, lower.half_extents, lower.angle},
                malloy::collide::Obb2{upper.position, upper.half_extents, upper.angle});
            MALLOY_CHECK_EQ(hits.size(), std::size_t{2});
            MALLOY_CHECK_TRUE(std::abs(hits[0].penetration - hits[1].penetration) > 0.03);
            const Real before = std::fmax(hits[0].penetration, hits[1].penetration);
            const Vec2 momentum = lower.velocity * lower.mass + upper.velocity * upper.mass;

            MALLOY_CHECK_TRUE(w.step().ok());
            const RigidBody2D& a = w.bodies()[0];
            const RigidBody2D& b = w.bodies()[1];
            const Real scale = std::fmax(1.0 + std::sqrt(1.0 / upper.inertia) * 0.4,
                                         std::sqrt(1.0 / lower.mass) +
                                             std::sqrt(1.0 / lower.inertia) * 1.0);
            const Real r = 1e-8 * scale;
            for (const malloy::collide::Contact& hit : hits)
            {
                const Vec2 arm_a = hit.point - lower.position;
                const Vec2 arm_b = hit.point - upper.position;
                const Vec2 va = a.velocity + Vec2{-a.angular_velocity * arm_a.y,
                                                  a.angular_velocity * arm_a.x};
                const Vec2 vb = b.velocity + Vec2{-b.angular_velocity * arm_b.y,
                                                  b.angular_velocity * arm_b.x};
                MALLOY_CHECK_NEAR(dot(vb - va, hit.normal), 0.0, r);
            }
            MALLOY_CHECK_VEC2_NEAR(a.velocity * a.mass + b.velocity * b.mass, momentum, 1e-15);
            const std::vector<malloy::collide::Contact> after = malloy::collide::contacts(
                malloy::collide::Obb2{a.position, a.half_extents, a.angle},
                malloy::collide::Obb2{b.position, b.half_extents, b.angle});
            for (const malloy::collide::Contact& hit : after)
            {
                MALLOY_CHECK_TRUE(hit.penetration < before);
            }
        }

        // --- 3D, box against box: a small box tilted about a generic axis on a
        //     wide flat one, four manifold points at four depths (about 0.0007
        //     to 0.040), both movable, no gravity, restitution 0, no friction.
        //     Each loaded point ends at a relative normal velocity of 0, bound
        //     as in 2D. To read the solved spins back the step is tiny, dt 1e-9:
        //     every point is in and nothing forces the bodies, so the contact
        //     solve does not depend on dt; the projection keeps each body's
        //     world angular momentum L through its turns; and one integration
        //     step moves L by a negligible |omega| dt. So the solved world spin
        //     is I_world^-1 L at the start pose. ---
        {
            RigidBody3D lower;
            lower.mass = 4.0;
            lower.half_extents = Vec3{1.2, 1.0, 0.3};
            lower.inertia = box_inertia(lower.half_extents, lower.mass);
            lower.orientation = from_axis_angle(Vec3{0.0, 0.0, 1.0}, 0.3);
            RigidBody3D upper;
            upper.mass = 1.0;
            upper.half_extents = Vec3{0.4, 0.3, 0.2};
            upper.inertia = box_inertia(upper.half_extents, upper.mass);
            upper.orientation = from_axis_angle(normalize(Vec3{1.0, 2.0, 0.0}), 0.04);
            upper.position = Vec3{0.2, -0.15, 0.479512};
            upper.velocity = Vec3{0.0, 0.0, -1.0};
            Rigid3DSettings settings;
            settings.restitution = 0.0;
            Rigid3DWorld w{SimulationSettings{1e-9}, {lower, upper}, settings};
            const std::vector<malloy::collide::Contact3> hits = malloy::collide::contacts(
                malloy::collide::Box3{lower.position, lower.half_extents, lower.orientation},
                malloy::collide::Box3{upper.position, upper.half_extents, upper.orientation});
            MALLOY_CHECK_EQ(hits.size(), std::size_t{4});
            Real deepest = 0.0;
            Real shallowest = 1.0;
            for (const malloy::collide::Contact3& hit : hits)
            {
                deepest = std::fmax(deepest, hit.penetration);
                shallowest = std::fmin(shallowest, hit.penetration);
            }
            MALLOY_CHECK_TRUE(deepest - shallowest > 0.03);
            const Vec3 momentum = lower.velocity * lower.mass + upper.velocity * upper.mass;

            MALLOY_CHECK_TRUE(w.step().ok());
            const RigidBody3D& a = w.bodies()[0];
            const RigidBody3D& b = w.bodies()[1];
            const auto angular_momentum = [](const RigidBody3D& body) {
                const Vec3& s = body.angular_velocity;
                return rotate(body.orientation, Vec3{body.inertia.x * s.x, body.inertia.y * s.y,
                                                     body.inertia.z * s.z});
            };
            const Vec3 spin_a = inverse_inertia_world(lower, angular_momentum(a));
            const Vec3 spin_b = inverse_inertia_world(upper, angular_momentum(b));
            const Real lightest_a = std::fmin(std::fmin(lower.inertia.x, lower.inertia.y),
                                              lower.inertia.z);
            const Real lightest_b = std::fmin(std::fmin(upper.inertia.x, upper.inertia.y),
                                              upper.inertia.z);
            const Real r = 1e-8 * std::fmax(1.0 + std::sqrt(1.0 / lightest_b) * 0.4,
                                            std::sqrt(1.0 / lower.mass) +
                                                std::sqrt(1.0 / lightest_a) * 1.2);
            for (const malloy::collide::Contact3& hit : hits)
            {
                const Vec3 va = a.velocity + cross(spin_a, hit.point - lower.position);
                const Vec3 vb = b.velocity + cross(spin_b, hit.point - upper.position);
                MALLOY_CHECK_NEAR(dot(vb - va, hit.normal), 0.0, r);
            }
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                a.velocity * a.mass + b.velocity * b.mass, momentum, 1e-15));
        }

        // --- A bare touch has no weights: every depth is exactly 0, so the
        //     centroid prepass falls back to equal weights instead of dividing
        //     zero by zero. Built from binary fractions so the first step lands
        //     the box's bottom face EXACTLY on the floor while it is still
        //     moving down at 0.25 (M40 finds the contact at the start of the
        //     next step); the second must stop it (restitution 0) with no spin,
        //     the centroid impulse through the centre of mass solving the level
        //     face exactly. Without the fallback the weights are NaN and the box
        //     keeps falling. ---
        {
            RigidBody2D square;
            square.half_extents = Vec2{0.5, 0.5};
            square.inertia = (0.25 + 0.25) / 3.0;
            square.position = Vec2{0.0, 0.625};
            square.velocity = Vec2{0.0, -0.25};
            RigidSettings flat;
            flat.restitution = 0.0;
            flat.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});
            RigidWorld w2{SimulationSettings{0.5}, {square}, flat};
            MALLOY_CHECK_TRUE(w2.step().ok());
            MALLOY_CHECK_EQ(w2.bodies()[0].position.y, 0.5); // exactly touching
            MALLOY_CHECK_EQ(w2.bodies()[0].velocity.y, -0.25);
            MALLOY_CHECK_TRUE(w2.step().ok());
            MALLOY_CHECK_EQ(w2.bodies()[0].position.y, 0.5);
            MALLOY_CHECK_EQ(w2.bodies()[0].velocity.y, 0.0);
            MALLOY_CHECK_EQ(w2.bodies()[0].angular_velocity, 0.0);

            RigidBody3D cube;
            cube.half_extents = Vec3{0.5, 0.5, 0.5};
            cube.inertia = box_inertia(cube.half_extents, cube.mass);
            cube.position = Vec3{0.0, 0.0, 0.625};
            cube.velocity = Vec3{0.0, 0.0, -0.25};
            Rigid3DSettings level;
            level.restitution = 0.0;
            level.ground.push_back(malloy::collide::Plane3{Vec3{0.0, 0.0, 1.0}, 0.0});
            Rigid3DWorld w3{SimulationSettings{0.5}, {cube}, level};
            MALLOY_CHECK_TRUE(w3.step().ok());
            MALLOY_CHECK_EQ(w3.bodies()[0].position.z, 0.5);
            MALLOY_CHECK_EQ(w3.bodies()[0].velocity.z, -0.25);
            MALLOY_CHECK_TRUE(w3.step().ok());
            MALLOY_CHECK_EQ(w3.bodies()[0].position.z, 0.5);
            MALLOY_CHECK_EQ(w3.bodies()[0].velocity.z, 0.0);
            MALLOY_CHECK_EQ(malloy::math::length(w3.bodies()[0].angular_velocity), 0.0);
        }

        // --- The issue's own reproduction, 2D: a box resting on its long face
        //     (half-widths 0.5 by 0.25, mass 1) nudged with spin 0.05, gravity
        //     9.8, restitution 0, no friction, dt 0.002, for 20 s. Every contact
        //     here only removes kinetic energy; the one thing that can ADD
        //     energy is the positional correction, which lifts the box by the
        //     deepest penetration.
        //
        //     The bound. Each step starts with every corner at or above the
        //     floor (the last correction lifted the deepest one exactly onto it),
        //     so a correction can lift the centre of mass by at most how far a
        //     corner sank in that step. From rest on the floor the first step
        //     sinks a corner by at most (g dt + a |w0|) dt, a the corner's
        //     horizontal lever: one step of gravity plus one step of free
        //     rotation, worth B = m g (g dt + a |w0|) dt = 8.74e-4 J. A contact
        //     that only dissipates after that never climbs past E0 + B; one
        //     that pumps does: the plain mean reached E0 + 1.49 B, and the
        //     weighted centroid stays at E0 + 0.65 B. (B bounds the lift of the
        //     FIRST step, not of every later one: a nudge well below one step's
        //     gravity, w0 = 0.02, briefly reaches 1.4 B even with the weighting,
        //     because each step-scale rock gets its own lift. The issue's nudge
        //     is the one pinned.) M40's per-point solver (issue #73) measures
        //     no climb at all: E never exceeds E0.
        //
        //     And it comes to REST: over the last 5 s its kinetic energy stays
        //     below 1e-4 of one step's gravity, m g^2 dt^2 (the weighted
        //     centroid measured 3e-10 J and the plain mean kept it rocking at up
        //     to 8.3e-4 J, still turning at 0.043 rad/s after 20 s; M40
        //     measures 8e-21 J). ---
        {
            const Real g = 9.8;
            const Real dt = 0.002;
            const Real w0 = 0.05;
            const Vec2 half{0.5, 0.25};
            RigidBody2D box;
            box.half_extents = half;
            box.inertia = (half.x * half.x + half.y * half.y) / 3.0; // mass 1
            box.position = Vec2{0.0, half.y};
            box.angular_velocity = w0;
            RigidSettings settings;
            settings.gravity = Vec2{0.0, -g};
            settings.restitution = 0.0;
            settings.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});
            RigidWorld w{SimulationSettings{dt}, {box}, settings};
            // Mechanical energy above the box resting flat, centre at y = h.
            const auto energy = [&] {
                return total_energy(w.bodies(), settings.gravity) - g * half.y;
            };
            const Real e0 = energy();
            const Real bound = g * (g * dt + half.x * w0) * dt;
            const Real step_scale = g * g * dt * dt;
            Real worst = e0;
            Real late_kinetic = 0.0;
            for (int i = 1; i <= 10000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                worst = std::fmax(worst, energy());
                if (i > 7500)
                {
                    late_kinetic = std::fmax(late_kinetic, total_kinetic_energy(w.bodies()));
                }
            }
            MALLOY_CHECK_TRUE(worst - e0 <= bound);
            MALLOY_CHECK_TRUE(late_kinetic < 1e-4 * step_scale);
            MALLOY_CHECK_TRUE(std::abs(w.bodies()[0].angular_velocity) < 1e-3);
        }

        // --- And 3D: a box (0.5, 0.5, 0.25) resting on its face, nudged with
        //     spin (0.3, 0.2, 0), the issue's reproduction. The same bound, with
        //     the lever now the corner's horizontal distance sqrt(a^2 + b^2):
        //     B = m g (g dt + sqrt(a^2 + b^2) |w0|) dt = 5.38e-3 J. The plain mean
        //     climbed to E0 + 1.95 B (2.5 times E0); the weighted centroid peaks
        //     at E0 + 0.06 B; M40 stays at E0.
        //
        //     Before M40 it did NOT come fully to rest: one contact cannot stop
        //     four corners at once, so the box settled into a step-scale rock
        //     (a four-step cycle around its corners, tilted about 1e-4 rad) at
        //     5.07e-4 J, 1.3 m g^2 dt^2 (the plain mean kept 7.52e-3 J). M40's
        //     per-point solver (issue #73) stops every corner, so the 2D rest
        //     bound now holds in 3D too: below 1e-4 m g^2 dt^2 over the last
        //     5 s (measured 6e-21 J). ---
        {
            const Real g = 9.8;
            const Real dt = 0.002;
            const Vec3 half{0.5, 0.5, 0.25};
            const Vec3 w0{0.3, 0.2, 0.0};
            RigidBody3D box;
            box.half_extents = half;
            box.inertia = box_inertia(half, 1.0);
            box.position = Vec3{0.0, 0.0, half.z};
            box.angular_velocity = w0;
            Rigid3DSettings settings;
            settings.gravity = Vec3{0.0, 0.0, -g};
            settings.restitution = 0.0;
            settings.ground.push_back(malloy::collide::Plane3{Vec3{0.0, 0.0, 1.0}, 0.0});
            Rigid3DWorld w{SimulationSettings{dt}, {box}, settings};
            const auto energy = [&] {
                return malloy::rigid::total_kinetic_energy3d(w.bodies()) +
                       malloy::rigid::total_potential_energy3d(w.bodies(), settings.gravity) -
                       g * half.z;
            };
            const Real e0 = energy();
            const Real lever = std::sqrt(half.x * half.x + half.y * half.y);
            const Real bound = g * (g * dt + lever * malloy::math::length(w0)) * dt;
            const Real step_scale = g * g * dt * dt;
            Real worst = e0;
            Real late_kinetic = 0.0;
            for (int i = 1; i <= 10000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                worst = std::fmax(worst, energy());
                if (i > 7500)
                {
                    late_kinetic = std::fmax(late_kinetic,
                                             malloy::rigid::total_kinetic_energy3d(w.bodies()));
                }
            }
            MALLOY_CHECK_TRUE(worst - e0 <= bound);
            MALLOY_CHECK_TRUE(late_kinetic < 1e-4 * step_scale);
        }
    }

    // --- M40 (issue #73): the warm-start impulses are part of the world's
    //     state. A failed step restores them with the bodies, and they are
    //     what lets a stack come to rest within the sweep cap. ---
    {
        using malloy::math::Vec3;
        using malloy::rigid::box_inertia;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;

        const auto same2 = [](const RigidWorld& x, const RigidWorld& y) {
            bool same = x.bodies().size() == y.bodies().size() &&
                        x.warm_contacts().size() == y.warm_contacts().size();
            for (std::size_t i = 0; same && i < x.bodies().size(); ++i)
            {
                const RigidBody2D& p = x.bodies()[i];
                const RigidBody2D& q = y.bodies()[i];
                same = p.position.x == q.position.x && p.position.y == q.position.y &&
                       p.velocity.x == q.velocity.x && p.velocity.y == q.velocity.y &&
                       p.angle == q.angle && p.angular_velocity == q.angular_velocity;
            }
            for (std::size_t i = 0; same && i < x.warm_contacts().size(); ++i)
            {
                const RigidWorld::WarmContact& p = x.warm_contacts()[i];
                const RigidWorld::WarmContact& q = y.warm_contacts()[i];
                same = p.a == q.a && p.b == q.b && p.ground == q.ground && p.key == q.key &&
                       p.local.x == q.local.x && p.local.y == q.local.y &&
                       p.normal_impulse == q.normal_impulse &&
                       p.tangent_impulse == q.tangent_impulse;
            }
            return same;
        };

        // 2D: two identical worlds hold a box tipping from a corner onto its
        // face (its contact impulses change every step, so one step's cache
        // differs from the next) and, far away, a body with no collider at
        // x = 1.34e154, just inside the squarable range (issue #16). At step
        // 20 one world gives that body an impulse of 1e154, so the next step
        // carries it past 1.3408e154, where its square overflows: the step
        // fails and must restore the bodies AND the cache. The opposite
        // impulse then puts the body back exactly (1e154 - 1e154 = 0), and
        // from there the two worlds must agree bit for bit, cache included. A
        // world that kept the failed step's cache would warm-start from
        // impulses one step ahead and differ at the rounding level.
        {
            const Real tilt = 0.1;
            const auto make = [tilt] {
                RigidBody2D box;
                box.half_extents = Vec2{0.5, 0.25};
                box.inertia = (0.25 + 0.0625) / 3.0;
                box.angle = tilt;
                box.position = Vec2{0.0, 0.5 * std::sin(tilt) + 0.25 * std::cos(tilt)};
                RigidBody2D far;
                far.inertia = 1.0;
                far.position = Vec2{1.34e154, 0.0};
                RigidSettings settings;
                settings.gravity = Vec2{0.0, -9.8};
                settings.restitution = 0.0;
                settings.friction = 0.5;
                settings.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});
                return RigidWorld{SimulationSettings{0.002}, {box, far}, settings};
            };
            RigidWorld failing = make();
            RigidWorld twin = make();
            for (int i = 0; i < 20; ++i)
            {
                MALLOY_CHECK_TRUE(failing.step().ok());
                MALLOY_CHECK_TRUE(twin.step().ok());
            }
            MALLOY_CHECK_FALSE(failing.warm_contacts().empty()); // tipping on its corner
            const Vec2 far_centre = failing.bodies()[1].position;
            MALLOY_CHECK_TRUE(failing.apply_impulse_at(1, Vec2{1e154, 0.0}, far_centre));
            MALLOY_CHECK_TRUE(failing.step().status == StepStatus::InvalidState);
            MALLOY_CHECK_EQ(failing.tick_count(), twin.tick_count());
            MALLOY_CHECK_TRUE(failing.apply_impulse_at(1, Vec2{-1e154, 0.0}, far_centre));
            MALLOY_CHECK_EQ(failing.bodies()[1].velocity.x, 0.0);
            MALLOY_CHECK_TRUE(same2(failing, twin));
            for (int i = 0; i < 200; ++i)
            {
                MALLOY_CHECK_TRUE(failing.step().ok());
                MALLOY_CHECK_TRUE(twin.step().ok());
                MALLOY_CHECK_TRUE(same2(failing, twin));
            }
        }

        // 3D: the same, with a box tipping from an edge. There is no impulse
        // to aim, so the far body starts at x = 1.3e154 moving at 1e154: it
        // crosses 1.3408e154 at step 21 (20.4 steps of 2e151), which fails.
        // That step must leave the bodies and the cache exactly as step 20
        // left them, and so must every retry.
        {
            const Real tilt = 0.1;
            const Vec3 half{0.5, 0.4, 0.25};
            RigidBody3D box;
            box.half_extents = half;
            box.inertia = box_inertia(half, box.mass);
            box.orientation = malloy::math::from_axis_angle(Vec3{0.0, 1.0, 0.0}, tilt);
            box.position = Vec3{0.0, 0.0, 0.5 * std::sin(tilt) + 0.25 * std::cos(tilt)};
            RigidBody3D far;
            far.position = Vec3{1.3e154, 0.0, 0.0};
            far.velocity = Vec3{1e154, 0.0, 0.0};
            Rigid3DSettings settings;
            settings.gravity = Vec3{0.0, 0.0, -9.8};
            settings.restitution = 0.0;
            settings.friction = 0.5;
            settings.ground.push_back(malloy::collide::Plane3{Vec3{0.0, 0.0, 1.0}, 0.0});
            Rigid3DWorld w{SimulationSettings{0.002}, {box, far}, settings};
            int ok_steps = 0;
            while (ok_steps < 40 && w.step().ok())
            {
                ++ok_steps;
            }
            MALLOY_CHECK_EQ(ok_steps, 20);
            MALLOY_CHECK_FALSE(w.warm_contacts().empty());
            const std::vector<RigidBody3D> bodies = w.bodies();
            const std::vector<Rigid3DWorld::WarmContact> cache = w.warm_contacts();
            for (int retry = 0; retry < 2; ++retry)
            {
                if (retry > 0)
                {
                    MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
                }
                MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{20});
                bool same = w.bodies().size() == bodies.size() &&
                            w.warm_contacts().size() == cache.size();
                for (std::size_t i = 0; same && i < bodies.size(); ++i)
                {
                    const RigidBody3D& p = w.bodies()[i];
                    const RigidBody3D& q = bodies[i];
                    same = p.position.x == q.position.x && p.position.y == q.position.y &&
                           p.position.z == q.position.z && p.velocity.x == q.velocity.x &&
                           p.velocity.y == q.velocity.y && p.velocity.z == q.velocity.z &&
                           p.orientation.w == q.orientation.w &&
                           p.orientation.v.x == q.orientation.v.x &&
                           p.orientation.v.y == q.orientation.v.y &&
                           p.orientation.v.z == q.orientation.v.z &&
                           p.angular_velocity.x == q.angular_velocity.x &&
                           p.angular_velocity.y == q.angular_velocity.y &&
                           p.angular_velocity.z == q.angular_velocity.z;
                }
                for (std::size_t i = 0; same && i < cache.size(); ++i)
                {
                    const Rigid3DWorld::WarmContact& p = w.warm_contacts()[i];
                    const Rigid3DWorld::WarmContact& q = cache[i];
                    same = p.a == q.a && p.b == q.b && p.ground == q.ground && p.key == q.key &&
                           p.local.x == q.local.x && p.local.y == q.local.y &&
                           p.local.z == q.local.z && p.normal_impulse == q.normal_impulse &&
                           p.tangent_impulse.x == q.tangent_impulse.x &&
                           p.tangent_impulse.y == q.tangent_impulse.y &&
                           p.tangent_impulse.z == q.tangent_impulse.z;
                }
                MALLOY_CHECK_TRUE(same);
            }
        }

        // A stack of three equal boxes resting on the floor stays at rest.
        // At rest the velocity scale is one step of gravity, g dt, so every
        // loaded point ends within r = 1e-8 g dt of 0. A box on face
        // contacts w = 1 apart (2D: the width; 3D: the narrower side) with
        // corners within R of its centre of mass then keeps |omega| <=
        // 2 sqrt(3) r / w and |v| <= 3 r (1 + 2 R / w) relative to the box
        // below it, so the k-th box up keeps at most k times that, and the
        // stack's kinetic energy is at most sum_k (1/2) k^2 (m v^2 + I w^2)
        // with those bounds: 1.6e-17 J in 3D, less in 2D. Over the last
        // second of 10 s it stays inside. The warm start is what makes that
        // possible within 100 sweeps: started from zero every step, the three
        // stacked manifolds hit the cap and kept 2.5e-11 J (3D).
        {
            const Real g = 9.8;
            const Real dt = 0.002;
            const Real r = 1e-8 * g * dt;
            std::vector<RigidBody2D> flat;
            std::vector<RigidBody3D> deep;
            for (int k = 0; k < 3; ++k)
            {
                RigidBody2D b2;
                b2.half_extents = Vec2{0.5, 0.25};
                b2.inertia = (0.25 + 0.0625) / 3.0;
                b2.position = Vec2{0.0, 0.25 + 0.5 * k};
                flat.push_back(b2);
                RigidBody3D b3;
                b3.half_extents = Vec3{0.5, 0.5, 0.25};
                b3.inertia = box_inertia(b3.half_extents, b3.mass);
                b3.position = Vec3{0.0, 0.0, 0.25 + 0.5 * k};
                deep.push_back(b3);
            }
            RigidSettings s2;
            s2.gravity = Vec2{0.0, -g};
            s2.restitution = 0.0;
            s2.friction = 0.5;
            s2.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});
            Rigid3DSettings s3;
            s3.gravity = Vec3{0.0, 0.0, -g};
            s3.restitution = 0.0;
            s3.friction = 0.5;
            s3.ground.push_back(malloy::collide::Plane3{Vec3{0.0, 0.0, 1.0}, 0.0});
            RigidWorld w2{SimulationSettings{dt}, flat, s2};
            Rigid3DWorld w3{SimulationSettings{dt}, deep, s3};
            Real late2 = 0.0;
            Real late3 = 0.0;
            for (int i = 1; i <= 5000; ++i)
            {
                MALLOY_CHECK_TRUE(w2.step().ok());
                MALLOY_CHECK_TRUE(w3.step().ok());
                if (i > 4500)
                {
                    late2 = std::fmax(late2, total_kinetic_energy(w2.bodies()));
                    late3 = std::fmax(late3, malloy::rigid::total_kinetic_energy3d(w3.bodies()));
                }
            }
            const auto bound = [r](Real reach, Real width, Real mass, Real moment) {
                const Real v = 3.0 * r * (1.0 + 2.0 * reach / width);
                const Real spin = 2.0 * std::sqrt(3.0) * r / width;
                Real total = 0.0;
                for (int k = 1; k <= 3; ++k)
                {
                    total += 0.5 * k * k * (mass * v * v + moment * spin * spin);
                }
                return total;
            };
            const Real moment3 = std::fmax(std::fmax(deep[0].inertia.x, deep[0].inertia.y),
                                           deep[0].inertia.z);
            MALLOY_CHECK_TRUE(late2 <= bound(std::sqrt(0.25 + 0.0625), 1.0, 1.0, flat[0].inertia));
            MALLOY_CHECK_TRUE(late3 <= bound(0.75, 1.0, 1.0, moment3));
        }
    }

    // --- Issue #51: a huge mass or inertia must not validate while every energy
    //     the world reports is inf. Since #16 the position and velocity must be
    //     squarable, but nothing bounded the mass or inertia, so a body whose
    //     own kinetic energy overflowed was accepted and every step reported Ok
    //     with E printed as inf. In 2D the body's own kinetic energy is now
    //     required to be finite, which bounds its momenta too; in 3D, where |p|
    //     and |L| are reported, its momenta must be squarable, which bounds its
    //     energy too. Infinite mass and inertia stay legal in 2D. ---
    {
        using malloy::math::Vec3;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;

        // 2D. The reviewer's body: mass and inertia 1e300, speed and spin 1e5.
        {
            RigidBody2D huge;
            huge.mass = 1e300;
            huge.inertia = 1e300;
            huge.radius = 0.5;
            huge.velocity = Vec2{1e5, 0.0};
            huge.angular_velocity = 1e5;
            MALLOY_CHECK_FALSE(huge.is_valid());
            RigidWorld w{SimulationSettings{0.001}, {huge}, RigidSettings{}};
            MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidState);
            MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
            MALLOY_CHECK_EQ(w.bodies()[0].position.x, 0.0); // left unchanged

            // Each half on its own: the translation, then the rotation.
            RigidBody2D heavy;
            heavy.mass = 1e300;
            heavy.velocity = Vec2{0.0, 1e5};
            MALLOY_CHECK_FALSE(heavy.is_valid());
            RigidBody2D stiff;
            stiff.inertia = 1e300;
            stiff.angular_velocity = 1e5;
            MALLOY_CHECK_FALSE(stiff.is_valid());

            // And their SUM, which is what total_kinetic_energy adds: each half
            // is 1e308, finite, and together they are not.
            RigidBody2D both;
            both.mass = 1e300;
            both.inertia = 1e300;
            both.velocity = Vec2{std::sqrt(2e8), 0.0};
            both.angular_velocity = std::sqrt(2e8);
            MALLOY_CHECK_TRUE(std::isfinite(0.5 * both.mass * both.velocity.x * both.velocity.x));
            MALLOY_CHECK_FALSE(both.is_valid());

            // A huge mass is not itself the problem: at rest, or slow enough,
            // it is a perfectly good body.
            RigidBody2D resting = huge;
            resting.velocity = Vec2{};
            resting.angular_velocity = 0.0;
            MALLOY_CHECK_TRUE(resting.is_valid());
            RigidBody2D creeping = huge;
            creeping.velocity = Vec2{1e-100, 0.0};
            creeping.angular_velocity = -1e-100;
            MALLOY_CHECK_TRUE(creeping.is_valid());

            // Infinite mass and inertia stay legal: a static body, a kinematic
            // one (immovable but moving), and each half-infinite kind. inf * 0
            // would be NaN and inf * v inf, so the infinite half is skipped.
            RigidBody2D wall;
            wall.mass = inf;
            wall.inertia = inf;
            wall.radius = 1.0;
            MALLOY_CHECK_TRUE(wall.is_valid());
            RigidBody2D conveyor = wall;
            conveyor.velocity = Vec2{2.0, 0.0};
            conveyor.angular_velocity = 0.5;
            MALLOY_CHECK_TRUE(conveyor.is_valid());
            RigidBody2D pinned; // cannot translate, can spin
            pinned.mass = inf;
            pinned.angular_velocity = 3.0;
            MALLOY_CHECK_TRUE(pinned.is_valid());
            RigidBody2D slider; // cannot spin, can translate
            slider.inertia = inf;
            slider.velocity = Vec2{1.0, 2.0};
            MALLOY_CHECK_TRUE(slider.is_valid());

            // A body that starts valid and is pushed past the limit by a step
            // is refused and rolled back, like any other invalid result. Mass
            // 1e300 at 18900 has energy 1.79e308 minus a little; one step of a
            // strong gravity takes it to 19900, and 1.98e308 overflows.
            RigidBody2D edge;
            edge.mass = 1e300;
            edge.velocity = Vec2{0.0, -18900.0};
            MALLOY_CHECK_TRUE(edge.is_valid());
            RigidSettings strong;
            strong.gravity = Vec2{0.0, -1000.0};
            RigidWorld fall{SimulationSettings{1.0}, {edge}, strong};
            MALLOY_CHECK_TRUE(fall.step().status == StepStatus::InvalidState);
            MALLOY_CHECK_EQ(fall.bodies()[0].velocity.y, -18900.0);
        }

        // 3D. The reviewer's body, then each term on its own.
        {
            RigidBody3D huge;
            huge.mass = 1e300;
            huge.inertia = Vec3{1e300, 1e300, 1e300};
            huge.radius = 0.5;
            huge.velocity = Vec3{1e5, 0.0, 0.0};
            huge.angular_velocity = Vec3{1e5, 0.0, 0.0};
            MALLOY_CHECK_FALSE(huge.is_valid());
            Rigid3DWorld w{SimulationSettings{0.001}, {huge}};
            MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidState);
            MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);

            RigidBody3D heavy;
            heavy.mass = 1e300;
            heavy.velocity = Vec3{0.0, 1e5, 0.0};
            MALLOY_CHECK_FALSE(heavy.is_valid());

            // Each axis pairs with its OWN moment: a huge moment about an axis
            // the body does not spin about is harmless, spinning about it is
            // not.
            for (int axis = 0; axis < 3; ++axis)
            {
                RigidBody3D stiff;
                Vec3 inertia{1.0, 1.0, 1.0};
                Vec3 spin{};
                Vec3 other_spin{};
                (axis == 0 ? inertia.x : axis == 1 ? inertia.y : inertia.z) = 1e300;
                (axis == 0 ? spin.x : axis == 1 ? spin.y : spin.z) = 1e5;
                (axis == 0 ? other_spin.y : axis == 1 ? other_spin.z : other_spin.x) = 1e5;
                stiff.inertia = inertia;
                stiff.angular_velocity = spin;
                MALLOY_CHECK_FALSE(stiff.is_valid());
                stiff.angular_velocity = other_spin; // a moment of 1, at 1e5
                MALLOY_CHECK_TRUE(stiff.is_valid());
            }

            // The momenta must be SQUARABLE, not just finite, because the
            // world reports |p| and |L|. Mass 1e300 at 1e-140 has energy 5e19
            // and momentum 1e160, whose square is inf; the same for a moment.
            RigidBody3D drifting;
            drifting.mass = 1e300;
            drifting.velocity = Vec3{1e-140, 0.0, 0.0};
            MALLOY_CHECK_TRUE(std::isfinite(0.5 * drifting.mass * 1e-140 * 1e-140));
            MALLOY_CHECK_FALSE(drifting.is_valid());
            RigidBody3D turning;
            turning.inertia = Vec3{1.0, 1e300, 1.0};
            turning.angular_velocity = Vec3{0.0, 1e-140, 0.0};
            MALLOY_CHECK_FALSE(turning.is_valid());

            // At rest, or slow enough, a huge mass is a good body.
            RigidBody3D resting = huge;
            resting.velocity = Vec3{};
            resting.angular_velocity = Vec3{};
            MALLOY_CHECK_TRUE(resting.is_valid());
            RigidBody3D creeping = huge;
            creeping.velocity = Vec3{0.0, 0.0, 1e-160};
            creeping.angular_velocity = Vec3{1e-160, 0.0, 0.0};
            MALLOY_CHECK_TRUE(creeping.is_valid());
        }
    }

    // =========================================================================
    // M40 acceptance (issue #73): the per-point iterative contact solver.
    //
    // Criteria AC-1 to AC-17, in 2D and 3D. Every tolerance is DERIVED from
    // the documented constants of the contact step (the world headers' step
    // contract), never fitted to a run:
    //
    //   rest phase    sweeps until no contact constraint is violated by
    //                 more than 1e-8 v_scale (checked before every sweep), or
    //                 100 sweeps; v_scale is the larger of one step of the
    //                 external acceleration (|g| dt) and the fastest material
    //                 point in contact, |v| + |w| size, taken after the warm
    //                 start and the centroid prepass;
    //   impact phase  the same with 1e-12 v_scale;
    //   correction    NGS on positions only: linear velocities untouched, and
    //                 through a turn the body's WORLD angular momentum kept;
    //                 slop 1e-6 x size;
    //   restitution   Newton on the mid-step normal velocity, only for
    //                 approach faster than v_th = 4 dt (|g| + |F| / m_min);
    //   margin        corners within 1e-4 x size join the manifold.
    //
    // A residual r bounds what a resting body can keep. Every contact-point
    // velocity (normal, and tangential where friction holds) lies within r
    // of its target 0, and the face contacts are at least w apart with every
    // corner within R of the centre of mass, so |omega| <= 2 sqrt(3) r / w
    // and |v| <= 3 r (1 + 2 R / w). Those two numbers, and the kinetic energy
    // they allow, are the rest bounds below. At rest the warm start and the
    // prepass leave every speed near zero, so the step's v_scale is at most
    // |g| dt plus the largest |v| + |w| size the tests measure around it.
    //
    // Each criterion prints its value and bound, and the block fails at the
    // end if any failed, so one run reports every number.
    // =========================================================================
    {
        using malloy::math::cross;
        using malloy::math::dot;
        using malloy::math::from_axis_angle;
        using malloy::math::length;
        using malloy::math::Quat;
        using malloy::math::rotate;
        using malloy::math::Vec3;
        using malloy::collide::Plane3;
        using malloy::rigid::box_inertia;
        using malloy::rigid::Rigid3DSettings;
        using malloy::rigid::Rigid3DWorld;
        using malloy::rigid::RigidBody3D;

        const Real m40_rest_tol = 1e-8;    // rest phase, times v_scale
        const Real m40_impact_tol = 1e-12; // impact phase, times v_scale
        const Real m40_slop_rel = 1e-6;    // NGS slop, times the size
        const Real m40_cap = 100.0;        // sweeps per phase, at most
        const Real m40_beta = 0.2;         // NGS: fraction of the excess depth per iteration
        const Real m40_ngs_iters = 3.0;    // NGS iterations per step
        const Real m40_eps = std::numeric_limits<Real>::epsilon();

        int m40_failed = 0;
        const auto m40_check = [&m40_failed](const char* name, Real value, Real bound) {
            const bool ok = value <= bound; // a NaN fails
            std::cout << (ok ? "  M40 ok    " : "  M40 FAIL  ") << name << ": " << value
                      << " (bound " << bound << ")\n";
            if (!ok)
            {
                ++m40_failed;
            }
        };
        std::cout.precision(4);

        // ---- builders ----
        const auto m40_box2 = [](Real mass, Vec2 half, Vec2 pos, Real angle, Vec2 vel,
                                 Real spin) {
            RigidBody2D b;
            b.mass = mass;
            b.inertia = mass * (half.x * half.x + half.y * half.y) / 3.0;
            b.half_extents = half;
            b.position = pos;
            b.angle = angle;
            b.velocity = vel;
            b.angular_velocity = spin;
            return b;
        };
        const auto m40_disc = [](Real mass, Real radius, Vec2 pos, Vec2 vel) {
            RigidBody2D b;
            b.mass = mass;
            b.inertia = 0.5 * mass * radius * radius;
            b.radius = radius;
            b.position = pos;
            b.velocity = vel;
            return b;
        };
        const auto m40_floor2 = [](Real g, Real e, Real mu) {
            RigidSettings s;
            s.gravity = Vec2{0.0, -g};
            s.restitution = e;
            s.friction = mu;
            s.ground.push_back(Halfplane{Vec2{0.0, 1.0}, 0.0});
            return s;
        };
        const auto m40_box3 = [](Real mass, Vec3 half, Vec3 pos, Quat q, Vec3 vel, Vec3 spin) {
            RigidBody3D b;
            b.mass = mass;
            b.inertia = box_inertia(half, mass);
            b.half_extents = half;
            b.position = pos;
            b.orientation = q;
            b.velocity = vel;
            b.angular_velocity = spin;
            return b;
        };
        const auto m40_sphere = [](Real mass, Real radius, Vec3 pos, Vec3 vel) {
            RigidBody3D b;
            b.mass = mass;
            const Real moment = 0.4 * mass * radius * radius;
            b.inertia = Vec3{moment, moment, moment};
            b.radius = radius;
            b.position = pos;
            b.velocity = vel;
            return b;
        };
        const auto m40_floor3 = [](Real g, Real e, Real mu) {
            Rigid3DSettings s;
            s.gravity = Vec3{0.0, 0.0, -g};
            s.restitution = e;
            s.friction = mu;
            s.ground.push_back(Plane3{Vec3{0.0, 0.0, 1.0}, 0.0});
            return s;
        };

        // ---- measurements ----
        const auto m40_vmax2 = [](const std::vector<RigidBody2D>& bodies) {
            Real v = 0.0;
            for (const RigidBody2D& b : bodies)
            {
                v = std::fmax(v, length(b.velocity));
            }
            return v;
        };
        const auto m40_vmax3 = [](const std::vector<RigidBody3D>& bodies) {
            Real v = 0.0;
            for (const RigidBody3D& b : bodies)
            {
                v = std::fmax(v, length(b.velocity));
            }
            return v;
        };
        // The contract's velocity scale: the fastest material point,
        // |v| + |w| times the collider size (a box's largest half-width).
        const auto m40_scale2 = [](const std::vector<RigidBody2D>& bodies) {
            Real v = 0.0;
            for (const RigidBody2D& b : bodies)
            {
                const Real size = b.half_extents.x > 0.0 && b.half_extents.y > 0.0
                                      ? std::fmax(b.half_extents.x, b.half_extents.y)
                                      : b.radius;
                v = std::fmax(v, length(b.velocity) + std::fabs(b.angular_velocity) * size);
            }
            return v;
        };
        const auto m40_scale3 = [](const std::vector<RigidBody3D>& bodies) {
            Real v = 0.0;
            for (const RigidBody3D& b : bodies)
            {
                const Real size = b.half_extents.x > 0.0 && b.half_extents.y > 0.0 && b.half_extents.z > 0.0
                                      ? std::fmax(std::fmax(b.half_extents.x, b.half_extents.y), b.half_extents.z)
                                      : b.radius;
                v = std::fmax(v, length(b.velocity) + length(b.angular_velocity) * size);
            }
            return v;
        };
        const auto m40_ke3 = [](const std::vector<RigidBody3D>& bodies) {
            return malloy::rigid::total_kinetic_energy3d(bodies);
        };
        const auto m40_energy3 = [](const Rigid3DWorld& w) {
            return malloy::rigid::total_kinetic_energy3d(w.bodies()) +
                   malloy::rigid::total_potential_energy3d(w.bodies(), w.settings().gravity);
        };
        // The exact pre-solve state of a step: the same bodies stepped with no
        // collider and no ground, since forcing and integration read neither.
        const auto m40_presolve2 = [](const RigidWorld& w) {
            std::vector<RigidBody2D> ghost = w.bodies();
            for (RigidBody2D& b : ghost)
            {
                b.half_extents = Vec2{};
                b.radius = 0.0;
            }
            RigidSettings s = w.settings();
            s.ground.clear();
            RigidWorld g{w.simulation_settings(), ghost, s};
            g.step();
            return g.bodies();
        };
        // In 3D the forcing before the solve is the velocity's forcing and the
        // spin's torque KICK alone, at the start pose (issue #88): the free
        // rotation now comes AFTER the solve, so a ghost step, which rotates
        // the bodies, is no longer the pre-solve state of an asymmetric one.
        // Written out instead, in the step's own arithmetic.
        const auto m40_presolve3 = [](const Rigid3DWorld& w) {
            std::vector<RigidBody3D> pre = w.bodies();
            const Rigid3DSettings& s = w.settings();
            const Real dt = w.simulation_settings().dt;
            for (RigidBody3D& b : pre)
            {
                b.velocity += s.gravity * dt;
                b.velocity += (s.force / b.mass) * dt;
                const Vec3 tb = rotate(malloy::math::conjugate(b.orientation), s.torque);
                b.angular_velocity += Vec3{tb.x / b.inertia.x, tb.y / b.inertia.y, tb.z / b.inertia.z} * dt;
            }
            return pre;
        };
        // The exact SOLVED state of a 3D step, from the bodies before it and
        // the world after it: the start pose, the end velocities, and the spin
        // that carries each body's end WORLD angular momentum in its start
        // pose. After the solve the step only rotates each body freely,
        // translates it and projects it, and each of those keeps every
        // velocity and every body's world angular momentum (rigid3d_world.hpp,
        // steps 2 to 4 and D), so this is the state the impulses left.
        const auto m40_solved3 = [](const std::vector<RigidBody3D>& before, const Rigid3DWorld& w) {
            std::vector<RigidBody3D> solved = before;
            for (std::size_t k = 0; k < solved.size(); ++k)
            {
                const RigidBody3D& after = w.bodies()[k];
                const Vec3& m = after.inertia;
                const Vec3& s = after.angular_velocity;
                const Vec3 l = rotate(after.orientation, Vec3{m.x * s.x, m.y * s.y, m.z * s.z});
                const Vec3 lb = rotate(malloy::math::conjugate(before[k].orientation), l);
                solved[k].angular_velocity = Vec3{lb.x / m.x, lb.y / m.y, lb.z / m.z};
                solved[k].velocity = after.velocity;
            }
            return solved;
        };
        // Tilt of a 2D box from its nearest face-down pose; of a 3D box, the
        // angle between `up` and its most nearly vertical face normal.
        const auto m40_tilt2 = [&pi](const RigidBody2D& b) {
            Real a = std::fmod(b.angle, pi / 2.0);
            if (a < 0.0)
            {
                a += pi / 2.0;
            }
            return std::fmin(a, pi / 2.0 - a);
        };
        const auto m40_tilt3 = [](const RigidBody3D& b, Vec3 up) {
            Real best = 0.0;
            for (const Vec3 axis : {Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}})
            {
                best = std::fmax(best, std::fabs(dot(rotate(b.orientation, axis), up)));
            }
            return std::acos(std::fmin(1.0, best));
        };
        // The rest bounds of the header comment.
        const auto m40_rest_speed = [](Real r, Real w, Real R) { return 3.0 * r * (1.0 + 2.0 * R / w); };
        const auto m40_rest_spin = [](Real r, Real w) { return 2.0 * std::sqrt(3.0) * r / w; };
        const auto m40_rest_ke = [&](Real mass, Real moments, Real r, Real w, Real R) {
            const Real v = m40_rest_speed(r, w, R);
            const Real s = m40_rest_spin(r, w);
            return 0.5 * mass * v * v + 0.5 * moments * s * s;
        };
        const auto m40_moments3 = [](const RigidBody3D& b) {
            return b.inertia.x + b.inertia.y + b.inertia.z;
        };
        const auto m40_corner2 = [](const RigidBody2D& b, Real sx, Real sy) {
            const Vec2 u0{std::cos(b.angle), std::sin(b.angle)};
            const Vec2 u1{-std::sin(b.angle), std::cos(b.angle)};
            return b.position + u0 * (sx * b.half_extents.x) + u1 * (sy * b.half_extents.y);
        };
        const auto m40_corner3 = [](const RigidBody3D& b, int i) {
            const Vec3 local{(i & 1) ? b.half_extents.x : -b.half_extents.x,
                             (i & 2) ? b.half_extents.y : -b.half_extents.y,
                             (i & 4) ? b.half_extents.z : -b.half_extents.z};
            return b.position + rotate(b.orientation, local);
        };
        const auto m40_same2 = [](const RigidBody2D& a, const RigidBody2D& b) {
            return a.position.x == b.position.x && a.position.y == b.position.y &&
                   a.angle == b.angle && a.velocity.x == b.velocity.x &&
                   a.velocity.y == b.velocity.y && a.angular_velocity == b.angular_velocity;
        };
        const auto m40_same3 = [](const RigidBody3D& a, const RigidBody3D& b) {
            return a.position.x == b.position.x && a.position.y == b.position.y &&
                   a.position.z == b.position.z && a.velocity.x == b.velocity.x &&
                   a.velocity.y == b.velocity.y && a.velocity.z == b.velocity.z &&
                   a.orientation.w == b.orientation.w && a.orientation.v.x == b.orientation.v.x &&
                   a.orientation.v.y == b.orientation.v.y &&
                   a.orientation.v.z == b.orientation.v.z &&
                   a.angular_velocity.x == b.angular_velocity.x &&
                   a.angular_velocity.y == b.angular_velocity.y &&
                   a.angular_velocity.z == b.angular_velocity.z;
        };

        const auto m40_same_warm2 = [](const RigidWorld& x, const RigidWorld& y) {
            const auto& p = x.warm_contacts();
            const auto& q = y.warm_contacts();
            if (p.size() != q.size())
            {
                return false;
            }
            for (std::size_t k = 0; k < p.size(); ++k)
            {
                if (p[k].a != q[k].a || p[k].b != q[k].b || p[k].ground != q[k].ground ||
                    p[k].key != q[k].key || p[k].local.x != q[k].local.x ||
                    p[k].local.y != q[k].local.y || p[k].normal_impulse != q[k].normal_impulse ||
                    p[k].tangent_impulse != q[k].tangent_impulse)
                {
                    return false;
                }
            }
            return true;
        };
        const auto m40_same_warm3 = [](const Rigid3DWorld& x, const Rigid3DWorld& y) {
            const auto& p = x.warm_contacts();
            const auto& q = y.warm_contacts();
            if (p.size() != q.size())
            {
                return false;
            }
            for (std::size_t k = 0; k < p.size(); ++k)
            {
                if (p[k].a != q[k].a || p[k].b != q[k].b || p[k].ground != q[k].ground ||
                    p[k].key != q[k].key || p[k].local.x != q[k].local.x ||
                    p[k].local.y != q[k].local.y || p[k].local.z != q[k].local.z ||
                    p[k].normal_impulse != q[k].normal_impulse ||
                    p[k].tangent_impulse.x != q[k].tangent_impulse.x ||
                    p[k].tangent_impulse.y != q[k].tangent_impulse.y ||
                    p[k].tangent_impulse.z != q[k].tangent_impulse.z)
                {
                    return false;
                }
            }
            return true;
        };

        const Real m40_g = 9.8;
        const Real m40_dt = 0.002;
        const Real m40_gdt = m40_g * m40_dt;
        const Vec2 m40_brick{0.5, 0.25};
        const Vec3 m40_block{0.5, 0.4, 0.25};
        const Real m40_brick_I = (0.25 + 0.0625) / 3.0; // mass 1
        const Real m40_brick_R = length(m40_brick);
        const Real m40_block_R = length(m40_block);

        // ---------------------------------------------------------------------
        // AC-1. A box landing tilted comes to REST (issue #41): a 2D brick and
        // a 3D box dropped from 1 m with a 20 degree tilt, e 0, mu 0 and 0.5.
        // Over 2 s to 30 s the kinetic energy stays within the rest bound for
        // r = 1e-8 v_scale. T = 2 s: 0.35 s of fall plus a tip-down about the
        // landing corner of at most 0.27 s (corner moment (4/3) m (a^2 + b^2),
        // angular acceleration about 9.4 rad/s^2 over 0.35 rad), doubled.
        // w is the smaller face (2 x 0.25), so the bound holds on any face.
        // Frictionless, nothing stops a spin about the vertical, and before
        // issue #88 the forward-Euler rotation of an ANISOTROPIC body tilted
        // its angular momentum while it rocked (the old M20 drift law), so the
        // frictionless 3D case uses a CUBE: no gyroscopic term, L_z stays 0
        // (AC-13b), and it must come fully to rest.
        // ---------------------------------------------------------------------
        for (const Real mu : {0.0, 0.5})
        {
            RigidWorld w2{SimulationSettings{m40_dt},
                          {m40_box2(1.0, m40_brick, Vec2{0.0, 1.0}, 20.0 * pi / 180.0, Vec2{}, 0.0)},
                          m40_floor2(m40_g, 0.0, mu)};
            const Vec3 half3 = mu > 0 ? m40_block : Vec3{0.5, 0.5, 0.5};
            Rigid3DWorld w3{SimulationSettings{m40_dt},
                            {m40_box3(1.0, half3, Vec3{0.0, 0.0, 1.0},
                                      from_axis_angle(Vec3{1.0, 2.0, 0.0}, 20.0 * pi / 180.0), Vec3{},
                                      Vec3{})},
                            m40_floor3(m40_g, 0.0, mu)};
            Real late2 = 0.0, late3 = 0.0, v2 = 0.0, v3 = 0.0;
            int failed_steps = 0;
            for (int i = 1; i <= 15000; ++i)
            {
                failed_steps += w2.step().ok() ? 0 : 1;
                failed_steps += w3.step().ok() ? 0 : 1;
                if (i > 1000)
                {
                    late2 = std::fmax(late2, total_kinetic_energy(w2.bodies()));
                    late3 = std::fmax(late3, m40_ke3(w3.bodies()));
                    v2 = std::fmax(v2, m40_scale2(w2.bodies()));
                    v3 = std::fmax(v3, m40_scale3(w3.bodies()));
                }
            }
            const RigidBody3D& b3 = w3.bodies()[0];
            m40_check(mu > 0 ? "AC-1 2D tilted landing rests, mu 0.5 (late KE)"
                             : "AC-1 2D tilted landing rests, mu 0 (late KE)",
                      late2,
                      m40_rest_ke(1.0, m40_brick_I, m40_rest_tol * (m40_gdt + v2), 0.5, m40_brick_R));
            m40_check(mu > 0 ? "AC-1 3D tilted landing rests, mu 0.5 (late KE)"
                             : "AC-1 3D tilted landing rests, mu 0 (late KE)",
                      late3,
                      m40_rest_ke(1.0, m40_moments3(b3), m40_rest_tol * (m40_gdt + v3), 0.5,
                                  length(half3)));
            m40_check("AC-1 every step valid", failed_steps, 0.0);
        }

        // ---------------------------------------------------------------------
        // AC-2 (issue #41, e 0) and AC-3 (issue #83, e 1). A box resting on its
        // face, nudged: 2D brick, spin 0.05; 3D box (0.5, 0.5, 0.25), spin
        // (0.3, 0.2, 0). Over 20 s the mechanical energy never climbs more than
        // B = m g (g dt + lever |w0|) dt above its start (the most one step's
        // undetected approach can sink a corner, the existing #41 bound), and
        // over the last 5 s the box is at rest (the AC-1 bound). At e 1 the
        // threshold makes the step-scale contacts inelastic, so it rests too.
        // The 3D box is a symmetric top about its vertical axis (Ix = Iy) and
        // its spin starts horizontal, so the gyroscopic term never makes a
        // vertical spin: frictionless, it must still come to rest.
        // ---------------------------------------------------------------------
        for (const Real e : {0.0, 1.0})
        {
            for (const Real mu : {0.0, 0.5})
            {
                const Real w0 = 0.05;
                RigidWorld w2{SimulationSettings{m40_dt},
                              {m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{}, w0)},
                              m40_floor2(m40_g, e, mu)};
                const Vec3 half3{0.5, 0.5, 0.25};
                const Vec3 w03{0.3, 0.2, 0.0};
                Rigid3DWorld w3{SimulationSettings{m40_dt},
                                {m40_box3(1.0, half3, Vec3{0.0, 0.0, 0.25}, Quat{}, Vec3{}, w03)},
                                m40_floor3(m40_g, e, mu)};
                const Real e02 = total_energy(w2.bodies(), w2.settings().gravity);
                const Real e03 = m40_energy3(w3);
                const Real B2 = m40_g * (m40_gdt + m40_brick.x * w0) * m40_dt;
                const Real B3 =
                    m40_g * (m40_gdt + std::sqrt(0.5) * length(w03)) * m40_dt;
                Real top2 = e02, top3 = e03, late2 = 0.0, late3 = 0.0, v2 = 0.0, v3 = 0.0;
                for (int i = 1; i <= 10000; ++i)
                {
                    w2.step();
                    w3.step();
                    top2 = std::fmax(top2, total_energy(w2.bodies(), w2.settings().gravity));
                    top3 = std::fmax(top3, m40_energy3(w3));
                    if (i > 7500)
                    {
                        late2 = std::fmax(late2, total_kinetic_energy(w2.bodies()));
                        late3 = std::fmax(late3, m40_ke3(w3.bodies()));
                        v2 = std::fmax(v2, m40_scale2(w2.bodies()));
                        v3 = std::fmax(v3, m40_scale3(w3.bodies()));
                    }
                }
                const bool elastic = e > 0.5;
                const bool sticky = mu > 0.0;
                m40_check(elastic ? (sticky ? "AC-3 2D nudge e 1 mu .5: (Emax - E0) / B"
                                            : "AC-3 2D nudge e 1 mu 0: (Emax - E0) / B")
                                  : (sticky ? "AC-2 2D nudge e 0 mu .5: (Emax - E0) / B"
                                            : "AC-2 2D nudge e 0 mu 0: (Emax - E0) / B"),
                          (top2 - e02) / B2, 1.0);
                m40_check(elastic ? (sticky ? "AC-3 3D nudge e 1 mu .5: (Emax - E0) / B"
                                            : "AC-3 3D nudge e 1 mu 0: (Emax - E0) / B")
                                  : (sticky ? "AC-2 3D nudge e 0 mu .5: (Emax - E0) / B"
                                            : "AC-2 3D nudge e 0 mu 0: (Emax - E0) / B"),
                          (top3 - e03) / B3, 1.0);
                m40_check(elastic ? "AC-3 2D nudge rests (late KE)" : "AC-2 2D nudge rests (late KE)",
                          late2,
                          m40_rest_ke(1.0, m40_brick_I, m40_rest_tol * (m40_gdt + v2), 0.5,
                                      m40_brick_R));
                m40_check(elastic ? "AC-3 3D nudge rests (late KE)" : "AC-2 3D nudge rests (late KE)",
                          late3,
                          m40_rest_ke(1.0, m40_moments3(w3.bodies()[0]),
                                      m40_rest_tol * (m40_gdt + v3), 0.5, length(half3)));
            }
        }

        // ---------------------------------------------------------------------
        // AC-4. Every step at e 0 is dissipative: the kinetic energy after the
        // step is at most that of the exact pre-solve state plus
        // (1 + mu) r Lambda. Why: at the solve's end every loaded normal point
        // is within r of its target (0, or -s/dt <= 0 for a speculative one)
        // and friction is within r of sticking or opposes the slide, so
        // KE_after - KE_pre = sum(lambda . v_after) - (1/2) lambda^T K lambda
        // <= r sum(lambda_n + |lambda_t|) <= (1 + mu) r Lambda. Lambda bounds
        // the total normal impulse: each manifold's is at most what crosses it,
        // sum m (|dv| + |g| dt) over all bodies, and there are at most
        // n (n + 1) / 2 manifolds. The NGS correction moves no velocity. r is
        // 1e-8 times the larger of the pre-solve and solved states' scale (the
        // warm start and the prepass act between the two) plus |g| dt.
        // Setups: a tilted landing, a nudged box, a staggered 3-stack.
        // In 3D the comparison is between the exact pre-solve and SOLVED
        // states (m40_presolve3, m40_solved3): since issue #88 the free
        // rotation follows the solve and changes an asymmetric body's energy
        // by its own second-order error, and the projection's turn re-expresses
        // the world momentum in a turned frame, neither of which is the solve's.
        // ---------------------------------------------------------------------
        {
            Real worst = 0.0;
            for (const Real mu : {0.0, 0.5})
            {
                std::vector<std::vector<RigidBody2D>> setups2{
                    {m40_box2(1.0, m40_brick, Vec2{0.0, 0.6}, 20.0 * pi / 180.0, Vec2{}, 0.0)},
                    {m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{0.3, 0.0}, 0.4)},
                    {m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{}, 0.0),
                     m40_box2(1.3, Vec2{0.4, 0.2}, Vec2{0.1, 0.7}, 0.0, Vec2{}, 0.0),
                     m40_box2(0.7, Vec2{0.3, 0.25}, Vec2{0.15, 1.15}, 0.0, Vec2{}, 0.0)}};
                for (const std::vector<RigidBody2D>& init : setups2)
                {
                    RigidWorld w{SimulationSettings{m40_dt}, init, m40_floor2(m40_g, 0.0, mu)};
                    const Real n = static_cast<Real>(init.size());
                    for (int i = 0; i < 1500; ++i)
                    {
                        const std::vector<RigidBody2D> pre = m40_presolve2(w);
                        w.step();
                        Real flow = 0.0;
                        for (std::size_t k = 0; k < pre.size(); ++k)
                        {
                            flow += pre[k].mass *
                                    (length(w.bodies()[k].velocity - pre[k].velocity) + m40_gdt);
                        }
                        const Real r = m40_rest_tol *
                                       (m40_gdt + std::fmax(m40_scale2(pre), m40_scale2(w.bodies())));
                        const Real allowed = (1.0 + mu) * r * 0.5 * n * (n + 1.0) * flow +
                                             64.0 * m40_eps * total_kinetic_energy(pre);
                        const Real gain = total_kinetic_energy(w.bodies()) - total_kinetic_energy(pre);
                        worst = std::fmax(worst, gain / allowed);
                    }
                }
                std::vector<std::vector<RigidBody3D>> setups3{
                    {m40_box3(1.0, m40_block, Vec3{0.0, 0.0, 0.6},
                              from_axis_angle(Vec3{1.0, 2.0, 0.0}, 20.0 * pi / 180.0), Vec3{}, Vec3{})},
                    {m40_box3(1.0, m40_block, Vec3{0.0, 0.0, 0.25}, Quat{}, Vec3{0.2, -0.1, 0.0},
                              Vec3{0.3, 0.2, 0.1})},
                    {m40_box3(1.0, m40_block, Vec3{0.0, 0.0, 0.25}, Quat{}, Vec3{}, Vec3{}),
                     m40_box3(1.3, Vec3{0.4, 0.3, 0.2}, Vec3{0.1, 0.05, 0.7}, Quat{}, Vec3{}, Vec3{}),
                     m40_box3(0.7, Vec3{0.3, 0.3, 0.25}, Vec3{0.15, -0.05, 1.15}, Quat{}, Vec3{},
                              Vec3{})}};
                for (const std::vector<RigidBody3D>& init : setups3)
                {
                    Rigid3DWorld w{SimulationSettings{m40_dt}, init, m40_floor3(m40_g, 0.0, mu)};
                    const Real n = static_cast<Real>(init.size());
                    for (int i = 0; i < 1500; ++i)
                    {
                        const std::vector<RigidBody3D> before = w.bodies();
                        const std::vector<RigidBody3D> pre = m40_presolve3(w);
                        w.step();
                        const std::vector<RigidBody3D> solved = m40_solved3(before, w);
                        Real flow = 0.0;
                        for (std::size_t k = 0; k < pre.size(); ++k)
                        {
                            flow += pre[k].mass *
                                    (length(solved[k].velocity - pre[k].velocity) + m40_gdt);
                        }
                        const Real r = m40_rest_tol *
                                       (m40_gdt + std::fmax(m40_scale3(pre), m40_scale3(solved)));
                        const Real allowed = (1.0 + mu) * r * 0.5 * n * (n + 1.0) * flow +
                                             64.0 * m40_eps * m40_ke3(pre);
                        const Real gain = m40_ke3(solved) - m40_ke3(pre);
                        worst = std::fmax(worst, gain / allowed);
                    }
                }
            }
            m40_check("AC-4 per-step dissipation at e 0: max gain / allowed", worst, 1.0);
        }

        // ---------------------------------------------------------------------
        // AC-5. Complementarity certificate after ONE step (g 0, e 0, mu 0):
        // a box whose bottom corners are already in the floor at the start of
        // the step (the solver detects at the start-of-step pose). From the
        // momentum and angular momentum each step's contact impulse is
        // recovered; then every corner leaves at normal speed >= -r, a corner
        // carrying load leaves at |speed| <= r, a corner leaving faster than r
        // carries none, and no impulse pulls. r = 1e-8 V plus rounding, V a
        // bound on the scale the sweeps use: it is taken after the centroid
        // prepass (no warm start on a first step), whose one impulse is at
        // most m (|v| + |w| R) (it stops the centroid's approach, K >= 1/m),
        // so V <= (|v| + |w| s) + (|v| + |w| R)(1 + m R s / I), s the size and
        // R the corner distance. The impulse tolerance is the recovery's
        // rounding, 1e-12 m V.
        // 2D: two corners (two unknowns, solved exactly); spins 0.7 and 5 rad/s
        // (at 5 one corner separates while both are in). 3D: a CUBE's four
        // corners; the centre of pressure must lie in the hull of the corners
        // that can carry load. A cube, because its world-frame spin is
        // unchanged by the step's rotation (isotropic inertia), so the
        // impulse's moment is I times the change of the world-frame spin
        // whatever the position phase does to the orientation: the correction
        // keeps the world angular momentum, which for a cube is the world
        // spin times I (AC-13b checks the same on a frictionless floor).
        // ---------------------------------------------------------------------
        for (const Real spin : {0.7, 5.0})
        {
            const Real tilt = 0.05;
            RigidBody2D b0 = m40_box2(1.0, m40_brick, Vec2{0.3, 0.0}, tilt, Vec2{0.4, -1.0}, spin);
            b0.position.y = 0.25 * std::cos(tilt) + 0.5 * std::sin(tilt) - 0.1;
            RigidWorld w{SimulationSettings{m40_dt}, {b0}, m40_floor2(0.0, 0.0, 0.0)};
            const RigidBody2D pre = m40_presolve2(w)[0];
            w.step();
            const RigidBody2D& post = w.bodies()[0];
            const Real reach = length(pre.velocity) + std::fabs(pre.angular_velocity) * m40_brick_R;
            const Real vs = m40_scale2({pre}) +
                            reach * (1.0 + pre.mass * m40_brick_R * 0.5 / pre.inertia);
            const Real tol_v = m40_rest_tol * vs + 64.0 * m40_eps * vs;
            const Real tol_l = 1e-12 * pre.mass * vs;
            // The two bottom corners of the START pose, midway to the surface.
            Vec2 arm[2];
            for (int k = 0; k < 2; ++k)
            {
                const Vec2 c = m40_corner2(b0, k == 0 ? -1.0 : 1.0, -1.0);
                arm[k] = Vec2{c.x, c.y / 2.0} - b0.position;
            }
            const Real Py = post.mass * (post.velocity.y - pre.velocity.y);
            const Real dL = post.inertia * (post.angular_velocity - pre.angular_velocity);
            const Real l1 = (dL - arm[0].x * Py) / (arm[1].x - arm[0].x);
            const Real l0 = Py - l1;
            Real worst_v = 0.0, worst_l = 0.0, worst_c = 0.0;
            for (int k = 0; k < 2; ++k)
            {
                const Real lam = k == 0 ? l0 : l1;
                const Real vn = post.velocity.y + post.angular_velocity * arm[k].x;
                worst_v = std::fmax(worst_v, -vn);  // approaching after the step
                worst_l = std::fmax(worst_l, -lam); // pulling
                if (lam > tol_l)
                {
                    worst_c = std::fmax(worst_c, std::fabs(vn)); // loaded: at target
                }
            }
            const bool fast = spin > 1.0;
            m40_check(fast ? "AC-5 2D certificate spin 5: no corner approaches"
                           : "AC-5 2D certificate spin .7: no corner approaches",
                      worst_v, tol_v);
            m40_check(fast ? "AC-5 2D certificate spin 5: no pull" : "AC-5 2D certificate spin .7: no pull",
                      worst_l, tol_l);
            m40_check(fast ? "AC-5 2D certificate spin 5: loaded corners at target"
                           : "AC-5 2D certificate spin .7: loaded corners at target",
                      worst_c, tol_v);
            m40_check("AC-5 2D frictionless: tangential velocity untouched",
                      std::fabs(post.velocity.x - pre.velocity.x), 64.0 * m40_eps * vs);
        }
        {
            const Quat q = from_axis_angle(Vec3{1.0, 2.0, 0.0}, 0.03);
            const Vec3 cube{0.5, 0.5, 0.5};
            RigidBody3D b0 = m40_box3(1.0, cube, Vec3{0.2, -0.1, 0.0}, q, Vec3{0.3, -0.2, -1.0},
                                      Vec3{0.5, -0.4, 0.3});
            Real lowest = 1e300;
            for (int i = 0; i < 8; ++i)
            {
                lowest = std::fmin(lowest, m40_corner3(b0, i).z);
            }
            b0.position.z = -lowest - 0.08; // every bottom corner is in (they span 0.06)
            Rigid3DWorld w{SimulationSettings{m40_dt}, {b0}, m40_floor3(0.0, 0.0, 0.0)};
            const RigidBody3D pre = m40_presolve3(w)[0];
            w.step();
            const RigidBody3D& post = w.bodies()[0];
            const Real reach = length(pre.velocity) + length(pre.angular_velocity) * length(cube);
            const Real vs = m40_scale3({pre}) +
                            reach * (1.0 + pre.mass * length(cube) * 0.5 / pre.inertia.x);
            const Real tol_v = m40_rest_tol * vs + 64.0 * m40_eps * vs;
            const Real tol_l = 1e-12 * pre.mass * vs;
            // The impulse's moment about the start-of-step centre: the moment
            // of inertia times the change of the world-frame spin.
            const Vec3 spin_world = rotate(post.orientation, post.angular_velocity);
            const Vec3 tau = (spin_world - rotate(pre.orientation, pre.angular_velocity)) * pre.inertia.x;
            const Real Lambda = post.mass * (post.velocity.z - pre.velocity.z);
            std::vector<Vec2> carriers;
            Real worst_v = 0.0;
            int in = 0;
            for (int i = 0; i < 8; ++i)
            {
                const Vec3 c = m40_corner3(b0, i);
                if (c.z > 0.0)
                {
                    continue;
                }
                ++in;
                const Vec3 arm = Vec3{c.x, c.y, c.z / 2.0} - b0.position;
                const Real vn = post.velocity.z + cross(spin_world, arm).z;
                worst_v = std::fmax(worst_v, -vn);
                if (vn <= tol_v)
                {
                    carriers.push_back(Vec2{arm.x, arm.y});
                }
            }
            // Centre of pressure: sum(lambda r) x z = tau, so its x is -tau.y
            // and its y is tau.x, divided by the total.
            Real outside = 0.0;
            if (Lambda > tol_l)
            {
                const Vec2 cop{-tau.y / Lambda, tau.x / Lambda};
                // Distance outside the hull of the carriers (0 when inside).
                Real best = 1e300;
                const std::size_t m = carriers.size();
                for (std::size_t a = 0; a < m; ++a)
                {
                    best = std::fmin(best, length(cop - carriers[a]));
                    for (std::size_t b = a + 1; b < m; ++b)
                    {
                        const Vec2 ab = carriers[b] - carriers[a];
                        const Real t = std::fmax(0.0, std::fmin(1.0, dot(cop - carriers[a], ab) / dot(ab, ab)));
                        best = std::fmin(best, length(cop - (carriers[a] + ab * t)));
                        for (std::size_t c = b + 1; c < m; ++c)
                        {
                            const Vec2 ac = carriers[c] - carriers[a];
                            const Real den = ab.x * ac.y - ab.y * ac.x;
                            if (std::fabs(den) > 0.0)
                            {
                                const Vec2 ap = cop - carriers[a];
                                const Real u = (ap.x * ac.y - ap.y * ac.x) / den;
                                const Real v = (ab.x * ap.y - ab.y * ap.x) / den;
                                if (u >= 0.0 && v >= 0.0 && u + v <= 1.0)
                                {
                                    best = 0.0;
                                }
                            }
                        }
                    }
                }
                outside = best;
            }
            m40_check("AC-5 3D certificate: four corners in at the start", std::fabs(in - 4.0), 0.0);
            m40_check("AC-5 3D certificate: no corner approaches", worst_v, tol_v);
            m40_check("AC-5 3D certificate: no pull (total impulse)", -Lambda, tol_l);
            m40_check("AC-5 3D certificate: centre of pressure in the carriers' hull", outside,
                      1e-9 * length(cube));
        }

        // ---------------------------------------------------------------------
        // AC-6. Restitution.
        // (a) A flat landing at the start of a step (bottom face exactly on
        //     the floor, approaching): the mid-step law, v_mid_out = -e v_mid_in
        //     with v_mid = v + (half a step of gravity), so the step ends at
        //     v_out = -e (v_in - g dt / 2) - g dt / 2, to the impact tolerance
        //     plus rounding, and with no spin. g 0 and g 9.8; |v_mid| is far
        //     above v_th = 4 g dt.
        // (b) Two corners of a tilted brick land at DIFFERENT speeds (g 0): each
        //     loaded corner leaves at exactly -e times its own approach speed.
        // (c) Below the threshold a resting box keeps no micro-bounce: placed at
        //     rest at e 0.5 and e 1, over 0.5 s to 1 s its speed stays within
        //     the rest bound (before M40: e g dt / (1 + e), and g dt at e 1).
        // (d) and (e), issue #88: the prediction under the kick-then-rotate
        //     order, and a fast anisotropic tumbler's bounces (below).
        // ---------------------------------------------------------------------
        for (const Real g : {0.0, m40_g})
        {
            const Real e = 0.5;
            const Real vin = -1.0;
            const Real expected = -e * (vin - g * m40_dt / 2.0) - g * m40_dt / 2.0;
            const Real tol = 4.0 * m40_impact_tol * 1.0 + 64.0 * m40_eps;
            RigidWorld w2{SimulationSettings{m40_dt},
                          {m40_box2(1.0, Vec2{0.5, 0.5}, Vec2{0.0, 0.5}, 0.0, Vec2{0.0, vin}, 0.0)},
                          m40_floor2(g, e, 0.0)};
            Rigid3DWorld w3{SimulationSettings{m40_dt},
                            {m40_box3(1.0, Vec3{0.5, 0.5, 0.5}, Vec3{0.0, 0.0, 0.5}, Quat{},
                                      Vec3{0.0, 0.0, vin}, Vec3{})},
                            m40_floor3(g, e, 0.0)};
            w2.step();
            w3.step();
            m40_check(g > 0 ? "AC-6a 2D flat bounce, g 9.8: mid-step law" : "AC-6a 2D flat bounce, g 0: v_out = -e v_in",
                      std::fabs(w2.bodies()[0].velocity.y - expected) +
                          std::fabs(w2.bodies()[0].angular_velocity),
                      tol);
            m40_check(g > 0 ? "AC-6a 3D flat bounce, g 9.8: mid-step law" : "AC-6a 3D flat bounce, g 0: v_out = -e v_in",
                      std::fabs(w3.bodies()[0].velocity.z - expected) +
                          length(w3.bodies()[0].angular_velocity),
                      tol);
        }
        {
            const Real e = 0.5;
            const Real tilt = 0.05;
            RigidBody2D b0 = m40_box2(1.0, m40_brick, Vec2{0.3, 0.0}, tilt, Vec2{0.0, -1.0}, 0.3);
            b0.position.y = 0.25 * std::cos(tilt) + 0.5 * std::sin(tilt) - 0.098;
            RigidWorld w{SimulationSettings{m40_dt}, {b0}, m40_floor2(0.0, e, 0.0)};
            w.step();
            const RigidBody2D& post = w.bodies()[0];
            Real worst = 0.0;
            int approaching = 0;
            for (int k = 0; k < 2; ++k)
            {
                const Vec2 c = m40_corner2(b0, k == 0 ? -1.0 : 1.0, -1.0);
                const Vec2 arm = Vec2{c.x, c.y / 2.0} - b0.position;
                const Real vin = b0.velocity.y + b0.angular_velocity * arm.x;
                const Real vout = post.velocity.y + post.angular_velocity * arm.x;
                approaching += vin < 0.0 ? 1 : 0;
                worst = std::fmax(worst, std::fabs(vout + e * vin));
            }
            m40_check("AC-6b both corners approach", std::fabs(approaching - 2.0), 0.0);
            m40_check("AC-6b each corner leaves at -e times its own approach", worst,
                      4.0 * m40_impact_tol * 1.15 + 64.0 * m40_eps);
        }
        for (const Real e : {0.5, 1.0})
        {
            RigidWorld w2{SimulationSettings{m40_dt},
                          {m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{}, 0.0)},
                          m40_floor2(m40_g, e, 0.0)};
            Rigid3DWorld w3{SimulationSettings{m40_dt},
                            {m40_box3(1.0, m40_block, Vec3{0.0, 0.0, 0.25}, Quat{}, Vec3{}, Vec3{})},
                            m40_floor3(m40_g, e, 0.0)};
            Real v2 = 0.0, v3 = 0.0, s23 = 0.0;
            for (int i = 1; i <= 500; ++i)
            {
                w2.step();
                w3.step();
                if (i > 250)
                {
                    v2 = std::fmax(v2, m40_vmax2(w2.bodies()));
                    v3 = std::fmax(v3, m40_vmax3(w3.bodies()));
                    s23 = std::fmax(s23, std::fmax(m40_scale2(w2.bodies()), m40_scale3(w3.bodies())));
                }
            }
            const Real r = m40_rest_tol * (m40_gdt + s23);
            m40_check(e > 0.75 ? "AC-6c 2D at rest, e 1: no micro-bounce" : "AC-6c 2D at rest, e .5: no micro-bounce",
                      v2, m40_rest_speed(r, 0.5, m40_brick_R));
            m40_check(e > 0.75 ? "AC-6c 3D at rest, e 1: no micro-bounce" : "AC-6c 3D at rest, e .5: no micro-bounce",
                      v3, m40_rest_speed(r, 0.5, m40_block_R));
        }
        // (d) Issue #88, the prediction. The impact phase predicts delta, the
        //     change the step makes to a point's normal velocity between the
        //     impact and rest phases. Since the fix that change is the forcing
        //     and the torque KICK alone, body b's minus body a's: the
        //     gyroscopic change of the spin happens in the free rotation, after
        //     the solve. A fast anisotropic box under a torque and a force lands
        //     on ONE corner at e 0.5; from the exact pre-solve and solved states
        //     the corner must leave at -e (v_n + delta / 2) + delta / 2, with
        //     delta read off the pre-solve state, to the impact tolerance.
        //     Forward Euler's gyroscopic term, kept in delta, would move that
        //     target by half of it, measured here too: it must exceed a million
        //     times the tolerance (it is 9e-4 m/s), so keeping it fails.
        {
            const Real e = 0.5;
            const Vec3 half{0.6, 0.35, 0.15};
            Rigid3DSettings s = m40_floor3(m40_g, e, 0.0);
            s.torque = Vec3{0.4, -0.7, 0.3};
            s.force = Vec3{0.3, 0.0, -0.5};
            RigidBody3D b0 = m40_box3(1.0, half, Vec3{}, from_axis_angle(Vec3{1.0, 1.0, 0.0}, 0.3),
                                      Vec3{0.2, -0.1, 0.0}, Vec3{2.0, 6.0, 4.0});
            int low = 0;
            for (int i = 1; i < 8; ++i)
            {
                if (m40_corner3(b0, i).z < m40_corner3(b0, low).z)
                {
                    low = i;
                }
            }
            Real clearance = 1e300; // of the next corner above the lowest
            for (int i = 0; i < 8; ++i)
            {
                if (i != low)
                {
                    clearance = std::fmin(clearance, m40_corner3(b0, i).z - m40_corner3(b0, low).z);
                }
            }
            const Vec3 arm = m40_corner3(b0, low) - b0.position;
            b0.position.z = -arm.z; // the corner exactly touching
            b0.velocity.z = -2.0 - cross(rotate(b0.orientation, b0.angular_velocity), arm).z;
            const auto normal_speed = [&arm](const RigidBody3D& b) {
                return b.velocity.z + cross(rotate(b.orientation, b.angular_velocity), arm).z;
            };

            Rigid3DWorld w{SimulationSettings{m40_dt}, {b0}, s};
            const std::vector<RigidBody3D> before = w.bodies();
            const RigidBody3D pre = m40_presolve3(w)[0];
            MALLOY_CHECK_TRUE(w.step().ok());
            const RigidBody3D solved = m40_solved3(before, w)[0];

            const Real vn = normal_speed(b0);
            const Real delta = normal_speed(pre) - vn;
            const Real middle = vn + delta / 2.0;
            const Real threshold = 4.0 * m40_dt * (m40_g + length(s.force) / b0.mass);
            const Real vs = std::fmax(m40_scale3({b0}), std::fmax(m40_scale3({pre}), m40_scale3({solved})));
            m40_check("AC-6d one corner in the manifold (clearance of the next, m)",
                      1e-4 * 0.6 / clearance, 1.0);
            m40_check("AC-6d the corner approaches past the threshold (v_mid + v_th)",
                      middle + threshold, 0.0);
            const Real tol = 4.0 * m40_impact_tol * vs + 64.0 * m40_eps * vs;
            m40_check("AC-6d corner leaves at -e (v_n + delta / 2) + delta / 2 (KD delta)",
                      std::fabs(normal_speed(solved) - (-e * middle + delta / 2.0)), tol);
            const Vec3& m = b0.inertia;
            const Vec3& sp = b0.angular_velocity;
            const Vec3 u = cross(Vec3{m.x * sp.x, m.y * sp.y, m.z * sp.z}, sp);
            const Vec3 gyro = rotate(b0.orientation, Vec3{u.x / m.x, u.y / m.y, u.z / m.z} * m40_dt);
            m40_check("AC-6d forward Euler's gyroscopic delta would move the target (1e6 tol / half of it)",
                      1e6 * tol / (std::fabs(cross(gyro, arm).z) / 2.0), 1.0);
        }
        // (e) Issue #88, the regression the scope asks for: a fast anisotropic
        //     tumbler bounces. The free rotation after the solve moves an
        //     asymmetric body's spin, and so a corner's normal velocity, by up
        //     to dt (I_max / I_min) |w|^2 |r| per step, unpredicted; here that
        //     is 132 dt, above the threshold's 8 |g| dt. The box (0.6, 0.35,
        //     0.15), spin (2, 6, 4), tilted 0.3 rad about (1, 1, 0), dropped from
        //     2 m at e 0.5, mu 0 and 0.5, dt 0.002, 0.001 and 0.0005, 4 s:
        //       - every impact step whose manifold is all impact points obeys
        //         the mid-step law at each of them, from the exact solved
        //         state, to the impact tolerance (a step that mixes in a point
        //         below the threshold is counted, not checked: the rest phase
        //         then acts on the same body after the impact). With friction,
        //         only a ONE-point impact is held to it: a frictional impact
        //         at several points can stop at the 100-sweep cap
        //         (kMaxSweeps, ADR 0011), short of the tolerance, as it
        //         already did on the forward-Euler step in this scenario
        //         (measured with an instrumented build: the impact phase
        //         stopped at the cap once there and 3 times here, the worst
        //         law error 2.2e-6 m/s), so those are counted and their worst
        //         printed;
        //       - the energy never rises above the start by more than this
        //         body's own rotational envelope, the largest change of its
        //         rotational energy rotating freely from the same start for the
        //         same time (translation in free flight only sheds energy, M22,
        //         and every impact at e < 1 dissipates), and no impact step
        //         raises it by more than the largest single-step change of that
        //         free rotation;
        //       - the first apex converges as dt halves. Its order is printed,
        //         not assumed: impacts are not smooth.
        {
            const Vec3 half{0.6, 0.35, 0.15};
            const Real g = 9.81;
            const Real e = 0.5;
            const Real margin = 1e-4 * 0.6;
            Real worst_law = 0.0;
            Real worst_rise = 0.0;
            Real worst_jump = 0.0;
            int impacts = 0;
            int mixed = 0;
            int checked = 0;
            int frictional = 0;
            Real worst_capped = 0.0;
            Real apex[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
            for (int f = 0; f < 2; ++f)
            {
                const Real mu = f == 0 ? 0.0 : 0.5;
                for (int k = 0; k < 3; ++k)
                {
                    const Real dt = 0.002 / static_cast<Real>(1 << k);
                    const int count = static_cast<int>(std::lround(4.0 / dt));
                    Rigid3DSettings s = m40_floor3(g, e, mu);
                    const RigidBody3D b0 = m40_box3(1.0, half, Vec3{0.0, 0.0, 2.0},
                                                    from_axis_angle(Vec3{1.0, 1.0, 0.0}, 0.3), Vec3{},
                                                    Vec3{2.0, 6.0, 4.0});
                    if (f == 0 && k == 0)
                    {
                        const Real ratio = std::fmax(std::fmax(b0.inertia.x, b0.inertia.y), b0.inertia.z) /
                                           std::fmin(std::fmin(b0.inertia.x, b0.inertia.y), b0.inertia.z);
                        m40_check("AC-6e the tumble is fast: 8 |g| / ((I_max / I_min) |w|^2 |r|)",
                                  8.0 * g / (ratio * dot(b0.angular_velocity, b0.angular_velocity) *
                                             length(half)),
                                  1.0);
                    }

                    // This body's own rotational envelope, rotating freely.
                    RigidBody3D lone = b0;
                    lone.half_extents = Vec3{};
                    Rigid3DWorld spinning{SimulationSettings{dt}, {lone}};
                    const Real t0 = malloy::rigid::rotational_energy(lone);
                    Real envelope = 0.0;
                    Real step_envelope = 0.0;
                    Real last = t0;
                    for (int i = 0; i < count; ++i)
                    {
                        MALLOY_CHECK_TRUE(spinning.step().ok());
                        const Real t = malloy::rigid::rotational_energy(spinning.bodies()[0]);
                        envelope = std::fmax(envelope, std::fabs(t - t0));
                        step_envelope = std::fmax(step_envelope, std::fabs(t - last));
                        last = t;
                    }

                    Rigid3DWorld w{SimulationSettings{dt}, {b0}, s};
                    const Real e0 = m40_energy3(w);
                    Real z1 = b0.position.z;
                    Real z2 = z1;
                    bool first = true;
                    for (int i = 1; i <= count; ++i)
                    {
                        const std::vector<RigidBody3D> before = w.bodies();
                        const RigidBody3D pre = m40_presolve3(w)[0];
                        const RigidBody3D& b = before[0];
                        const Real threshold = 4.0 * dt * g;
                        int in = 0;
                        int active = 0;
                        Vec3 arms[8];
                        Real middles[8];
                        Real deltas[8];
                        for (int c = 0; c < 8; ++c)
                        {
                            const Vec3 corner = m40_corner3(b, c);
                            if (corner.z > margin)
                            {
                                continue;
                            }
                            const Vec3 arm = corner - b.position;
                            const Real vn = b.velocity.z + cross(rotate(b.orientation, b.angular_velocity), arm).z;
                            const Real delta =
                                pre.velocity.z + cross(rotate(pre.orientation, pre.angular_velocity), arm).z - vn;
                            const Real middle = vn + delta / 2.0;
                            const Real gap = corner.z;
                            ++in;
                            if (middle < -threshold && (!(gap > 0.0) || -middle * dt >= gap))
                            {
                                arms[active] = arm;
                                middles[active] = middle;
                                deltas[active] = delta;
                                ++active;
                            }
                        }
                        const Real energy_before = m40_energy3(w);
                        MALLOY_CHECK_TRUE(w.step().ok());
                        if (active > 0)
                        {
                            ++impacts;
                            worst_jump = std::fmax(worst_jump, (m40_energy3(w) - energy_before) / step_envelope);
                            if (active < in)
                            {
                                ++mixed;
                            }
                            else
                            {
                                const RigidBody3D solved = m40_solved3(before, w)[0];
                                const Real vs = std::fmax(m40_scale3({b}), m40_scale3({solved}));
                                const bool held = mu == 0.0 || active == 1;
                                for (int a = 0; a < active; ++a)
                                {
                                    const Real out =
                                        solved.velocity.z +
                                        cross(rotate(solved.orientation, solved.angular_velocity), arms[a]).z;
                                    const Real error = std::fabs(out - (-e * middles[a] + deltas[a] / 2.0));
                                    if (held)
                                    {
                                        worst_law = std::fmax(
                                            worst_law, error / (4.0 * m40_impact_tol * vs + 64.0 * m40_eps * vs));
                                    }
                                    else
                                    {
                                        worst_capped = std::fmax(worst_capped, error);
                                    }
                                }
                                if (held)
                                {
                                    ++checked;
                                }
                                else
                                {
                                    ++frictional;
                                }
                            }
                        }
                        worst_rise = std::fmax(worst_rise, (m40_energy3(w) - e0) / envelope);
                        const Real z = w.bodies()[0].position.z;
                        if (first && i >= 2 && z1 > z2 && z1 >= z)
                        {
                            apex[f][k] = z1;
                            first = false;
                        }
                        z2 = z1;
                        z1 = z;
                    }
                }
            }
            m40_check("AC-6e every held impact: mid-step law, |error| / tolerance", worst_law, 1.0);
            m40_check("AC-6e impacts held to the law (1 / count)", 1.0 / static_cast<Real>(checked), 1.0);
            std::cout << "  M40 info  AC-6e impact steps " << impacts << ": " << checked << " held, " << mixed
                      << " mixed with a resting point, " << frictional
                      << " frictional at several points (worst " << worst_capped << " m/s)\n";
            m40_check("AC-6e energy rise above the start / the rotational envelope", worst_rise, 1.0);
            m40_check("AC-6e energy rise across an impact step / the single-step envelope", worst_jump, 1.0);
            for (int f = 0; f < 2; ++f)
            {
                const Real coarse = std::fabs(apex[f][0] - apex[f][1]);
                const Real fine = std::fabs(apex[f][1] - apex[f][2]);
                m40_check(f == 0 ? "AC-6e first apex converges, mu 0 (fine / coarse difference)"
                                 : "AC-6e first apex converges, mu .5 (fine / coarse difference)",
                          fine / coarse, 1.0);
                std::cout << "  M40 info  AC-6e first apex, mu " << (f == 0 ? "0" : ".5") << ": "
                          << apex[f][0] << ", " << apex[f][1] << ", " << apex[f][2] << " m; observed order "
                          << std::log2(coarse / fine) << "\n";
            }
        }

        // ---------------------------------------------------------------------
        // AC-7. Symmetric drops stay symmetric (box_drop setups: from 2 m in
        // g 4, e 0.5, dt 0.002, 2000 steps). The level manifold's prepass is
        // symmetric, so spin can come only from the solve stopping inside its
        // tolerance band: |omega| <= 2 sqrt(3) r / w with r = 1e-8 v_scale_max
        // (v_scale_max = the impact speed plus g dt). Frictionless, the normals
        // are exactly vertical, so the horizontal velocity stays exactly 0;
        // with friction it stays within the rest speed bound. And the box ends
        // AT REST (no micro-bounce) at its resting HEIGHT, within the slop,
        // 1e-6 x size, plus what the rest speed can sink while the projection
        // works, 4 v_rest dt (it removes 1 - 0.8^3 of the excess per step). The
        // first landing is found a step late, up to u dt = 7e-3 deep, but at
        // e 0.5 the bounce lifts the box out within that same step, so these
        // drops pass with the projection switched off; AC-17 tests it directly.
        // ---------------------------------------------------------------------
        {
            const Real g = 4.0;
            for (const Real mu : {0.0, 0.5})
            {
                const std::vector<Vec2> halves2{Vec2{0.5, 0.5}, Vec2{0.5, 0.25}};
                for (const Vec2 h : halves2)
                {
                    RigidWorld w{SimulationSettings{m40_dt},
                                 {m40_box2(1.0, h, Vec2{0.0, 2.0}, 0.0, Vec2{}, 0.0)},
                                 m40_floor2(g, 0.5, mu)};
                    Real spin = 0.0, side = 0.0, vmax = 0.0;
                    for (int i = 0; i < 2000; ++i)
                    {
                        w.step();
                        const RigidBody2D& b = w.bodies()[0];
                        spin = std::fmax(spin, std::fabs(b.angular_velocity));
                        side = std::fmax(side, std::fabs(b.velocity.x));
                        vmax = std::fmax(vmax, m40_scale2(w.bodies()));
                    }
                    const Real r = m40_rest_tol * (vmax + g * m40_dt);
                    const Real wmin = 2.0 * std::fmin(h.x, h.y);
                    m40_check(mu > 0 ? "AC-7 2D symmetric drop mu .5: max |w|" : "AC-7 2D symmetric drop mu 0: max |w|",
                              spin, m40_rest_spin(r, wmin));
                    m40_check(mu > 0 ? "AC-7 2D symmetric drop mu .5: max |vx|" : "AC-7 2D symmetric drop mu 0: max |vx|",
                              side, mu > 0 ? m40_rest_speed(r, wmin, length(h)) : 0.0);
                    const Real rr = m40_rest_tol * (g * m40_dt + m40_scale2(w.bodies()));
                    m40_check("AC-7 2D symmetric drop ends at rest: |v|", m40_vmax2(w.bodies()),
                              m40_rest_speed(rr, wmin, length(h)));
                    m40_check("AC-7 2D symmetric drop rests at its height: |y - h|",
                              std::fabs(w.bodies()[0].position.y - h.y),
                              m40_slop_rel * std::fmax(h.x, h.y) +
                                  4.0 * m40_rest_speed(rr, wmin, length(h)) * m40_dt);
                }
                const std::vector<Vec3> halves3{Vec3{0.5, 0.5, 0.5}, Vec3{0.7, 0.3, 0.2}};
                for (const Vec3 h : halves3)
                {
                    Rigid3DWorld w{SimulationSettings{m40_dt},
                                   {m40_box3(1.0, h, Vec3{0.0, 0.0, 2.0}, Quat{}, Vec3{}, Vec3{})},
                                   m40_floor3(g, 0.5, mu)};
                    Real spin = 0.0, side = 0.0, vmax = 0.0;
                    for (int i = 0; i < 2000; ++i)
                    {
                        w.step();
                        const RigidBody3D& b = w.bodies()[0];
                        spin = std::fmax(spin, length(b.angular_velocity));
                        side = std::fmax(side, std::fabs(b.velocity.x) + std::fabs(b.velocity.y));
                        vmax = std::fmax(vmax, m40_scale3(w.bodies()));
                    }
                    const Real r = m40_rest_tol * (vmax + g * m40_dt);
                    const Real wmin = 2.0 * std::fmin(h.x, h.y);
                    m40_check(mu > 0 ? "AC-7 3D symmetric drop mu .5: max |w|" : "AC-7 3D symmetric drop mu 0: max |w|",
                              spin, m40_rest_spin(r, wmin));
                    m40_check(mu > 0 ? "AC-7 3D symmetric drop mu .5: max |vx| + |vy|" : "AC-7 3D symmetric drop mu 0: max |vx| + |vy|",
                              side, mu > 0 ? 2.0 * m40_rest_speed(r, wmin, length(h)) : 0.0);
                    const Real rr = m40_rest_tol * (g * m40_dt + m40_scale3(w.bodies()));
                    m40_check("AC-7 3D symmetric drop ends at rest: |v|", m40_vmax3(w.bodies()),
                              m40_rest_speed(rr, wmin, length(h)));
                    m40_check("AC-7 3D symmetric drop rests at its height: |z - h|",
                              std::fabs(w.bodies()[0].position.z - h.z),
                              m40_slop_rel * std::fmax(std::fmax(h.x, h.y), h.z) +
                                  4.0 * m40_rest_speed(rr, wmin, length(h)) * m40_dt);
                }
            }
        }

        // ---------------------------------------------------------------------
        // AC-8. A rounding-level tilt lands like a flat box: tilt 1e-15 and
        // 1e-12 rad, landing at 1.872 m/s (no gravity, e 0.5), at 32 phases
        // of the landing within one step. The corners are within the margin of
        // each other, so they join one level manifold: the spin left is at
        // most 2 sqrt(3) (1e-8 + 1e-12) v_scale / w plus the tilt's own share,
        // 10 delta v_scale / h. A lone-corner impulse (before M40) gave about 2.4.
        // ---------------------------------------------------------------------
        {
            const Real u = 1.872;
            Real worst2 = 0.0, worst3 = 0.0;
            for (const Real delta : {1e-15, 1e-12})
            {
                for (int p = 0; p < 32; ++p)
                {
                    const Real lift = u * m40_dt * (0.5 + p / 32.0);
                    RigidWorld w2{SimulationSettings{m40_dt},
                                  {m40_box2(1.0, m40_brick,
                                            Vec2{0.0, 0.25 * std::cos(delta) + 0.5 * std::sin(delta) + lift},
                                            delta, Vec2{0.0, -u}, 0.0)},
                                  m40_floor2(0.0, 0.5, 0.0)};
                    RigidBody3D b3 = m40_box3(1.0, m40_block, Vec3{}, from_axis_angle(Vec3{1.0, 2.0, 0.0}, delta),
                                              Vec3{0.0, 0.0, -u}, Vec3{});
                    Real lowest = 1e300;
                    for (int i = 0; i < 8; ++i)
                    {
                        lowest = std::fmin(lowest, m40_corner3(b3, i).z);
                    }
                    b3.position.z = -lowest + lift;
                    Rigid3DWorld w3{SimulationSettings{m40_dt}, {b3}, m40_floor3(0.0, 0.5, 0.0)};
                    for (int i = 0; i < 4; ++i)
                    {
                        w2.step();
                        w3.step();
                    }
                    worst2 = std::fmax(worst2, std::fabs(w2.bodies()[0].angular_velocity) /
                                                   (m40_rest_spin((m40_rest_tol + m40_impact_tol) * u, 0.5) +
                                                    10.0 * delta * u / 0.25));
                    worst3 = std::fmax(worst3, length(w3.bodies()[0].angular_velocity) /
                                                   (m40_rest_spin((m40_rest_tol + m40_impact_tol) * u, 0.5) +
                                                    10.0 * delta * u / 0.25));
                }
            }
            m40_check("AC-8 2D near-flat landing: spin / bound", worst2, 1.0);
            m40_check("AC-8 3D near-flat landing: spin / bound", worst3, 1.0);
        }

        // ---------------------------------------------------------------------
        // AC-9. Stacks stand (e 0, mu 0.5, 10 s, placed exactly touching):
        // staggered 3- and 5-stacks, an aligned 5-stack, a 3D stack yawed 30
        // degrees per level (eight-point manifolds) and a sphere resting off
        // centre on a box.
        //   top-box horizontal drift <= 1 mm: a REQUIREMENT (1 per cent of the
        //     0.1 stagger), not a derivation;
        //   every box's tilt <= 1e-3 rad: a requirement as well;
        //   sag of the top box <= n times one interface's allowance. Each
        //     interface sits within the NGS slop, plus what the solve leaves
        //     sinking. A converged solve leaves the rest speed, sinking at
        //     most 4 v_rest dt over the steps NGS takes. A solve stopped by
        //     the cap removes at least 1 - rho^cap of the approach it starts
        //     with (rho, the Gauss-Seidel rate of the stack's normal rows),
        //     so it leaves at most rho^cap / (1 - rho^cap) g dt approaching,
        //     sinking that times dt per step, which NGS removes at the rate
        //     1 - (1 - beta)^3 per step. rho is the 1D chain's: cos^2(pi/2n)
        //     for n equal boxes, M / (M + 1) for M on a unit box (0.99 at
        //     M 100, whose cold first steps do hit the cap);
        //   kinetic energy over the last 2 s <= n^2 times the sum of the rest
        //     bounds: box k moves at most k times one interface's allowance.
        // ---------------------------------------------------------------------
        {
            // One interface's sag allowance, as derived above (size 0.5).
            const auto m40_interface = [&](int n, Real mass_top, Real v_rest) {
                const Real rho = n == 2 ? mass_top / (mass_top + 1.0)
                                        : std::pow(std::cos(pi / (2.0 * n)), 2.0);
                const Real left = std::pow(rho, m40_cap);
                const Real removal = 1.0 - std::pow(1.0 - m40_beta, m40_ngs_iters);
                return m40_slop_rel * 0.5 + 4.0 * v_rest * m40_dt +
                       left / ((1.0 - left) * removal) * m40_gdt * m40_dt;
            };
            const auto stack2 = [&](const char* name, int n, Real stagger, Real mass_top) {
                std::vector<RigidBody2D> bodies;
                for (int k = 0; k < n; ++k)
                {
                    bodies.push_back(m40_box2(k == n - 1 ? mass_top : 1.0, m40_brick,
                                              Vec2{stagger * k, 0.25 + 0.5 * k}, 0.0, Vec2{}, 0.0));
                }
                const Vec2 top0 = bodies.back().position;
                RigidWorld w{SimulationSettings{m40_dt}, bodies, m40_floor2(m40_g, 0.0, 0.5)};
                Real drift = 0.0, tilt = 0.0, sag = 0.0, late = 0.0, vlate = 0.0;
                int failed_steps = 0;
                for (int i = 1; i <= 5000; ++i)
                {
                    failed_steps += w.step().ok() ? 0 : 1;
                    const RigidBody2D& top = w.bodies().back();
                    drift = std::fmax(drift, std::fabs(top.position.x - top0.x));
                    sag = std::fmax(sag, top0.y - top.position.y);
                    for (const RigidBody2D& b : w.bodies())
                    {
                        tilt = std::fmax(tilt, m40_tilt2(b));
                    }
                    if (i > 4000)
                    {
                        late = std::fmax(late, total_kinetic_energy(w.bodies()));
                        vlate = std::fmax(vlate, m40_scale2(w.bodies()));
                    }
                }
                const Real r = m40_rest_tol * (m40_gdt + vlate);
                const Real v_rest = m40_rest_speed(r, 0.5, m40_brick_R);
                Real ke = 0.0;
                for (const RigidBody2D& b : w.bodies())
                {
                    ke += m40_rest_ke(b.mass, b.inertia, r, 0.5, m40_brick_R);
                }
                const Real nn = static_cast<Real>(n);
                std::cout << "  M40 ---- " << name << '\n';
                m40_check("AC-9/10 2D stack: top drift (m)", drift, 1e-3);
                m40_check("AC-9/10 2D stack: max tilt (rad)", tilt, 1e-3);
                m40_check("AC-9/10 2D stack: top sag (m)", sag, nn * m40_interface(n, mass_top, v_rest));
                m40_check("AC-9/10 2D stack: late KE", late, nn * nn * ke);
                m40_check("AC-9/10 2D stack: every step valid", failed_steps, 0.0);
            };
            stack2("2D staggered 3-stack", 3, 0.1, 1.0);
            stack2("2D staggered 5-stack", 5, 0.1, 1.0);
            stack2("2D aligned 5-stack", 5, 0.0, 1.0);
            stack2("AC-10 2D 100 kg on 1 kg", 2, 0.1, 100.0);

            const auto stack3 = [&](const char* name, int n, Vec3 stagger, Real yaw, Real mass_top) {
                std::vector<RigidBody3D> bodies;
                for (int k = 0; k < n; ++k)
                {
                    bodies.push_back(m40_box3(k == n - 1 ? mass_top : 1.0, m40_block,
                                              Vec3{stagger.x * k, stagger.y * k, 0.25 + 0.5 * k},
                                              from_axis_angle(Vec3{0.0, 0.0, 1.0}, yaw * k), Vec3{}, Vec3{}));
                }
                const Vec3 top0 = bodies.back().position;
                Rigid3DWorld w{SimulationSettings{m40_dt}, bodies, m40_floor3(m40_g, 0.0, 0.5)};
                Real drift = 0.0, tilt = 0.0, sag = 0.0, late = 0.0, vlate = 0.0;
                int failed_steps = 0;
                for (int i = 1; i <= 5000; ++i)
                {
                    failed_steps += w.step().ok() ? 0 : 1;
                    const RigidBody3D& top = w.bodies().back();
                    const Vec3 d = top.position - top0;
                    drift = std::fmax(drift, std::sqrt(d.x * d.x + d.y * d.y));
                    sag = std::fmax(sag, -d.z);
                    for (const RigidBody3D& b : w.bodies())
                    {
                        tilt = std::fmax(tilt, m40_tilt3(b, Vec3{0.0, 0.0, 1.0}));
                    }
                    if (i > 4000)
                    {
                        late = std::fmax(late, m40_ke3(w.bodies()));
                        vlate = std::fmax(vlate, m40_scale3(w.bodies()));
                    }
                }
                const Real r = m40_rest_tol * (m40_gdt + vlate);
                const Real v_rest = m40_rest_speed(r, 0.5, m40_block_R);
                Real ke = 0.0;
                for (const RigidBody3D& b : w.bodies())
                {
                    ke += m40_rest_ke(b.mass, m40_moments3(b), r, 0.5, m40_block_R);
                }
                const Real nn = static_cast<Real>(n);
                std::cout << "  M40 ---- " << name << '\n';
                m40_check("AC-9/10 3D stack: top drift (m)", drift, 1e-3);
                m40_check("AC-9/10 3D stack: max tilt (rad)", tilt, 1e-3);
                m40_check("AC-9/10 3D stack: top sag (m)", sag, nn * m40_interface(n, mass_top, v_rest));
                m40_check("AC-9/10 3D stack: late KE", late, nn * nn * ke);
                m40_check("AC-9/10 3D stack: every step valid", failed_steps, 0.0);
            };
            stack3("3D staggered 3-stack", 3, Vec3{0.1, 0.05, 0.0}, 0.0, 1.0);
            stack3("3D staggered 5-stack", 5, Vec3{0.1, 0.05, 0.0}, 0.0, 1.0);
            stack3("3D aligned 5-stack", 5, Vec3{}, 0.0, 1.0);
            stack3("3D 3-stack yawed 30 degrees per level", 3, Vec3{0.05, 0.0, 0.0}, pi / 6.0, 1.0);
            stack3("AC-10 3D 100 kg on 1 kg", 2, Vec3{0.1, 0.05, 0.0}, 0.0, 100.0);

            // A sphere (r 0.3, m 0.5) resting off centre on a cube on the floor.
            {
                const Vec3 s0{0.2, 0.1, 1.3};
                Rigid3DWorld w{SimulationSettings{m40_dt},
                               {m40_sphere(0.5, 0.3, s0, Vec3{}),
                                m40_box3(1.0, Vec3{0.5, 0.5, 0.5}, Vec3{0.0, 0.0, 0.5}, Quat{}, Vec3{}, Vec3{})},
                               m40_floor3(m40_g, 0.0, 0.5)};
                Real drift = 0.0, sag = 0.0, late = 0.0, vlate = 0.0;
                for (int i = 1; i <= 5000; ++i)
                {
                    w.step();
                    const Vec3 d = w.bodies()[0].position - s0;
                    drift = std::fmax(drift, std::sqrt(d.x * d.x + d.y * d.y));
                    sag = std::fmax(sag, -d.z);
                    if (i > 4000)
                    {
                        late = std::fmax(late, m40_ke3(w.bodies()));
                        vlate = std::fmax(vlate, m40_scale3(w.bodies()));
                    }
                }
                const Real r = m40_rest_tol * (m40_gdt + vlate);
                // The sphere touches at one point, so its speed is bounded
                // through the cube's: the cube's allowance twice (the point
                // moves with the cube) plus r at the point. Its rolling mode
                // has no residual, but nothing excites it: every normal here
                // is vertical and passes through its centre. n = 2 bodies, so
                // the KE bound is n^2 = 4 times the sum.
                const Real v_rest = 2.0 * m40_rest_speed(r, 0.5, length(Vec3{0.5, 0.5, 0.5})) + r;
                std::cout << "  M40 ---- sphere resting off centre on a cube\n";
                m40_check("AC-9 sphere on cube: drift (m)", drift, 1e-3);
                m40_check("AC-9 sphere on cube: sag (m)", sag, 2.0 * (m40_slop_rel * 0.5 + 4.0 * v_rest * m40_dt));
                m40_check("AC-9 sphere on cube: late KE", late,
                          4.0 * (m40_rest_ke(1.0, 0.5, r, 0.5, length(Vec3{0.5, 0.5, 0.5})) +
                                 0.5 * 0.5 * v_rest * v_rest + 0.5 * 0.054 * m40_rest_spin(r, 0.6) *
                                                                    m40_rest_spin(r, 0.6)));
            }

            // AC-10, mass ratio 1e4: its normal rows converge at 1e4 / (1e4 + 1)
            // per sweep, so every solve stops at the cap and none of the
            // converged-solve bounds apply. What must still hold, as
            // REQUIREMENTS: every step valid, the heavy box never past the
            // light one's middle (sag < 0.25), and the mechanical energy never
            // above its start (a passive system at e 0 has no source), to
            // rounding.
            {
                Rigid3DWorld w{SimulationSettings{m40_dt},
                               {m40_box3(1.0, m40_block, Vec3{0.0, 0.0, 0.25}, Quat{}, Vec3{}, Vec3{}),
                                m40_box3(1e4, m40_block, Vec3{0.1, 0.05, 0.75}, Quat{}, Vec3{}, Vec3{})},
                               m40_floor3(m40_g, 0.0, 0.5)};
                const Real e0 = m40_energy3(w);
                Real climb = 0.0, sag = 0.0;
                int failed_steps = 0;
                for (int i = 0; i < 2000; ++i)
                {
                    failed_steps += w.step().ok() ? 0 : 1;
                    climb = std::fmax(climb, m40_energy3(w) - e0);
                    sag = std::fmax(sag, 0.75 - w.bodies()[1].position.z);
                }
                m40_check("AC-10 3D 1e4 kg on 1 kg: every step valid", failed_steps, 0.0);
                m40_check("AC-10 3D 1e4 kg on 1 kg: sag (m)", sag, 0.25);
                m40_check("AC-10 3D 1e4 kg on 1 kg: E - E0", climb, 64.0 * m40_eps * std::fabs(e0));
            }
        }

        // ---------------------------------------------------------------------
        // AC-11. Static friction holds on a slope (issue #44): a box at rest on
        // a 20 degree slope, mu 0.6 > tan 20 deg, e 0, 10 s, at dt 0.002 and
        // 0.01. Contacts are solved before positions move, so nothing creeps
        // but the residual: moved <= T v_rest + slop. Before M40 it crept at
        // g sin(theta) dt per unit time (0.067 m in 10 s at dt 0.002).
        // ---------------------------------------------------------------------
        for (const Real dt : {0.002, 0.01})
        {
            const Real th = 20.0 * pi / 180.0;
            const int steps = static_cast<int>(10.0 / dt + 0.5);
            RigidSettings s2;
            s2.gravity = Vec2{0.0, -m40_g};
            s2.restitution = 0.0;
            s2.friction = 0.6;
            const Vec2 n2{std::sin(th), std::cos(th)};
            s2.ground.push_back(Halfplane{n2, 0.0});
            RigidWorld w2{SimulationSettings{dt}, {m40_box2(1.0, m40_brick, n2 * 0.25, -th, Vec2{}, 0.0)}, s2};
            Rigid3DSettings s3;
            s3.gravity = Vec3{0.0, 0.0, -m40_g};
            s3.restitution = 0.0;
            s3.friction = 0.6;
            const Vec3 n3{std::sin(th), 0.0, std::cos(th)};
            s3.ground.push_back(Plane3{n3, 0.0});
            Rigid3DWorld w3{SimulationSettings{dt},
                            {m40_box3(1.0, m40_block, n3 * 0.25, from_axis_angle(Vec3{0.0, 1.0, 0.0}, th), Vec3{},
                                      Vec3{})},
                            s3};
            Real v2 = 0.0, v3 = 0.0;
            for (int i = 0; i < steps; ++i)
            {
                w2.step();
                w3.step();
                v2 = std::fmax(v2, m40_scale2(w2.bodies()));
                v3 = std::fmax(v3, m40_scale3(w3.bodies()));
            }
            const Real moved2 = dot(w2.bodies()[0].position - n2 * 0.25, Vec2{std::cos(th), -std::sin(th)});
            const Real moved3 = dot(w3.bodies()[0].position - n3 * 0.25, Vec3{std::cos(th), 0.0, -std::sin(th)});
            const Real r2 = m40_rest_tol * (m40_g * dt + v2);
            const Real r3 = m40_rest_tol * (m40_g * dt + v3);
            m40_check(dt > 0.005 ? "AC-11 2D slope dt .01: moved (m)" : "AC-11 2D slope dt .002: moved (m)",
                      std::fabs(moved2), 10.0 * m40_rest_speed(r2, 0.5, m40_brick_R) + m40_slop_rel * 0.5);
            m40_check(dt > 0.005 ? "AC-11 3D slope dt .01: moved (m)" : "AC-11 3D slope dt .002: moved (m)",
                      std::fabs(moved3), 10.0 * m40_rest_speed(r3, 0.5, m40_block_R) + m40_slop_rel * 0.5);
        }

        // ---------------------------------------------------------------------
        // AC-12. Coulomb sliding: a box sliding at 2 m/s on the floor, mu 0.3,
        // decelerates at exactly mu g. Each step the friction impulse is mu
        // times the normal impulse, m (g dt + dv_n) with |dv_n| <= 2 v_rest,
        // so after k steps the speed is within
        // v0 - mu g t: 2 mu k v_rest, with v_rest the AC-1 rest speed for
        // r = 1e-8 v_scale, v_scale about the speed. In 3D the
        // slide is DIAGONAL (30 degrees): a friction pyramid would decelerate
        // it at mu g (cos + sin) and turn it, so the direction is checked too.
        // ---------------------------------------------------------------------
        {
            const Real mu = 0.3;
            const int k = 150;
            const Real v0 = 2.0;
            RigidWorld w2{SimulationSettings{m40_dt},
                          {m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{v0, 0.0}, 0.0)},
                          m40_floor2(m40_g, 0.0, mu)};
            const Vec3 dir{std::cos(pi / 6.0), std::sin(pi / 6.0), 0.0};
            Rigid3DWorld w3{SimulationSettings{m40_dt},
                            {m40_box3(1.0, m40_block, Vec3{0.0, 0.0, 0.25}, Quat{}, dir * v0, Vec3{})},
                            m40_floor3(m40_g, 0.0, mu)};
            for (int i = 0; i < k; ++i)
            {
                w2.step();
                w3.step();
            }
            const Real expected = v0 - mu * m40_g * m40_dt * k;
            const Real r = m40_rest_tol * (v0 + m40_gdt);
            const Real tol2 = 2.0 * mu * k * m40_rest_speed(r, 1.0, m40_brick_R) + 64.0 * m40_eps * v0;
            const Real tol3 = 2.0 * mu * k * m40_rest_speed(r, 0.8, m40_block_R) + 64.0 * m40_eps * v0;
            const Vec3 v3 = w3.bodies()[0].velocity;
            const Vec3 side{-dir.y, dir.x, 0.0};
            m40_check("AC-12 2D sliding: |v - (v0 - mu g t)|", std::fabs(w2.bodies()[0].velocity.x - expected), tol2);
            m40_check("AC-12 3D diagonal slide: |v . dir - (v0 - mu g t)|", std::fabs(dot(v3, dir) - expected), tol3);
            m40_check("AC-12 3D diagonal slide: sideways speed", std::fabs(dot(v3, side)), tol3);
        }

        // ---------------------------------------------------------------------
        // AC-13. Conservation.
        // (a) Body-body contacts, no gravity, no ground (the box_collision
        //     setups, aligned and offset 0.3, e 0.5, mu 0 and 0.5): linear
        //     momentum within 4 eps |p| per impulse application, at most
        //     2 phases x (1 + 2 x 100 x points) per contact step.
        // (b) A frictionless level floor exerts only vertical impulses: a CUBE
        //     (isotropic inertia, so free rotation cannot move L either) dropped
        //     tilted keeps p_x, p_y and the vertical angular momentum L_z, to
        //     rounding (4 eps per step of their scale). A correction that turns
        //     the body but keeps its body-frame spin would rotate L.
        // (c) A contact at ZERO depth (no correction) between two boxes, off
        //     centre: the impulse acts at one shared point, so the angular
        //     momentum about the origin is unchanged to rounding.
        // (d) Aligned head-on boxes give the exact 1D answer (heavy box stops,
        //     light one leaves at 1, no spin) to the impact tolerance.
        // ---------------------------------------------------------------------
        {
            Real worst = 0.0;
            for (const Real off : {0.0, 0.3})
            {
                for (const Real mu : {0.0, 0.5})
                {
                    RigidSettings s2;
                    s2.restitution = 0.5;
                    s2.friction = mu;
                    RigidWorld w2{SimulationSettings{m40_dt},
                                  {m40_box2(2.0, Vec2{0.5, 0.5}, Vec2{-1.2, 0.0}, 0.0, Vec2{1.0, 0.0}, 0.0),
                                   m40_box2(1.0, Vec2{0.5, 0.5}, Vec2{1.2, off}, 0.0, Vec2{-1.0, 0.0}, 0.0)},
                                  s2};
                    Rigid3DSettings s3;
                    s3.restitution = 0.5;
                    s3.friction = mu;
                    Rigid3DWorld w3{SimulationSettings{m40_dt},
                                    {m40_box3(2.0, Vec3{0.5, 0.5, 0.5}, Vec3{-1.2, 0.0, 0.0}, Quat{},
                                              Vec3{1.0, 0.0, 0.0}, Vec3{}),
                                     m40_box3(1.0, Vec3{0.5, 0.5, 0.5}, Vec3{1.2, off, 0.0}, Quat{},
                                              Vec3{-1.0, 0.0, 0.0}, Vec3{})},
                                    s3};
                    const Vec2 p2 = total_linear_momentum(w2.bodies());
                    const Vec3 p3 = malloy::rigid::total_linear_momentum3d(w3.bodies());
                    int contact2 = 0, contact3 = 0;
                    for (int i = 0; i < 2000; ++i)
                    {
                        const Vec2 before2 = w2.bodies()[0].velocity;
                        const Vec3 before3 = w3.bodies()[0].velocity;
                        w2.step();
                        w3.step();
                        contact2 += (w2.bodies()[0].velocity.x != before2.x || w2.bodies()[0].velocity.y != before2.y) ? 1 : 0;
                        contact3 += (w3.bodies()[0].velocity.x != before3.x || w3.bodies()[0].velocity.y != before3.y ||
                                     w3.bodies()[0].velocity.z != before3.z) ? 1 : 0;
                        const Real tol2 = 4.0 * m40_eps * (contact2 + 1) * 2.0 * (1.0 + 2.0 * m40_cap * 2.0) * 3.0;
                        const Real tol3 = 4.0 * m40_eps * (contact3 + 1) * 2.0 * (1.0 + 2.0 * m40_cap * 8.0) * 3.0;
                        worst = std::fmax(worst, length(total_linear_momentum(w2.bodies()) - p2) / tol2);
                        worst = std::fmax(worst, length(malloy::rigid::total_linear_momentum3d(w3.bodies()) - p3) / tol3);
                    }
                }
            }
            m40_check("AC-13a body-body momentum: max |dp| / allowed", worst, 1.0);
        }
        {
            Rigid3DWorld w{SimulationSettings{m40_dt},
                           {m40_box3(1.0, Vec3{0.5, 0.5, 0.5}, Vec3{0.3, -0.2, 1.0},
                                     from_axis_angle(Vec3{1.0, 2.0, 0.0}, 20.0 * pi / 180.0), Vec3{}, Vec3{})},
                           m40_floor3(m40_g, 0.0, 0.0)};
            Real pxy = 0.0, lz = 0.0, scale = 0.0;
            for (int i = 1; i <= 15000; ++i)
            {
                w.step();
                const Vec3 p = malloy::rigid::total_linear_momentum3d(w.bodies());
                const Vec3 L = malloy::rigid::total_angular_momentum3d(w.bodies());
                pxy = std::fmax(pxy, std::fabs(p.x) + std::fabs(p.y));
                lz = std::fmax(lz, std::fabs(L.z));
                scale = std::fmax(scale, length(malloy::rigid::spin_angular_momentum(w.bodies()[0])) +
                                             length(w.bodies()[0].position) * length(p));
            }
            m40_check("AC-13b frictionless floor: |p_x| + |p_y|", pxy, 64.0 * m40_eps * 3.0);
            m40_check("AC-13b frictionless floor: |L_z|", lz, 4.0 * m40_eps * 15000.0 * scale);
        }
        {
            Real worst = 0.0;
            for (const Real mu : {0.0, 0.5})
            {
                RigidSettings s2;
                s2.restitution = 0.5;
                s2.friction = mu;
                RigidWorld w2{SimulationSettings{m40_dt},
                              {m40_box2(2.0, Vec2{0.5, 0.5}, Vec2{-0.5, 0.0}, 0.0, Vec2{1.0, 0.2}, 0.3),
                               m40_box2(1.0, Vec2{0.5, 0.5}, Vec2{0.5, 0.3}, 0.0, Vec2{-1.0, 0.0}, -0.2)},
                              s2};
                const Real L2 = total_angular_momentum(w2.bodies());
                w2.step();
                worst = std::fmax(worst, std::fabs(total_angular_momentum(w2.bodies()) - L2) /
                                             (4.0 * m40_eps * 2.0 * (1.0 + 2.0 * m40_cap * 2.0) * 3.1 * 1.3));
                Rigid3DSettings s3;
                s3.restitution = 0.5;
                s3.friction = mu;
                Rigid3DWorld w3{SimulationSettings{m40_dt},
                                {m40_box3(2.0, Vec3{0.5, 0.5, 0.5}, Vec3{-0.5, 0.0, 0.0}, Quat{}, Vec3{1.0, 0.2, -0.1},
                                          Vec3{0.1, 0.2, 0.3}),
                                 m40_box3(1.0, Vec3{0.5, 0.5, 0.5}, Vec3{0.5, 0.3, 0.2}, Quat{}, Vec3{-1.0, 0.0, 0.0},
                                          Vec3{-0.2, 0.0, 0.1})},
                                s3};
                const Vec3 L3 = malloy::rigid::total_angular_momentum3d(w3.bodies());
                w3.step();
                worst = std::fmax(worst, length(malloy::rigid::total_angular_momentum3d(w3.bodies()) - L3) /
                                             (4.0 * m40_eps * 2.0 * (1.0 + 2.0 * m40_cap * 8.0) * 3.1 * 1.5));
            }
            // Allowed: 4 eps per impulse application, at most 2 phases of
            // (1 + 2 x 100 x points) applications (2 points in 2D, 8 in 3D),
            // times the impulse scale sum m |v| = 3.1 and the arm, 1.3 (2D) and
            // 1.5 (3D).
            m40_check("AC-13c zero-depth contact: |dL| / allowed", worst, 1.0);
        }
        {
            RigidWorld w2{SimulationSettings{m40_dt},
                          {m40_box2(2.0, Vec2{0.5, 0.5}, Vec2{-1.2, 0.0}, 0.0, Vec2{1.0, 0.0}, 0.0),
                           m40_box2(1.0, Vec2{0.5, 0.5}, Vec2{1.2, 0.0}, 0.0, Vec2{-1.0, 0.0}, 0.0)},
                          RigidSettings{0.5}};
            Rigid3DSettings s3;
            s3.restitution = 0.5;
            Rigid3DWorld w3{SimulationSettings{m40_dt},
                            {m40_box3(2.0, Vec3{0.5, 0.5, 0.5}, Vec3{-1.2, 0.0, 0.0}, Quat{}, Vec3{1.0, 0.0, 0.0}, Vec3{}),
                             m40_box3(1.0, Vec3{0.5, 0.5, 0.5}, Vec3{1.2, 0.0, 0.0}, Quat{}, Vec3{-1.0, 0.0, 0.0},
                                      Vec3{})},
                            s3};
            for (int i = 0; i < 2000; ++i)
            {
                w2.step();
                w3.step();
            }
            const Real err2 = length(w2.bodies()[0].velocity) + length(w2.bodies()[1].velocity - Vec2{1.0, 0.0}) +
                              std::fabs(w2.bodies()[0].angular_velocity) + std::fabs(w2.bodies()[1].angular_velocity);
            const Real err3 = length(w3.bodies()[0].velocity) + length(w3.bodies()[1].velocity - Vec3{1.0, 0.0, 0.0}) +
                              length(w3.bodies()[0].angular_velocity) + length(w3.bodies()[1].angular_velocity);
            // Four velocity errors, each within 2 sqrt(3) x 1e-12 v_scale / w
            // (w 1, v_scale 1) of the exact answer, plus rounding.
            m40_check("AC-13d 2D aligned head-on: exact 1D answer", err2, 4.0 * 2.0 * std::sqrt(3.0) * m40_impact_tol + 1e-14);
            m40_check("AC-13d 3D aligned head-on: exact 1D answer", err3, 4.0 * 2.0 * std::sqrt(3.0) * m40_impact_tol + 1e-14);
        }

        // ---------------------------------------------------------------------
        // AC-14. Determinism and the documented order.
        // (a) The same inputs give the same output bit for bit, and a world
        //     copied mid-run (its warm-start cache included) steps exactly as
        //     its original.
        // (b) The order, pinned by REPLICAS of one rest-phase solve written
        //     from the documented contract, not by golden literals. Two discs
        //     (spheres in 3D) stacked on the floor, exactly touching, sliding
        //     together at 10 m/s so v_scale is 10 and the solve stops early
        //     enough for its path to show: the residual is checked before every
        //     sweep, and the pair manifold comes before the ground. The reversed order
        //     lands 5e-8 m/s away; the replica agrees to rounding, bounded by a
        //     separation of 4 ulp turned into a speculative target, 4 ulp / dt
        //     (5.5e-13), times 4.
        //     And a brick sliding on the floor (mu 0.3): the manifold prepass at
        //     its weighted centroid (normal, then friction from the manifold
        //     total), then per sweep all normal rows, then all friction rows.
        //     Interleaving, friction first, no prepass or reversed points each
        //     land 7.8e-10 to 1.1e-8 away; the replica agrees within 1e-13.
        // ---------------------------------------------------------------------
        {
            const auto scene3 = [&]() {
                return std::vector<RigidBody3D>{
                    m40_box3(1.0, m40_block, Vec3{0.0, 0.0, 0.25}, Quat{}, Vec3{}, Vec3{}),
                    m40_box3(1.3, Vec3{0.4, 0.3, 0.2}, Vec3{0.1, 0.05, 0.7}, from_axis_angle(Vec3{0.0, 0.0, 1.0}, 0.3),
                             Vec3{}, Vec3{}),
                    m40_sphere(0.5, 0.2, Vec3{0.15, 0.0, 1.1}, Vec3{0.1, 0.0, 0.0}),
                    m40_box3(0.8, Vec3{0.3, 0.25, 0.2}, Vec3{2.0, 0.0, 1.0},
                             from_axis_angle(Vec3{1.0, 2.0, 0.0}, 0.4), Vec3{-0.5, 0.0, 0.0}, Vec3{0.3, -0.2, 0.5})};
            };
            Rigid3DWorld a{SimulationSettings{m40_dt}, scene3(), m40_floor3(m40_g, 0.3, 0.5)};
            Rigid3DWorld b{SimulationSettings{m40_dt}, scene3(), m40_floor3(m40_g, 0.3, 0.5)};
            for (int i = 0; i < 500; ++i)
            {
                a.step();
                b.step();
            }
            Rigid3DWorld c = a;
            for (int i = 0; i < 500; ++i)
            {
                a.step();
                b.step();
                c.step();
            }
            int differ = 0;
            for (std::size_t k = 0; k < a.bodies().size(); ++k)
            {
                differ += m40_same3(a.bodies()[k], b.bodies()[k]) ? 0 : 1;
                differ += m40_same3(a.bodies()[k], c.bodies()[k]) ? 0 : 1;
            }
            differ += m40_same_warm3(a, b) ? 0 : 1;
            differ += m40_same_warm3(a, c) ? 0 : 1;
            differ += a.warm_contacts().empty() ? 1 : 0; // the cache was in play
            const auto scene2 = [&]() {
                return std::vector<RigidBody2D>{
                    m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{}, 0.0),
                    m40_box2(1.3, Vec2{0.4, 0.2}, Vec2{0.1, 0.7}, 0.1, Vec2{}, 0.0),
                    m40_box2(0.8, Vec2{0.3, 0.2}, Vec2{1.5, 1.0}, 0.4, Vec2{-0.5, 0.0}, 0.5)};
            };
            RigidWorld a2{SimulationSettings{m40_dt}, scene2(), m40_floor2(m40_g, 0.3, 0.5)};
            RigidWorld b2{SimulationSettings{m40_dt}, scene2(), m40_floor2(m40_g, 0.3, 0.5)};
            for (int i = 0; i < 500; ++i)
            {
                a2.step();
                b2.step();
            }
            RigidWorld c2 = a2;
            for (int i = 0; i < 500; ++i)
            {
                a2.step();
                b2.step();
                c2.step();
            }
            for (std::size_t k = 0; k < a2.bodies().size(); ++k)
            {
                differ += m40_same2(a2.bodies()[k], b2.bodies()[k]) ? 0 : 1;
                differ += m40_same2(a2.bodies()[k], c2.bodies()[k]) ? 0 : 1;
            }
            differ += m40_same_warm2(a2, b2) ? 0 : 1;
            differ += m40_same_warm2(a2, c2) ? 0 : 1;
            differ += a2.warm_contacts().empty() ? 1 : 0;
            m40_check("AC-14a rerun and mid-run copy are bit-identical, caches too (differences)", differ,
                      0.0);
        }
        {
            // The disc replica: masses 2 (radius 0.5) and 0.5 (radius 0.25),
            // binary-exact so the touch is exact. Rows: pair (n up, K = 1/m0 +
            // 1/m1, centred, so no rotation), then ground (n down, K = 1/m0).
            const Real m0 = 2.0, m1 = 0.5, vx = 10.0;
            RigidWorld w2{SimulationSettings{m40_dt},
                          {m40_disc(m0, 0.5, Vec2{0.0, 0.5}, Vec2{vx, 0.0}), m40_disc(m1, 0.25, Vec2{0.0, 1.25}, Vec2{vx, 0.0})},
                          m40_floor2(m40_g, 0.0, 0.0)};
            Rigid3DWorld w3{SimulationSettings{m40_dt},
                            {m40_sphere(m0, 0.5, Vec3{0.0, 0.0, 0.5}, Vec3{vx, 0.0, 0.0}),
                             m40_sphere(m1, 0.25, Vec3{0.0, 0.0, 1.25}, Vec3{vx, 0.0, 0.0})},
                            m40_floor3(m40_g, 0.0, 0.0)};
            w2.step();
            w3.step();
            // Replica.
            Real v0 = -m40_gdt, v1 = -m40_gdt; // vertical, after forcing
            const Real vs = std::sqrt(vx * vx + m40_gdt * m40_gdt);
            Real lp = 0.0, lg = 0.0;
            for (int it = 0; it < 100; ++it)
            {
                const Real res = std::fmax(lp > 0.0 ? std::fabs(v1 - v0) : std::fmax(0.0, v0 - v1),
                                           lg > 0.0 ? std::fabs(v0) : std::fmax(0.0, -v0));
                if (res <= m40_rest_tol * vs)
                {
                    break;
                }
                {
                    const Real vn = v1 - v0;
                    Real d = -vn / (1.0 / m0 + 1.0 / m1);
                    const Real nl = std::fmax(0.0, lp + d);
                    d = nl - lp;
                    lp = nl;
                    v0 -= d / m0;
                    v1 += d / m1;
                }
                {
                    const Real vn = v0;
                    Real d = -vn / (1.0 / m0);
                    const Real nl = std::fmax(0.0, lg + d);
                    d = nl - lg;
                    lg = nl;
                    v0 += d / m0;
                }
            }
            const Real tol = 4.0 * 4.0 * 1.25 * m40_eps / m40_dt;
            const Real err2 = std::fabs(w2.bodies()[0].velocity.y - v0) + std::fabs(w2.bodies()[1].velocity.y - v1) +
                              std::fabs(w2.bodies()[0].velocity.x - vx) + std::fabs(w2.bodies()[1].velocity.x - vx);
            const Real err3 = std::fabs(w3.bodies()[0].velocity.z - v0) + std::fabs(w3.bodies()[1].velocity.z - v1) +
                              std::fabs(w3.bodies()[0].velocity.x - vx) + std::fabs(w3.bodies()[1].velocity.x - vx);
            m40_check("AC-14b 2D disc-stack replica (pair, then ground)", err2, tol);
            m40_check("AC-14b 3D sphere-stack replica (pair, then ground)", err3, tol);
        }
        {
            // The sliding-brick replica: half (0.5, 0.25), mass 1, on y = 0
            // exactly, sliding at 2 m/s, mu 0.3, e 0. Corners in collide order
            // (left, right), both exactly on the floor: a LEVEL manifold, so
            // equal weights and the centroid arm (0, -0.25).
            const Real mu = 0.3;
            RigidWorld w{SimulationSettings{m40_dt},
                         {m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{2.0, 0.0}, 0.0)},
                         m40_floor2(m40_g, 0.0, mu)};
            w.step();
            const Real I = m40_brick_I;
            Real vx = 2.0, vy = -m40_gdt, om = 0.0;
            const Vec2 arms[2]{Vec2{-0.5, -0.25}, Vec2{0.5, -0.25}};
            Real ln[2]{0.0, 0.0}, lt[2]{0.0, 0.0};
            // The box is body a, the floor b; n = (0, -1) from box to floor,
            // t = (1, 0). Velocity of the box at arm r: (vx - om r.y, vy + om r.x).
            // An impulse J on the floor puts -J on the box at r.
            const auto push = [&](Vec2 r, Vec2 J) {
                vx -= J.x;
                vy -= J.y;
                om -= (r.x * J.y - r.y * J.x) / I;
            };
            {
                const Vec2 rc{0.0, -0.25};
                const Real vn = vy + om * rc.x; // (0 - v_box) . n
                const Real kc = 1.0 + (rc.x * rc.x) / I;
                const Real lam = std::fmax(0.0, (0.0 - vn) / kc);
                push(rc, Vec2{0.0, -lam});
                ln[0] = ln[1] = lam * 0.5;
                const Real kt = 1.0 + (rc.y * rc.y) / I;
                const Real vt = -(vx - om * rc.y);
                const Real total = ln[0] + ln[1];
                const Real P = std::fmax(-mu * total, std::fmin(mu * total, -vt / kt));
                push(rc, Vec2{P, 0.0});
                lt[0] = lt[1] = P * 0.5;
            }
            // v_scale after the prepass: |v| + |w| size (0.5), or g dt.
            const Real vs = std::fmax(m40_gdt, std::sqrt(vx * vx + vy * vy) + std::fabs(om) * 0.5);
            const auto residual = [&]() {
                Real res = 0.0;
                for (int k = 0; k < 2; ++k)
                {
                    const Vec2 r = arms[k];
                    const Real vn = vy + om * r.x;
                    const Real vt = -(vx - om * r.y);
                    res = std::fmax(res, ln[k] > 0.0 ? std::fabs(vn) : std::fmax(0.0, -vn));
                    const Real lim = mu * ln[k];
                    if (std::fabs(lt[k]) < lim * (1.0 - 1e-12))
                    {
                        res = std::fmax(res, std::fabs(vt));
                    }
                    else if (lt[k] != 0.0)
                    {
                        res = std::fmax(res, std::fmax(0.0, vt * (lt[k] > 0.0 ? 1.0 : -1.0)));
                    }
                }
                return res;
            };
            for (int it = 0; it < 100; ++it)
            {
                if (residual() <= m40_rest_tol * vs)
                {
                    break;
                }
                for (int k = 0; k < 2; ++k)
                {
                    const Vec2 r = arms[k];
                    const Real vn = vy + om * r.x;
                    Real d = (0.0 - vn) / (1.0 + r.x * r.x / I);
                    const Real nl = std::fmax(0.0, ln[k] + d);
                    d = nl - ln[k];
                    ln[k] = nl;
                    push(r, Vec2{0.0, -d});
                }
                for (int k = 0; k < 2; ++k)
                {
                    const Vec2 r = arms[k];
                    const Real vt = -(vx - om * r.y);
                    Real d = -vt / (1.0 + r.y * r.y / I);
                    const Real lim = mu * ln[k];
                    const Real nl = std::fmax(-lim, std::fmin(lim, lt[k] + d));
                    d = nl - lt[k];
                    lt[k] = nl;
                    push(r, Vec2{d, 0.0});
                }
            }
            const RigidBody2D& b = w.bodies()[0];
            const Real err = std::fabs(b.velocity.x - vx) + std::fabs(b.velocity.y - vy) +
                             std::fabs(b.angular_velocity - om);
            m40_check("AC-14b 2D sliding-brick replica (prepass, normals, frictions)", err, 1e-13);
        }

        // ---------------------------------------------------------------------
        // AC-15. A rolled-back step restores the WHOLE state, the warm-start
        // cache included. A far body with no collider sits at x = 1.339e154,
        // squarable, but not after one step at 1e154 m/s. A 2D stack is mid
        // landing (a transient, so the cache changes every step). An impulse
        // of 1e154 makes the next step fail AFTER the solve (whose v_scale is
        // then 1e154, so it stops after one sweep and leaves a different
        // cache), and the opposite impulse undoes it exactly (0 + 1e154 -
        // 1e154). A copy taken before must then step bit for bit the same.
        // ---------------------------------------------------------------------
        {
            RigidBody2D far;
            far.mass = 1.0;
            far.inertia = 1.0;
            far.position = Vec2{1.339e154, 0.0};
            RigidWorld w{SimulationSettings{m40_dt},
                         {m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{}, 0.0),
                          m40_box2(0.8, Vec2{0.4, 0.2}, Vec2{0.1, 0.75}, 0.05, Vec2{0.3, -0.5}, 0.2), far},
                         m40_floor2(m40_g, 0.2, 0.5)};
            for (int i = 0; i < 30; ++i)
            {
                w.step(); // the top box lands within these steps
            }
            RigidWorld copy = w;
            const Vec2 J{1e154, 0.0};
            w.apply_impulse_at(2, J, w.bodies()[2].position);
            const bool failed = !w.step().ok();
            int restored = m40_same_warm2(w, copy) ? 0 : 1; // the cache, right after
            restored += copy.warm_contacts().empty() ? 1 : 0; // and it held something
            w.apply_impulse_at(2, J * -1.0, w.bodies()[2].position);
            copy.apply_impulse_at(2, J, copy.bodies()[2].position);
            copy.apply_impulse_at(2, J * -1.0, copy.bodies()[2].position);
            int differ = 0;
            for (int i = 0; i < 300; ++i)
            {
                w.step();
                copy.step();
            }
            for (std::size_t k = 0; k < w.bodies().size(); ++k)
            {
                differ += m40_same2(w.bodies()[k], copy.bodies()[k]) ? 0 : 1;
            }
            differ += m40_same_warm2(w, copy) ? 0 : 1;
            m40_check("AC-15 the trigger step failed", failed ? 0.0 : 1.0, 0.0);
            m40_check("AC-15 the failed step left the warm-start cache as it was", restored, 0.0);
            m40_check("AC-15 after a rolled-back step, differences from a copy (bodies, cache)", differ, 0.0);
        }

        // ---------------------------------------------------------------------
        // AC-16. Bodies that touch nothing are untouched: a free tumbling body
        // beside a resting stack steps bit for bit as the same body alone, and
        // a world whose colliders never meet steps bit for bit as the same
        // world with the colliders removed.
        // ---------------------------------------------------------------------
        {
            const RigidBody3D tumbler = m40_box3(0.8, Vec3{0.3, 0.25, 0.2}, Vec3{5.0, 0.0, 30.0},
                                                 from_axis_angle(Vec3{1.0, 2.0, 0.0}, 0.4), Vec3{0.5, 0.0, 2.0},
                                                 Vec3{0.3, -2.0, 0.5});
            Rigid3DWorld with{SimulationSettings{m40_dt},
                              {m40_box3(1.0, m40_block, Vec3{0.0, 0.0, 0.25}, Quat{}, Vec3{}, Vec3{}),
                               m40_box3(1.0, m40_block, Vec3{0.1, 0.0, 0.75}, Quat{}, Vec3{}, Vec3{}), tumbler},
                              m40_floor3(m40_g, 0.0, 0.5)};
            Rigid3DWorld alone{SimulationSettings{m40_dt}, {tumbler}, m40_floor3(m40_g, 0.0, 0.5)};
            RigidBody3D ghost = tumbler;
            ghost.half_extents = Vec3{};
            Rigid3DSettings bare = m40_floor3(m40_g, 0.0, 0.5);
            bare.ground.clear();
            Rigid3DWorld none{SimulationSettings{m40_dt}, {ghost}, bare};
            for (int i = 0; i < 1000; ++i)
            {
                with.step();
                alone.step();
                none.step();
            }
            int differ = m40_same3(with.bodies()[2], alone.bodies()[0]) ? 0 : 1;
            RigidBody3D a = alone.bodies()[0];
            a.half_extents = Vec3{};
            differ += m40_same3(a, none.bodies()[0]) ? 0 : 1;
            const RigidBody2D tumbler2 =
                m40_box2(0.8, Vec2{0.3, 0.2}, Vec2{5.0, 30.0}, 0.4, Vec2{0.5, 2.0}, -2.0);
            RigidWorld with2{SimulationSettings{m40_dt},
                             {m40_box2(1.0, m40_brick, Vec2{0.0, 0.25}, 0.0, Vec2{}, 0.0),
                              m40_box2(1.0, m40_brick, Vec2{0.1, 0.75}, 0.0, Vec2{}, 0.0), tumbler2},
                             m40_floor2(m40_g, 0.0, 0.5)};
            RigidWorld alone2{SimulationSettings{m40_dt}, {tumbler2}, m40_floor2(m40_g, 0.0, 0.5)};
            for (int i = 0; i < 1000; ++i)
            {
                with2.step();
                alone2.step();
            }
            differ += m40_same2(with2.bodies()[2], alone2.bodies()[0]) ? 0 : 1;
            m40_check("AC-16 untouched bodies are bit-identical (differing)", differ, 0.0);
        }

        // ---------------------------------------------------------------------
        // AC-17. An overlap is projected out by positions alone. A box placed
        // level and at rest 0.01 into the floor (deeper than the slop, and
        // within one capped push, 0.2 x size) is lifted in ONE step to the
        // slop's depth, 1e-6 x size: the level manifold's centroid push takes
        // the whole excess out at once, and the per-point pushes after it find
        // none. The push is exact to rounding (1e-12, about 1e4 ulps of the
        // height), and the step's own motion is the solved velocity times dt,
        // at most v_rest dt. The push must not touch a velocity, so the box
        // leaves the step within the rest bounds, |v| <= v_rest and
        // |omega| <= omega_rest for r = 1e-8 (|g| dt + |v| + |w| size), not at
        // the 0.01 / dt = 5 m/s a velocity-level correction would give it. Then
        // it stays: 1000 steps on, it is at the slop's depth within 4 v_rest dt.
        // With the projection off it would stay 0.01 deep, and with a larger
        // slop it would stop short of 1e-6 x size.
        // ---------------------------------------------------------------------
        {
            const Real d0 = 0.01;
            const std::vector<Vec2> halves2{Vec2{0.5, 0.5}, Vec2{0.5, 0.25}};
            for (const Vec2 h : halves2)
            {
                RigidWorld w{SimulationSettings{m40_dt},
                             {m40_box2(1.0, h, Vec2{0.0, h.y - d0}, 0.0, Vec2{}, 0.0)},
                             m40_floor2(m40_g, 0.5, 0.5)};
                const Real slop = m40_slop_rel * std::fmax(h.x, h.y);
                const Real wmin = 2.0 * std::fmin(h.x, h.y);
                w.step();
                const Real r = m40_rest_tol * (m40_gdt + m40_scale2(w.bodies()));
                const Real v_rest = m40_rest_speed(r, wmin, length(h));
                m40_check("AC-17 2D overlap: one step leaves it at the slop's depth, |y - (h - slop)|",
                          std::fabs(w.bodies()[0].position.y - (h.y - slop)), v_rest * m40_dt + 1e-12);
                m40_check("AC-17 2D overlap: the push left the velocity at rest, |v|", m40_vmax2(w.bodies()),
                          v_rest);
                m40_check("AC-17 2D overlap: and the spin, |omega|", std::fabs(w.bodies()[0].angular_velocity),
                          m40_rest_spin(r, wmin));
                for (int i = 0; i < 1000; ++i)
                {
                    w.step();
                }
                const Real rr = m40_rest_tol * (m40_gdt + m40_scale2(w.bodies()));
                m40_check("AC-17 2D overlap: it stays there, |y - (h - slop)|",
                          std::fabs(w.bodies()[0].position.y - (h.y - slop)),
                          4.0 * m40_rest_speed(rr, wmin, length(h)) * m40_dt + 1e-12);
            }
            const std::vector<Vec3> halves3{Vec3{0.5, 0.5, 0.5}, Vec3{0.7, 0.3, 0.2}};
            for (const Vec3 h : halves3)
            {
                Rigid3DWorld w{SimulationSettings{m40_dt},
                               {m40_box3(1.0, h, Vec3{0.0, 0.0, h.z - d0}, Quat{}, Vec3{}, Vec3{})},
                               m40_floor3(m40_g, 0.5, 0.5)};
                const Real slop = m40_slop_rel * std::fmax(std::fmax(h.x, h.y), h.z);
                const Real wmin = 2.0 * std::fmin(h.x, h.y);
                w.step();
                const Real r = m40_rest_tol * (m40_gdt + m40_scale3(w.bodies()));
                const Real v_rest = m40_rest_speed(r, wmin, length(h));
                m40_check("AC-17 3D overlap: one step leaves it at the slop's depth, |z - (h - slop)|",
                          std::fabs(w.bodies()[0].position.z - (h.z - slop)), v_rest * m40_dt + 1e-12);
                m40_check("AC-17 3D overlap: the push left the velocity at rest, |v|", m40_vmax3(w.bodies()),
                          v_rest);
                m40_check("AC-17 3D overlap: and the spin, |omega|", length(w.bodies()[0].angular_velocity),
                          m40_rest_spin(r, wmin));
                for (int i = 0; i < 1000; ++i)
                {
                    w.step();
                }
                const Real rr = m40_rest_tol * (m40_gdt + m40_scale3(w.bodies()));
                m40_check("AC-17 3D overlap: it stays there, |z - (h - slop)|",
                          std::fabs(w.bodies()[0].position.z - (h.z - slop)),
                          4.0 * m40_rest_speed(rr, wmin, length(h)) * m40_dt + 1e-12);
            }
        }

        std::cout << "  M40 acceptance: " << m40_failed << " criteria failed\n";
        MALLOY_CHECK_EQ(m40_failed, 0);
    }
    std::cout << "malloy_rigid_tests passed\n";
    return 0;
}
