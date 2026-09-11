#include <malloy/math/math.hpp>

#include <test_check.hpp>

#include <iostream>
#include <limits>

// Quaternion tests, split out of vec2_tests.cpp once M20 gave them enough
// coverage to be their own translation unit. `Quat` is what a scalar angle
// becomes in three dimensions, and the properties worth pinning are the ones a
// scalar angle does not have: composition does not commute, and a quaternion
// carries a magnitude the rotation formulas depend on.
int main()
{
    using malloy::math::Real;

    // --- M20: quaternions. This is what a scalar angle becomes in three
    //     dimensions, and the properties worth pinning are the ones a scalar
    //     angle does not have. ---
    {
        using malloy::math::Quat;
        using malloy::math::Vec3;
        const Real inf = std::numeric_limits<Real>::infinity();

        const Vec3 ex{1.0, 0.0, 0.0};
        const Vec3 ey{0.0, 1.0, 0.0};
        const Vec3 ez{0.0, 0.0, 1.0};
        const Real quarter = 1.57079632679489661923;

        // The default is the IDENTITY rotation, not the zero quaternion, and
        // rotating by it is exact rather than nearly so.
        MALLOY_CHECK_TRUE(malloy::math::is_unit(Quat{}));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::rotate(Quat{}, Vec3{3.0, -4.0, 5.0}),
            Vec3{3.0, -4.0, 5.0}, 0.0));

        // A quarter turn about z takes x to y, right-handed.
        const Quat qz = malloy::math::from_axis_angle(ez, quarter);
        MALLOY_CHECK_TRUE(malloy::math::is_unit(qz));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::rotate(qz, ex), ey, 1e-15));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::rotate(qz, ez), ez, 1e-15)); // the axis is fixed

        // Rotation preserves length, which is the property every use here
        // depends on: it is why |L| does not depend on the orientation.
        const Vec3 p{1.0, -2.0, 3.0};
        const Quat tilted = malloy::math::from_axis_angle(Vec3{1.0, 2.0, -3.0}, 0.7);
        MALLOY_CHECK_NEAR(malloy::math::length(malloy::math::rotate(tilted, p)),
                          malloy::math::length(p), 1e-14);

        // --- Composition does NOT commute, which is the whole content of 3D
        //     rotation and has no 2D analogue: turning about z then about x
        //     does not land where the reverse order does. ---
        const Quat qx = malloy::math::from_axis_angle(ex, quarter);
        const Vec3 zx = malloy::math::rotate(qz * qx, ex);
        const Vec3 xz = malloy::math::rotate(qx * qz, ex);
        MALLOY_CHECK_FALSE(malloy::math::approx_equal(zx, xz, 1e-6));

        // And it agrees with applying them one at a time, in the matching
        // order. q1*q2 applied to p means q1 applied to (q2 applied to p).
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::rotate(qz * qx, p),
            malloy::math::rotate(qz, malloy::math::rotate(qx, p)), 1e-14));

        // The conjugate undoes the rotation.
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::rotate(malloy::math::conjugate(tilted),
                                 malloy::math::rotate(tilted, p)),
            p, 1e-14));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            tilted * malloy::math::conjugate(tilted), Quat{}, 1e-15));

        // Two quarter turns about the same axis make a half turn.
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::rotate(qz * qz, ex), -ex, 1e-15));

        // --- is_unit is stricter than is_finite, and the rotation formulas
        //     assume the strict one. ---
        const Quat scaled{2.0, Vec3{0.0, 0.0, 0.0}};
        MALLOY_CHECK_TRUE(malloy::math::is_finite(scaled));
        MALLOY_CHECK_FALSE(malloy::math::is_unit(scaled));
        MALLOY_CHECK_TRUE(malloy::math::is_unit(malloy::math::normalize(scaled)));
        MALLOY_CHECK_FALSE(malloy::math::is_unit(Quat{inf, Vec3{}}));

        // A quaternion with no norm has no direction, so normalize returns the
        // identity rather than NaN, matching the Vec2 and Vec3 rule.
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::normalize(Quat{0.0, Vec3{}}), Quat{}, 0.0));

        // A zero axis gives the identity rather than NaN.
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::from_axis_angle(Vec3{}, 1.0), Quat{}, 0.0));

        // --- Why normalization is required: `rotate` uses the optimized
        //     unit-only form, which is NOT scale invariant. Fed a quaternion of
        //     norm k it returns p + k^2 (R p - p), scaling the displacement
        //     rather than performing the rotation. A unit q gives k = 1 and the
        //     exact rotation; a q of norm 2 overshoots by a factor of four on
        //     the part that moves. ---
        {
            const Quat doubled = qz * Real{2};
            const Vec3 turned = malloy::math::rotate(qz, ex); // R ex = ey
            const Vec3 scaled_rot = malloy::math::rotate(doubled, ex);
            const Vec3 predicted = ex + (turned - ex) * Real{4};
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(scaled_rot, predicted, 1e-14));
            // Not the rotation, so this is a defect the normalize call prevents
            // rather than a harmless scaling.
            MALLOY_CHECK_FALSE(malloy::math::approx_equal(scaled_rot, turned, 1e-6));
            // Normalizing first restores it exactly.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::math::rotate(malloy::math::normalize(doubled), ex), turned, 1e-15));
        }
    }

    std::cout << "malloy_quat_tests passed\n";
    return 0;
}
