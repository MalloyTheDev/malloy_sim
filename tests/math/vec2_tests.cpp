#include <malloy/math/math.hpp>

#include <test_check.hpp>

#include <iostream>
#include <limits>

int main()
{
    using malloy::math::approx_equal;
    using malloy::math::distance;
    using malloy::math::distance_squared;
    using malloy::math::dot;
    using malloy::math::is_finite;
    using malloy::math::length;
    using malloy::math::length_squared;
    using malloy::math::normalize;
    using malloy::math::Real;
    using malloy::math::Vec2;

    const Real eps = 1e-12;

    // Constructors: default is the zero vector; two-arg stores components.
    {
        const Vec2 zero;
        MALLOY_CHECK_NEAR(zero.x, 0.0, eps);
        MALLOY_CHECK_NEAR(zero.y, 0.0, eps);

        const Vec2 v{3.0, 4.0};
        MALLOY_CHECK_NEAR(v.x, 3.0, eps);
        MALLOY_CHECK_NEAR(v.y, 4.0, eps);
    }

    // Arithmetic operators, including symmetric scalar multiply and unary minus.
    {
        const Vec2 a{1.0, 2.0};
        const Vec2 b{3.0, 5.0};

        MALLOY_CHECK_VEC2_NEAR(a + b, Vec2(4.0, 7.0), eps);
        MALLOY_CHECK_VEC2_NEAR(b - a, Vec2(2.0, 3.0), eps);
        MALLOY_CHECK_VEC2_NEAR(a * 2.0, Vec2(2.0, 4.0), eps);
        MALLOY_CHECK_VEC2_NEAR(2.0 * a, Vec2(2.0, 4.0), eps);
        MALLOY_CHECK_VEC2_NEAR(b / 2.0, Vec2(1.5, 2.5), eps);
        MALLOY_CHECK_VEC2_NEAR(-a, Vec2(-1.0, -2.0), eps);
    }

    // Compound assignment operators mutate in place.
    {
        Vec2 v{1.0, 1.0};
        v += Vec2{2.0, 3.0};
        MALLOY_CHECK_VEC2_NEAR(v, Vec2(3.0, 4.0), eps);
        v -= Vec2{1.0, 1.0};
        MALLOY_CHECK_VEC2_NEAR(v, Vec2(2.0, 3.0), eps);
        v *= 2.0;
        MALLOY_CHECK_VEC2_NEAR(v, Vec2(4.0, 6.0), eps);
        v /= 2.0;
        MALLOY_CHECK_VEC2_NEAR(v, Vec2(2.0, 3.0), eps);
    }

    // dot, lengths, distances (3-4-5 triangle).
    {
        const Vec2 a{3.0, 4.0};
        MALLOY_CHECK_NEAR(dot(a, a), 25.0, eps);
        MALLOY_CHECK_NEAR(length_squared(a), 25.0, eps);
        MALLOY_CHECK_NEAR(length(a), 5.0, eps);
        MALLOY_CHECK_NEAR(dot(Vec2(1.0, 0.0), Vec2(0.0, 1.0)), 0.0, eps);

        MALLOY_CHECK_NEAR(distance_squared(a, Vec2(0.0, 0.0)), 25.0, eps);
        MALLOY_CHECK_NEAR(distance(a, Vec2(0.0, 0.0)), 5.0, eps);
    }

    // normalize: unit length, and a defined, finite answer for the zero vector.
    {
        const Vec2 n = normalize(Vec2{3.0, 4.0});
        MALLOY_CHECK_NEAR(length(n), 1.0, eps);
        MALLOY_CHECK_VEC2_NEAR(n, Vec2(0.6, 0.8), eps);

        const Vec2 z = normalize(Vec2{0.0, 0.0});
        MALLOY_CHECK_TRUE(is_finite(z));
        MALLOY_CHECK_VEC2_NEAR(z, Vec2(0.0, 0.0), eps);
    }

    // approx_equal (scalar): explicit tolerance, never ==.
    {
        MALLOY_CHECK_TRUE(approx_equal(1.0, 1.0 + 1e-13, eps));
        MALLOY_CHECK_FALSE(approx_equal(1.0, 1.1, eps));
    }

    // approx_equal (vector): component-wise.
    {
        MALLOY_CHECK_TRUE(approx_equal(Vec2(1.0, 2.0), Vec2(1.0 + 1e-13, 2.0), eps));
        MALLOY_CHECK_FALSE(approx_equal(Vec2(1.0, 2.0), Vec2(1.0, 2.5), eps));
    }

    // is_finite: false for infinity and NaN, on scalars and vectors.
    {
        const Real inf = std::numeric_limits<Real>::infinity();
        const Real nan = std::numeric_limits<Real>::quiet_NaN();

        MALLOY_CHECK_TRUE(is_finite(0.0));
        MALLOY_CHECK_TRUE(is_finite(-123.456));
        MALLOY_CHECK_FALSE(is_finite(inf));
        MALLOY_CHECK_FALSE(is_finite(nan));

        MALLOY_CHECK_TRUE(is_finite(Vec2(1.0, 2.0)));
        MALLOY_CHECK_FALSE(is_finite(Vec2(inf, 0.0)));
        MALLOY_CHECK_FALSE(is_finite(Vec2(0.0, nan)));
    }

    std::cout << "malloy_math_tests passed\n";
    return 0;
}
