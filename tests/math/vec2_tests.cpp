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

    // --- distance between two points that are BOTH away from the origin.
    //     Every existing case passed Vec2(0,0) as the second argument, where
    //     a - b and a + b are identical, so a sign error was invisible. ---
    {
        const Vec2 a(1.0, 2.0);
        const Vec2 b(4.0, 6.0);
        MALLOY_CHECK_NEAR(distance_squared(a, b), 25.0, eps); // a+b would give 89
        MALLOY_CHECK_NEAR(distance(a, b), 5.0, eps);
        MALLOY_CHECK_NEAR(distance(b, a), 5.0, eps); // symmetric
    }

    // --- dot of two distinct vectors with a nonzero result. The existing cases
    //     are dot(a, a) and a perpendicular pair, so neither pins the cross
    //     terms of a general product. ---
    {
        MALLOY_CHECK_NEAR(dot(Vec2(1.0, 2.0), Vec2(3.0, 5.0)), 13.0, eps);
        MALLOY_CHECK_NEAR(dot(Vec2(3.0, 5.0), Vec2(1.0, 2.0)), 13.0, eps); // symmetric
        MALLOY_CHECK_NEAR(dot(Vec2(1.0, 2.0), Vec2(-3.0, 1.0)), -1.0, eps);
    }

    // --- approx_equal(Vec2) with x as the SOLE difference. The existing true
    //     case differs only in x within tolerance and the false case only in y,
    //     so an implementation checking y alone passes both. ---
    {
        MALLOY_CHECK_FALSE(approx_equal(Vec2(1.0, 2.0), Vec2(1.5, 2.0), eps));
        MALLOY_CHECK_FALSE(approx_equal(Vec2(1.0, 2.0), Vec2(1.5, 2.5), eps));
        MALLOY_CHECK_TRUE(approx_equal(Vec2(1.0, 2.0), Vec2(1.0, 2.0), eps));
    }

    // --- approx_equal is documented as <=, so a difference of exactly epsilon
    //     compares equal. The boundary is exactly representable here. ---
    {
        MALLOY_CHECK_TRUE(approx_equal(1.0, 1.25, 0.25));
        MALLOY_CHECK_FALSE(approx_equal(1.0, 1.26, 0.25));
        MALLOY_CHECK_TRUE(approx_equal(Vec2(0.0, 0.0), Vec2(0.25, -0.25), 0.25));
    }

    // --- Issue #16: is_squarable is not is_finite. length_squared squares each
    //     component, so it reaches infinity above sqrt(DBL_MAX), about
    //     1.34e154, while the vector itself is still finite and every
    //     component is a normal number.
    //
    //     That gap is why validation testing only is_finite accepted state up
    //     to 1.8e308 while every energy computed from it was inf. ---
    {
        const Real inf = std::numeric_limits<Real>::infinity();
        const Real nan = std::numeric_limits<Real>::quiet_NaN();

        MALLOY_CHECK_TRUE(malloy::math::is_squarable(Vec2{}));
        MALLOY_CHECK_TRUE(malloy::math::is_squarable(Vec2{3.0, -4.0}));
        MALLOY_CHECK_TRUE(malloy::math::is_squarable(Vec2{1.0e154, 0.0}));

        // Finite, and not squarable. Both halves matter.
        const Vec2 huge{1.4e154, 0.0};
        MALLOY_CHECK_TRUE(malloy::math::is_finite(huge));
        MALLOY_CHECK_FALSE(malloy::math::is_squarable(huge));

        // A pair each of which is squarable alone but not together.
        const Vec2 both{1.0e154, 1.0e154};
        MALLOY_CHECK_TRUE(malloy::math::is_finite(both));
        MALLOY_CHECK_FALSE(malloy::math::is_squarable(both));

        MALLOY_CHECK_FALSE(malloy::math::is_squarable(Vec2{inf, 0.0}));
        MALLOY_CHECK_FALSE(malloy::math::is_squarable(Vec2{nan, 0.0}));

        // Underflow is not an error: a squared magnitude below the smallest
        // normal is zero, which is a real answer rather than a broken one.
        MALLOY_CHECK_TRUE(malloy::math::is_squarable(Vec2{1.0e-200, 0.0}));
    }

    // --- M19: Vec3. Only the parts that differ from Vec2 are worth dwelling
    //     on, and the cross product is all of them. ---
    {
        using malloy::math::Vec3;
        const Real inf = std::numeric_limits<Real>::infinity();
        const Real nan = std::numeric_limits<Real>::quiet_NaN();

        // Component arithmetic, with every component distinct so a dropped or
        // transposed one shows.
        const Vec3 a{1.0, 2.0, 3.0};
        const Vec3 b{-4.0, 5.0, -6.0};
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(a + b, Vec3{-3.0, 7.0, -3.0}, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(a - b, Vec3{5.0, -3.0, 9.0}, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(a * 2.0, Vec3{2.0, 4.0, 6.0}, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(-a, Vec3{-1.0, -2.0, -3.0}, eps));
        MALLOY_CHECK_NEAR(malloy::math::dot(a, b), -4.0 + 10.0 - 18.0, eps);

        // 3-4-5 in three dimensions, so length is exact rather than nearly so.
        MALLOY_CHECK_NEAR(malloy::math::length(Vec3{2.0, 3.0, 6.0}), 7.0, eps);
        MALLOY_CHECK_NEAR(malloy::math::distance(Vec3{1.0, 1.0, 1.0},
                                                 Vec3{3.0, 4.0, 7.0}), 7.0, eps);

        // --- The cross product, which is where 3D stops being 2D with an
        //     extra component. ---

        // Right-handed, checked on all three axis pairs. A sign error here is
        // invisible to every magnitude, so it needs its own assertion.
        const Vec3 ex{1.0, 0.0, 0.0};
        const Vec3 ey{0.0, 1.0, 0.0};
        const Vec3 ez{0.0, 0.0, 1.0};
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(malloy::math::cross(ex, ey), ez, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(malloy::math::cross(ey, ez), ex, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(malloy::math::cross(ez, ex), ey, eps));

        // Anticommutative, and zero on a parallel pair.
        const Vec3 c = malloy::math::cross(a, b);
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(malloy::math::cross(b, a), -c, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::cross(a, a * 3.0), Vec3{}, eps));

        // Perpendicular to both inputs. This is the property the whole of 3D
        // angular momentum rests on.
        MALLOY_CHECK_NEAR(malloy::math::dot(c, a), 0.0, 1e-12);
        MALLOY_CHECK_NEAR(malloy::math::dot(c, b), 0.0, 1e-12);

        // |a x b|^2 + (a . b)^2 = |a|^2 |b|^2, which pins the magnitude without
        // needing a sine.
        const Real dot_ab = malloy::math::dot(a, b);
        MALLOY_CHECK_NEAR(malloy::math::length_squared(c) + dot_ab * dot_ab,
                          malloy::math::length_squared(a) *
                              malloy::math::length_squared(b),
                          1e-12);

        // Worked by hand, off every axis, so no component can be dropped:
        // (1,2,3) x (-4,5,-6) = (2*-6 - 3*5, 3*-4 - 1*-6, 1*5 - 2*-4)
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(c, Vec3{-27.0, -6.0, 13.0}, eps));

        // --- The shared vocabulary means the same thing it does for Vec2. ---
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::normalize(Vec3{0.0, 0.0, -5.0}), Vec3{0.0, 0.0, -1.0}, eps));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            malloy::math::normalize(Vec3{}), Vec3{}, eps)); // no direction to give

        MALLOY_CHECK_TRUE(malloy::math::is_finite(a));
        MALLOY_CHECK_FALSE(malloy::math::is_finite(Vec3{0.0, inf, 0.0}));
        MALLOY_CHECK_FALSE(malloy::math::is_finite(Vec3{0.0, 0.0, nan}));

        MALLOY_CHECK_TRUE(malloy::math::is_squarable(a));
        MALLOY_CHECK_TRUE(malloy::math::is_squarable(Vec3{1.0e154, 0.0, 0.0}));
        // Finite in every component, and its square is not.
        const Vec3 huge{0.0, 0.0, 1.4e154};
        MALLOY_CHECK_TRUE(malloy::math::is_finite(huge));
        MALLOY_CHECK_FALSE(malloy::math::is_squarable(huge));
    }

    // --- M25: Mat3, and the diagonalization of a symmetric matrix into
    //     principal values and axes. This is the numerical core the 3D mass
    //     properties rest on, so it is tested against reconstruction and
    //     round-trips before anything is built on it. ---
    {
        using malloy::math::diagonal3;
        using malloy::math::eigen_symmetric;
        using malloy::math::identity3;
        using malloy::math::Mat3;
        using malloy::math::outer;
        using malloy::math::Quat;
        using malloy::math::SymmetricEigen;
        using malloy::math::to_quat;
        using malloy::math::transpose;
        using malloy::math::Vec3;

        // Basics: matrix times vector is a combination of the columns.
        {
            const Mat3 m{Vec3{1.0, 2.0, 3.0}, Vec3{4.0, 5.0, 6.0}, Vec3{7.0, 8.0, 9.0}};
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                m * Vec3{1.0, 0.0, 0.0}, Vec3{1.0, 2.0, 3.0}, 0.0)); // first column
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                m * Vec3{1.0, 1.0, 1.0}, Vec3{12.0, 15.0, 18.0}, 1e-14));
            // Transpose swaps rows and columns.
            MALLOY_CHECK_NEAR(malloy::math::element(transpose(m), 0, 1),
                              malloy::math::element(m, 1, 0), 0.0);
            // Outer product a b^T has (i,j) = a_i b_j.
            const Mat3 op = outer(Vec3{1.0, 2.0, 3.0}, Vec3{4.0, 5.0, 6.0});
            MALLOY_CHECK_NEAR(malloy::math::element(op, 0, 0), 4.0, 0.0);  // 1*4
            MALLOY_CHECK_NEAR(malloy::math::element(op, 2, 1), 15.0, 0.0); // 3*5
        }

        // A quaternion turned into its rotation matrix, columns = rotated basis.
        const auto matrix_of = [](const Quat& q) {
            return Mat3{malloy::math::rotate(q, Vec3{1.0, 0.0, 0.0}),
                        malloy::math::rotate(q, Vec3{0.0, 1.0, 0.0}),
                        malloy::math::rotate(q, Vec3{0.0, 0.0, 1.0})};
        };

        // Already diagonal, out of order: the values come back SORTED ascending,
        // and the eigenvectors are the canonical axes.
        {
            const SymmetricEigen e = eigen_symmetric(diagonal3(Vec3{3.0, 1.0, 2.0}));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(e.values, Vec3{1.0, 2.0, 3.0}, 1e-14));
            // Reconstruct: V diag(values) V^T is the original.
            const Mat3 rebuilt =
                e.vectors * diagonal3(e.values) * transpose(e.vectors);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                rebuilt, diagonal3(Vec3{3.0, 1.0, 2.0}), 1e-13));
            // The eigenvectors are a right-handed rotation (det +1).
            MALLOY_CHECK_NEAR(malloy::math::dot(malloy::math::cross(e.vectors.col0,
                                                                    e.vectors.col1),
                                                e.vectors.col2),
                              1.0, 1e-13);
        }

        // A full symmetric matrix with known eigenvalues, built as
        // R diag(lambda) R^T for a tilted R. Diagonalization must recover the
        // sorted eigenvalues and reconstruct the matrix, whatever basis it picks.
        {
            const Vec3 lambda{7.0, 1.0, 3.0};
            const Quat r = malloy::math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.9);
            const Mat3 rmat = matrix_of(r);
            const Mat3 symmetric = rmat * diagonal3(lambda) * transpose(rmat);

            const SymmetricEigen e = eigen_symmetric(symmetric);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(e.values, Vec3{1.0, 3.0, 7.0}, 1e-12));
            const Mat3 rebuilt = e.vectors * diagonal3(e.values) * transpose(e.vectors);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(rebuilt, symmetric, 1e-12));
            MALLOY_CHECK_NEAR(malloy::math::dot(malloy::math::cross(e.vectors.col0,
                                                                    e.vectors.col1),
                                                e.vectors.col2),
                              1.0, 1e-12);

            // to_quat turns the eigenvector rotation back into a quaternion that
            // rotates the basis onto those same columns.
            const Quat q = to_quat(e.vectors);
            MALLOY_CHECK_TRUE(malloy::math::is_unit(q));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::math::rotate(q, Vec3{1.0, 0.0, 0.0}), e.vectors.col0, 1e-13));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::math::rotate(q, Vec3{0.0, 0.0, 1.0}), e.vectors.col2, 1e-13));
        }

        // A repeated eigenvalue (an axisymmetric tensor): still reconstructs,
        // and the repeated value appears twice.
        {
            const Mat3 axisymmetric = diagonal3(Vec3{2.0, 2.0, 5.0});
            const SymmetricEigen e = eigen_symmetric(axisymmetric);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(e.values, Vec3{2.0, 2.0, 5.0}, 1e-14));
            const Mat3 rebuilt = e.vectors * diagonal3(e.values) * transpose(e.vectors);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(rebuilt, axisymmetric, 1e-14));
        }

        // A 180-degree rotation has a negative trace, so to_quat must take its
        // other branch; the well-conditioned root would be zero on the trace
        // branch. Checked by action, since q and -q are the same rotation.
        {
            const Vec3 axis = malloy::math::normalize(Vec3{1.0, 2.0, 2.0});
            const Quat half_turn = malloy::math::from_axis_angle(axis, 3.14159265358979323846);
            const Quat back = to_quat(matrix_of(half_turn));
            const Vec3 probe{0.7, -1.3, 0.4};
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::math::rotate(back, probe),
                malloy::math::rotate(half_turn, probe), 1e-12));
        }

        // to_quat of the identity is the identity rotation, and of a quarter
        // turn about z is that quarter turn.
        {
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(to_quat(identity3()), Quat{}, 1e-15));
            const Quat qz = malloy::math::from_axis_angle(Vec3{0.0, 0.0, 1.0}, 1.5707963267948966);
            const Quat back = to_quat(matrix_of(qz));
            // A quaternion and its negation are the same rotation, so compare by
            // action, not components.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                malloy::math::rotate(back, Vec3{1.0, 0.0, 0.0}), Vec3{0.0, 1.0, 0.0}, 1e-14));
        }
    }

    // --- M26: to_mat3, a unit quaternion as its rotation matrix. The columns
    //     must be the rotated basis vectors, so to_mat3(q) * v == rotate(q, v)
    //     for every v; that defining contract is checked column by column
    //     against `rotate` (M20), so the two cannot drift. It is the inverse of
    //     to_quat, and the piece that lets `combine` rebuild a part's inertia
    //     tensor from its stored orientation. ---
    {
        using malloy::math::from_axis_angle;
        using malloy::math::Mat3;
        using malloy::math::Quat;
        using malloy::math::rotate;
        using malloy::math::to_mat3;
        using malloy::math::to_quat;
        using malloy::math::Vec3;

        const Quat q = from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.9);
        const Mat3 m = to_mat3(q);
        // Each column is the corresponding basis vector, rotated: this pins the
        // columns independently, so a swapped or wrong one shows.
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            m.col0, rotate(q, Vec3{1.0, 0.0, 0.0}), 1e-15));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            m.col1, rotate(q, Vec3{0.0, 1.0, 0.0}), 1e-15));
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            m.col2, rotate(q, Vec3{0.0, 0.0, 1.0}), 1e-15));
        // And so, for a generic vector, to_mat3(q) * v agrees with rotate(q, v).
        const Vec3 probe{0.7, -1.3, 0.4};
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(m * probe, rotate(q, probe), 1e-14));
        // Inverse of to_quat: matrix -> quat -> action returns the same rotation.
        const Quat round = to_quat(m);
        MALLOY_CHECK_TRUE(malloy::math::approx_equal(
            rotate(round, probe), rotate(q, probe), 1e-13));
        // The identity quaternion is the identity matrix.
        MALLOY_CHECK_TRUE(
            malloy::math::approx_equal(to_mat3(Quat{}), malloy::math::identity3(), 1e-15));
    }

    // --- M27: trace, the sum of the diagonal. It is what turns a covariance
    //     tensor C = integral x x^T into an inertia tensor, I = trace(C) I - C,
    //     which is how the mesh mass properties get their inertia. ---
    {
        using malloy::math::diagonal3;
        using malloy::math::identity3;
        using malloy::math::Mat3;
        using malloy::math::trace;
        using malloy::math::Vec3;
        // Columns (1,2,3), (4,5,6), (7,8,9): the diagonal is 1, 5, 9.
        const Mat3 m{Vec3{1.0, 2.0, 3.0}, Vec3{4.0, 5.0, 6.0}, Vec3{7.0, 8.0, 9.0}};
        MALLOY_CHECK_NEAR(trace(m), 15.0, 0.0);
        MALLOY_CHECK_NEAR(trace(identity3()), 3.0, 0.0);
        MALLOY_CHECK_NEAR(trace(diagonal3(Vec3{2.0, -3.0, 7.0})), 6.0, 0.0);
    }

    std::cout << "malloy_math_tests passed\n";
    return 0;
}
