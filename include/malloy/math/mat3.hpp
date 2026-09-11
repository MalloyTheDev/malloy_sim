#pragma once

#include <cmath>

#include <malloy/math/quat.hpp>
#include <malloy/math/real.hpp>
#include <malloy/math/vec3.hpp>

namespace malloy::math
{
// A 3x3 matrix, stored as its three COLUMNS. Added in M25 for one job: an
// inertia tensor and its diagonalization into principal moments. It is a
// concrete type like Vec3 and Quat, not a general linear-algebra layer, and it
// carries only the operations mass properties needs.
//
// Column storage makes `Mat3 * Vec3` read as a combination of the columns, and
// makes a rotation's columns its basis vectors, which is exactly what
// diagonalization returns.
struct Mat3
{
    Vec3 col0{};
    Vec3 col1{};
    Vec3 col2{};

    constexpr Mat3() = default;
    constexpr Mat3(const Vec3& c0, const Vec3& c1, const Vec3& c2)
        : col0{c0}, col1{c1}, col2{c2}
    {
    }
};

// Element at row r, column c (both 0..2). Reads, does not write.
inline Real element(const Mat3& m, int r, int c)
{
    const Vec3& column = (c == 0) ? m.col0 : (c == 1) ? m.col1 : m.col2;
    return (r == 0) ? column.x : (r == 1) ? column.y : column.z;
}

inline Mat3 identity3()
{
    return Mat3{Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
}

inline Mat3 diagonal3(const Vec3& d)
{
    return Mat3{Vec3{d.x, 0.0, 0.0}, Vec3{0.0, d.y, 0.0}, Vec3{0.0, 0.0, d.z}};
}

// Matrix times vector: a combination of the columns.
inline Vec3 operator*(const Mat3& m, const Vec3& v)
{
    return m.col0 * v.x + m.col1 * v.y + m.col2 * v.z;
}

inline Mat3 operator*(const Mat3& a, const Mat3& b)
{
    return Mat3{a * b.col0, a * b.col1, a * b.col2};
}

inline Mat3 operator+(const Mat3& a, const Mat3& b)
{
    return Mat3{a.col0 + b.col0, a.col1 + b.col1, a.col2 + b.col2};
}

inline Mat3 operator-(const Mat3& a, const Mat3& b)
{
    return Mat3{a.col0 - b.col0, a.col1 - b.col1, a.col2 - b.col2};
}

inline Mat3 operator*(const Mat3& m, Real s)
{
    return Mat3{m.col0 * s, m.col1 * s, m.col2 * s};
}

inline Mat3 transpose(const Mat3& m)
{
    return Mat3{Vec3{m.col0.x, m.col1.x, m.col2.x},
                Vec3{m.col0.y, m.col1.y, m.col2.y},
                Vec3{m.col0.z, m.col1.z, m.col2.z}};
}

// Outer product a b^T, a symmetric-looking 3x3 whose (i, j) entry is a_i b_j.
// The piece the parallel-axis theorem subtracts: m (|d|^2 I - d d^T).
inline Mat3 outer(const Vec3& a, const Vec3& b)
{
    // (a b^T)(i, j) = a_i b_j, so column j is a scaled by b_j.
    return Mat3{a * b.x, a * b.y, a * b.z};
}

inline bool is_finite(const Mat3& m)
{
    return is_finite(m.col0) && is_finite(m.col1) && is_finite(m.col2);
}

inline bool approx_equal(const Mat3& a, const Mat3& b, Real epsilon)
{
    return approx_equal(a.col0, b.col0, epsilon) &&
           approx_equal(a.col1, b.col1, epsilon) &&
           approx_equal(a.col2, b.col2, epsilon);
}

// A unit quaternion as its rotation matrix: the columns are the rotated basis
// vectors, so `to_mat3(q) * v` equals `rotate(q, v)`. The inverse of `to_quat`
// for a proper rotation. Used to turn a stored orientation back into the tensor
// R diag(moments) R^T when combining mass properties.
inline Mat3 to_mat3(const Quat& q)
{
    return Mat3{rotate(q, Vec3{1.0, 0.0, 0.0}), rotate(q, Vec3{0.0, 1.0, 0.0}),
                rotate(q, Vec3{0.0, 0.0, 1.0})};
}

// The eigen-decomposition of a SYMMETRIC 3x3: three real eigenvalues (the
// principal values) and an orthonormal set of eigenvectors (the principal
// axes), returned as the columns of a proper rotation. Reconstructing the
// input is `vectors * diagonal3(values) * transpose(vectors)`.
//
// `values` is sorted ascending and `vectors` is a right-handed rotation
// (determinant +1), so the output is canonical: the same symmetric matrix
// always gives the same decomposition, which the tests and the physics both
// rely on.
struct SymmetricEigen
{
    Vec3 values{};
    Mat3 vectors{};
};

// Classical cyclic Jacobi rotation, the standard method for a small symmetric
// matrix: it repeatedly zeroes the largest off-diagonal entry with an orthogonal
// rotation, and for a 3x3 converges to machine precision in a handful of sweeps.
// The input is assumed symmetric; only its lower triangle would differ, and it
// is not read.
inline SymmetricEigen eigen_symmetric(const Mat3& matrix)
{
    // Work in a plain array for the rotations; the result goes back into Mat3.
    Real a[3][3];
    Real v[3][3];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            a[i][j] = element(matrix, i, j);
            v[i][j] = (i == j) ? Real{1} : Real{0};
        }
    }

    // A generous cap: 3x3 Jacobi needs only a few sweeps, but the loop exits
    // early once the off-diagonal is negligible, so the cap is just a guard.
    for (int sweep = 0; sweep < 50; ++sweep)
    {
        const Real off =
            std::abs(a[0][1]) + std::abs(a[0][2]) + std::abs(a[1][2]);
        if (!(off > Real{0}))
        {
            break; // already diagonal (also the exit for a zero matrix)
        }

        // Rotate the three (p, q) pairs in turn.
        for (int p = 0; p < 2; ++p)
        {
            for (int q = p + 1; q < 3; ++q)
            {
                if (!(std::abs(a[p][q]) > Real{0}))
                {
                    continue;
                }
                // Angle that zeroes a[p][q]: cot(2 theta) = (a_qq - a_pp)/(2 a_pq).
                const Real difference = a[q][q] - a[p][p];
                Real t; // tan(theta)
                if (std::abs(a[p][q]) < std::abs(difference) * Real{1e-38})
                {
                    t = a[p][q] / difference;
                }
                else
                {
                    const Real theta = difference / (Real{2} * a[p][q]);
                    t = (theta >= Real{0} ? Real{1} : Real{-1}) /
                        (std::abs(theta) + std::sqrt(theta * theta + Real{1}));
                }
                const Real c = Real{1} / std::sqrt(t * t + Real{1});
                const Real s = t * c;

                // Apply the rotation to A on both sides (it stays symmetric).
                const Real app = a[p][p];
                const Real aqq = a[q][q];
                const Real apq = a[p][q];
                a[p][p] = c * c * app - Real{2} * s * c * apq + s * s * aqq;
                a[q][q] = s * s * app + Real{2} * s * c * apq + c * c * aqq;
                a[p][q] = Real{0};
                a[q][p] = Real{0};
                const int r = 3 - p - q; // the third index
                const Real arp = a[r][p];
                const Real arq = a[r][q];
                a[r][p] = c * arp - s * arq;
                a[p][r] = a[r][p];
                a[r][q] = s * arp + c * arq;
                a[q][r] = a[r][q];

                // Accumulate the rotation into the eigenvector columns.
                for (int i = 0; i < 3; ++i)
                {
                    const Real vip = v[i][p];
                    const Real viq = v[i][q];
                    v[i][p] = c * vip - s * viq;
                    v[i][q] = s * vip + c * viq;
                }
            }
        }
    }

    // Collect (eigenvalue, eigenvector column) and sort ascending by value.
    int order[3] = {0, 1, 2};
    for (int i = 0; i < 2; ++i)
    {
        for (int j = i + 1; j < 3; ++j)
        {
            if (a[order[j]][order[j]] < a[order[i]][order[i]])
            {
                const int tmp = order[i];
                order[i] = order[j];
                order[j] = tmp;
            }
        }
    }

    const auto column = [&](int k) {
        return Vec3{v[0][k], v[1][k], v[2][k]};
    };
    SymmetricEigen result;
    result.values = Vec3{a[order[0]][order[0]], a[order[1]][order[1]],
                         a[order[2]][order[2]]};
    Vec3 e0 = column(order[0]);
    Vec3 e1 = column(order[1]);
    Vec3 e2 = column(order[2]);
    // Make it a RIGHT-handed rotation (determinant +1), so it converts cleanly
    // to a quaternion: flip the third axis if the triple is left-handed.
    if (dot(cross(e0, e1), e2) < Real{0})
    {
        e2 = -e2;
    }
    result.vectors = Mat3{e0, e1, e2};
    return result;
}

// A proper rotation matrix (orthonormal columns, determinant +1) as a unit
// quaternion. Shepperd's method: pick the largest of the four components from
// the trace and diagonal to stay well-conditioned, then read the rest off.
inline Quat to_quat(const Mat3& r)
{
    const Real m00 = element(r, 0, 0);
    const Real m11 = element(r, 1, 1);
    const Real m22 = element(r, 2, 2);
    const Real trace = m00 + m11 + m22;
    Quat q;
    if (trace > Real{0})
    {
        Real root = std::sqrt(trace + Real{1});
        q.w = Real{0.5} * root;
        root = Real{0.5} / root;
        q.v.x = (element(r, 2, 1) - element(r, 1, 2)) * root;
        q.v.y = (element(r, 0, 2) - element(r, 2, 0)) * root;
        q.v.z = (element(r, 1, 0) - element(r, 0, 1)) * root;
    }
    else
    {
        // Largest diagonal element decides which axis component leads.
        int i = 0;
        if (m11 > m00) i = 1;
        if (m22 > element(r, i, i)) i = 2;
        const int j = (i + 1) % 3;
        const int k = (i + 2) % 3;
        Real root =
            std::sqrt(element(r, i, i) - element(r, j, j) - element(r, k, k) + Real{1});
        Real axis[3] = {Real{0}, Real{0}, Real{0}};
        axis[i] = Real{0.5} * root;
        root = Real{0.5} / root;
        q.w = (element(r, k, j) - element(r, j, k)) * root;
        axis[j] = (element(r, j, i) + element(r, i, j)) * root;
        axis[k] = (element(r, k, i) + element(r, i, k)) * root;
        q.v = Vec3{axis[0], axis[1], axis[2]};
    }
    return normalize(q);
}
} // namespace malloy::math
