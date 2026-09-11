#pragma once

#include <cmath>

#include <malloy/math/real.hpp>

namespace malloy::math
{
// A plain 3D vector of Real, concrete and non-templated, for the same reason
// Vec2 is: this project prefers a concrete type over a general one until the
// general one earns itself (rule 11).
//
// It is a SEPARATE type rather than `Vector<N>` with N = 3, and that is a
// decision rather than an omission. ADR 0009 left the choice open; it is
// settled here because building it made the answer obvious.
//
// The two are not the same algebra with a different count. `cross` returns a
// SCALAR in two dimensions and a VECTOR in three, so a template would need a
// specialization for the one operation that matters most to rigid-body
// dynamics. `perp` and rotation by a scalar angle exist only in 2D. A template
// would unify the easy half, the component-wise arithmetic, and then need
// specializing for every part that is actually interesting.
//
// What is shared is spelled the same way on purpose: dot, length,
// length_squared, distance, normalize, is_finite, is_squarable and
// approx_equal all mean here exactly what they mean for Vec2, so code that
// reads one reads the other.
struct Vec3
{
    Real x{0.0};
    Real y{0.0};
    Real z{0.0};

    constexpr Vec3() = default;
    constexpr Vec3(Real x_value, Real y_value, Real z_value)
        : x{x_value}, y{y_value}, z{z_value}
    {
    }

    constexpr Vec3& operator+=(const Vec3& other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    constexpr Vec3& operator-=(const Vec3& other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    constexpr Vec3& operator*=(Real scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    constexpr Vec3& operator/=(Real scalar)
    {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }
};

constexpr Vec3 operator+(const Vec3& a, const Vec3& b)
{
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

constexpr Vec3 operator-(const Vec3& a, const Vec3& b)
{
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

constexpr Vec3 operator-(const Vec3& v)
{
    return Vec3{-v.x, -v.y, -v.z};
}

constexpr Vec3 operator*(const Vec3& v, Real scalar)
{
    return Vec3{v.x * scalar, v.y * scalar, v.z * scalar};
}

constexpr Vec3 operator*(Real scalar, const Vec3& v)
{
    return v * scalar;
}

constexpr Vec3 operator/(const Vec3& v, Real scalar)
{
    return Vec3{v.x / scalar, v.y / scalar, v.z / scalar};
}

constexpr Real dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// The cross product, which is where 3D stops being 2D with an extra component.
//
// In two dimensions this is a scalar, the signed area, and the project keeps it
// privately in the domains that need it. In three it is a vector perpendicular
// to both inputs, and it is the reason angular momentum and torque become
// vectors rather than scalars.
constexpr Vec3 cross(const Vec3& a, const Vec3& b)
{
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x};
}

constexpr Real length_squared(const Vec3& v)
{
    return dot(v, v);
}

inline Real length(const Vec3& v)
{
    return std::sqrt(length_squared(v));
}

constexpr Real distance_squared(const Vec3& a, const Vec3& b)
{
    return length_squared(a - b);
}

inline Real distance(const Vec3& a, const Vec3& b)
{
    return std::sqrt(distance_squared(a, b));
}

// The zero vector has no direction, so it is returned unchanged rather than
// producing NaN. Same rule as Vec2.
inline Vec3 normalize(const Vec3& v)
{
    const Real len = length(v);
    if (!(len > Real{0}) || !is_finite(len))
    {
        return Vec3{};
    }
    return v / len;
}

inline bool is_finite(const Vec3& v)
{
    return is_finite(v.x) && is_finite(v.y) && is_finite(v.z);
}

// True when |v|^2 is representable, not merely when v is. See the Vec2 version:
// everything that squares a vector reaches infinity well below DBL_MAX, and
// validating only finiteness leaves a window in which every reported energy is
// infinite while the state looks sound.
inline bool is_squarable(const Vec3& v)
{
    return is_finite(length_squared(v));
}

inline bool approx_equal(const Vec3& a, const Vec3& b, Real epsilon)
{
    return approx_equal(a.x, b.x, epsilon) && approx_equal(a.y, b.y, epsilon) &&
           approx_equal(a.z, b.z, epsilon);
}
} // namespace malloy::math
