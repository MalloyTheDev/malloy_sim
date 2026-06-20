#pragma once

#include <cmath>

#include <malloy/math/real.hpp>

namespace malloy::math
{
// A plain 2D vector of Real. Concrete and non-templated by design
// (docs/04): no Vec3/Quat/templates in v1.
struct Vec2
{
    Real x{0.0};
    Real y{0.0};

    constexpr Vec2() = default;
    constexpr Vec2(Real x_value, Real y_value) : x{x_value}, y{y_value} {}

    constexpr Vec2& operator+=(const Vec2& other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    constexpr Vec2& operator-=(const Vec2& other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    constexpr Vec2& operator*=(Real scalar)
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    constexpr Vec2& operator/=(Real scalar)
    {
        x /= scalar;
        y /= scalar;
        return *this;
    }
};

constexpr Vec2 operator+(const Vec2& a, const Vec2& b)
{
    return Vec2{a.x + b.x, a.y + b.y};
}

constexpr Vec2 operator-(const Vec2& a, const Vec2& b)
{
    return Vec2{a.x - b.x, a.y - b.y};
}

constexpr Vec2 operator-(const Vec2& v)
{
    return Vec2{-v.x, -v.y};
}

constexpr Vec2 operator*(const Vec2& v, Real scalar)
{
    return Vec2{v.x * scalar, v.y * scalar};
}

constexpr Vec2 operator*(Real scalar, const Vec2& v)
{
    return Vec2{v.x * scalar, v.y * scalar};
}

constexpr Vec2 operator/(const Vec2& v, Real scalar)
{
    return Vec2{v.x / scalar, v.y / scalar};
}

constexpr Real dot(const Vec2& a, const Vec2& b)
{
    return a.x * b.x + a.y * b.y;
}

constexpr Real length_squared(const Vec2& v)
{
    return dot(v, v);
}

inline Real length(const Vec2& v)
{
    return std::sqrt(length_squared(v));
}

constexpr Real distance_squared(const Vec2& a, const Vec2& b)
{
    return length_squared(a - b);
}

inline Real distance(const Vec2& a, const Vec2& b)
{
    return std::sqrt(distance_squared(a, b));
}

// Unit vector in the direction of v. A zero vector has no direction, so it
// is returned unchanged instead of dividing by zero and yielding NaN/Inf.
inline Vec2 normalize(const Vec2& v)
{
    const Real len = length(v);
    if (len > Real{0})
    {
        return v / len;
    }
    return Vec2{};
}

inline bool is_finite(const Vec2& v)
{
    return is_finite(v.x) && is_finite(v.y);
}

// Component-wise approximate equality, reusing the scalar tolerance check
// so there is one authoritative definition of "approximately equal".
inline bool approx_equal(const Vec2& a, const Vec2& b, Real epsilon)
{
    return approx_equal(a.x, b.x, epsilon) && approx_equal(a.y, b.y, epsilon);
}
} // namespace malloy::math
