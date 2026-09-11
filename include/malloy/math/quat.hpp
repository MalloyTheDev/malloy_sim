#pragma once

#include <cmath>

#include <malloy/math/real.hpp>
#include <malloy/math/vec3.hpp>

namespace malloy::math
{
// A quaternion, used here only ever as a UNIT quaternion representing an
// orientation.
//
// This is what a scalar angle becomes in three dimensions, and ADR 0007 said so
// when it chose the scalar: the 2D form is this one with two components forced
// to zero, not a different model.
//
// Stored as a scalar part w and a vector part v, which is the spelling that
// makes the algebra readable: the product is
//
//     (w1 w2 - v1 . v2,  w1 v2 + w2 v1 + v1 x v2)
//
// and the cross product in that last term is exactly why composing rotations
// does not commute in 3D while it does in 2D.
struct Quat
{
    Real w{1.0}; // the identity rotation, not the zero quaternion
    Vec3 v{};

    constexpr Quat() = default;
    constexpr Quat(Real w_value, const Vec3& v_value) : w{w_value}, v{v_value} {}
    constexpr Quat(Real w_value, Real x, Real y, Real z)
        : w{w_value}, v{Vec3{x, y, z}}
    {
    }
};

// Hamilton product. NOT commutative, which is the whole content of 3D rotation:
// turning right then up does not land where turning up then right does.
constexpr Quat operator*(const Quat& a, const Quat& b)
{
    return Quat{a.w * b.w - dot(a.v, b.v),
                b.v * a.w + a.v * b.w + cross(a.v, b.v)};
}

constexpr Quat operator+(const Quat& a, const Quat& b)
{
    return Quat{a.w + b.w, a.v + b.v};
}

constexpr Quat operator*(const Quat& q, Real scalar)
{
    return Quat{q.w * scalar, q.v * scalar};
}

constexpr Quat operator*(Real scalar, const Quat& q)
{
    return q * scalar;
}

// The inverse of a UNIT quaternion. For a non-unit one this is the conjugate
// rather than the inverse, hence the name.
constexpr Quat conjugate(const Quat& q)
{
    return Quat{q.w, -q.v};
}

constexpr Real norm_squared(const Quat& q)
{
    return q.w * q.w + dot(q.v, q.v);
}

inline Real norm(const Quat& q)
{
    return std::sqrt(norm_squared(q));
}

// Returns the identity for a quaternion with no norm, matching the rule Vec2
// and Vec3 use for normalize: a thing with no magnitude has no direction, so
// invent nothing.
inline Quat normalize(const Quat& q)
{
    const Real n = norm(q);
    if (!(n > Real{0}) || !is_finite(n))
    {
        return Quat{};
    }
    return q * (Real{1} / n);
}

// Rotate a vector by a unit quaternion: q v q*.
//
// Written in the expanded form rather than as two quaternion products, because
// it is both faster and, more importantly, exact for the identity: with
// v_q = 0 the terms vanish cleanly instead of accumulating rounding through
// two multiplications.
inline Vec3 rotate(const Quat& q, const Vec3& p)
{
    const Vec3 t = cross(q.v, p) * Real{2};
    return p + t * q.w + cross(q.v, t);
}

// A rotation of `angle` radians about `axis`, right-handed.
//
// The axis is normalized here, so a caller need not, and a zero axis gives the
// identity rather than NaN.
inline Quat from_axis_angle(const Vec3& axis, Real angle)
{
    const Vec3 unit = normalize(axis);
    if (!(length_squared(unit) > Real{0}))
    {
        return Quat{};
    }
    const Real half = angle * Real{0.5};
    return Quat{std::cos(half), unit * std::sin(half)};
}

inline bool is_finite(const Quat& q)
{
    return is_finite(q.w) && is_finite(q.v);
}

// Valid as an ORIENTATION, which is stricter than being finite: it must be a
// unit quaternion, because every rotation formula above assumes it.
//
// The tolerance is on the SQUARE of the norm, so it needs no square root, and
// it is loose enough for a hand-written orientation whose components are not
// exact in binary.
inline bool is_unit(const Quat& q)
{
    return is_finite(q) && std::abs(norm_squared(q) - Real{1}) <= Real{1e-12};
}

inline bool approx_equal(const Quat& a, const Quat& b, Real epsilon)
{
    return approx_equal(a.w, b.w, epsilon) && approx_equal(a.v, b.v, epsilon);
}
} // namespace malloy::math
