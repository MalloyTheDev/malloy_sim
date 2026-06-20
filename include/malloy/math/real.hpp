#pragma once

#include <cmath>

namespace malloy::math
{
// The single scalar type for all simulation state. See
// docs/04_NUMERIC_AND_PHYSICS_CONVENTIONS.md: double only, never float.
using Real = double;

// True only for finite values; false for +/-infinity and NaN.
inline bool is_finite(Real value)
{
    return std::isfinite(value);
}

// Absolute-tolerance approximate equality. The tolerance is required and
// explicit on purpose: simulation reals must never be compared with ==.
inline bool approx_equal(Real a, Real b, Real epsilon)
{
    return std::abs(a - b) <= epsilon;
}
} // namespace malloy::math
