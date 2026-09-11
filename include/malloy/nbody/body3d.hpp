#pragma once

#include <malloy/math/vec3.hpp>

namespace malloy::nbody
{
// A point mass in three dimensions.
//
// This lives in `malloy_nbody` beside `Body2D` rather than in a library of its
// own, and that is the second question ADR 0009 left open, answered here.
//
// Gravity is ONE domain. Two and three dimensions are the same physics with a
// different component count, not two domains that happen to share a name, so
// splitting them across libraries would put the same equations in two places
// and then require them to be kept in step by hand. The unit that ADR 0006
// makes a library is the domain, and the dimension is not one.
//
// The 2D types are untouched, which answers the third open question too: 2D
// stays supported. Both worlds are selected by the scenario `type` key,
// dispatched by the same plain switch as every other domain.
struct Body3D
{
    math::Vec3 position{};
    math::Vec3 velocity{};
    math::Real mass{1.0};

    // Valid when mass is strictly positive and finite, and when position and
    // velocity are SQUARABLE rather than merely finite (docs/04).
    bool is_valid() const;
};
} // namespace malloy::nbody
