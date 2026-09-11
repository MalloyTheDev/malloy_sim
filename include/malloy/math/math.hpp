#pragma once

// Umbrella header for malloy_math.
//
// Vec2 is the 2D vocabulary every domain shipped so far is written in. Vec3
// arrived with M19 and is a separate concrete type rather than a template, for
// the reason spelled out in vec3.hpp: the cross product is a scalar in two
// dimensions and a vector in three, so the two are not one algebra with a
// different component count.
#include <malloy/math/real.hpp>
#include <malloy/math/vec2.hpp>
#include <malloy/math/vec3.hpp>
