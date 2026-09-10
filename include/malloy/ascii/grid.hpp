#pragma once

#include <string>
#include <vector>

#include <malloy/math/vec2.hpp>

namespace malloy::ascii
{
// An axis-aligned 2D view rectangle in world coordinates.
struct Viewport
{
    math::Real min_x{-1.0};
    math::Real max_x{1.0};
    math::Real min_y{-1.0};
    math::Real max_y{1.0};
};

// A viewport that contains every finite point with a margin on each side.
// For empty/all-non-finite input, or a zero-width/height extent, it falls back
// to a non-degenerate rectangle so points still render at the center.
//
// It also falls back when the extent or the margin OVERFLOWS to infinity even
// though every input was finite, which happens near the limits of double. An
// infinite viewport would make render's normalized coordinates inf/inf = NaN.
// The consequence is worth knowing: the fallback box does not contain the
// points it was built from, so they are not drawn at all.
Viewport fit_viewport(const std::vector<math::Vec2>& points,
                      math::Real margin_fraction = 0.1);

// Render points into a width x height character grid framed by a border.
// A cell holding at least one point shows point_char, others empty_char. Row 0
// is the top (y = max_y); y increases upward. Points outside the viewport (or
// non-finite) are clipped. Returns the framed grid as a multi-line string.
std::string render(const std::vector<math::Vec2>& points, const Viewport& view,
                   int width, int height, char point_char = '*', char empty_char = ' ');
} // namespace malloy::ascii
