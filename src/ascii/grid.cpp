#include <malloy/ascii/grid.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <malloy/math/vec2.hpp>

namespace malloy::ascii
{
Viewport fit_viewport(const std::vector<math::Vec2>& points, math::Real margin_fraction)
{
    math::Real min_x = std::numeric_limits<math::Real>::infinity();
    math::Real max_x = -std::numeric_limits<math::Real>::infinity();
    math::Real min_y = std::numeric_limits<math::Real>::infinity();
    math::Real max_y = -std::numeric_limits<math::Real>::infinity();
    bool any = false;

    for (const math::Vec2& point : points)
    {
        if (!math::is_finite(point))
        {
            continue;
        }
        any = true;
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y);
        max_y = std::max(max_y, point.y);
    }

    if (!any)
    {
        return Viewport{};
    }

    // Expand a zero-size extent so a single point (or a colinear set) still has
    // a visible, non-degenerate view.
    if (max_x - min_x <= math::Real{0})
    {
        const math::Real center = (min_x + max_x) / math::Real{2};
        min_x = center - math::Real{1};
        max_x = center + math::Real{1};
    }
    if (max_y - min_y <= math::Real{0})
    {
        const math::Real center = (min_y + max_y) / math::Real{2};
        min_y = center - math::Real{1};
        max_y = center + math::Real{1};
    }

    const math::Real margin_x = (max_x - min_x) * margin_fraction;
    const math::Real margin_y = (max_y - min_y) * margin_fraction;
    return Viewport{min_x - margin_x, max_x + margin_x, min_y - margin_y, max_y + margin_y};
}

std::string render(const std::vector<math::Vec2>& points, const Viewport& view, int width,
                   int height, char point_char, char empty_char)
{
    if (width < 1)
    {
        width = 1;
    }
    if (height < 1)
    {
        height = 1;
    }

    std::vector<std::string> grid(static_cast<std::size_t>(height),
                                  std::string(static_cast<std::size_t>(width), empty_char));

    const math::Real extent_x = view.max_x - view.min_x;
    const math::Real extent_y = view.max_y - view.min_y;

    if (extent_x > math::Real{0} && extent_y > math::Real{0})
    {
        for (const math::Vec2& point : points)
        {
            if (!math::is_finite(point))
            {
                continue;
            }
            const math::Real tx = (point.x - view.min_x) / extent_x;
            const math::Real ty = (point.y - view.min_y) / extent_y;
            if (tx < math::Real{0} || tx > math::Real{1} || ty < math::Real{0} ||
                ty > math::Real{1})
            {
                continue; // outside the view; clip
            }
            const long col = std::lround(tx * (width - 1));
            const long row = std::lround((math::Real{1} - ty) * (height - 1)); // flip y
            grid[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] = point_char;
        }
    }

    const std::string border = "+" + std::string(static_cast<std::size_t>(width), '-') + "+";
    std::string out = border + "\n";
    for (const std::string& line : grid)
    {
        out += "|" + line + "|\n";
    }
    out += border + "\n";
    return out;
}
} // namespace malloy::ascii
