#include <malloy/ascii/grid.hpp>

#include <test_check.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

using malloy::ascii::fit_viewport;
using malloy::ascii::render;
using malloy::ascii::Viewport;
using malloy::math::Real;
using malloy::math::Vec2;

int main()
{
    const Real eps = 1e-12;

    // A single centered point renders at the middle of a 3x3 grid.
    {
        const Viewport view{0.0, 2.0, 0.0, 2.0};
        const std::string grid = render({Vec2{1.0, 1.0}}, view, 3, 3);
        const std::string expected = "+---+\n"
                                     "|   |\n"
                                     "| * |\n"
                                     "|   |\n"
                                     "+---+\n";
        MALLOY_CHECK_TRUE(grid == expected);
    }

    // y points up: max_y maps to the top row, min_x to the left column.
    {
        const Viewport view{0.0, 2.0, 0.0, 2.0};
        const std::string grid = render({Vec2{0.0, 2.0}, Vec2{2.0, 0.0}}, view, 3, 3);
        const std::string expected = "+---+\n"
                                     "|*  |\n"
                                     "|   |\n"
                                     "|  *|\n"
                                     "+---+\n";
        MALLOY_CHECK_TRUE(grid == expected);
    }

    // Points outside the viewport are clipped, leaving an empty grid.
    {
        const Viewport view{0.0, 1.0, 0.0, 1.0};
        const std::string grid = render({Vec2{5.0, 5.0}, Vec2{-3.0, 0.5}}, view, 3, 3);
        const std::string expected = "+---+\n"
                                     "|   |\n"
                                     "|   |\n"
                                     "|   |\n"
                                     "+---+\n";
        MALLOY_CHECK_TRUE(grid == expected);
    }

    // Empty input renders an empty framed grid.
    {
        const std::string grid = render({}, Viewport{}, 2, 2);
        const std::string expected = "+--+\n"
                                     "|  |\n"
                                     "|  |\n"
                                     "+--+\n";
        MALLOY_CHECK_TRUE(grid == expected);
    }

    // fit_viewport contains the points with margin.
    {
        const Viewport view = fit_viewport({Vec2{-1.0, -2.0}, Vec2{3.0, 4.0}}, 0.1);
        MALLOY_CHECK_NEAR(view.min_x, -1.4, eps); // extent 4 -> margin 0.4
        MALLOY_CHECK_NEAR(view.max_x, 3.4, eps);
        MALLOY_CHECK_NEAR(view.min_y, -2.6, eps); // extent 6 -> margin 0.6
        MALLOY_CHECK_NEAR(view.max_y, 4.6, eps);
    }

    // fit_viewport on a single point yields a non-degenerate view.
    {
        const Viewport view = fit_viewport({Vec2{5.0, 5.0}}, 0.1);
        MALLOY_CHECK_TRUE(view.max_x > view.min_x);
        MALLOY_CHECK_TRUE(view.max_y > view.min_y);
    }

    // Non-finite points are skipped rather than indexing outside the grid.
    {
        const Viewport view{0.0, 2.0, 0.0, 2.0};
        const Real nan_value = std::numeric_limits<Real>::quiet_NaN();
        const Real inf_value = std::numeric_limits<Real>::infinity();
        const std::string grid =
            render({Vec2{nan_value, 1.0}, Vec2{1.0, inf_value}}, view, 3, 3);
        const std::string expected = "+---+\n"
                                     "|   |\n"
                                     "|   |\n"
                                     "|   |\n"
                                     "+---+\n";
        MALLOY_CHECK_TRUE(grid == expected);
    }

    // A zero-extent viewport renders empty instead of dividing by zero.
    {
        const Viewport view{1.0, 1.0, 1.0, 1.0};
        const std::string grid = render({Vec2{1.0, 1.0}}, view, 2, 2);
        const std::string expected = "+--+\n"
                                     "|  |\n"
                                     "|  |\n"
                                     "+--+\n";
        MALLOY_CHECK_TRUE(grid == expected);
    }

    // Non-positive grid sizes are clamped to 1x1 rather than producing a
    // negatively sized grid.
    {
        const std::string grid = render({Vec2{0.0, 0.0}}, Viewport{}, 0, -4);
        const std::string expected = "+-+\n"
                                     "|*|\n"
                                     "+-+\n";
        MALLOY_CHECK_TRUE(grid == expected);
    }

    // Custom point/empty characters fill the cells and the background.
    {
        const Viewport view{0.0, 2.0, 0.0, 2.0};
        const std::string grid = render({Vec2{1.0, 1.0}}, view, 3, 3, 'O', '.');
        const std::string expected = "+---+\n"
                                     "|...|\n"
                                     "|.O.|\n"
                                     "|...|\n"
                                     "+---+\n";
        MALLOY_CHECK_TRUE(grid == expected);
    }

    // fit_viewport ignores non-finite points and fits the finite ones.
    {
        const Real nan_value = std::numeric_limits<Real>::quiet_NaN();
        const Viewport view = fit_viewport(
            {Vec2{nan_value, nan_value}, Vec2{0.0, 0.0}, Vec2{10.0, 10.0}}, 0.0);
        MALLOY_CHECK_NEAR(view.min_x, 0.0, eps);
        MALLOY_CHECK_NEAR(view.max_x, 10.0, eps);
        MALLOY_CHECK_NEAR(view.min_y, 0.0, eps);
        MALLOY_CHECK_NEAR(view.max_y, 10.0, eps);
    }

    // All-non-finite input falls back to the default viewport.
    {
        const Real inf_value = std::numeric_limits<Real>::infinity();
        const Viewport view = fit_viewport({Vec2{inf_value, 0.0}, Vec2{0.0, -inf_value}});
        MALLOY_CHECK_NEAR(view.min_x, -1.0, eps);
        MALLOY_CHECK_NEAR(view.max_x, 1.0, eps);
        MALLOY_CHECK_NEAR(view.min_y, -1.0, eps);
        MALLOY_CHECK_NEAR(view.max_y, 1.0, eps);
    }

    // --- Issue #5: a viewport that overflowed to infinity used to make the
    //     normalized coordinates inf/inf = NaN. Every comparison against NaN is
    //     false, so the point escaped the clip test and reached std::lround,
    //     whose NaN result is unspecified: 0 on MSVC, LONG_MIN on GCC, which
    //     then casts to a huge index. ---
    {
        const Real infinity_value = std::numeric_limits<Real>::infinity();
        const Viewport infinite{-infinity_value, infinity_value, -1.0, 1.0};
        const std::string grid = render({Vec2{0.0, 0.0}}, infinite, 3, 3);
        const std::string expected = "+---+\n"
                                     "|   |\n"
                                     "|   |\n"
                                     "|   |\n"
                                     "+---+\n";
        MALLOY_CHECK_TRUE(grid == expected); // skipped, not drawn at column 0
    }

    // --- Issue #5: fit_viewport never returns a non-finite rectangle, even
    //     when every input point is finite. Coordinates near the limits of
    //     double make the extent, and then the margin, overflow. ---
    {
        const Viewport view = fit_viewport({Vec2{-1e308, 0.0}, Vec2{1e308, 0.0}}, 0.1);
        MALLOY_CHECK_TRUE(std::isfinite(view.min_x));
        MALLOY_CHECK_TRUE(std::isfinite(view.max_x));
        MALLOY_CHECK_TRUE(std::isfinite(view.min_y));
        MALLOY_CHECK_TRUE(std::isfinite(view.max_y));

        // And the pair renders without indexing out of range.
        const std::string grid =
            render({Vec2{-1e308, 0.0}, Vec2{1e308, 0.0}}, view, 5, 3);
        MALLOY_CHECK_TRUE(grid.size() > 0);
    }

    std::cout << "malloy_ascii_tests passed\n";
    return 0;
}
