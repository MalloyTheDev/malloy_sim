#include "self_test.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace viewer
{
namespace view3d = malloy::view3d;
using view3d::Real;
using view3d::Vec3;

namespace
{
// What the checks know of one subject: its bounding sphere, and the radius of
// a sphere it certainly covers (for a box, its smallest half-width).
struct Subject
{
    Vec3 center;
    Real outer{0};
    Real inner{0};
};

std::vector<Subject> subjects_of(const view3d::Frame& frame, const Look& look)
{
    std::vector<Subject> out;
    for (const view3d::Sphere& s : frame.spheres)
    {
        out.push_back(Subject{s.center, s.radius, s.radius});
    }
    for (const view3d::Box& b : frame.boxes)
    {
        const Vec3& h = b.half_extents;
        out.push_back(Subject{b.center, malloy::math::length(h), std::min({h.x, h.y, h.z})});
    }
    for (const view3d::Point& p : frame.points)
    {
        const Real r = point_radius(look, p.position);
        if (r > 0)
        {
            out.push_back(Subject{p.position, r, r});
        }
    }
    return out;
}

bool is_background(const std::vector<std::uint8_t>& rgba, int width, int x, int y, Rgb bg)
{
    const std::size_t i = 4U * (static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                                static_cast<std::size_t>(x));
    return rgba[i] == bg.r && rgba[i + 1] == bg.g && rgba[i + 2] == bg.b;
}
} // namespace

PixelReport check_pixels(const std::vector<std::uint8_t>& rgba, int width, int height,
                         Rgb background, const view3d::Frame& frame, const Look& look)
{
    PixelReport report;
    std::ostringstream why;
    if (width <= 0 || height <= 0 ||
        rgba.size() != 4U * static_cast<std::size_t>(width) * static_cast<std::size_t>(height))
    {
        report.failure = "no pixels were read back";
        return report;
    }

    // (a) Not blank.
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            report.drawn_pixels += is_background(rgba, width, x, y, background) ? 0U : 1U;
        }
    }
    if (report.drawn_pixels == 0)
    {
        report.failure = "the picture is blank";
        return report;
    }

    const view3d::Camera& cam = look.camera;
    const Real t = std::tan(cam.fovy * 0.5);
    const Real pixels_per_unit = 0.5 * look.viewport_height / t; // at depth 1
    const std::vector<Subject> subjects = subjects_of(frame, look);

    // (b) Every subject at least 2 px across is drawn where the model puts it.
    for (const Subject& s : subjects)
    {
        const Vec3 ndc = view3d::project(cam, look.aspect, s.center);
        if (!(ndc.z > 0))
        {
            continue;
        }
        if (s.inner / ndc.z * pixels_per_unit < 1.0)
        {
            continue; // under 2 px across: too small to probe
        }
        const Real px = (ndc.x + 1.0) * 0.5 * width;
        const Real py = (1.0 - ndc.y) * 0.5 * height;
        const int x = static_cast<int>(std::floor(px));
        const int y = static_cast<int>(std::floor(py));
        ++report.subject_probes;
        if (x < 0 || y < 0 || x >= width || y >= height)
        {
            why << "a subject projects off screen at (" << px << ", " << py << ")";
            report.failure = why.str();
            return report;
        }
        if (is_background(rgba, width, x, y, background))
        {
            why << "the subject at pixel (" << x << ", " << y << ") is not drawn there";
            report.failure = why.str();
            return report;
        }
    }

    // (c) Clear pixels are background. A pixel is clear of a subject when its
    // view ray misses the subject's bounding sphere by an angle of more than
    // 3 px + 0.1 of the sphere's radius, measured at the centre of the view
    // (where a pixel subtends the largest angle, so the margin only grows
    // toward the edges). The ray is built in the model's own camera convention
    // (view3d::project's: looking down -z, x right, y up), which the drawing
    // must agree with for the check to pass.
    const Vec3 f = malloy::math::normalize(cam.target - cam.position);
    const Vec3 side = malloy::math::normalize(malloy::math::cross(f, cam.up));
    const Vec3 u = malloy::math::cross(side, f);
    const Real pixel_angle = 2.0 * t / look.viewport_height;
    for (int y = 8; y < height; y += 16)
    {
        for (int x = 8; x < width; x += 16)
        {
            const Real nx = (2.0 * (x + 0.5) / width) - 1.0;
            const Real ny = 1.0 - (2.0 * (y + 0.5) / height);
            const Vec3 ray = malloy::math::normalize(f + side * (nx * look.aspect * t) + u * (ny * t));
            bool clear = true;
            for (const Subject& s : subjects)
            {
                const Vec3 to = s.center - cam.position;
                const Real d = malloy::math::length(to);
                if (!(d > s.outer))
                {
                    clear = false; // the eye is inside it
                    break;
                }
                const Real cosine = std::clamp(malloy::math::dot(ray, to) / d, -1.0, 1.0);
                const Real angle = std::acos(cosine);
                const Real radius = std::asin(s.outer / d);
                const Real radius_pixels = radius / pixel_angle;
                if (angle <= radius + (3.0 + 0.1 * radius_pixels) * pixel_angle)
                {
                    clear = false;
                    break;
                }
            }
            if (!clear)
            {
                continue;
            }
            ++report.clear_probes;
            if (!is_background(rgba, width, x, y, background))
            {
                why << "pixel (" << x << ", " << y
                    << ") is drawn where the model has nothing";
                report.failure = why.str();
                return report;
            }
        }
    }
    return report;
}
} // namespace viewer
