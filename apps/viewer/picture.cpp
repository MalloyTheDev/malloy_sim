#include "picture.hpp"

#include <algorithm>
#include <cmath>

namespace viewer
{
namespace view3d = malloy::view3d;
using view3d::Real;
using view3d::Role;
using view3d::Vec3;

namespace
{
constexpr Real kAmbient = 0.3;
constexpr Real kPi = 3.14159265358979323846;

constexpr Rgb kBody{232, 142, 62};
constexpr Rgb kProxy{150, 132, 210};
constexpr Rgb kMass{240, 214, 120};
constexpr Rgb kPositive{232, 72, 60};
constexpr Rgb kNegative{72, 124, 236};
constexpr Rgb kNeutral{170, 172, 180};
constexpr Rgb kNode{190, 200, 216};
constexpr Rgb kGround{78, 84, 96};
constexpr Rgb kGrid{112, 118, 130};
constexpr Rgb kWall{150, 154, 166};
constexpr Rgb kAxisX{226, 70, 70};
constexpr Rgb kAxisY{80, 200, 96};
constexpr Rgb kAxisZ{82, 132, 240};
constexpr Rgb kElectric{240, 208, 60};
constexpr Rgb kMagnetic{70, 206, 232};
constexpr Rgb kSpringRest{170, 172, 180};

Rgb scaled(const Rgb& c, Real k)
{
    const auto channel = [k](std::uint8_t v) {
        return static_cast<std::uint8_t>(std::clamp(std::lround(v * k), 0L, 255L));
    };
    return Rgb{channel(c.r), channel(c.g), channel(c.b)};
}

Rgb mix(const Rgb& a, const Rgb& b, Real t)
{
    const auto channel = [t](std::uint8_t x, std::uint8_t y) {
        return static_cast<std::uint8_t>(std::clamp(std::lround(x + (y - x) * t), 0L, 255L));
    };
    return Rgb{channel(a.r, b.r), channel(a.g, b.g), channel(a.b, b.b)};
}

Rgb point_color(const view3d::Point& p)
{
    switch (p.role)
    {
    case Role::Charge:
        return p.weight > 0 ? kPositive : (p.weight < 0 ? kNegative : kNeutral);
    case Role::Node:
        return kNode;
    default:
        return kMass;
    }
}

Rgb line_color(Role role)
{
    switch (role)
    {
    case Role::AxisX:
        return kAxisX;
    case Role::AxisY:
        return kAxisY;
    case Role::AxisZ:
        return kAxisZ;
    case Role::Ground:
        return kGrid;
    case Role::Wall:
        return kWall;
    case Role::Electric:
        return kElectric;
    case Role::Magnetic:
        return kMagnetic;
    default:
        return kSpringRest;
    }
}

// Appends positions relative to `origin`, in float, last.
class Sink
{
public:
    Sink(DrawList& out, const Vec3& origin) : out_{out}, origin_{origin} {}

    void triangle(const Vec3& a, const Vec3& b, const Vec3& c, const Rgb& color)
    {
        for (const Vec3* p : {&a, &b, &c})
        {
            push(out_.triangles, *p);
        }
        out_.triangle_colors.push_back(color);
    }

    void line(const Vec3& a, const Vec3& b, const Rgb& color)
    {
        push(out_.lines, a);
        push(out_.lines, b);
        out_.line_colors.push_back(color);
    }

private:
    void push(std::vector<float>& to, const Vec3& p)
    {
        const Vec3 r = p - origin_;
        to.push_back(static_cast<float>(r.x));
        to.push_back(static_cast<float>(r.y));
        to.push_back(static_cast<float>(r.z));
    }

    DrawList& out_;
    Vec3 origin_;
};
} // namespace

double point_radius(const Look& look, const Vec3& position)
{
    const Vec3 ndc = view3d::project(look.camera, look.aspect, position);
    if (!(ndc.z > 0))
    {
        return 0.0;
    }
    return 0.5 * view3d::world_size_of_pixels(look.camera, look.viewport_height, ndc.z,
                                              look.point_pixels);
}

void Trails::record(const view3d::Frame& frame)
{
    for (const view3d::Point& p : frame.points)
    {
        if (p.role != Role::Mass && p.role != Role::Charge)
        {
            continue;
        }
        Path& path = paths_[p.source];
        path.role = p.role;
        path.positions.push_back(p.position);
        // Kept between capacity and twice it, so dropping the oldest is cheap.
        if (path.positions.size() >= 2 * capacity_)
        {
            path.positions.erase(path.positions.begin(),
                                 path.positions.begin() +
                                     static_cast<std::ptrdiff_t>(path.positions.size() - capacity_));
        }
    }
}

const std::vector<view3d::Triangle>& Painter::unit_sphere(const view3d::SphereDetail& d)
{
    auto& unit = units_[{d.rings, d.slices}];
    if (unit.empty())
    {
        unit = view3d::sphere_triangles(Vec3{}, 1.0, d);
    }
    return unit;
}

void Painter::paint(const view3d::Frame& frame, const Look& look, const Trails* trails,
                    bool subjects_only, DrawList& out)
{
    out.clear();
    const view3d::Camera& cam = look.camera;
    Sink sink{out, cam.target};

    // A headlight from above and to the left of the eye, fixed to the camera,
    // so the faces turned toward the viewer are always the lit ones.
    const Vec3 forward = malloy::math::normalize(cam.target - cam.position);
    const Vec3 right = malloy::math::normalize(malloy::math::cross(forward, cam.up));
    const Vec3 up = malloy::math::cross(right, forward);
    const Vec3 to_light = malloy::math::normalize(-forward + up * 0.55 - right * 0.35);
    const auto lit = [&](const Rgb& base, const Vec3& normal) {
        return scaled(base, view3d::shade(normal, to_light, kAmbient));
    };

    // Ground patches first, both faces, so the floor is seen from below too.
    if (!subjects_only)
    {
        for (const view3d::Patch& patch : frame.patches)
        {
            for (const view3d::Triangle& t : view3d::patch_triangles(patch))
            {
                sink.triangle(t.a, t.b, t.c, lit(kGround, t.normal));
                sink.triangle(t.a, t.c, t.b, scaled(lit(kGround, -t.normal), 0.6));
            }
            for (const view3d::Segment& s : view3d::patch_grid(patch))
            {
                sink.line(s.a, s.b, kGrid);
            }
        }
    }

    // Spheres, from one unit mesh per level of detail.
    const view3d::SphereDetail detail = view3d::sphere_detail(frame.spheres.size());
    const std::vector<view3d::Triangle>& unit = unit_sphere(detail);
    for (const view3d::Sphere& s : frame.spheres)
    {
        placed_.clear();
        view3d::place_sphere(unit, s.center, s.radius, placed_);
        for (const view3d::Triangle& t : placed_)
        {
            sink.triangle(t.a, t.b, t.c, lit(kBody, t.normal));
        }
    }

    // Boxes: twelve faces each, and their edges a shade darker.
    for (const view3d::Box& b : frame.boxes)
    {
        const Rgb base = b.role == Role::Proxy ? kProxy : kBody;
        for (const view3d::Triangle& t : view3d::box_triangles(b))
        {
            sink.triangle(t.a, t.b, t.c, lit(base, t.normal));
        }
        if (!subjects_only)
        {
            for (const auto& e : view3d::box_edges)
            {
                sink.line(b.corners[static_cast<std::size_t>(e[0])],
                          b.corners[static_cast<std::size_t>(e[1])], scaled(base, 0.45));
            }
        }
    }

    // Points have no size of their own: a small sphere of a fixed size on screen.
    const std::vector<view3d::Triangle>& small = unit_sphere(view3d::SphereDetail{6, 10});
    for (const view3d::Point& p : frame.points)
    {
        const double radius = point_radius(look, p.position);
        if (!(radius > 0))
        {
            continue;
        }
        placed_.clear();
        view3d::place_sphere(small, p.position, radius, placed_);
        const Rgb base = point_color(p);
        for (const view3d::Triangle& t : placed_)
        {
            sink.triangle(t.a, t.b, t.c, lit(base, t.normal));
        }
    }

    if (subjects_only)
    {
        return;
    }

    // Lines: axis triads, springs coloured by how far each is from its rest
    // length (red stretched, blue compressed), and the context.
    Real largest_extension = 0;
    for (const view3d::Segment& s : frame.segments)
    {
        if (s.role == Role::Spring)
        {
            largest_extension = std::max(largest_extension, std::abs(s.weight));
        }
    }
    for (const view3d::Segment& s : frame.segments)
    {
        Rgb color = line_color(s.role);
        if (s.role == Role::Spring && largest_extension > 0)
        {
            const Real t = std::abs(s.weight) / largest_extension;
            color = mix(kSpringRest, s.weight > 0 ? kPositive : kNegative, t);
        }
        sink.line(s.a, s.b, color);
    }
    for (const view3d::Arrow& a : frame.arrows)
    {
        for (const view3d::Segment& s : view3d::arrow_lines(a))
        {
            sink.line(s.a, s.b, line_color(a.role));
        }
    }
    if (trails != nullptr)
    {
        for (const auto& [source, path] : trails->paths())
        {
            const Rgb base = path.role == Role::Charge ? kNeutral : kMass;
            const Rgb color = scaled(base, 0.7);
            for (std::size_t i = 1; i < path.positions.size(); ++i)
            {
                sink.line(path.positions[i - 1], path.positions[i], color);
            }
        }
    }
}

CameraSetup camera_setup(const view3d::Camera& camera)
{
    CameraSetup setup;
    const Vec3 eye = camera.position - camera.target;
    setup.position[0] = static_cast<float>(eye.x);
    setup.position[1] = static_cast<float>(eye.y);
    setup.position[2] = static_cast<float>(eye.z);
    setup.up[0] = static_cast<float>(camera.up.x);
    setup.up[1] = static_cast<float>(camera.up.y);
    setup.up[2] = static_cast<float>(camera.up.z);
    setup.fovy_degrees = camera.fovy * 180.0 / kPi;
    setup.near_plane = camera.near_plane;
    setup.far_plane = camera.far_plane;
    return setup;
}
} // namespace viewer
