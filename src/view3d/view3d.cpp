#include <malloy/view3d/view3d.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <vector>

#include <malloy/charges/charge3d_world.hpp>
#include <malloy/collide/shapes3d.hpp>
#include <malloy/math/quat.hpp>
#include <malloy/math/vec3.hpp>
#include <malloy/nbody/nbody3d_world.hpp>
#include <malloy/particles/particle3d_world.hpp>
#include <malloy/rigid/rigid3d_world.hpp>
#include <malloy/springs/spring_world3d.hpp>

namespace malloy::view3d
{
namespace
{
constexpr Real pi = 3.14159265358979323846;

bool finite(const Vec3& v)
{
    return math::is_finite(v);
}

template <std::size_t N>
bool finite(const std::array<Vec3, N>& points)
{
    return std::all_of(points.begin(), points.end(), [](const Vec3& p) { return finite(p); });
}

Real largest(const Vec3& v)
{
    return std::max({v.x, v.y, v.z});
}

// The direction from a target to an eye at this azimuth and elevation.
Vec3 eye_direction(Real yaw, Real pitch)
{
    return Vec3{std::cos(pitch) * std::cos(yaw), std::cos(pitch) * std::sin(yaw),
                std::sin(pitch)};
}

bool sane(const Bounds& bounds)
{
    return finite(bounds.center) && math::is_finite(bounds.radius) &&
           bounds.radius > Real{0};
}

ViewAngles sanitized(const ViewAngles& in)
{
    const ViewAngles defaults{};
    ViewAngles a = in;
    if (!math::is_finite(a.fovy) || !(a.fovy > Real{0}) || !(a.fovy < pi))
    {
        a.fovy = defaults.fovy;
    }
    if (!math::is_finite(a.aspect) || !(a.aspect > Real{0}))
    {
        a.aspect = defaults.aspect;
    }
    if (!math::is_finite(a.margin) || !(a.margin >= Real{0}))
    {
        a.margin = defaults.margin;
    }
    if (!math::is_finite(a.yaw))
    {
        a.yaw = defaults.yaw;
    }
    if (!math::is_finite(a.pitch))
    {
        a.pitch = defaults.pitch;
    }
    a.pitch = std::clamp(a.pitch, -pitch_limit, pitch_limit);
    return a;
}

// Three equal arms along the body's own axes. Equal, so the axis normal to a
// flat plate (a gyroscope's spin axis) shows as plainly as the others.
std::array<Segment, 3> triad(const Vec3& position, const Quat& q, Real arm,
                             std::size_t source)
{
    return {Segment{position, position + math::rotate(q, Vec3{arm, 0.0, 0.0}), 0.0,
                    source, Role::AxisX},
            Segment{position, position + math::rotate(q, Vec3{0.0, arm, 0.0}), 0.0,
                    source, Role::AxisY},
            Segment{position, position + math::rotate(q, Vec3{0.0, 0.0, arm}), 0.0,
                    source, Role::AxisZ}};
}

bool finite_segment(const Segment& s)
{
    return finite(s.a) && finite(s.b) && math::is_finite(s.weight);
}

Patch make_patch(const collide::Plane3& plane, const Bounds& view, std::size_t source)
{
    Patch patch;
    patch.normal = plane.normal;
    patch.offset = plane.offset;
    patch.source = source;
    patch.center =
        view.center - plane.normal * (math::dot(plane.normal, view.center) - plane.offset);
    tangents(plane.normal, patch.u, patch.v);
    patch.half_size = Real{2} * view.radius;
    const Real h = patch.half_size;
    patch.corners = {patch.center - patch.u * h - patch.v * h,
                     patch.center + patch.u * h - patch.v * h,
                     patch.center + patch.u * h + patch.v * h,
                     patch.center - patch.u * h + patch.v * h};
    return patch;
}

bool finite_patch(const Patch& p)
{
    return finite(p.center) && finite(p.normal) && math::is_finite(p.offset) &&
           finite(p.u) && finite(p.v) && math::is_finite(p.half_size) && finite(p.corners);
}

// A field arrow beside the view centre (side +1 or -1 keeps E and B apart),
// half the view radius long, pointing along the field.
Arrow make_arrow(const Vec3& field, const Bounds& view, Real side, Role role)
{
    const Vec3 direction = math::normalize(field);
    Vec3 u;
    Vec3 v;
    tangents(direction, u, v);
    const Vec3 origin = view.center + u * (side * Real{0.75} * view.radius);
    return Arrow{origin, origin + direction * (Real{0.5} * view.radius), role};
}

void include(Vec3& lo, Vec3& hi, const Vec3& p)
{
    lo = Vec3{std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = Vec3{std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
}
} // namespace

bool is_subject(Role role)
{
    return role != Role::Ground && role != Role::Wall && role != Role::Electric &&
           role != Role::Magnetic;
}

std::array<Vec3, 8> box_corners(const Vec3& center, const Quat& orientation,
                                const Vec3& half_extents)
{
    std::array<Vec3, 8> out{};
    for (std::size_t i = 0; i < out.size(); ++i)
    {
        const Vec3 local{(i & 1U) ? half_extents.x : -half_extents.x,
                         (i & 2U) ? half_extents.y : -half_extents.y,
                         (i & 4U) ? half_extents.z : -half_extents.z};
        out[i] = center + math::rotate(orientation, local);
    }
    return out;
}

Vec3 inertia_box(const Vec3& inertia, Real mass, Real min_fraction)
{
    const Real k = Real{3} / (Real{2} * mass);
    const Real a2 = k * (inertia.y + inertia.z - inertia.x);
    const Real b2 = k * (inertia.x + inertia.z - inertia.y);
    const Real c2 = k * (inertia.x + inertia.y - inertia.z);
    const Vec3 h{std::sqrt(std::max(a2, Real{0})), std::sqrt(std::max(b2, Real{0})),
                 std::sqrt(std::max(c2, Real{0}))};
    const Real floor = largest(h) * min_fraction;
    return Vec3{std::max(h.x, floor), std::max(h.y, floor), std::max(h.z, floor)};
}

void tangents(const Vec3& normal, Vec3& u, Vec3& v)
{
    const Real ax = std::abs(normal.x);
    const Real ay = std::abs(normal.y);
    const Real az = std::abs(normal.z);
    Vec3 axis{1.0, 0.0, 0.0};
    if (ay < ax && ay <= az)
    {
        axis = Vec3{0.0, 1.0, 0.0};
    }
    else if (az < ax && az < ay)
    {
        axis = Vec3{0.0, 0.0, 1.0};
    }
    u = math::normalize(math::cross(axis, normal));
    v = math::cross(normal, u);
}

Frame build_frame(const rigid::Rigid3DWorld& world, const Bounds& view)
{
    Frame frame;
    frame.tick = world.tick_count();
    frame.time = world.elapsed_time();
    const std::vector<rigid::RigidBody3D>& bodies = world.bodies();
    for (std::size_t i = 0; i < bodies.size(); ++i)
    {
        const rigid::RigidBody3D& body = bodies[i];
        if (!finite(body.position) || !math::is_finite(body.orientation) ||
            !finite(body.half_extents) || !math::is_finite(body.radius) ||
            !finite(body.inertia) || !math::is_finite(body.mass))
        {
            ++frame.dropped;
            continue;
        }
        // The collider rule: a box when all three half-widths are positive,
        // else a sphere when the radius is, else no collider at all, drawn as
        // the box with the body's inertia.
        const Vec3& h = body.half_extents;
        const bool is_box = h.x > Real{0} && h.y > Real{0} && h.z > Real{0};
        const bool is_sphere = !is_box && body.radius > Real{0};
        const Vec3 extents = is_box      ? h
                             : is_sphere ? Vec3{body.radius, body.radius, body.radius}
                                         : inertia_box(body.inertia, body.mass);
        const Box box{body.position, body.orientation, extents,
                      box_corners(body.position, body.orientation, extents), i,
                      is_box ? Role::Body : Role::Proxy};
        const std::array<Segment, 3> axes =
            triad(body.position, body.orientation, Real{1.25} * largest(extents), i);

        // A body is drawn whole or not at all.
        if (!finite(extents) || !(largest(extents) > Real{0}) || !finite(box.corners) ||
            !std::all_of(axes.begin(), axes.end(), finite_segment))
        {
            ++frame.dropped;
            continue;
        }
        if (is_sphere)
        {
            frame.spheres.push_back(Sphere{body.position, body.radius, i, Role::Body});
        }
        else
        {
            frame.boxes.push_back(box);
        }
        frame.segments.insert(frame.segments.end(), axes.begin(), axes.end());
    }

    const std::vector<collide::Plane3>& ground = world.settings().ground;
    for (std::size_t j = 0; j < ground.size(); ++j)
    {
        const Patch patch = make_patch(ground[j], view, j);
        if (!ground[j].is_valid() || !finite_patch(patch))
        {
            ++frame.dropped;
            continue;
        }
        frame.patches.push_back(patch);
    }
    return frame;
}

Frame build_frame(const particles::ParticleWorld3D& world, const Bounds& /*view*/)
{
    Frame frame;
    frame.tick = world.tick_count();
    frame.time = world.elapsed_time();
    const std::vector<particles::Particle3D>& list = world.particles();
    for (std::size_t i = 0; i < list.size(); ++i)
    {
        if (!finite(list[i].position) || !math::is_finite(list[i].radius))
        {
            ++frame.dropped;
            continue;
        }
        frame.spheres.push_back(Sphere{list[i].position, list[i].radius, i, Role::Body});
    }

    // The walls' corners are taken from min and max directly, by sign bit: a
    // centre plus or minus a half-width does not land back on max exactly.
    const collide::Aabb3& b = world.particle_settings().bounds;
    if (!finite(b.min) || !finite(b.max))
    {
        ++frame.dropped;
        return frame;
    }
    std::array<Vec3, 8> c{};
    for (std::size_t i = 0; i < c.size(); ++i)
    {
        c[i] = Vec3{(i & 1U) ? b.max.x : b.min.x, (i & 2U) ? b.max.y : b.min.y,
                    (i & 4U) ? b.max.z : b.min.z};
    }
    for (const std::array<int, 2>& e : box_edges)
    {
        frame.segments.push_back(Segment{c[static_cast<std::size_t>(e[0])],
                                         c[static_cast<std::size_t>(e[1])], 0.0, 0,
                                         Role::Wall});
    }
    return frame;
}

Frame build_frame(const nbody::NBody3DWorld& world, const Bounds& /*view*/)
{
    Frame frame;
    frame.tick = world.tick_count();
    frame.time = world.elapsed_time();
    const std::vector<nbody::Body3D>& list = world.bodies();
    for (std::size_t i = 0; i < list.size(); ++i)
    {
        if (!finite(list[i].position) || !math::is_finite(list[i].mass))
        {
            ++frame.dropped;
            continue;
        }
        frame.points.push_back(Point{list[i].position, list[i].mass, i, Role::Mass});
    }
    return frame;
}

Frame build_frame(const springs::SpringWorld3D& world, const Bounds& /*view*/)
{
    Frame frame;
    frame.tick = world.tick_count();
    frame.time = world.elapsed_time();
    const std::vector<springs::SpringBody3D>& list = world.bodies();
    for (std::size_t i = 0; i < list.size(); ++i)
    {
        if (!finite(list[i].position) || !math::is_finite(list[i].mass))
        {
            ++frame.dropped;
            continue;
        }
        frame.points.push_back(Point{list[i].position, list[i].mass, i, Role::Node});
    }

    // Extension rather than strain: a rest length of zero is valid, and the
    // sign (stretched or compressed) is what a colour shows.
    const std::vector<springs::Spring>& springs = world.network().springs();
    for (std::size_t k = 0; k < springs.size(); ++k)
    {
        const springs::Spring& s = springs[k];
        if (s.a >= list.size() || s.b >= list.size())
        {
            ++frame.dropped;
            continue;
        }
        const Vec3 a = list[s.a].position;
        const Vec3 b = list[s.b].position;
        const Segment segment{a, b, math::distance(a, b) - s.rest_length, k, Role::Spring};
        if (!finite_segment(segment))
        {
            ++frame.dropped;
            continue;
        }
        frame.segments.push_back(segment);
    }
    return frame;
}

Frame build_frame(const charges::Charge3DWorld& world, const Bounds& view)
{
    Frame frame;
    frame.tick = world.tick_count();
    frame.time = world.elapsed_time();
    const std::vector<charges::ChargedParticle3D>& list = world.particles();
    for (std::size_t i = 0; i < list.size(); ++i)
    {
        if (!finite(list[i].position) || !math::is_finite(list[i].charge))
        {
            ++frame.dropped;
            continue;
        }
        frame.points.push_back(Point{list[i].position, list[i].charge, i, Role::Charge});
    }

    const charges::Charge3DSettings& settings = world.charge_settings();
    const Real sides[2] = {Real{1}, Real{-1}};
    const Vec3 fields[2] = {settings.electric, settings.magnetic};
    const Role roles[2] = {Role::Electric, Role::Magnetic};
    for (std::size_t k = 0; k < 2; ++k)
    {
        if (finite(fields[k]) && math::length_squared(fields[k]) == Real{0})
        {
            continue; // no field, no arrow
        }
        const Arrow arrow = make_arrow(fields[k], view, sides[k], roles[k]);
        if (!finite(fields[k]) || !finite(arrow.origin) || !finite(arrow.tip))
        {
            ++frame.dropped;
            continue;
        }
        frame.arrows.push_back(arrow);
    }
    return frame;
}

bool is_finite(const Frame& f)
{
    const bool spheres = std::all_of(f.spheres.begin(), f.spheres.end(), [](const Sphere& s) {
        return finite(s.center) && math::is_finite(s.radius);
    });
    const bool boxes = std::all_of(f.boxes.begin(), f.boxes.end(), [](const Box& b) {
        return finite(b.center) && math::is_finite(b.orientation) && finite(b.half_extents) &&
               finite(b.corners);
    });
    const bool patches = std::all_of(f.patches.begin(), f.patches.end(), finite_patch);
    const bool segments = std::all_of(f.segments.begin(), f.segments.end(), finite_segment);
    const bool points = std::all_of(f.points.begin(), f.points.end(), [](const Point& p) {
        return finite(p.position) && math::is_finite(p.weight);
    });
    const bool arrows = std::all_of(f.arrows.begin(), f.arrows.end(), [](const Arrow& a) {
        return finite(a.origin) && finite(a.tip);
    });
    return spheres && boxes && patches && segments && points && arrows &&
           math::is_finite(f.time);
}

Bounds fit_bounds(const Frame& frame)
{
    const Real inf = std::numeric_limits<Real>::infinity();
    Vec3 lo{inf, inf, inf};
    Vec3 hi{-inf, -inf, -inf};
    bool any = false;

    for (const Sphere& s : frame.spheres)
    {
        if (!is_subject(s.role) || !finite(s.center) || !math::is_finite(s.radius))
        {
            continue;
        }
        const Vec3 r{s.radius, s.radius, s.radius};
        include(lo, hi, s.center - r);
        include(lo, hi, s.center + r);
        any = true;
    }
    for (const Box& b : frame.boxes)
    {
        if (!is_subject(b.role))
        {
            continue;
        }
        for (const Vec3& c : b.corners)
        {
            if (finite(c))
            {
                include(lo, hi, c);
                any = true;
            }
        }
    }
    for (const Segment& s : frame.segments)
    {
        if (!is_subject(s.role))
        {
            continue;
        }
        for (const Vec3& p : {s.a, s.b})
        {
            if (finite(p))
            {
                include(lo, hi, p);
                any = true;
            }
        }
    }
    for (const Point& p : frame.points)
    {
        if (is_subject(p.role) && finite(p.position))
        {
            include(lo, hi, p.position);
            any = true;
        }
    }
    if (!any)
    {
        return Bounds{};
    }

    const Vec3 extent = hi - lo;
    Bounds out{(lo + hi) * Real{0.5},
               std::hypot(extent.x, extent.y, extent.z) * Real{0.5}};
    // Everything at one point: show a unit neighbourhood, as fit_viewport
    // widens a zero extent by one either side.
    if (!(out.radius > Real{0}))
    {
        out.radius = Real{1};
    }
    // Finite subjects whose extent overflows: the documented default.
    if (!finite(out.center) || !math::is_finite(out.radius))
    {
        return Bounds{};
    }
    return out;
}

Bounds merge(const Bounds& a, const Bounds& b)
{
    if (!sane(b))
    {
        return a;
    }
    if (!sane(a))
    {
        return b;
    }
    const Real d = math::distance(a.center, b.center);
    if (d + b.radius <= a.radius)
    {
        return a;
    }
    if (d + a.radius <= b.radius)
    {
        return b;
    }
    const Real radius = (d + a.radius + b.radius) * Real{0.5};
    const Vec3 center = a.center + (b.center - a.center) * ((radius - a.radius) / d);
    const Bounds out{center, radius};
    return sane(out) ? out : a;
}

OrbitView fit_view(const Bounds& bounds_in, const ViewAngles& angles_in)
{
    const ViewAngles a = sanitized(angles_in);
    const Bounds bounds = sane(bounds_in) ? bounds_in : Bounds{};

    const Real half_y = a.fovy * Real{0.5};
    const Real half_x = std::atan(a.aspect * std::tan(half_y));
    const Real theta = std::min(half_y, half_x);
    const Real framed = bounds.radius * (Real{1} + a.margin);
    const OrbitView view{bounds.center, a.yaw, a.pitch, framed / std::sin(theta), a.fovy};
    if (!math::is_finite(view.distance) || !(view.distance > Real{0}) ||
        !finite(view.target + eye_direction(view.yaw, view.pitch) * view.distance))
    {
        return fit_view(Bounds{}, a);
    }
    return view;
}

Camera to_camera(const OrbitView& view_in, const Bounds& scene_in)
{
    const Bounds scene = sane(scene_in) ? scene_in : Bounds{};
    OrbitView view = view_in;
    if (!finite(view.target) || !math::is_finite(view.distance) ||
        !(view.distance > Real{0}) || !math::is_finite(view.yaw) ||
        !math::is_finite(view.pitch) || !math::is_finite(view.fovy) ||
        !(view.fovy > Real{0}) || !(view.fovy < pi))
    {
        view = fit_view(scene, ViewAngles{});
    }
    view.pitch = std::clamp(view.pitch, -pitch_limit, pitch_limit);

    const Vec3 direction = eye_direction(view.yaw, view.pitch);
    Camera camera;
    camera.target = view.target;
    camera.position = view.target + direction * view.distance;
    camera.up = Vec3{0.0, 0.0, 1.0};
    camera.fovy = view.fovy;

    // Where the scene sphere lies along the view direction decides the clip
    // planes: twice the distance to its back, and half the gap to its front,
    // but never nearer than a thousandth of far, which keeps the depth
    // buffer's precision; inside the sphere that floor is the near plane.
    const Real middle = math::dot(-direction, scene.center - camera.position);
    const Real front = middle - scene.radius;
    const Real back = middle + scene.radius;
    camera.far_plane =
        back > Real{0} ? Real{2} * back : Real{2} * (std::abs(middle) + scene.radius);
    camera.near_plane = std::max(Real{0.5} * front, Real{1e-3} * camera.far_plane);

    if (!finite(camera.position) || !math::is_finite(camera.near_plane) ||
        !math::is_finite(camera.far_plane) || !(camera.near_plane > Real{0}) ||
        !(camera.far_plane > camera.near_plane))
    {
        return to_camera(fit_view(Bounds{}, ViewAngles{}), Bounds{});
    }
    return camera;
}

Camera fit_camera(const Bounds& bounds, const ViewAngles& angles)
{
    return to_camera(fit_view(bounds, angles), bounds);
}

void orbit(OrbitView& view, Real dx_pixels, Real dy_pixels, Real radians_per_pixel)
{
    const Real yaw = view.yaw - dx_pixels * radians_per_pixel;
    const Real pitch = view.pitch + dy_pixels * radians_per_pixel;
    if (!math::is_finite(yaw) || !math::is_finite(pitch))
    {
        return;
    }
    view.yaw = yaw;
    view.pitch = std::clamp(pitch, -pitch_limit, pitch_limit);
}

void pan(OrbitView& view, Real dx_pixels, Real dy_pixels, Real viewport_height_pixels)
{
    if (!(viewport_height_pixels > Real{0}))
    {
        return;
    }
    const Vec3 forward = -eye_direction(view.yaw, view.pitch);
    const Vec3 right = math::normalize(math::cross(forward, Vec3{0.0, 0.0, 1.0}));
    const Vec3 up = math::cross(right, forward);
    const Real per_pixel =
        Real{2} * view.distance * std::tan(view.fovy * Real{0.5}) / viewport_height_pixels;
    // The scene follows the cursor: dragging right moves the target left, and
    // dragging down (screen y grows down) moves it up.
    const Vec3 target =
        view.target - right * (dx_pixels * per_pixel) + up * (dy_pixels * per_pixel);
    if (finite(target))
    {
        view.target = target;
    }
}

void zoom(OrbitView& view, Real wheel_steps, Real min_distance, Real max_distance)
{
    const Real distance = view.distance * std::pow(Real{0.85}, wheel_steps);
    if (!math::is_finite(distance) || !(min_distance <= max_distance))
    {
        return;
    }
    view.distance = std::clamp(distance, min_distance, max_distance);
}

Vec3 project(const Camera& camera, Real aspect, const Vec3& p)
{
    const Vec3 f = math::normalize(camera.target - camera.position);
    const Vec3 s = math::normalize(math::cross(f, camera.up));
    const Vec3 u = math::cross(s, f);
    const Vec3 rel = p - camera.position;
    const Real depth = math::dot(f, rel);
    const Real t = std::tan(camera.fovy * Real{0.5});
    return Vec3{math::dot(s, rel) / (depth * aspect * t), math::dot(u, rel) / (depth * t),
                depth};
}

Real world_size_of_pixels(const Camera& camera, Real viewport_height_pixels, Real depth,
                          Real pixels)
{
    return pixels * Real{2} * depth * std::tan(camera.fovy * Real{0.5}) /
           viewport_height_pixels;
}

Pacer::Pacer(Real dt, std::uint64_t max_steps, std::uint64_t end_tick)
    : dt_{dt}, max_steps_{max_steps}, end_tick_{end_tick}
{
}

void Pacer::rebase(Real wall_now, std::uint64_t tick_now)
{
    base_wall_ = wall_now;
    base_tick_ = tick_now;
}

std::uint64_t Pacer::due(Real wall_now, std::uint64_t tick_now)
{
    if (paused_ || tick_now >= end_tick_)
    {
        // Nothing is owed while paused or once the run has ended, and nothing
        // accrues: the clock restarts from now.
        rebase(wall_now, tick_now);
        behind_ = false;
        return 0;
    }
    if (!math::is_finite(dt_) || !(dt_ > Real{0}))
    {
        return 0;
    }
    const Real owed = std::floor((wall_now - base_wall_) * speed_ / dt_);
    if (!math::is_finite(owed) || !(owed >= Real{0}))
    {
        return 0;
    }
    // In Real, so an absurd wall time cannot overflow a cast: the ticks stay
    // far below 2^53, where the sum is exact.
    const Real wanted = static_cast<Real>(base_tick_) + owed - static_cast<Real>(tick_now);
    if (!(wanted >= Real{1}))
    {
        behind_ = false;
        return 0;
    }
    const std::uint64_t left = end_tick_ - tick_now;
    if (wanted > static_cast<Real>(max_steps_) && max_steps_ < left)
    {
        behind_ = true;
        rebase(wall_now, tick_now + max_steps_);
        return max_steps_;
    }
    behind_ = false;
    const std::uint64_t n = wanted >= static_cast<Real>(left)
                                ? left
                                : static_cast<std::uint64_t>(wanted);
    return std::min(n, max_steps_);
}

void Pacer::set_speed(Real speed, Real wall_now, std::uint64_t tick_now)
{
    if (!math::is_finite(speed) || !(speed > Real{0}))
    {
        return;
    }
    speed_ = std::clamp(speed, Real{1} / Real{64}, Real{64});
    rebase(wall_now, tick_now);
}

void Pacer::set_paused(bool paused, Real wall_now, std::uint64_t tick_now)
{
    paused_ = paused;
    rebase(wall_now, tick_now);
}
} // namespace malloy::view3d
