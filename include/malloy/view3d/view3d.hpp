#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
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
// The render model (M41): the 3D sibling of malloy_ascii. It turns a 3D
// world's state into a FRAME, a list of primitives in world coordinates whose
// geometry is already final, and it frames them with a camera. A viewer draws
// a frame with no rotation, scaling or projection maths of its own, so every
// number a viewer could get wrong is computed here, where it is tested without
// a window. mesh.hpp turns the primitives into triangles.
//
// It knows no graphics library and never will: it is a client of the domain
// libraries like any other, and it only ever READS a world (every builder
// takes it by const reference). Building frames does not change a simulation,
// and the tests hold that to bit identity on every 3D template.
//
// Coordinates are the simulation's: z is up (gravity3 0 0 -g, plane3 0 0 1),
// right-handed. Everything stays in Real; a viewer converts to float last.

using math::Quat;
using math::Real;
using math::Vec3;

// What a primitive stands for, so a viewer can colour it without knowing which
// domain made it. SUBJECT roles are what the camera frames; CONTEXT roles
// (Ground, Wall, Electric, Magnetic) are drawn but never enlarge the fit, so a
// particle box a hundred times wider than its particles does not shrink them
// to dots.
enum class Role : std::uint8_t
{
    Body,     // a rigid body or particle with a collider
    Proxy,    // a rigid body with no collider, drawn as the box with its inertia
    Mass,     // an N-body point mass; weight = mass
    Charge,   // a charged particle; weight = charge, whose sign is its colour
    Node,     // a spring body; weight = mass
    Spring,   // a spring; weight = extension (length - rest), positive stretched
    AxisX,    // a rigid body's own x axis
    AxisY,    // ... its y axis
    AxisZ,    // ... its z axis
    Ground,   // a ground plane patch, or a line of its grid   (context)
    Wall,     // an edge of a particle box                    (context)
    Electric, // the uniform electric field direction         (context)
    Magnetic, // the uniform magnetic field direction         (context)
};

bool is_subject(Role role);

struct Sphere
{
    Vec3 center{};
    Real radius{0.0};
    std::size_t source{0}; // index in the world's body or particle list
    Role role{Role::Body};
};

// An oriented box. `corners` is final: corner i has the local sign bits
// (i & 1 -> +x, i & 2 -> +y, i & 4 -> +z), the order collide uses, so
// `box_edges` and `box_faces` below index it directly.
struct Box
{
    Vec3 center{};
    Quat orientation{};
    Vec3 half_extents{};
    std::array<Vec3, 8> corners{};
    std::size_t source{0};
    Role role{Role::Body};
};

// A finite square of an infinite plane dot(normal, p) == offset, centred on the
// point of the plane nearest the view centre. (u, v, normal) is right-handed,
// and the corners run counter-clockwise seen from the FREE side, the side the
// normal points into.
struct Patch
{
    Vec3 center{};
    Vec3 normal{0.0, 0.0, 1.0};
    Real offset{0.0};
    Vec3 u{1.0, 0.0, 0.0};
    Vec3 v{0.0, 1.0, 0.0};
    Real half_size{0.0};
    std::array<Vec3, 4> corners{};
    std::size_t source{0}; // index in the world's ground planes
    Role role{Role::Ground};
};

struct Segment
{
    Vec3 a{};
    Vec3 b{};
    Real weight{0.0};
    std::size_t source{0};
    Role role{Role::Spring};
};

// Drawn at a fixed size on screen (`world_size_of_pixels`): a point mass, a
// charge or a spring body has no radius of its own.
struct Point
{
    Vec3 position{};
    Real weight{0.0};
    std::size_t source{0};
    Role role{Role::Mass};
};

struct Arrow
{
    Vec3 origin{};
    Vec3 tip{};
    Role role{Role::Electric};
};

struct Frame
{
    std::vector<Sphere> spheres;
    std::vector<Box> boxes;
    std::vector<Patch> patches;
    std::vector<Segment> segments;
    std::vector<Point> points;
    std::vector<Arrow> arrows;
    std::uint64_t tick{0};
    Real time{0.0};
    // Entries left out because their input was not finite (or made a
    // non-finite shape). A valid world never produces one; the count is there
    // so a drop is visible rather than silent.
    std::size_t dropped{0};
};

// The sphere the camera frames.
struct Bounds
{
    Vec3 center{};
    Real radius{1.0};
};

// The twelve edges and six faces of a box, as indices into Box::corners. The
// faces are -x, +x, -y, +y, -z, +z in that order, each wound counter-clockwise
// seen from OUTSIDE, so a viewer culling back faces shows the near ones.
inline constexpr std::array<std::array<int, 2>, 12> box_edges{{
    {0, 1}, {2, 3}, {4, 5}, {6, 7}, // along x
    {0, 2}, {1, 3}, {4, 6}, {5, 7}, // along y
    {0, 4}, {1, 5}, {2, 6}, {3, 7}, // along z
}};
inline constexpr std::array<std::array<int, 4>, 6> box_faces{{
    {0, 4, 6, 2}, // -x
    {1, 3, 7, 5}, // +x
    {0, 1, 5, 4}, // -y
    {2, 6, 7, 3}, // +y
    {0, 2, 3, 1}, // -z
    {4, 5, 7, 6}, // +z
}};

// The eight world-space corners of a box: each local corner (signs times half
// widths) turned by the orientation, then moved to the centre. Scale, then
// rotate, then translate; any other order is a different box.
std::array<Vec3, 8> box_corners(const Vec3& center, const Quat& orientation,
                                const Vec3& half_extents);

// Half-widths of the uniform solid box whose principal moments are `inertia`:
// from I_x = m (b^2 + c^2) / 3 and its two siblings, a^2 = 3 (I_y + I_z - I_x)
// / (2 m). Used to draw a rigid body that has no collider, so its rotation is
// visible and the drawn shape has its inertia. A moment triple that breaks the
// triangle inequality (gyroscope.scn: I_z > I_x + I_y) has no such box, and a
// flat plate has a zero half-width, so every half-width is held to at least
// `min_fraction` of the largest.
Vec3 inertia_box(const Vec3& inertia, Real mass, Real min_fraction = 0.05);

// Two unit tangents with (u, v, normal) right-handed, for a unit normal, by a
// fixed rule (the world axis least aligned with the normal, ties to the
// earlier axis), so a plane is always drawn the same way.
void tangents(const Vec3& normal, Vec3& u, Vec3& v);

// Builders: one overload per 3D world, no template and no base class
// (ADR 0006). `view` is the camera's current bounding sphere; only context
// (ground patches, field arrows) is placed from it.
//
//   rigid3d      a Box when all three half_extents are positive (the collider
//                rule), else a Sphere when the radius is, else a Proxy box
//                (`inertia_box`); every body gets an axis triad of three equal
//                arms, 1.25 times its largest half-width; each ground plane a
//                Patch of half-size 2 R
//   particles3d  a Sphere per particle and the twelve Wall edges of the bounds
//   nbody3d      a Point per body, weight = mass
//   springs3d    a Node Point per body and a Spring Segment per spring,
//                weight = length - rest length
//   charges3d    a Point per particle, weight = charge, and an Arrow per
//                nonzero uniform field, half the view radius long
Frame build_frame(const rigid::Rigid3DWorld& world, const Bounds& view);
Frame build_frame(const particles::ParticleWorld3D& world, const Bounds& view);
Frame build_frame(const nbody::NBody3DWorld& world, const Bounds& view);
Frame build_frame(const springs::SpringWorld3D& world, const Bounds& view);
Frame build_frame(const charges::Charge3DWorld& world, const Bounds& view);

// Every value in the frame is finite.
bool is_finite(const Frame& frame);

// The bounding sphere of the frame's SUBJECTS: the sphere around their
// axis-aligned bounding box. The guards are fit_viewport's: no finite subject
// gives the default {0, 1}; all subjects at one point give radius 1; an extent
// that overflows to infinity falls back to the default, which then does not
// contain the subjects (as fit_viewport, issue #5).
Bounds fit_bounds(const Frame& frame);

// The smallest sphere containing both, for a grow-only auto-fit: grown after
// every step, the view never loses a body.
Bounds merge(const Bounds& a, const Bounds& b);

// How the default view looks at a scene. The default direction lies along no
// world axis and in no world plane, so motion along any axis shows.
struct ViewAngles
{
    Real yaw{-1.0471975511965976};  // -60 degrees: the eye's azimuth about +z, from +x
    Real pitch{0.4363323129985824}; // 25 degrees above the horizon
    Real fovy{0.7853981633974483};  // 45 degrees, vertical, in RADIANS
    Real aspect{16.0 / 9.0};        // width / height
    Real margin{0.1};               // fraction of the radius left clear
};

// The pitch is held within this of the poles, so the camera basis never
// degenerates (z up and a vertical view direction have no cross product).
inline constexpr Real pitch_limit = 1.5533430342749532; // 89 degrees

// An orbit camera: the state a viewer's mouse changes.
struct OrbitView
{
    Vec3 target{};
    Real yaw{-1.0471975511965976};
    Real pitch{0.4363323129985824};
    Real distance{3.0};
    Real fovy{0.7853981633974483};
};

// A perspective camera, z up, OpenGL conventions: it looks down its own -z,
// NDC x points right and NDC y up. `fovy` is in RADIANS (raylib wants degrees).
struct Camera
{
    Vec3 target{};
    Vec3 position{0.0, -3.0, 0.0};
    Vec3 up{0.0, 0.0, 1.0};
    Real fovy{0.7853981633974483};
    Real near_plane{0.01};
    Real far_plane{100.0};
};

// Frames `bounds`: the sphere of radius R (1 + margin) exactly fills the
// narrower field of view, dist = R (1 + margin) / sin(theta) with
// theta = min(fovy / 2, atan(aspect tan(fovy / 2))). A sphere fit is the same
// in every direction, so orbiting or resizing never clips a body. Non-finite
// or out-of-range angles, aspect or margin are replaced by the defaults, and a
// non-finite result by the fit of the default bounds.
OrbitView fit_view(const Bounds& bounds, const ViewAngles& angles);

// The drawing camera for an orbit view. The clip planes come from where the
// scene sphere lies along the view direction, never from fixed constants:
// near is half the gap to its front (or a thousandth of far when the camera is
// inside it), far is twice the distance to its back. Both are always positive
// and finite, with near < far.
Camera to_camera(const OrbitView& view, const Bounds& scene);

// fit_view, then to_camera on the same bounds.
Camera fit_camera(const Bounds& bounds, const ViewAngles& angles);

// Mouse controls. Orbit turns the eye about the target (pitch held to
// pitch_limit); pan moves the target so a point at the target's depth follows
// the cursor exactly (screen y grows DOWN, as mouse coordinates do); zoom
// scales the distance by 0.85 per wheel step, clamped. None of them changes
// anything else.
void orbit(OrbitView& view, Real dx_pixels, Real dy_pixels, Real radians_per_pixel);
void pan(OrbitView& view, Real dx_pixels, Real dy_pixels, Real viewport_height_pixels);
void zoom(OrbitView& view, Real wheel_steps, Real min_distance, Real max_distance);

// NDC x, NDC y, and view depth (distance along the view direction) of p. The
// tests and the viewer's pixel self-test use it; drawing does not.
Vec3 project(const Camera& camera, Real aspect, const Vec3& p);

// The world length that spans `pixels` on screen at `depth`, for a viewport
// `viewport_height_pixels` tall: how a viewer sizes a Point.
Real world_size_of_pixels(const Camera& camera, Real viewport_height_pixels, Real depth,
                          Real pixels);

// Fixed-timestep pacing for a real-time loop. Given an ABSOLUTE wall clock in
// seconds (a double, not summed float frame times), it says how many steps are
// due this frame: target = base_tick + floor((wall - base_wall) speed / dt).
// Steps are never resized or skipped, so the states are the terminal app's
// states whatever the frame rate. More than `max_steps` owed in one frame
// rebases the clock: the run falls behind (behind() is true) rather than
// repaying a stall later. Pause owes nothing; pausing, a speed change and
// rebase() restart the clock from the given moment, so nothing jumps. Speeds
// are clamped to [1/64, 64]; a non-positive or non-finite one is ignored. No
// step is ever due at or past `end_tick` (the scenario's `steps`).
class Pacer
{
public:
    Pacer(Real dt, std::uint64_t max_steps,
          std::uint64_t end_tick = std::numeric_limits<std::uint64_t>::max());
    std::uint64_t due(Real wall_now, std::uint64_t tick_now);
    void rebase(Real wall_now, std::uint64_t tick_now);
    void set_speed(Real speed, Real wall_now, std::uint64_t tick_now);
    void set_paused(bool paused, Real wall_now, std::uint64_t tick_now);
    Real speed() const { return speed_; }
    bool paused() const { return paused_; }
    bool behind() const { return behind_; }

private:
    Real dt_;
    std::uint64_t max_steps_;
    std::uint64_t end_tick_;
    Real speed_{1.0};
    bool paused_{false};
    bool behind_{false};
    Real base_wall_{0.0};
    std::uint64_t base_tick_{0};
};
} // namespace malloy::view3d
