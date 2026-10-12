#pragma once

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include <malloy/view3d/mesh.hpp>
#include <malloy/view3d/view3d.hpp>

#include "window.hpp"

// Turning a render-model frame into the triangles and lines the window draws:
// which colour each role gets, the headlight shade of each face, and the float
// conversion relative to the camera target. Every position, mesh and size comes
// from malloy_view3d; nothing here computes geometry.
namespace viewer
{
// How the frame is being looked at.
struct Look
{
    malloy::view3d::Camera camera;
    double aspect{16.0 / 9.0};
    double viewport_height{720.0}; // render pixels
    double point_pixels{10.0};     // on-screen diameter of a Point
};

// The world radius a Point is drawn with, so it spans look.point_pixels on
// screen at its depth (zero when it is behind the camera). The self-test uses
// the same number as the drawing.
double point_radius(const Look& look, const malloy::view3d::Vec3& position);

// The past positions of the point masses and charges, so an orbit or a helix
// shows its path. Sampled by the stepping loop at fixed step intervals, never
// per drawn frame, so a trail does not depend on the frame rate.
class Trails
{
public:
    explicit Trails(std::size_t capacity = 2048) : capacity_{capacity} {}
    void clear() { paths_.clear(); }
    void record(const malloy::view3d::Frame& frame);
    // The kept positions of each traced point, oldest first, with its role.
    struct Path
    {
        malloy::view3d::Role role{malloy::view3d::Role::Mass};
        std::vector<malloy::view3d::Vec3> positions;
    };
    const std::map<std::size_t, Path>& paths() const { return paths_; }

private:
    std::size_t capacity_;
    std::map<std::size_t, Path> paths_;
};

class Painter
{
public:
    // Fills `out` with `frame` as seen through `look`. With `subjects_only`,
    // only spheres, boxes and points are drawn, flat on the background, which
    // is what the pixel self-test reads.
    void paint(const malloy::view3d::Frame& frame, const Look& look, const Trails* trails,
               bool subjects_only, DrawList& out);

private:
    const std::vector<malloy::view3d::Triangle>& unit_sphere(const malloy::view3d::SphereDetail& d);

    std::map<std::pair<int, int>, std::vector<malloy::view3d::Triangle>> units_;
    std::vector<malloy::view3d::Triangle> placed_;
};

// The window's camera for a render-model camera: positions relative to the
// target, fovy in degrees.
CameraSetup camera_setup(const malloy::view3d::Camera& camera);

inline constexpr Rgb background_color{22, 24, 30};
} // namespace viewer
