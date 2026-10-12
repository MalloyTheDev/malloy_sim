#pragma once

#include <cstdint>
#include <string>
#include <vector>

// The window: the only part of the viewer that talks to raylib (window.cpp is
// the one file that includes raylib.h and rlgl.h). Everything here is plain
// data, so the rest of the viewer, and the render model behind it, never sees a
// graphics type. Geometry arrives already final, in float, relative to the
// camera target; this layer draws it and reads the input, and computes nothing.
namespace viewer
{
struct Rgb
{
    std::uint8_t r{0};
    std::uint8_t g{0};
    std::uint8_t b{0};
};

// One frame's worth of flat-shaded triangles and coloured lines, positions
// relative to the camera target (so float keeps its precision far from the
// world origin).
struct DrawList
{
    std::vector<float> triangles; // 9 floats per triangle, counter-clockwise from its front
    std::vector<Rgb> triangle_colors;
    std::vector<float> lines; // 6 floats per line
    std::vector<Rgb> line_colors;

    void clear()
    {
        triangles.clear();
        triangle_colors.clear();
        lines.clear();
        line_colors.clear();
    }
};

// The drawing camera, relative to its target, which is therefore the origin.
// fovy is in DEGREES here, which is what raylib takes; the render model's
// camera is in radians.
struct CameraSetup
{
    float position[3]{0.0F, -3.0F, 0.0F};
    float up[3]{0.0F, 0.0F, 1.0F};
    double fovy_degrees{45.0};
    double near_plane{0.01};
    double far_plane{100.0};
};

// What the user did since the last frame. Mouse motion is in screen pixels,
// with y growing downward, as the render model's pan expects.
struct Input
{
    bool quit{false};
    bool toggle_pause{false};
    bool single_step{false};
    bool reset{false};
    bool faster{false};
    bool slower{false};
    bool refit{false};
    bool toggle_auto_fit{false};
    bool toggle_trails{false};
    double orbit_dx{0.0};
    double orbit_dy{0.0};
    double pan_dx{0.0};
    double pan_dy{0.0};
    double wheel{0.0};
};

struct WindowOptions
{
    int width{1280};
    int height{720};
    bool msaa{true};
    bool vsync{true};
    bool resizable{true};
};

// Opens the window and its OpenGL context. False when no context can be made
// (no display, or no OpenGL 3.3 driver); nothing else is then valid to call.
bool open_window(const std::string& title, const WindowOptions& options);
void close_window();

// Seconds since the window opened, from the platform's monotonic clock.
double wall_seconds();

Input poll_input();

// The framebuffer size (which HiDPI scaling may make larger than the window),
// and the window height in screen pixels, which mouse motion is measured in.
int render_width();
int render_height();
int screen_height();

// Draws one frame and presents it. When `pixels` is given, the back buffer is
// read into it first (RGBA, top row first, render_width x render_height), so
// what is read is exactly what is about to be shown.
void present(const CameraSetup& camera, const DrawList& list, const std::string& hud,
             Rgb background, std::vector<std::uint8_t>* pixels);
} // namespace viewer
