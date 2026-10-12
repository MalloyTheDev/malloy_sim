#include "window.hpp"

#include <cstddef>
#include <cstring>

// The one file in the repository that includes a graphics header. raylib's
// include directory is a SYSTEM one (CMakeLists.txt), so its headers are
// compiled under /external:W0 and -isystem and the project's own warnings,
// errors included, stay on for this file.
#include <raylib.h>
#include <rlgl.h>

namespace viewer
{
namespace
{
Color color_of(const Rgb& c)
{
    return Color{c.r, c.g, c.b, 255};
}
} // namespace

bool open_window(const std::string& title, const WindowOptions& options)
{
    // raylib logs at INFO level to stdout by default; the viewer's stdout
    // carries the report lines, so only warnings and errors are let through.
    SetTraceLogLevel(LOG_WARNING);
    unsigned int flags = 0;
    if (options.msaa)
    {
        flags |= FLAG_MSAA_4X_HINT;
    }
    if (options.vsync)
    {
        flags |= FLAG_VSYNC_HINT;
    }
    if (options.resizable)
    {
        flags |= FLAG_WINDOW_RESIZABLE;
    }
    SetConfigFlags(flags);
    InitWindow(options.width, options.height, title.c_str());
    // InitWindow returns early, leaving the window not ready, when the
    // platform layer cannot make a window or an OpenGL 3.3 context.
    if (!IsWindowReady())
    {
        return false;
    }
    SetExitKey(KEY_ESCAPE);
    return true;
}

void close_window()
{
    CloseWindow();
}

double wall_seconds()
{
    return GetTime();
}

Input poll_input()
{
    Input in;
    in.quit = WindowShouldClose(); // Esc or the window's close button
    in.toggle_pause = IsKeyPressed(KEY_SPACE);
    in.single_step = IsKeyPressed(KEY_N) || IsKeyPressedRepeat(KEY_N) ||
                     IsKeyPressed(KEY_RIGHT) || IsKeyPressedRepeat(KEY_RIGHT);
    in.reset = IsKeyPressed(KEY_R);
    in.faster = IsKeyPressed(KEY_RIGHT_BRACKET);
    in.slower = IsKeyPressed(KEY_LEFT_BRACKET);
    in.refit = IsKeyPressed(KEY_F);
    in.toggle_auto_fit = IsKeyPressed(KEY_A);
    in.toggle_trails = IsKeyPressed(KEY_T);

    const bool shift = IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT);
    const Vector2 delta = GetMouseDelta();
    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) || (shift && IsMouseButtonDown(MOUSE_BUTTON_LEFT)))
    {
        in.pan_dx = delta.x;
        in.pan_dy = delta.y;
    }
    else if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
    {
        in.orbit_dx = delta.x;
        in.orbit_dy = delta.y;
    }
    in.wheel = GetMouseWheelMove();
    return in;
}

int render_width()
{
    return GetRenderWidth();
}

int render_height()
{
    return GetRenderHeight();
}

int screen_height()
{
    return GetScreenHeight();
}

void present(const CameraSetup& setup, const DrawList& list, const std::string& hud,
             Rgb background, std::vector<std::uint8_t>* pixels)
{
    BeginDrawing();
    ClearBackground(color_of(background));

    // raylib's perspective takes its clip planes from rlgl's cull distances,
    // which default to fixed constants; the render model's planes replace them
    // every frame, before BeginMode3D builds the projection.
    rlSetClipPlanes(setup.near_plane, setup.far_plane);
    Camera3D camera{};
    camera.position = Vector3{setup.position[0], setup.position[1], setup.position[2]};
    camera.target = Vector3{0.0F, 0.0F, 0.0F};
    camera.up = Vector3{setup.up[0], setup.up[1], setup.up[2]};
    camera.fovy = static_cast<float>(setup.fovy_degrees);
    camera.projection = CAMERA_PERSPECTIVE;
    BeginMode3D(camera);

    // One batch per primitive kind: rlVertex3f flushes the batch itself when
    // it fills, at a whole triangle or line.
    const std::size_t triangle_count = list.triangle_colors.size();
    if (triangle_count > 0)
    {
        rlBegin(RL_TRIANGLES);
        for (std::size_t i = 0; i < triangle_count; ++i)
        {
            const Rgb& c = list.triangle_colors[i];
            rlColor4ub(c.r, c.g, c.b, 255);
            const float* v = &list.triangles[9 * i];
            rlVertex3f(v[0], v[1], v[2]);
            rlVertex3f(v[3], v[4], v[5]);
            rlVertex3f(v[6], v[7], v[8]);
        }
        rlEnd();
    }
    const std::size_t line_count = list.line_colors.size();
    if (line_count > 0)
    {
        rlBegin(RL_LINES);
        for (std::size_t i = 0; i < line_count; ++i)
        {
            const Rgb& c = list.line_colors[i];
            rlColor4ub(c.r, c.g, c.b, 255);
            const float* v = &list.lines[6 * i];
            rlVertex3f(v[0], v[1], v[2]);
            rlVertex3f(v[3], v[4], v[5]);
        }
        rlEnd();
    }
    EndMode3D();

    if (!hud.empty())
    {
        DrawText(hud.c_str(), 10, 10, 20, Color{230, 232, 238, 255});
    }

    if (pixels != nullptr)
    {
        // Everything queued so far is drawn first, then the back buffer is read
        // before EndDrawing swaps it away.
        rlDrawRenderBatchActive();
        Image image = LoadImageFromScreen();
        const std::size_t bytes = static_cast<std::size_t>(image.width) *
                                  static_cast<std::size_t>(image.height) * 4U;
        pixels->resize(bytes);
        if (image.data != nullptr && bytes > 0)
        {
            std::memcpy(pixels->data(), image.data, bytes);
        }
        UnloadImage(image);
    }
    EndDrawing();
}
} // namespace viewer
