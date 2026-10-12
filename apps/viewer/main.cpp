// malloy_viewer (M41, issue #80): the five 3D scenario types, drawn in real time.
//
// A client of the libraries and nothing more (rules 14 and 15). The run layer
// builds the world and formats the report line, malloy_view3d turns the world
// into a frame, frames it with a camera and paces the steps, and window.cpp
// draws what it is handed. The physics advances only through world.step(), one
// whole fixed step at a time, so every state the viewer shows is a state the
// terminal app reaches, and its report lines are the terminal app's lines.
//
//   malloy_viewer <scenario.scn> [--paused] [--stop-at N] [--frames N]
//                 [--frame-time S] [--headless | --self-test]
//
//   --headless     no window and no OpenGL: the same loop on a synthetic clock
//                  of --frame-time seconds a frame (default 1/60), printing the
//                  report line at every report step, then exiting at the end
//   --frame-time S the clock advances S seconds a frame instead of following
//                  the wall clock (S may be written 1/60)
//   --stop-at N    end the run at tick N (or at the scenario's steps, if fewer)
//   --frames N     exit after N frames
//   --paused       start paused
//   --self-test    draw, read the pixels back and check them against the render
//                  model, run, reset, and check the reset picture is the first
//                  one; exit 0 on a pass, 1 on a failure, 77 with no OpenGL
//
// Keys: Space pause, N or Right one step, R restart, ] and [ speed x2 and /2,
// F refit, A auto-fit on or off, T trails, Esc quit. Mouse: left drag orbits,
// right drag or shift and left drag pans, the wheel zooms.

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <vector>

#include <malloy/scenario/run.hpp>
#include <malloy/scenario/scenario.hpp>
#include <malloy/view3d/view3d.hpp>

#include "picture.hpp"
#include "self_test.hpp"
#include "window.hpp"

namespace
{
namespace view3d = malloy::view3d;
namespace run = malloy::scenario;
using malloy::scenario::Scenario;
using malloy::scenario::ScenarioType;

constexpr int kExitNoGl = 77; // the automake/CTest "skipped" convention

struct Options
{
    std::string path;
    bool headless{false};
    bool self_test{false};
    bool paused{false};
    std::optional<double> frame_time;
    std::optional<std::uint64_t> stop_at;
    std::optional<std::uint64_t> frames;
};

void print_usage(std::ostream& out)
{
    out << "usage: malloy_viewer <scenario.scn> [--paused] [--stop-at N] [--frames N]\n"
           "                     [--frame-time S] [--headless | --self-test]\n";
}

std::optional<std::uint64_t> parse_count(std::string_view text)
{
    std::uint64_t value = 0;
    const auto [end, ec] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (ec != std::errc{} || end != text.data() + text.size())
    {
        return std::nullopt;
    }
    return value;
}

std::optional<double> parse_real(std::string_view text)
{
    double value = 0;
    const auto [end, ec] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (ec != std::errc{} || end != text.data() + text.size())
    {
        return std::nullopt;
    }
    return value;
}

// A positive, finite number of seconds, written as a decimal or as a/b.
std::optional<double> parse_seconds(std::string_view text)
{
    const std::size_t slash = text.find('/');
    std::optional<double> value;
    if (slash == std::string_view::npos)
    {
        value = parse_real(text);
    }
    else
    {
        const auto top = parse_real(text.substr(0, slash));
        const auto bottom = parse_real(text.substr(slash + 1));
        if (top && bottom && *bottom != 0)
        {
            value = *top / *bottom;
        }
    }
    if (!value || !std::isfinite(*value) || !(*value > 0))
    {
        return std::nullopt;
    }
    return value;
}

// The hand-rolled parser the design allows: one positional path and a few flags.
std::optional<Options> parse_options(int argc, char** argv)
{
    Options options;
    bool have_path = false;
    for (int i = 1; i < argc; ++i)
    {
        const std::string_view arg{argv[i]};
        const auto value = [&]() -> std::optional<std::string_view> {
            if (i + 1 >= argc)
            {
                std::cerr << "malloy_viewer: " << arg << " needs a value\n";
                return std::nullopt;
            }
            return std::string_view{argv[++i]};
        };
        if (arg == "--headless")
        {
            options.headless = true;
        }
        else if (arg == "--self-test")
        {
            options.self_test = true;
        }
        else if (arg == "--paused")
        {
            options.paused = true;
        }
        else if (arg == "--frame-time" || arg == "--stop-at" || arg == "--frames")
        {
            const auto text = value();
            if (!text)
            {
                return std::nullopt;
            }
            if (arg == "--frame-time")
            {
                options.frame_time = parse_seconds(*text);
                if (!options.frame_time)
                {
                    std::cerr << "malloy_viewer: --frame-time needs a positive number of "
                                 "seconds, such as 0.02 or 1/60\n";
                    return std::nullopt;
                }
            }
            else
            {
                const auto count = parse_count(*text);
                if (!count)
                {
                    std::cerr << "malloy_viewer: " << arg << " needs a whole number\n";
                    return std::nullopt;
                }
                (arg == "--stop-at" ? options.stop_at : options.frames) = count;
            }
        }
        else if (arg.size() >= 2 && arg.substr(0, 2) == "--")
        {
            std::cerr << "malloy_viewer: unknown option " << arg << '\n';
            return std::nullopt;
        }
        else if (have_path)
        {
            std::cerr << "malloy_viewer: one scenario path only\n";
            return std::nullopt;
        }
        else
        {
            options.path = std::string{arg};
            have_path = true;
        }
    }
    if (!have_path)
    {
        std::cerr << "malloy_viewer: no scenario file given\n";
        return std::nullopt;
    }
    if (options.headless && options.self_test)
    {
        std::cerr << "malloy_viewer: --headless and --self-test exclude each other "
                     "(the self-test needs a window)\n";
        return std::nullopt;
    }
    return options;
}

const char* type_name(ScenarioType type)
{
    switch (type)
    {
    case ScenarioType::NBody:
        return "nbody";
    case ScenarioType::Particles:
        return "particles";
    case ScenarioType::Rigid:
        return "rigid";
    case ScenarioType::Springs:
        return "springs";
    case ScenarioType::Charges:
        return "charges";
    case ScenarioType::NBody3D:
        return "nbody3d";
    case ScenarioType::Rigid3D:
        return "rigid3d";
    case ScenarioType::Charges3D:
        return "charges3d";
    case ScenarioType::Springs3D:
        return "springs3d";
    case ScenarioType::Particles3D:
        return "particles3d";
    }
    return "unknown";
}

std::string file_name(const std::string& path)
{
    const std::size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

// The bounding sphere of what a frame draws but does not frame (ground,
// walls, field arrows), so the clip planes never cut the context off. It goes
// through the render model's own fit_bounds by presenting the context as
// points; the camera FIT still ignores it.
std::optional<view3d::Bounds> context_bounds(const view3d::Frame& frame)
{
    view3d::Frame context;
    const auto add = [&context](const view3d::Vec3& p) {
        context.points.push_back(view3d::Point{p, 0.0, 0, view3d::Role::Mass});
    };
    for (const view3d::Segment& s : frame.segments)
    {
        if (!view3d::is_subject(s.role))
        {
            add(s.a);
            add(s.b);
        }
    }
    for (const view3d::Patch& p : frame.patches)
    {
        for (const view3d::Vec3& c : p.corners)
        {
            add(c);
        }
    }
    for (const view3d::Arrow& a : frame.arrows)
    {
        add(a.origin);
        add(a.tip);
    }
    if (context.points.empty())
    {
        return std::nullopt;
    }
    return view3d::fit_bounds(context);
}

std::string speed_label(int exponent)
{
    std::ostringstream out;
    if (exponent >= 0)
    {
        out << 'x' << (1 << exponent);
    }
    else
    {
        out << "x1/" << (1 << -exponent);
    }
    return out.str();
}

// One run of one concrete world type. A template over the world, not a base
// class: the five 3D worlds share nothing but step(), validate() and the
// run-layer and render-model overloads (ADR 0006, rule 12).
template <typename World, typename Make>
class Session
{
public:
    Session(const Options& options, const Scenario& scenario, Make make)
        : options_{options},
          scenario_{scenario},
          make_{make},
          end_tick_{std::min<std::uint64_t>(
              static_cast<std::uint64_t>(std::max(scenario.steps, 0)),
              options.stop_at.value_or(UINT64_MAX))},
          dt_{scenario.simulation.dt},
          pacer_{dt_, frame_cap(64), end_tick_},
          trail_every_{std::max<std::uint64_t>(
              1, (static_cast<std::uint64_t>(std::max(scenario.steps, 0)) + 2047) / 2048)}
    {
    }

    int run()
    {
        world_.emplace(make_(scenario_));
        if (world_->validate() != malloy::sim_core::StepStatus::Ok)
        {
            std::cerr << options_.path << ": invalid configuration\n";
            return 1;
        }
        return options_.headless ? run_headless() : run_window();
    }

private:
    // Points and charges leave trails; the other types have nothing to trace.
    static constexpr bool kTraced =
        std::is_same_v<World, malloy::nbody::NBody3DWorld> ||
        std::is_same_v<World, malloy::charges::Charge3DWorld>;

    // At most a quarter of a second of simulated time at `speed` in one frame,
    // so a stall (a dragged window, a breakpoint) is not repaid in a burst.
    std::uint64_t frame_cap(double speed) const
    {
        const double cap = std::floor(0.25 * speed / dt_);
        if (!std::isfinite(cap) || cap < 1.0)
        {
            return 1;
        }
        return cap > 1e15 ? static_cast<std::uint64_t>(1e15) : static_cast<std::uint64_t>(cap);
    }

    void report(std::uint64_t tick) const
    {
        std::cout << run::report_line(static_cast<std::int64_t>(tick), run::readout(*world_))
                  << '\n';
    }

    // Back to tick 0 by building the world again from the parsed scenario,
    // never from bodies(): the rigid warm-start impulses (ADR 0011) and the
    // particle generator (ADR 0010) are world state too.
    void restart(double now)
    {
        world_.emplace(make_(scenario_));
        tick_ = 0;
        failed_ = false;
        trails_.clear();
        pacer_.rebase(now, 0);
        report(0);
        bounds_ = view3d::fit_bounds(view3d::build_frame(*world_, view3d::Bounds{}));
        if constexpr (kTraced)
        {
            trails_.record(view3d::build_frame(*world_, bounds_));
        }
    }

    bool step_once()
    {
        if (!world_->step().ok())
        {
            std::cerr << options_.path << ": step " << (tick_ + 1) << " failed\n";
            failed_ = true;
            return false;
        }
        ++tick_;
        if (scenario_.output_every > 0 &&
            tick_ % static_cast<std::uint64_t>(scenario_.output_every) == 0)
        {
            report(tick_);
        }
        if constexpr (kTraced)
        {
            if (tick_ % trail_every_ == 0)
            {
                trails_.record(view3d::build_frame(*world_, bounds_));
            }
        }
        return true;
    }

    // The steps due at `now`, run one at a time. `budget` is a cap on the real
    // seconds spent stepping in this frame (none on a synthetic clock, so a
    // headless run is the same on any machine). Whatever cannot be run now is
    // dropped by rebasing the pacer: the run falls behind instead of bursting.
    void advance(double now, std::optional<double> budget)
    {
        std::uint64_t due = pacer_.due(now, tick_);
        bool cut = false;
        const std::uint64_t cap = frame_cap(pacer_.speed());
        if (due > cap)
        {
            due = cap;
            cut = true;
        }
        const auto start = std::chrono::steady_clock::now();
        for (std::uint64_t k = 0; k < due; ++k)
        {
            if (budget && k > 0 &&
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() >
                    *budget)
            {
                cut = true;
                break;
            }
            if (!step_once())
            {
                return;
            }
        }
        if (cut)
        {
            pacer_.rebase(now, tick_);
        }
        if (cut || pacer_.behind())
        {
            behind_until_ = now + 0.5;
        }
    }

    int run_headless()
    {
        const double frame_time = options_.frame_time.value_or(1.0 / 60.0);
        pacer_.set_paused(options_.paused, 0.0, 0);
        restart(0.0);
        std::uint64_t frames = 0;
        while (!failed_ && tick_ < end_tick_ && !pacer_.paused())
        {
            if (options_.frames && frames >= *options_.frames)
            {
                break;
            }
            ++frames;
            advance(static_cast<double>(frames) * frame_time, std::nullopt);
            // The same frame work as a windowed run, minus the drawing, so a
            // headless line also says that building frames changes nothing.
            const view3d::Frame frame = view3d::build_frame(*world_, bounds_);
            bounds_ = view3d::merge(bounds_, view3d::fit_bounds(frame));
        }
        return failed_ ? 1 : 0;
    }

    view3d::ViewAngles angles(double aspect) const
    {
        view3d::ViewAngles a;
        a.yaw = view_.yaw;
        a.pitch = view_.pitch;
        a.aspect = aspect;
        return a;
    }

    std::string hud(double fps) const
    {
        std::ostringstream out;
        out << file_name(options_.path) << "  " << type_name(scenario_.type) << "  tick "
            << tick_ << " / " << end_tick_ << "  t " << std::fixed << std::setprecision(4)
            << static_cast<double>(tick_) * dt_ << "  " << speed_label(speed_exponent_) << "  ";
        if (failed_)
        {
            out << "FAILED at step " << (tick_ + 1);
        }
        else if (tick_ >= end_tick_)
        {
            out << "FINISHED";
        }
        else
        {
            out << (pacer_.paused() ? "PAUSED" : "RUNNING");
        }
        out << "  auto-fit " << (auto_fit_ ? "on" : "off");
        if (behind_)
        {
            out << "  BEHIND";
        }
        out << "  " << std::setprecision(0) << fps << " fps";
        return out.str();
    }

    int run_window()
    {
        viewer::WindowOptions window;
        if (options_.self_test)
        {
            // Exact colours for the pixel checks: no multisampling, and no
            // waiting for the display.
            window.width = 960;
            window.height = 540;
            window.msaa = false;
            window.vsync = false;
            window.resizable = false;
        }
        if (!viewer::open_window("MalloySim viewer - " + file_name(options_.path), window))
        {
            std::cerr << "malloy_viewer: no OpenGL 3.3 context could be created "
                         "(no display, or no OpenGL 3.3 driver)\n";
            return options_.self_test ? kExitNoGl : 1;
        }
        const int code = options_.self_test ? self_test() : interactive();
        viewer::close_window();
        return code;
    }

    double aspect() const
    {
        const int h = viewer::render_height();
        return h > 0 ? static_cast<double>(viewer::render_width()) / h : 16.0 / 9.0;
    }

    // Tick 0 in the default view, as at start and after R.
    void start_view(double now)
    {
        restart(now);
        view_ = view3d::OrbitView{};
        view_ = view3d::fit_view(bounds_, angles(aspect()));
    }

    // The frame for the current state, with the camera following it under
    // auto-fit, drawn into the window. `pixels` asks for the picture back.
    void draw(double fps, bool subjects_only, std::vector<std::uint8_t>* pixels,
              view3d::Frame* shown, viewer::Look* look_out)
    {
        const double a = aspect();
        view3d::Frame frame = view3d::build_frame(*world_, bounds_);
        bounds_ = view3d::merge(bounds_, view3d::fit_bounds(frame));
        if (auto_fit_)
        {
            view_ = view3d::fit_view(bounds_, angles(a));
        }
        view3d::Bounds scene = bounds_;
        if (const auto context = context_bounds(frame))
        {
            scene = view3d::merge(scene, *context);
        }
        viewer::Look look;
        look.camera = view3d::to_camera(view_, scene);
        look.aspect = a;
        look.viewport_height = static_cast<double>(viewer::render_height());
        painter_.paint(frame, look, trails_on_ ? &trails_ : nullptr, subjects_only, list_);
        viewer::present(viewer::camera_setup(look.camera), list_,
                        subjects_only ? std::string{} : hud(fps), viewer::background_color,
                        pixels);
        if (shown != nullptr)
        {
            *shown = std::move(frame);
        }
        if (look_out != nullptr)
        {
            *look_out = look;
        }
    }

    double clock(std::uint64_t frame) const
    {
        return options_.frame_time ? static_cast<double>(frame) * *options_.frame_time
                                   : viewer::wall_seconds();
    }

    int interactive()
    {
        std::cout << "keys: Space pause, N or Right step, R restart, ] [ speed, F refit, "
                     "A auto-fit, T trails, Esc quit; mouse: left drag orbit, right or "
                     "shift drag pan, wheel zoom\n";
        pacer_.set_paused(options_.paused, clock(0), 0);
        start_view(clock(0));
        std::uint64_t frames = 0;
        double fps = 0.0;
        double last = viewer::wall_seconds();
        while (true)
        {
            if (options_.frames && frames >= *options_.frames)
            {
                break;
            }
            const viewer::Input in = viewer::poll_input();
            if (in.quit)
            {
                break;
            }
            ++frames;
            const double now = clock(frames);

            if (in.reset)
            {
                start_view(now);
            }
            if (in.toggle_pause)
            {
                pacer_.set_paused(!pacer_.paused(), now, tick_);
            }
            if (in.faster || in.slower)
            {
                speed_exponent_ = std::clamp(speed_exponent_ + (in.faster ? 1 : -1), -6, 6);
                pacer_.set_speed(std::ldexp(1.0, speed_exponent_), now, tick_);
            }
            if (in.single_step && pacer_.paused() && !failed_ && tick_ < end_tick_)
            {
                step_once();
            }
            if (in.refit)
            {
                bounds_ = view3d::fit_bounds(view3d::build_frame(*world_, bounds_));
                view_ = view3d::fit_view(bounds_, angles(aspect()));
            }
            if (in.toggle_auto_fit)
            {
                auto_fit_ = !auto_fit_;
            }
            if (in.toggle_trails)
            {
                trails_on_ = !trails_on_;
            }
            view3d::orbit(view_, in.orbit_dx, in.orbit_dy, 0.005);
            if (in.pan_dx != 0.0 || in.pan_dy != 0.0)
            {
                auto_fit_ = false; // the user has taken the camera
                view3d::pan(view_, in.pan_dx, in.pan_dy, viewer::screen_height());
            }
            if (in.wheel != 0.0)
            {
                auto_fit_ = false;
                view3d::zoom(view_, in.wheel, 1e-3 * bounds_.radius, 1e3 * bounds_.radius);
            }

            if (!failed_)
            {
                advance(now, options_.frame_time ? std::nullopt : std::optional<double>{0.012});
            }
            behind_ = now < behind_until_;

            const double wall = viewer::wall_seconds();
            const double elapsed = wall - last;
            last = wall;
            if (elapsed > 0)
            {
                fps = fps == 0.0 ? 1.0 / elapsed : 0.9 * fps + 0.1 / elapsed;
            }
            draw(fps, false, nullptr, nullptr, nullptr);
        }
        return failed_ ? 1 : 0;
    }

    // Draws the subjects, reads them back and checks them; reports and
    // returns false on a failure.
    bool checked_frame(const char* when, std::vector<std::uint8_t>& pixels,
                       viewer::PixelReport& total)
    {
        view3d::Frame frame;
        viewer::Look look;
        draw(0.0, true, &pixels, &frame, &look);
        const viewer::PixelReport r =
            viewer::check_pixels(pixels, viewer::render_width(), viewer::render_height(),
                                 viewer::background_color, frame, look);
        if (!r.ok())
        {
            std::cerr << "self-test: " << when << " (tick " << tick_ << "): " << r.failure
                      << '\n';
            return false;
        }
        total.drawn_pixels += r.drawn_pixels;
        total.subject_probes += r.subject_probes;
        total.clear_probes += r.clear_probes;
        return true;
    }

    int self_test()
    {
        if (!options_.frame_time)
        {
            options_.frame_time = 1.0 / 60.0; // a synthetic clock: the same on any machine
        }
        viewer::PixelReport total;
        std::vector<std::uint8_t> first;
        std::vector<std::uint8_t> pixels;
        pacer_.set_paused(options_.paused, clock(0), 0);
        start_view(clock(0));
        if (!checked_frame("first frame", first, total))
        {
            return 1;
        }
        std::uint64_t frames = 0;
        std::uint64_t checked = 1;
        while (!failed_ && tick_ < end_tick_ && !pacer_.paused() &&
               !(options_.frames && frames >= *options_.frames))
        {
            ++frames;
            advance(clock(frames), std::nullopt);
            if (frames % 15 == 0)
            {
                if (!checked_frame("running", pixels, total))
                {
                    return 1;
                }
                ++checked;
            }
            else
            {
                draw(0.0, true, nullptr, nullptr, nullptr);
            }
        }
        if (failed_)
        {
            return 1;
        }
        if (!checked_frame("last frame", pixels, total))
        {
            return 1;
        }
        ++checked;
        const std::uint64_t last_tick = tick_;

        // Restart: tick 0 again, and the very same picture as the first frame.
        start_view(clock(frames + 1));
        if (!checked_frame("after restart", pixels, total))
        {
            return 1;
        }
        ++checked;
        if (pixels != first)
        {
            std::size_t differing = 0;
            for (std::size_t i = 0; i < std::min(pixels.size(), first.size()); ++i)
            {
                differing += pixels[i] != first[i] ? 1U : 0U;
            }
            std::cerr << "self-test: the picture after restart differs from the first one ("
                      << differing << " bytes)\n";
            return 1;
        }
        std::cout << "self-test passed: " << checked << " frames checked, " << frames
                  << " frames run to tick " << last_tick << ", " << total.subject_probes
                  << " subject probes, " << total.clear_probes << " clear probes, "
                  << viewer::render_width() << "x" << viewer::render_height()
                  << ", restart picture identical\n";
        return 0;
    }

    Options options_;
    const Scenario& scenario_;
    Make make_;
    std::uint64_t end_tick_;
    double dt_;
    view3d::Pacer pacer_;
    std::uint64_t trail_every_;
    std::optional<World> world_;
    std::uint64_t tick_{0};
    bool failed_{false};
    view3d::Bounds bounds_{};
    view3d::OrbitView view_{};
    bool auto_fit_{true};
    bool trails_on_{true};
    int speed_exponent_{0};
    double behind_until_{-1.0};
    bool behind_{false};
    viewer::Trails trails_;
    viewer::Painter painter_;
    viewer::DrawList list_;
};

template <typename Make>
int run_session(const Options& options, const Scenario& scenario, Make make)
{
    using World = decltype(make(scenario));
    Session<World, Make> session{options, scenario, make};
    return session.run();
}
} // namespace

int main(int argc, char** argv)
{
    std::cout << "MalloySim viewer\n";
    const std::optional<Options> options = parse_options(argc, argv);
    if (!options)
    {
        print_usage(std::cerr);
        return 1;
    }

    const run::ScenarioParseResult parsed = run::parse_scenario_file(options->path);
    if (!parsed.ok)
    {
        std::cerr << "scenario error: " << parsed.error << '\n';
        return 1;
    }
    const Scenario& s = parsed.scenario;
    // As the terminal app: a file that defines nothing is refused, not run
    // (issue #31). Checked before the type, since an empty file parses as the
    // default nbody type.
    if (run::entity_count(s) == 0)
    {
        std::cerr << "scenario error: the file defines nothing to simulate "
                     "(no bodies, particles, charges or springs)\n";
        return 1;
    }

    // One switch, one concrete world per type (ADR 0006). The 2D types are not
    // drawn in M41; the terminal app runs them.
    switch (s.type)
    {
    case ScenarioType::NBody3D:
        return run_session(*options, s, run::make_nbody3d_world);
    case ScenarioType::Rigid3D:
        return run_session(*options, s, run::make_rigid3d_world);
    case ScenarioType::Charges3D:
        return run_session(*options, s, run::make_charge3d_world);
    case ScenarioType::Springs3D:
        return run_session(*options, s, run::make_spring3d_world);
    case ScenarioType::Particles3D:
        return run_session(*options, s, run::make_particle3d_world);
    case ScenarioType::NBody:
    case ScenarioType::Particles:
    case ScenarioType::Rigid:
    case ScenarioType::Springs:
    case ScenarioType::Charges:
        std::cerr << "malloy_viewer: " << options->path << " is a 2D scenario (type "
                  << type_name(s.type)
                  << "); the viewer draws the 3D types only (nbody3d, particles3d, rigid3d, "
                     "springs3d, charges3d). Run it with malloy_nbody_terminal.\n";
        return 1;
    }
    return 1;
}
