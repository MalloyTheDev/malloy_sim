#include <malloy/ascii/grid.hpp>
#include <malloy/math/math.hpp>
#include <malloy/charges/charges.hpp>
#include <malloy/nbody/nbody.hpp>
#include <malloy/particles/particles.hpp>
#include <malloy/rigid/rigid.hpp>
#include <malloy/springs/springs.hpp>
#include <malloy/scenario/scenario.hpp>
#include <malloy/sim_core/sim_core.hpp>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

// M8 terminal demo: runs an N-body scenario and reports the conserved system
// diagnostics (separation, total energy, total angular momentum) together with
// an ASCII view of the body positions. The scenario is either loaded from a
// file given as the single argument, or one of the built-in scenarios when no
// argument is given. The app stays dumb -- it owns no physics, no rendering
// math, and has no CLI parser (docs/03, rules 14-15).

using malloy::ascii::fit_viewport;
using malloy::ascii::render;
using malloy::ascii::Viewport;
using malloy::math::distance;
using malloy::math::Real;
using malloy::math::Vec2;
using malloy::nbody::Body2D;
using malloy::nbody::NBodySettings;
using malloy::charges::Charge3DSettings;
using malloy::charges::Charge3DWorld;
using malloy::charges::ChargedParticle2D;
using malloy::charges::ChargedParticle3D;
using malloy::charges::ChargeSettings;
using malloy::charges::ChargeWorld;
using malloy::nbody::NBodyWorld;
using malloy::nbody::total_angular_momentum;
using malloy::particles::Particle2D;
using malloy::rigid::RigidBody2D;
using malloy::rigid::RigidSettings;
using malloy::rigid::RigidWorld;
using malloy::springs::SpringBody2D;
using malloy::springs::SpringBody3D;
using malloy::springs::SpringNetwork;
using malloy::springs::SpringWorld;
using malloy::springs::SpringWorld3D;
using malloy::nbody::total_energy;
using malloy::particles::total_energy;
using malloy::particles::total_momentum;
using malloy::particles::Particle3D;
using malloy::particles::ParticleSettings;
using malloy::particles::ParticleSettings3D;
using malloy::particles::ParticleWorld;
using malloy::particles::ParticleWorld3D;
using malloy::scenario::parse_scenario_file;
using malloy::scenario::ScenarioType;

namespace
{
// Three domains each define total_momentum / total_kinetic_energy in their
// own namespace. Aliasing here keeps the runners readable without pulling
// three same-named overload sets into one scope.
constexpr auto& rigid_linear_momentum = malloy::rigid::total_linear_momentum;
constexpr auto& rigid_angular_momentum = malloy::rigid::total_angular_momentum;
constexpr auto& rigid_total_energy = malloy::rigid::total_energy;
constexpr auto& charge_momentum = malloy::charges::total_momentum;
constexpr auto& charge_total_energy = malloy::charges::total_energy;
constexpr auto& charge_total_energy3d = malloy::charges::total_energy3d;
constexpr auto& spring_momentum = malloy::springs::total_momentum;
constexpr auto& spring_kinetic_energy = malloy::springs::total_kinetic_energy;
constexpr auto& spring_elastic_energy = malloy::springs::total_elastic_energy;
constexpr auto& spring_kinetic_energy3d = malloy::springs::total_kinetic_energy3d;
constexpr auto& spring_elastic_energy3d = malloy::springs::total_elastic_energy3d;
constexpr auto& particle_total_energy3d = malloy::particles::total_energy3d;
constexpr auto& particle_momentum3d = malloy::particles::total_momentum3d;
} // namespace
using malloy::scenario::ScenarioParseResult;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;

namespace
{
// Size of the ASCII debug view. Roughly 2:1 because terminal cells are about
// twice as tall as they are wide, so a square region of world space reads as
// square on screen.
constexpr int view_width = 61;
constexpr int view_height = 25;

// The positions of any bodies that have one, as the plain point list the ASCII
// view works with.
//
// At M10 this was two copies and a comment saying the duplication was three
// lines and the abstraction would be a template or a base class, so the copies
// won. It is now three copies, which is the threshold this project uses as
// evidence rather than coincidence (ADR 0008). A template over "has a position"
// is not a base class and invents no vocabulary, so rule 12 is untouched.
//
// malloy_rigid keeps its own extractor below, because it emits two points per
// body rather than one and is genuinely a different function.
template <typename Body>
std::vector<Vec2> positions_of(const std::vector<Body>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size());
    for (const Body& body : bodies)
    {
        points.push_back(body.position);
    }
    return points;
}

// Draws the current body positions, first expanding `view` so that every body
// is inside it. The view only ever grows: refitting it from scratch each frame
// would rescale the picture every report and make a stationary body appear to
// drift, while growing keeps the framing steady and still guarantees no body is
// silently clipped. The extents are printed so a change of scale is visible.
void print_view(Viewport& view, const std::vector<Vec2>& points)
{
    const Viewport frame = fit_viewport(points);

    view.min_x = std::min(view.min_x, frame.min_x);
    view.max_x = std::max(view.max_x, frame.max_x);
    view.min_y = std::min(view.min_y, frame.min_y);
    view.max_y = std::max(view.max_y, frame.max_y);

    std::cout << render(points, view, view_width, view_height) << "  view x ["
              << view.min_x << ", " << view.max_x << "]  y [" << view.min_y << ", "
              << view.max_y << "]\n";
}

// Drives one world to completion, reporting every output_every steps. Returns 0
// on success, 1 on step failure.
//
// One copy rather than four. It is a template over the world type, not a base
// class: RigidWorld, ParticleWorld, NBodyWorld and SpringWorld share nothing but
// a step() returning StepResult, and giving them a common base is exactly what
// CLAUDE.md rule 12 forbids. Duck typing here costs nothing and invents no
// vocabulary.
//
// The int64 counter is the reason this is worth factoring at all. With
// steps == INT_MAX an int counter reaches INT_MAX, passes the loop test, and
// overflows on ++, which is undefined behaviour and in practice never
// terminates. That fix was replicated in all four runners but explained in only
// one, so three copies carried a correctness-critical detail with no note
// saying why it mattered.
template <typename World, typename Report>
int drive(const char* title, World& world, const Report& report, int steps,
          int output_every)
{
    report(0);
    for (std::int64_t step = 1; step <= steps; ++step)
    {
        if (!world.step().ok())
        {
            std::cerr << title << ": step " << step << " failed\n";
            return 1;
        }
        if (output_every > 0 && step % output_every == 0)
        {
            report(step);
        }
    }
    return 0;
}

// Runs one particle scenario to completion, reporting the conserved quantities
// this domain actually has. Deliberately a separate function from
// run_scenario: ParticleWorld and NBodyWorld share no base class, and ADR 0006
// says a dispatch switch is the whole mechanism.
int run_particles(const char* title, const SimulationSettings& sim,
                  const ParticleSettings& settings, std::vector<Particle2D> particles,
                  int steps, int output_every)
{
    ParticleWorld world{sim, settings, std::move(particles)};

    std::cout << "\n== " << title << " ==  particles=" << world.particles().size()
              << "  dt=" << sim.dt << "  steps=" << steps
              << "  restitution=" << settings.restitution << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(positions_of(world.particles()));

    const auto report = [&world, &view, &settings](std::int64_t step_index) {
        const auto& now = world.particles();
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        std::cout << "   E " << std::setw(16) << total_energy(now, settings.gravity)
                  << "   p " << std::setw(16) << total_momentum(now).x << std::setw(16)
                  << total_momentum(now).y << '\n';
        std::cout << std::fixed;
        print_view(view, positions_of(now));
    };

    return drive(title, world, report, steps, output_every);
}

// The orthographic flattening for 3D particles: plot (x, y), drop z, the same
// honest limit as the other 3D runners until rendering has its own milestone.
std::vector<Vec2> projected(const std::vector<malloy::particles::Particle3D>& particles)
{
    std::vector<Vec2> points;
    points.reserve(particles.size());
    for (const auto& p : particles)
    {
        points.push_back(Vec2{p.position.x, p.position.y});
    }
    return points;
}

// Runs one 3D colliding-particle scenario. A tenth concrete runner and the
// tenth branch of the dispatch switch: ADR 0006 once more, still no base class.
int run_particles3d(const char* title, const SimulationSettings& sim,
                    const ParticleSettings3D& settings,
                    std::vector<Particle3D> particles, int steps, int output_every)
{
    ParticleWorld3D world{sim, settings, std::move(particles)};

    std::cout << "\n== " << title << " ==  particles=" << world.particles().size()
              << "  dt=" << sim.dt << "  steps=" << steps
              << "  restitution=" << settings.restitution << "  (3D)" << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(projected(world.particles()));

    const auto report = [&world, &view, &settings](std::int64_t step_index) {
        const auto& now = world.particles();
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        std::cout << "   E " << std::setw(16) << particle_total_energy3d(now, settings.gravity)
                  << "   |p| " << std::setw(16)
                  << malloy::math::length(particle_momentum3d(now)) << '\n';
        std::cout << std::fixed;
        print_view(view, projected(now));
    };

    return drive(title, world, report, steps, output_every);
}

// The points the ASCII view shows for a rigid body: its centre of mass and its
// body origin. Both are real state rather than an invented marker, and when the
// two are offset the origin visibly orbits the centre, which is how rotation
// becomes legible in a character grid.
std::vector<Vec2> rigid_points(const std::vector<RigidBody2D>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size() * 2);
    for (const RigidBody2D& body : bodies)
    {
        points.push_back(malloy::rigid::center_of_mass(body));
        points.push_back(body.position);
    }
    return points;
}

// Runs one rigid-body scenario. A third concrete runner: RigidWorld shares no
// base class with the other two, and the dispatch switch is the whole
// mechanism (ADR 0006).
int run_rigid(const char* title, const SimulationSettings& sim,
              std::vector<RigidBody2D> bodies, RigidSettings settings, int steps,
              int output_every)
{
    RigidWorld world{sim, std::move(bodies), settings};

    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << sim.dt << "  steps=" << steps << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(rigid_points(world.bodies()));

    const auto report = [&world, &view, &settings](std::int64_t step_index) {
        const auto& now = world.bodies();
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        std::cout << "   E " << std::setw(16) << rigid_total_energy(now, settings.gravity)
                  << "   L " << std::setw(16) << rigid_angular_momentum(now)
                  << "   p " << std::setw(16) << rigid_linear_momentum(now).x
                  << std::setw(16) << rigid_linear_momentum(now).y << '\n';
        std::cout << std::fixed;
        print_view(view, rigid_points(now));
    };

    return drive(title, world, report, steps, output_every);
}

// Flattens 3D positions onto the xy plane for the debug view.
//
// An orthographic projection along z, and nothing cleverer. malloy_ascii plots
// points that are already in the view plane, and drawing a 3D scene properly is
// a different problem rather than a wider one (ADR 0009). Turning this into a
// camera would make it a graphics project, which rule 2 exists to prevent until
// rendering has its own milestone.
//
// The consequence is worth knowing while reading a frame: motion along z is
// invisible here, so an orbit tilted out of the xy plane looks like an ellipse
// rather than a circle, and two bodies at different depths can overlap.
std::vector<Vec2> projected(const std::vector<malloy::nbody::Body3D>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size());
    for (const auto& body : bodies)
    {
        points.push_back(Vec2{body.position.x, body.position.y});
    }
    return points;
}

// Runs one 3D N-body scenario. A sixth concrete runner and a sixth branch of
// the dispatch switch, with no base class: ADR 0006 holds in three dimensions
// exactly as it does in two.
int run_nbody3d(const char* title, const SimulationSettings& sim,
                NBodySettings settings, std::vector<malloy::nbody::Body3D> bodies,
                int steps, int output_every)
{
    malloy::nbody::NBody3DWorld world{sim, settings, std::move(bodies)};

    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << sim.dt << "  steps=" << steps << "  (3D)" << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(projected(world.bodies()));

    const auto report = [&world, &view, &settings](std::int64_t step_index) {
        const auto& now = world.bodies();
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        // Angular momentum is a vector in three dimensions, and its magnitude
        // is the single number worth watching: it is exactly conserved, so any
        // movement in this column is rounding.
        std::cout << "   E " << std::setw(16)
                  << malloy::nbody::total_energy(now, settings.g, settings.softening)
                  << "   |L| " << std::setw(16)
                  << malloy::math::length(malloy::nbody::total_angular_momentum(now))
                  << "   |p| " << std::setw(16)
                  << malloy::math::length(malloy::nbody::total_momentum(now)) << '\n';
        std::cout << std::fixed;
        print_view(view, projected(now));
    };

    return drive(title, world, report, steps, output_every);
}

// The same orthographic flattening for rigid bodies. A second overload rather
// than a template, for the same reason the worlds are separate types: Body3D
// and RigidBody3D share a field name and nothing else.
std::vector<Vec2> projected(const std::vector<malloy::rigid::RigidBody3D>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size());
    for (const auto& body : bodies)
    {
        points.push_back(Vec2{body.position.x, body.position.y});
    }
    return points;
}

// Runs one 3D rigid-body scenario. A seventh concrete runner and a seventh
// branch of the dispatch switch: ADR 0006 again, and still no base class.
//
// Worth knowing while reading a frame: the view plots the CENTRES of the
// bodies, so a body tumbling in place does not move on screen at all. The
// rotation is in the diagnostics columns, not the picture, which is the honest
// limit of an ASCII view of a 3D scene until rendering has its own milestone.
int run_rigid3d(const char* title, const SimulationSettings& sim,
                std::vector<malloy::rigid::RigidBody3D> bodies,
                malloy::rigid::Rigid3DSettings settings, int steps,
                int output_every)
{
    malloy::rigid::Rigid3DWorld world{sim, std::move(bodies), settings};

    std::cout << "\n\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << sim.dt << "  steps=" << steps << "  (3D)" << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(projected(world.bodies()));

    const malloy::math::Vec3 gravity = world.settings().gravity;
    const auto report = [&world, &view, gravity](std::int64_t step_index) {
        const auto& now = world.bodies();
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        // E is kinetic plus gravitational potential, the quantity conserved in
        // free flight (semi-implicit Euler sheds (1/2)(sum m)|g|^2 dt^2 per
        // step, and a bounce with e<1 removes energy too). |L| is the total
        // angular momentum, not conserved once gravity or a contact acts.
        std::cout << "   E " << std::setw(16)
                  << (malloy::rigid::total_kinetic_energy3d(now) +
                      malloy::rigid::total_potential_energy3d(now, gravity))
                  << "   |L| " << std::setw(16)
                  << malloy::math::length(malloy::rigid::total_angular_momentum3d(now))
                  << "   |p| " << std::setw(16)
                  << malloy::math::length(malloy::rigid::total_linear_momentum3d(now))
                  << '\n';
        std::cout << std::fixed;
        print_view(view, projected(now));
    };

    return drive(title, world, report, steps, output_every);
}

// Runs one charged-particle scenario. A fifth concrete runner, and the fifth
// branch of the dispatch switch: no base class, exactly as ADR 0006 intends.
int run_charges(const char* title, const SimulationSettings& sim,
                ChargeSettings settings, std::vector<ChargedParticle2D> particles,
                int steps, int output_every)
{
    ChargeWorld world{sim, settings, std::move(particles)};

    std::cout << "\n== " << title << " ==  charges=" << world.particles().size()
              << "  dt=" << sim.dt << "  steps=" << steps << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(positions_of(world.particles()));

    const auto report = [&world, &view, &settings](std::int64_t step_index) {
        const auto& now = world.particles();
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        // Speed is reported because it is the quantity a magnetic field must
        // leave alone: it is the column that shows the rotation is exact.
        std::cout << "   E " << std::setw(16) << charge_total_energy(now, settings)
                  << "   |v| " << std::setw(16)
                  << malloy::math::length(now.front().velocity) << "   p "
                  << std::setw(16) << charge_momentum(now).x << std::setw(16)
                  << charge_momentum(now).y << '\n';
        std::cout << std::fixed;
        print_view(view, positions_of(now));
    };

    return drive(title, world, report, steps, output_every);
}

// The same orthographic flattening for 3D charges: plot (x, y) and drop z, the
// honest limit of an ASCII view of a 3D scene until rendering has its own
// milestone. A helical orbit reads as a circle from this angle, with the drift
// along the field carrying it off screen only if the field is not along z.
std::vector<Vec2> projected(const std::vector<malloy::charges::ChargedParticle3D>& particles)
{
    std::vector<Vec2> points;
    points.reserve(particles.size());
    for (const auto& particle : particles)
    {
        points.push_back(Vec2{particle.position.x, particle.position.y});
    }
    return points;
}

// Runs one 3D charged-particle scenario. An eighth concrete runner and an
// eighth branch of the dispatch switch: ADR 0006 once more, still no base class.
int run_charges3d(const char* title, const SimulationSettings& sim,
                  Charge3DSettings settings,
                  std::vector<ChargedParticle3D> particles, int steps,
                  int output_every)
{
    Charge3DWorld world{sim, settings, std::move(particles)};

    std::cout << "\n== " << title << " ==  charges=" << world.particles().size()
              << "  dt=" << sim.dt << "  steps=" << steps << "  (3D)" << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(projected(world.particles()));

    const auto report = [&world, &view, &settings](std::int64_t step_index) {
        const auto& now = world.particles();
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        // |v| is the column a magnetic field must leave alone, exactly as in 2D:
        // the rotation is exact, so speed does not move even as the particle
        // spirals along the field.
        std::cout << "   E " << std::setw(16) << charge_total_energy3d(now, settings)
                  << "   |v| " << std::setw(16)
                  << malloy::math::length(now.front().velocity) << "   |p| "
                  << std::setw(16)
                  << malloy::math::length(malloy::charges::total_momentum3d(now)) << '\n';
        std::cout << std::fixed;
        print_view(view, projected(now));
    };

    return drive(title, world, report, steps, output_every);
}

// Runs one spring scenario. A fourth concrete runner: SpringWorld shares no
// base class with the other three (ADR 0006, ADR 0008).
int run_springs(const char* title, const SimulationSettings& sim, SpringNetwork network,
                std::vector<SpringBody2D> bodies, int steps, int output_every)
{
    SpringWorld world{sim, std::move(network), std::move(bodies)};

    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  springs=" << world.network().size() << "  dt=" << sim.dt
              << "  steps=" << steps << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(positions_of(world.bodies()));

    const auto report = [&world, &view](std::int64_t step_index) {
        const auto& now = world.bodies();
        const auto kinetic = spring_kinetic_energy(now);
        const auto elastic = spring_elastic_energy(world.network(), now);
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        std::cout << "   KE " << std::setw(15) << kinetic << "   PE " << std::setw(15)
                  << elastic << "   E " << std::setw(15) << kinetic + elastic
                  << "   p " << std::setw(15) << spring_momentum(now).x << std::setw(15)
                  << spring_momentum(now).y << '\n';
        std::cout << std::fixed;
        print_view(view, positions_of(now));
    };

    return drive(title, world, report, steps, output_every);
}

// The orthographic flattening for 3D spring bodies: plot (x, y), drop z, the
// same honest limit as the other 3D runners until rendering has its milestone.
std::vector<Vec2> projected(const std::vector<malloy::springs::SpringBody3D>& bodies)
{
    std::vector<Vec2> points;
    points.reserve(bodies.size());
    for (const auto& body : bodies)
    {
        points.push_back(Vec2{body.position.x, body.position.y});
    }
    return points;
}

// Runs one 3D spring scenario. A ninth concrete runner and a ninth branch of
// the dispatch switch: ADR 0006 once more, still no base class.
int run_springs3d(const char* title, const SimulationSettings& sim,
                  SpringNetwork network, std::vector<SpringBody3D> bodies, int steps,
                  int output_every)
{
    SpringWorld3D world{sim, std::move(network), std::move(bodies)};

    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  springs=" << world.network().size() << "  dt=" << sim.dt
              << "  steps=" << steps << "  (3D)" << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    Viewport view = fit_viewport(projected(world.bodies()));

    const auto report = [&world, &view](std::int64_t step_index) {
        const auto& now = world.bodies();
        const auto kinetic = spring_kinetic_energy3d(now);
        const auto elastic = spring_elastic_energy3d(world.network(), now);
        std::cout << "step " << std::setw(6) << step_index;
        std::cout << std::scientific;
        std::cout << "   KE " << std::setw(15) << kinetic << "   PE " << std::setw(15)
                  << elastic << "   E " << std::setw(15) << kinetic + elastic
                  << "   |p| " << std::setw(15)
                  << malloy::math::length(malloy::springs::total_momentum3d(now)) << '\n';
        std::cout << std::fixed;
        print_view(view, projected(now));
    };

    return drive(title, world, report, steps, output_every);
}

// Runs one scenario to completion, printing system diagnostics every
// output_every steps. Returns 0 on success, 1 on validation/step failure.
int run_scenario(const char* title, const SimulationSettings& sim,
                 const NBodySettings& nbody, std::vector<Body2D> bodies, int steps,
                 int output_every)
{
    NBodyWorld world{sim, nbody, std::move(bodies)};

    std::cout << "\n== " << title << " ==  bodies=" << world.bodies().size()
              << "  dt=" << sim.dt << "  steps=" << steps << '\n';

    if (world.validate() != StepStatus::Ok)
    {
        std::cerr << title << ": invalid configuration\n";
        return 1;
    }

    // One view for the whole run, seeded from the initial state and grown as
    // needed, so successive frames share a scale and can be compared.
    Viewport view = fit_viewport(positions_of(world.bodies()));

    const auto report = [&world, &nbody, &view](std::int64_t step_index) {
        const auto& bodies_now = world.bodies();
        std::cout << "step " << std::setw(6) << step_index;
        if (bodies_now.size() >= 2)
        {
            std::cout << "   sep " << std::setw(11)
                      << distance(bodies_now[0].position, bodies_now[1].position);
        }
        // The conserved quantities are the whole point of this output, and for
        // the normalized two-body demo they are of order 1e-6. Under fixed(8)
        // every step printed the identical -0.00000050, which would look
        // perfectly conserved even if it were not. Scientific notation keeps
        // significant digits at any magnitude; separation stays fixed because
        // it is order 1 and reads better that way.
        std::cout << std::scientific;
        std::cout << "   E_total " << std::setw(16)
                  << total_energy(bodies_now, nbody.g, nbody.softening)
                  << "   L_total " << std::setw(16) << total_angular_momentum(bodies_now)
                  << '\n';
        std::cout << std::fixed;
        print_view(view, positions_of(bodies_now));
    };

    return drive(title, world, report, steps, output_every);
}

// The built-in scenarios, run when no scenario file is given.
int run_builtin_scenarios()
{
    const int rc = run_scenario(
        "two-body orbit (normalized units)", SimulationSettings{0.001},
        NBodySettings{1.0, 0.000001},
        {Body2D{Vec2{0.0, 0.0}, Vec2{0.0, 0.0}, 1.0},
         Body2D{Vec2{1.0, 0.0}, Vec2{0.0, 1.0}, 0.000001}},
        10000, 1000);
    if (rc != 0)
    {
        return rc;
    }

    const Real omega = 0.75984; // circular angular velocity for the unit triangle
    return run_scenario(
        "equilateral three-body (Lagrange)", SimulationSettings{0.001},
        NBodySettings{1.0, 0.001},
        {Body2D{Vec2{1.0, 0.0}, Vec2{0.0, omega}, 1.0},
         Body2D{Vec2{-0.5, 0.8660254}, Vec2{-0.8660254 * omega, -0.5 * omega}, 1.0},
         Body2D{Vec2{-0.5, -0.8660254}, Vec2{0.8660254 * omega, -0.5 * omega}, 1.0}},
        10000, 1000);
}
} // namespace

int main(int argc, char** argv)
{
    std::cout << "MalloySim nbody_terminal\n";
    std::cout << std::fixed << std::setprecision(8);

    // With a path argument, load and run that scenario file; otherwise run the
    // built-in scenarios. A single positional path -- no flags, no CLI parser.
    if (argc >= 2)
    {
        const ScenarioParseResult parsed = parse_scenario_file(argv[1]);
        if (!parsed.ok)
        {
            std::cerr << "scenario error: " << parsed.error << '\n';
            return 1;
        }
        const auto& s = parsed.scenario;
        // The whole multi-domain mechanism: one switch, one concrete runner per
        // domain, no base class (ADR 0006).
        switch (s.type)
        {
        case ScenarioType::NBody:
            return run_scenario(argv[1], s.simulation, s.nbody_settings, s.bodies,
                                s.steps, s.output_every);
        case ScenarioType::Particles:
            return run_particles(argv[1], s.simulation, s.particle_settings,
                                 s.particle_list, s.steps, s.output_every);
        case ScenarioType::Rigid:
            return run_rigid(argv[1], s.simulation, s.rigid_bodies,
                             s.rigid_settings, s.steps, s.output_every);
        case ScenarioType::Springs:
            return run_springs(argv[1], s.simulation, s.spring_network,
                               s.spring_bodies, s.steps, s.output_every);
        case ScenarioType::Charges:
            return run_charges(argv[1], s.simulation, s.charge_settings,
                               s.charge_list, s.steps, s.output_every);
        case ScenarioType::NBody3D:
            return run_nbody3d(argv[1], s.simulation, s.nbody_settings, s.bodies3d,
                               s.steps, s.output_every);
        case ScenarioType::Rigid3D:
            return run_rigid3d(argv[1], s.simulation, s.rigid_bodies3d,
                               s.rigid3d_settings, s.steps, s.output_every);
        case ScenarioType::Charges3D:
            return run_charges3d(argv[1], s.simulation, s.charge3d_settings,
                                 s.charge_list3d, s.steps, s.output_every);
        case ScenarioType::Springs3D:
            return run_springs3d(argv[1], s.simulation, s.spring_network3d,
                                 s.spring_bodies3d, s.steps, s.output_every);
        case ScenarioType::Particles3D:
            return run_particles3d(argv[1], s.simulation, s.particle3d_settings,
                                   s.particle_list3d, s.steps, s.output_every);
        }
        return 1;
    }

    return run_builtin_scenarios();
}
