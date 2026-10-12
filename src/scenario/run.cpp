#include <malloy/scenario/run.hpp>

#include <iomanip>
#include <ios>
#include <locale>
#include <sstream>

// Every value below is the expression the terminal app evaluated inline before
// M41, in the same order and association, and every line is the one it
// printed, so the app's output is unchanged to the last byte.
namespace malloy::scenario
{
namespace
{
// A wall's heat is reported only for a thermal wall; a specular wall's entry is
// 0 by definition and was never printed (M39).
std::optional<math::Real> heat_if_thermal(math::Real temperature, math::Real heat)
{
    if (temperature > math::Real{0})
    {
        return heat;
    }
    return std::nullopt;
}

// The first particle's speed, or empty for a world with none: reading front()
// of an empty list was undefined behaviour (issue #24).
template <typename Particles>
std::optional<math::Real> first_speed(const Particles& particles)
{
    if (particles.empty())
    {
        return std::nullopt;
    }
    return math::length(particles.front().velocity);
}

// The stream every report line is written to: the terminal app's own number
// format (fixed, 8 digits) in the classic locale, so a global locale with a
// different decimal point or digit grouping cannot change a byte.
std::ostringstream line_stream()
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::fixed << std::setprecision(8);
    return out;
}

void step_column(std::ostringstream& out, std::int64_t step)
{
    out << "step " << std::setw(6) << step;
}

// A column whose value may not exist yet: the value, or a dash in its place.
void optional_column(std::ostringstream& out, const std::optional<math::Real>& value)
{
    if (value)
    {
        out << std::setw(16) << *value;
    }
    else
    {
        out << std::setw(16) << "-";
    }
}

// A wall's heat column, printed only when the wall is thermal.
void heat_column(std::ostringstream& out, const char* wall,
                 const std::optional<math::Real>& heat)
{
    if (heat)
    {
        out << "   Q " << wall << std::setw(16) << *heat;
    }
}
} // namespace

// --- Building a world ---

nbody::NBodyWorld make_nbody_world(const Scenario& s)
{
    return nbody::NBodyWorld{s.simulation, s.nbody_settings, s.bodies};
}

nbody::NBody3DWorld make_nbody3d_world(const Scenario& s)
{
    return nbody::NBody3DWorld{s.simulation, s.nbody_settings, s.bodies3d};
}

particles::ParticleWorld make_particle_world(const Scenario& s)
{
    return particles::ParticleWorld{s.simulation, s.particle_settings, s.particle_list};
}

particles::ParticleWorld3D make_particle3d_world(const Scenario& s)
{
    return particles::ParticleWorld3D{s.simulation, s.particle3d_settings,
                                      s.particle_list3d};
}

rigid::RigidWorld make_rigid_world(const Scenario& s)
{
    return rigid::RigidWorld{s.simulation, s.rigid_bodies, s.rigid_settings};
}

rigid::Rigid3DWorld make_rigid3d_world(const Scenario& s)
{
    return rigid::Rigid3DWorld{s.simulation, s.rigid_bodies3d, s.rigid3d_settings};
}

springs::SpringWorld make_spring_world(const Scenario& s)
{
    return springs::SpringWorld{s.simulation, s.spring_network, s.spring_bodies};
}

springs::SpringWorld3D make_spring3d_world(const Scenario& s)
{
    return springs::SpringWorld3D{s.simulation, s.spring_network3d, s.spring_bodies3d};
}

charges::ChargeWorld make_charge_world(const Scenario& s)
{
    return charges::ChargeWorld{s.simulation, s.charge_settings, s.charge_list};
}

charges::Charge3DWorld make_charge3d_world(const Scenario& s)
{
    return charges::Charge3DWorld{s.simulation, s.charge3d_settings, s.charge_list3d};
}

bool has_thermal_wall(const particles::ParticleSettings& settings)
{
    const auto& walls = settings.wall_temperature;
    return walls.min_x > 0.0 || walls.max_x > 0.0 || walls.min_y > 0.0 ||
           walls.max_y > 0.0;
}

bool has_thermal_wall(const particles::ParticleSettings3D& settings)
{
    const auto& walls = settings.wall_temperature;
    return walls.min_x > 0.0 || walls.max_x > 0.0 || walls.min_y > 0.0 ||
           walls.max_y > 0.0 || walls.min_z > 0.0 || walls.max_z > 0.0;
}

// --- Reading back ---

NBodyReadout readout(const nbody::NBodyWorld& world)
{
    const auto& bodies = world.bodies();
    const auto& settings = world.nbody_settings();
    NBodyReadout r;
    if (bodies.size() >= 2)
    {
        r.separation = math::distance(bodies[0].position, bodies[1].position);
    }
    r.energy = nbody::total_energy(bodies, settings.g, settings.softening);
    r.angular_momentum = nbody::total_angular_momentum(bodies);
    return r;
}

NBody3DReadout readout(const nbody::NBody3DWorld& world)
{
    const auto& bodies = world.bodies();
    const auto& settings = world.nbody_settings();
    NBody3DReadout r;
    r.energy = nbody::total_energy(bodies, settings.g, settings.softening);
    r.angular_momentum = nbody::total_angular_momentum(bodies);
    r.momentum = nbody::total_momentum(bodies);
    return r;
}

ParticleReadout readout(const particles::ParticleWorld& world)
{
    const auto& particles = world.particles();
    const auto& settings = world.particle_settings();
    const auto& walls = settings.wall_temperature;
    const auto& heat = world.wall_heat();
    ParticleReadout r;
    r.energy = particles::total_energy(particles, settings.gravity);
    r.momentum = particles::total_momentum(particles);
    r.temperature = particles::temperature(particles);
    r.pressure =
        particles::pressure(world.wall_impulse(), settings.bounds, world.elapsed_time());
    r.heat.min_x = heat_if_thermal(walls.min_x, heat.min_x);
    r.heat.max_x = heat_if_thermal(walls.max_x, heat.max_x);
    r.heat.min_y = heat_if_thermal(walls.min_y, heat.min_y);
    r.heat.max_y = heat_if_thermal(walls.max_y, heat.max_y);
    return r;
}

Particle3DReadout readout(const particles::ParticleWorld3D& world)
{
    const auto& particles = world.particles();
    const auto& settings = world.particle_settings();
    const auto& walls = settings.wall_temperature;
    const auto& heat = world.wall_heat();
    Particle3DReadout r;
    r.energy = particles::total_energy3d(particles, settings.gravity);
    r.momentum = particles::total_momentum3d(particles);
    r.temperature = particles::temperature3d(particles);
    r.pressure = particles::pressure3d(world.wall_impulse(), settings.bounds,
                                       world.elapsed_time());
    r.heat.min_x = heat_if_thermal(walls.min_x, heat.min_x);
    r.heat.max_x = heat_if_thermal(walls.max_x, heat.max_x);
    r.heat.min_y = heat_if_thermal(walls.min_y, heat.min_y);
    r.heat.max_y = heat_if_thermal(walls.max_y, heat.max_y);
    r.heat.min_z = heat_if_thermal(walls.min_z, heat.min_z);
    r.heat.max_z = heat_if_thermal(walls.max_z, heat.max_z);
    return r;
}

RigidReadout readout(const rigid::RigidWorld& world)
{
    const auto& bodies = world.bodies();
    RigidReadout r;
    r.energy = rigid::total_energy(bodies, world.settings().gravity);
    r.angular_momentum = rigid::total_angular_momentum(bodies);
    r.momentum = rigid::total_linear_momentum(bodies);
    return r;
}

Rigid3DReadout readout(const rigid::Rigid3DWorld& world)
{
    const auto& bodies = world.bodies();
    Rigid3DReadout r;
    // (kinetic + gravity) + force: the association the app has always used.
    r.energy = rigid::total_kinetic_energy3d(bodies) +
               rigid::total_potential_energy3d(bodies, world.settings().gravity) +
               rigid::total_force_potential3d(bodies, world.settings().force);
    r.angular_momentum = rigid::total_angular_momentum3d(bodies);
    r.momentum = rigid::total_linear_momentum3d(bodies);
    return r;
}

SpringReadout readout(const springs::SpringWorld& world)
{
    const auto& bodies = world.bodies();
    SpringReadout r;
    r.kinetic = springs::total_kinetic_energy(bodies);
    r.elastic = springs::total_elastic_energy(world.network(), bodies);
    r.energy = r.kinetic + r.elastic;
    r.momentum = springs::total_momentum(bodies);
    return r;
}

Spring3DReadout readout(const springs::SpringWorld3D& world)
{
    const auto& bodies = world.bodies();
    Spring3DReadout r;
    r.kinetic = springs::total_kinetic_energy3d(bodies);
    r.elastic = springs::total_elastic_energy3d(world.network(), bodies);
    r.energy = r.kinetic + r.elastic;
    r.momentum = springs::total_momentum3d(bodies);
    return r;
}

ChargeReadout readout(const charges::ChargeWorld& world)
{
    const auto& particles = world.particles();
    ChargeReadout r;
    r.energy = charges::total_energy(particles, world.charge_settings());
    r.first_speed = first_speed(particles);
    r.momentum = charges::total_momentum(particles);
    return r;
}

Charge3DReadout readout(const charges::Charge3DWorld& world)
{
    const auto& particles = world.particles();
    Charge3DReadout r;
    r.energy = charges::total_energy3d(particles, world.charge_settings());
    r.first_speed = first_speed(particles);
    r.momentum = charges::total_momentum3d(particles);
    return r;
}

// --- The report line ---

std::string report_line(std::int64_t step, const NBodyReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    if (r.separation)
    {
        // Fixed, because it is order 1 and reads better that way. The
        // conserved quantities after it are of order 1e-6 in the normalized
        // demo, where fixed(8) would print the same digits every step and look
        // perfectly conserved even if they were not, so they are scientific.
        out << "   sep " << std::setw(11) << *r.separation;
    }
    out << std::scientific;
    out << "   E_total " << std::setw(16) << r.energy << "   L_total " << std::setw(16)
        << r.angular_momentum;
    return out.str();
}

std::string report_line(std::int64_t step, const NBody3DReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    // Angular momentum is a vector in three dimensions, and its magnitude is
    // the single number worth watching: it is exactly conserved, so any
    // movement in this column is rounding.
    out << "   E " << std::setw(16) << r.energy << "   |L| " << std::setw(16)
        << math::length(r.angular_momentum) << "   |p| " << std::setw(16)
        << math::length(r.momentum);
    return out.str();
}

std::string report_line(std::int64_t step, const ParticleReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    out << "   E " << std::setw(16) << r.energy << "   p " << std::setw(16) << r.momentum.x
        << std::setw(16) << r.momentum.y << "   T " << std::setw(16) << r.temperature
        << "   P ";
    optional_column(out, r.pressure);
    heat_column(out, "min_x", r.heat.min_x);
    heat_column(out, "max_x", r.heat.max_x);
    heat_column(out, "min_y", r.heat.min_y);
    heat_column(out, "max_y", r.heat.max_y);
    return out.str();
}

std::string report_line(std::int64_t step, const Particle3DReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    out << "   E " << std::setw(16) << r.energy << "   |p| " << std::setw(16)
        << math::length(r.momentum) << "   T " << std::setw(16) << r.temperature
        << "   P ";
    optional_column(out, r.pressure);
    heat_column(out, "min_x", r.heat.min_x);
    heat_column(out, "max_x", r.heat.max_x);
    heat_column(out, "min_y", r.heat.min_y);
    heat_column(out, "max_y", r.heat.max_y);
    heat_column(out, "min_z", r.heat.min_z);
    heat_column(out, "max_z", r.heat.max_z);
    return out.str();
}

std::string report_line(std::int64_t step, const RigidReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    out << "   E " << std::setw(16) << r.energy << "   L " << std::setw(16)
        << r.angular_momentum << "   p " << std::setw(16) << r.momentum.x << std::setw(16)
        << r.momentum.y;
    return out.str();
}

std::string report_line(std::int64_t step, const Rigid3DReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    // |L| is the total angular momentum, not conserved once gravity, a force or
    // a contact acts.
    out << "   E " << std::setw(16) << r.energy << "   |L| " << std::setw(16)
        << math::length(r.angular_momentum) << "   |p| " << std::setw(16)
        << math::length(r.momentum);
    return out.str();
}

std::string report_line(std::int64_t step, const SpringReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    out << "   KE " << std::setw(15) << r.kinetic << "   PE " << std::setw(15) << r.elastic
        << "   E " << std::setw(15) << r.energy << "   p " << std::setw(15) << r.momentum.x
        << std::setw(15) << r.momentum.y;
    return out.str();
}

std::string report_line(std::int64_t step, const Spring3DReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    out << "   KE " << std::setw(15) << r.kinetic << "   PE " << std::setw(15) << r.elastic
        << "   E " << std::setw(15) << r.energy << "   |p| " << std::setw(15)
        << math::length(r.momentum);
    return out.str();
}

std::string report_line(std::int64_t step, const ChargeReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    // Speed is reported because it is the quantity a magnetic field must leave
    // alone: it is the column that shows the rotation is exact.
    out << "   E " << std::setw(16) << r.energy << "   |v| ";
    optional_column(out, r.first_speed);
    out << "   p " << std::setw(16) << r.momentum.x << std::setw(16) << r.momentum.y;
    return out.str();
}

std::string report_line(std::int64_t step, const Charge3DReadout& r)
{
    std::ostringstream out = line_stream();
    step_column(out, step);
    out << std::scientific;
    out << "   E " << std::setw(16) << r.energy << "   |v| ";
    optional_column(out, r.first_speed);
    out << "   |p| " << std::setw(16) << math::length(r.momentum);
    return out.str();
}
} // namespace malloy::scenario
