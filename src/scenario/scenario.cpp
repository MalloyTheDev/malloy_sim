#include <malloy/scenario/scenario.hpp>

#include <cstddef>
#include <fstream>
#include <limits>
#include <istream>
#include <sstream>
#include <string>

#include <malloy/collide/shapes.hpp>
#include <malloy/collide/shapes3d.hpp>
#include <malloy/math/vec2.hpp>
#include <malloy/rigid/mass_properties3d.hpp>
#include <malloy/rigid/rigid_body2d.hpp>
#include <malloy/springs/springs.hpp>

namespace malloy::scenario
{
namespace
{
ScenarioParseResult make_error(std::size_t line_number, const std::string& message)
{
    return ScenarioParseResult{false,
                               "line " + std::to_string(line_number) + ": " + message,
                               {}};
}

// True when nothing but whitespace remains. A directive that parsed its fields
// and left tokens behind is a malformed line, not a valid one: silently
// dropping the remainder would accept a body line written with 3D fields and
// run a different simulation than the file describes.
bool at_end_of_directive(std::istringstream& tokens)
{
    tokens >> std::ws;
    return tokens.eof();
}
} // namespace

ScenarioParseResult parse_scenario(std::istream& input)
{
    Scenario scenario;
    std::string line;
    std::size_t line_number = 0;
    // `type` selects which domain keys are legal, so it has to be settled
    // before any of them appear.
    bool saw_domain_key = false;

    while (std::getline(input, line))
    {
        ++line_number;

        const std::size_t comment = line.find('#');
        if (comment != std::string::npos)
        {
            line.erase(comment);
        }

        std::istringstream tokens(line);
        std::string key;
        if (!(tokens >> key))
        {
            continue; // blank or comment-only line
        }

        if (key == "type")
        {
            std::string value;
            if (!(tokens >> value))
            {
                return make_error(line_number,
                                  "type requires nbody, particles, rigid, springs, "
                                  "charges, nbody3d, rigid3d or charges3d");
            }
            if (value == "nbody")
            {
                scenario.type = ScenarioType::NBody;
            }
            else if (value == "particles")
            {
                scenario.type = ScenarioType::Particles;
            }
            else if (value == "rigid")
            {
                scenario.type = ScenarioType::Rigid;
            }
            else if (value == "springs")
            {
                scenario.type = ScenarioType::Springs;
            }
            else if (value == "charges")
            {
                scenario.type = ScenarioType::Charges;
            }
            else if (value == "nbody3d")
            {
                scenario.type = ScenarioType::NBody3D;
            }
            else if (value == "rigid3d")
            {
                scenario.type = ScenarioType::Rigid3D;
            }
            else if (value == "charges3d")
            {
                scenario.type = ScenarioType::Charges3D;
            }
            else
            {
                return make_error(line_number, "unknown type '" + value + "'");
            }
            if (saw_domain_key)
            {
                return make_error(line_number,
                                  "type must come before any domain-specific key");
            }
        }
        else if (key == "body")
        {
            math::Real mass{};
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (scenario.type != ScenarioType::NBody)
            {
                return make_error(line_number, "body belongs to type nbody");
            }
            saw_domain_key = true;
            if (!(tokens >> mass >> px >> py >> vx >> vy))
            {
                return make_error(line_number, "body requires: mass px py vx vy");
            }
            scenario.bodies.push_back(
                nbody::Body2D{math::Vec2{px, py}, math::Vec2{vx, vy}, mass});
        }
        else if (key == "dt")
        {
            if (!(tokens >> scenario.simulation.dt))
            {
                return make_error(line_number, "dt requires a number");
            }
        }
        else if (key == "g")
        {
            // Two dimensions or three: G means the same thing in either, and
            // both use NBodySettings.
            if (scenario.type != ScenarioType::NBody &&
                scenario.type != ScenarioType::NBody3D)
            {
                return make_error(line_number, "g belongs to type nbody or nbody3d");
            }
            saw_domain_key = true;
            if (!(tokens >> scenario.nbody_settings.g))
            {
                return make_error(line_number, "g requires a number");
            }
        }
        else if (key == "softening")
        {
            // Two domains now, for the same reason: both have a 1/r^2 pair
            // term that goes to infinity for a coincident pair.
            if (scenario.type != ScenarioType::NBody &&
                scenario.type != ScenarioType::NBody3D &&
                scenario.type != ScenarioType::Charges &&
                scenario.type != ScenarioType::Charges3D)
            {
                return make_error(
                    line_number,
                    "softening belongs to type nbody, nbody3d, charges or charges3d");
            }
            saw_domain_key = true;
            math::Real value{};
            if (!(tokens >> value))
            {
                return make_error(line_number, "softening requires a number");
            }
            if (scenario.type == ScenarioType::NBody ||
                scenario.type == ScenarioType::NBody3D)
            {
                scenario.nbody_settings.softening = value;
            }
            else if (scenario.type == ScenarioType::Charges3D)
            {
                scenario.charge3d_settings.softening = value;
            }
            else
            {
                scenario.charge_settings.softening = value;
            }
        }
        else if (key == "steps")
        {
            if (!(tokens >> scenario.steps))
            {
                return make_error(line_number, "steps requires an integer");
            }
            if (scenario.steps < 0)
            {
                return make_error(line_number, "steps must be >= 0");
            }
        }
        else if (key == "output_every")
        {
            if (!(tokens >> scenario.output_every))
            {
                return make_error(line_number, "output_every requires an integer");
            }
            if (scenario.output_every < 0)
            {
                return make_error(line_number, "output_every must be >= 0");
            }
        }
        else if (key == "restitution")
        {
            if (scenario.type != ScenarioType::Particles &&
                scenario.type != ScenarioType::Rigid &&
                scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(
                    line_number,
                    "restitution belongs to type particles, rigid or rigid3d");
            }
            saw_domain_key = true;
            math::Real value{};
            if (!(tokens >> value))
            {
                return make_error(line_number, "restitution requires a number");
            }
            // Write only the field the declared domain actually reads. Setting
            // both would leave a rigid scenario carrying a particle setting it
            // never uses, which reads as a mistake to anyone inspecting it.
            if (scenario.type == ScenarioType::Rigid3D)
            {
                scenario.rigid3d_settings.restitution = value;
            }
            else if (scenario.type == ScenarioType::Particles)
            {
                scenario.particle_settings.restitution = value;
            }
            else
            {
                scenario.rigid_settings.restitution = value;
            }
        }
        else if (key == "gravity")
        {
            if (scenario.type != ScenarioType::Particles &&
                scenario.type != ScenarioType::Rigid)
            {
                return make_error(line_number,
                                  "gravity belongs to type particles or rigid");
            }
            saw_domain_key = true;
            math::Real gx{};
            math::Real gy{};
            if (!(tokens >> gx >> gy))
            {
                return make_error(line_number, "gravity requires: gx gy");
            }
            if (scenario.type == ScenarioType::Particles)
            {
                scenario.particle_settings.gravity = math::Vec2{gx, gy};
            }
            else
            {
                scenario.rigid_settings.gravity = math::Vec2{gx, gy};
            }
        }
        else if (key == "body3")
        {
            if (scenario.type != ScenarioType::NBody3D)
            {
                return make_error(line_number, "body3 belongs to type nbody3d");
            }
            saw_domain_key = true;
            nbody::Body3D body;
            math::Real px{};
            math::Real py{};
            math::Real pz{};
            math::Real vx{};
            math::Real vy{};
            math::Real vz{};
            if (!(tokens >> body.mass >> px >> py >> pz >> vx >> vy >> vz))
            {
                return make_error(line_number,
                                  "body3 requires: mass px py pz vx vy vz");
            }
            body.position = math::Vec3{px, py, pz};
            body.velocity = math::Vec3{vx, vy, vz};
            scenario.bodies3d.push_back(body);
        }
        else if (key == "rigid_body3d")
        {
            if (scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number,
                                  "rigid_body3d belongs to type rigid3d");
            }
            saw_domain_key = true;
            rigid::RigidBody3D body;
            math::Real ix{};
            math::Real iy{};
            math::Real iz{};
            math::Real radius{};
            math::Real px{};
            math::Real py{};
            math::Real pz{};
            math::Real ax{};
            math::Real ay{};
            math::Real az{};
            math::Real angle{};
            math::Real vx{};
            math::Real vy{};
            math::Real vz{};
            math::Real wx{};
            math::Real wy{};
            math::Real wz{};
            if (!(tokens >> body.mass >> ix >> iy >> iz >> radius >> px >> py >>
                  pz >> ax >> ay >> az >> angle >> vx >> vy >> vz >> wx >> wy >>
                  wz))
            {
                return make_error(line_number,
                                  "rigid_body3d requires: mass Ix Iy Iz radius "
                                  "px py pz axisx axisy axisz angle vx vy vz "
                                  "wx wy wz");
            }
            body.inertia = math::Vec3{ix, iy, iz};
            body.radius = radius;
            body.position = math::Vec3{px, py, pz};
            body.velocity = math::Vec3{vx, vy, vz};
            body.angular_velocity = math::Vec3{wx, wy, wz};
            // Always a unit quaternion, which is why the file gives an axis and
            // an angle rather than four components.
            body.orientation = math::from_axis_angle(math::Vec3{ax, ay, az}, angle);
            scenario.rigid_bodies3d.push_back(body);
        }
        else if (key == "rigid_box3d")
        {
            if (scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number, "rigid_box3d belongs to type rigid3d");
            }
            saw_domain_key = true;
            math::Real mass{};
            math::Real hx{};
            math::Real hy{};
            math::Real hz{};
            math::Real px{};
            math::Real py{};
            math::Real pz{};
            math::Real ax{};
            math::Real ay{};
            math::Real az{};
            math::Real angle{};
            math::Real vx{};
            math::Real vy{};
            math::Real vz{};
            math::Real wx{};
            math::Real wy{};
            math::Real wz{};
            if (!(tokens >> mass >> hx >> hy >> hz >> px >> py >> pz >> ax >> ay >>
                  az >> angle >> vx >> vy >> vz >> wx >> wy >> wz))
            {
                return make_error(line_number,
                                  "rigid_box3d requires: mass hx hy hz px py pz "
                                  "axisx axisy axisz angle vx vy vz wx wy wz");
            }
            // The inertia is COMPUTED from the box's geometry and mass (M26),
            // not typed in, so the collision shape and the inertia always
            // describe the one box. Density follows from mass and volume; an
            // axis-aligned box's principal frame is the identity, so the moments
            // drop straight in. A degenerate box gives a zero-mass result whose
            // invalid body the domain's validate() then rejects.
            const math::Vec3 half_extents{hx, hy, hz};
            const math::Real volume = math::Real{8} * hx * hy * hz;
            const math::Real density = mass / volume;
            const rigid::MassProperties3D properties = rigid::mass_properties_3d(
                rigid::SolidBox{math::Vec3{}, half_extents, math::Quat{}, density});
            rigid::RigidBody3D body;
            body.mass = mass;
            body.inertia = properties.inertia;
            body.half_extents = half_extents;
            body.position = math::Vec3{px, py, pz};
            body.velocity = math::Vec3{vx, vy, vz};
            body.angular_velocity = math::Vec3{wx, wy, wz};
            body.orientation = math::from_axis_angle(math::Vec3{ax, ay, az}, angle);
            scenario.rigid_bodies3d.push_back(body);
        }
        else if (key == "torque")
        {
            if (scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number, "torque belongs to type rigid3d");
            }
            saw_domain_key = true;
            math::Real tx{};
            math::Real ty{};
            math::Real tz{};
            if (!(tokens >> tx >> ty >> tz))
            {
                return make_error(line_number, "torque requires: tx ty tz");
            }
            scenario.rigid3d_settings.torque = math::Vec3{tx, ty, tz};
        }
        else if (key == "gravity3")
        {
            if (scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number, "gravity3 belongs to type rigid3d");
            }
            saw_domain_key = true;
            math::Real gx{};
            math::Real gy{};
            math::Real gz{};
            if (!(tokens >> gx >> gy >> gz))
            {
                return make_error(line_number, "gravity3 requires: gx gy gz");
            }
            scenario.rigid3d_settings.gravity = math::Vec3{gx, gy, gz};
        }
        else if (key == "plane3")
        {
            if (scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number, "plane3 belongs to type rigid3d");
            }
            saw_domain_key = true;
            math::Real nx{};
            math::Real ny{};
            math::Real nz{};
            math::Real offset{};
            if (!(tokens >> nx >> ny >> nz >> offset))
            {
                return make_error(line_number, "plane3 requires: nx ny nz offset");
            }
            // Normalized here, at the input boundary, so collide::Plane3 keeps a
            // strict unit-normal invariant rather than normalizing on every
            // query. A direction that cannot be normalized is refused rather
            // than quietly turned into one.
            const math::Vec3 direction{nx, ny, nz};
            const math::Real magnitude = math::length(direction);
            if (!(magnitude > math::Real{0}) || !math::is_finite(magnitude) ||
                !math::is_finite(offset))
            {
                return make_error(line_number,
                                  "plane3 normal must be non-zero and finite");
            }
            scenario.rigid3d_settings.ground.push_back(
                collide::Plane3{direction / magnitude, offset});
        }
        else if (key == "coulomb")
        {
            if (scenario.type != ScenarioType::Charges &&
                scenario.type != ScenarioType::Charges3D)
            {
                return make_error(line_number,
                                  "coulomb belongs to type charges or charges3d");
            }
            saw_domain_key = true;
            math::Real& k = scenario.type == ScenarioType::Charges3D
                                ? scenario.charge3d_settings.k
                                : scenario.charge_settings.k;
            if (!(tokens >> k))
            {
                return make_error(line_number, "coulomb requires a value");
            }
        }
        else if (key == "efield")
        {
            if (scenario.type != ScenarioType::Charges)
            {
                return make_error(line_number, "efield belongs to type charges");
            }
            saw_domain_key = true;
            math::Real ex{};
            math::Real ey{};
            if (!(tokens >> ex >> ey))
            {
                return make_error(line_number, "efield requires: ex ey");
            }
            scenario.charge_settings.electric = math::Vec2{ex, ey};
        }
        else if (key == "bfield")
        {
            if (scenario.type != ScenarioType::Charges)
            {
                return make_error(line_number, "bfield belongs to type charges");
            }
            saw_domain_key = true;
            if (!(tokens >> scenario.charge_settings.magnetic))
            {
                return make_error(line_number, "bfield requires a value");
            }
        }
        else if (key == "charge")
        {
            if (scenario.type != ScenarioType::Charges)
            {
                return make_error(line_number, "charge belongs to type charges");
            }
            saw_domain_key = true;
            charges::ChargedParticle2D particle;
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (!(tokens >> particle.mass >> particle.charge >> px >> py >> vx >> vy))
            {
                return make_error(line_number,
                                  "charge requires: mass q px py vx vy");
            }
            particle.position = math::Vec2{px, py};
            particle.velocity = math::Vec2{vx, vy};
            scenario.charge_list.push_back(particle);
        }
        else if (key == "efield3")
        {
            if (scenario.type != ScenarioType::Charges3D)
            {
                return make_error(line_number, "efield3 belongs to type charges3d");
            }
            saw_domain_key = true;
            math::Real ex{};
            math::Real ey{};
            math::Real ez{};
            if (!(tokens >> ex >> ey >> ez))
            {
                return make_error(line_number, "efield3 requires: ex ey ez");
            }
            scenario.charge3d_settings.electric = math::Vec3{ex, ey, ez};
        }
        else if (key == "bfield3")
        {
            if (scenario.type != ScenarioType::Charges3D)
            {
                return make_error(line_number, "bfield3 belongs to type charges3d");
            }
            saw_domain_key = true;
            math::Real bx{};
            math::Real by{};
            math::Real bz{};
            if (!(tokens >> bx >> by >> bz))
            {
                return make_error(line_number, "bfield3 requires: bx by bz");
            }
            scenario.charge3d_settings.magnetic = math::Vec3{bx, by, bz};
        }
        else if (key == "charge3")
        {
            if (scenario.type != ScenarioType::Charges3D)
            {
                return make_error(line_number, "charge3 belongs to type charges3d");
            }
            saw_domain_key = true;
            charges::ChargedParticle3D particle;
            math::Real px{};
            math::Real py{};
            math::Real pz{};
            math::Real vx{};
            math::Real vy{};
            math::Real vz{};
            if (!(tokens >> particle.mass >> particle.charge >> px >> py >> pz >> vx >>
                  vy >> vz))
            {
                return make_error(line_number,
                                  "charge3 requires: mass q px py pz vx vy vz");
            }
            particle.position = math::Vec3{px, py, pz};
            particle.velocity = math::Vec3{vx, vy, vz};
            scenario.charge_list3d.push_back(particle);
        }
        else if (key == "friction")
        {
            if (scenario.type != ScenarioType::Rigid &&
                scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number,
                                  "friction belongs to type rigid or rigid3d");
            }
            saw_domain_key = true;
            math::Real value{};
            if (!(tokens >> value))
            {
                return make_error(line_number, "friction requires a value");
            }
            if (scenario.type == ScenarioType::Rigid3D)
            {
                scenario.rigid3d_settings.friction = value;
            }
            else
            {
                scenario.rigid_settings.friction = value;
            }
        }
        else if (key == "ground")
        {
            if (scenario.type != ScenarioType::Rigid)
            {
                return make_error(line_number, "ground belongs to type rigid");
            }
            saw_domain_key = true;
            math::Real nx{};
            math::Real ny{};
            math::Real offset{};
            if (!(tokens >> nx >> ny >> offset))
            {
                return make_error(line_number, "ground requires: nx ny offset");
            }
            // Normalized here, at the input boundary, so collide::Halfplane can
            // keep a strict unit-normal invariant rather than normalizing on
            // every query. A direction that cannot be normalized is refused
            // rather than quietly turned into one.
            const math::Vec2 direction{nx, ny};
            const math::Real magnitude = math::length(direction);
            if (!(magnitude > math::Real{0}) || !math::is_finite(magnitude) ||
                !math::is_finite(offset))
            {
                return make_error(line_number,
                                  "ground normal must be non-zero and finite");
            }
            scenario.rigid_settings.ground.push_back(
                collide::Halfplane{direction / magnitude, offset});
        }
        else if (key == "bounds")
        {
            if (scenario.type != ScenarioType::Particles)
            {
                return make_error(line_number, "bounds belongs to type particles");
            }
            saw_domain_key = true;
            math::Real min_x{};
            math::Real min_y{};
            math::Real max_x{};
            math::Real max_y{};
            if (!(tokens >> min_x >> min_y >> max_x >> max_y))
            {
                return make_error(line_number, "bounds requires: minx miny maxx maxy");
            }
            scenario.particle_settings.bounds =
                collide::Aabb{math::Vec2{min_x, min_y}, math::Vec2{max_x, max_y}};
        }
        else if (key == "particle")
        {
            if (scenario.type != ScenarioType::Particles)
            {
                return make_error(line_number, "particle belongs to type particles");
            }
            saw_domain_key = true;
            math::Real mass{};
            math::Real radius{};
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (!(tokens >> mass >> radius >> px >> py >> vx >> vy))
            {
                return make_error(line_number,
                                  "particle requires: mass radius px py vx vy");
            }
            scenario.particle_list.push_back(particles::Particle2D{
                math::Vec2{px, py}, math::Vec2{vx, vy}, mass, radius});
        }
        else if (key == "rigid_body")
        {
            if (scenario.type != ScenarioType::Rigid)
            {
                return make_error(line_number, "rigid_body belongs to type rigid");
            }
            saw_domain_key = true;
            rigid::RigidBody2D body;
            math::Real comx{};
            math::Real comy{};
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (!(tokens >> body.mass >> body.inertia >> body.radius >> comx >> comy >>
                  px >> py >> body.angle >> vx >> vy >> body.angular_velocity))
            {
                return make_error(line_number,
                                  "rigid_body requires: mass inertia radius comx comy "
                                  "px py angle vx vy omega");
            }
            body.local_center_of_mass = math::Vec2{comx, comy};
            body.position = math::Vec2{px, py};
            body.velocity = math::Vec2{vx, vy};
            scenario.rigid_bodies.push_back(body);
        }
        else if (key == "rigid_static")
        {
            if (scenario.type != ScenarioType::Rigid)
            {
                return make_error(line_number, "rigid_static belongs to type rigid");
            }
            saw_domain_key = true;
            // A separate key rather than writing "inf inf" on a rigid_body line,
            // because MSVC's num_get does not accept the token "inf", and
            // "static" reads better in a scenario than a pair of infinities.
            rigid::RigidBody2D body;
            math::Real px{};
            math::Real py{};
            if (!(tokens >> body.radius >> px >> py >> body.angle))
            {
                return make_error(line_number,
                                  "rigid_static requires: radius px py angle");
            }
            body.mass = std::numeric_limits<math::Real>::infinity();
            body.inertia = std::numeric_limits<math::Real>::infinity();
            body.position = math::Vec2{px, py};
            scenario.rigid_bodies.push_back(body);
        }
        else if (key == "spring_body")
        {
            if (scenario.type != ScenarioType::Springs)
            {
                return make_error(line_number, "spring_body belongs to type springs");
            }
            saw_domain_key = true;
            springs::SpringBody2D body;
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (!(tokens >> body.mass >> px >> py >> vx >> vy))
            {
                return make_error(line_number,
                                  "spring_body requires: mass px py vx vy");
            }
            body.position = math::Vec2{px, py};
            body.velocity = math::Vec2{vx, vy};
            scenario.spring_bodies.push_back(body);
        }
        else if (key == "spring")
        {
            if (scenario.type != ScenarioType::Springs)
            {
                return make_error(line_number, "spring belongs to type springs");
            }
            saw_domain_key = true;
            springs::Spring spring;
            if (!(tokens >> spring.a >> spring.b >> spring.rest_length >>
                  spring.stiffness >> spring.damping))
            {
                return make_error(
                    line_number, "spring requires: a b rest_length stiffness damping");
            }
            scenario.spring_network.add(spring);
        }
        else
        {
            return make_error(line_number, "unknown key '" + key + "'");
        }

        if (!at_end_of_directive(tokens))
        {
            std::string extra;
            tokens >> extra;
            return make_error(line_number, "unexpected extra token '" + extra + "'");
        }
    }

    // --- #7: getline stops on a read error as well as on end of input, so a
    //     truncated read would otherwise be reported as a successful parse.
    if (input.bad())
    {
        return make_error(line_number, "read error");
    }

    return ScenarioParseResult{true, "", std::move(scenario)};
}

ScenarioParseResult parse_scenario_file(const std::string& path)
{
    // Binary mode on purpose: in text mode the MSVC CRT treats 0x1A as end of
    // file, silently truncating a scenario mid-parse and reporting success.
    // CRLF still parses, because a carriage return is whitespace to the
    // stream extractor.
    std::ifstream file(path, std::ios::binary);
    if (!file)
    {
        return ScenarioParseResult{false, "cannot open scenario file: " + path, {}};
    }
    return parse_scenario(file);
}
} // namespace malloy::scenario
