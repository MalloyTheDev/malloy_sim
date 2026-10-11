#include <malloy/scenario/scenario.hpp>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <limits>
#include <istream>
#include <optional>
#include <sstream>
#include <string>

#include <malloy/collide/shapes.hpp>
#include <malloy/collide/shapes3d.hpp>
#include <malloy/math/vec2.hpp>
#include <malloy/particles/thermo.hpp>
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

// Upper bounds on the WORK and SIZE a scenario may request. A scenario file is
// the one untrusted input this project takes, and while the parser already
// rejects non-finite and out-of-range values, it did not bound sheer magnitude:
// a single `steps 2000000000` line makes the app spin for billions of steps, and
// an enormous body list allocates memory in proportion. These caps close that
// (a local denial of service, per the 2026-09-13 security review) without
// affecting any real scenario: the shipped templates use at most a few thousand
// steps and a handful of bodies, orders of magnitude below these limits.
constexpr int max_steps = 100000000;         // 1e8 steps
constexpr std::size_t max_entities = 100000; // total bodies plus springs

// The caps above bound steps and entities SEPARATELY, but a step of the
// pairwise domains costs O(N^2), so together they still allowed about 1e8 x 5e9
// pair tests: a 41-byte `lattice` file ran for hours (issue #25). This bounds
// their product, the work a run will do, counted in pair tests (plus a body
// against each ground plane, and a body or spring per step in the linear
// springs domain). At the roughly 7 ns a pair test costs in Release, 1e11 is
// about a quarter of an hour; the largest shipped template uses 3e7.
constexpr double max_work = 1e11;

// How many diagnostic reports a run may print. Output volume is bounded only by
// this, and a report with its ASCII frame is about 2 KB, so 1e5 reports is about
// 200 MB (issue #25). Raise `output_every` to report less often.
constexpr long long max_reports = 100000;

// Ground planes (`ground`, `plane3`) are iterated by every body every step, and
// were counted by no cap at all (issue #26). A container needs a handful.
constexpr std::size_t max_planes = 1000;

// The longest line read. Lines were read whole before any check, so one
// enormous line could allocate without limit (issue #29); the longest shipped
// template line is about 100 characters.
constexpr std::size_t max_line_length = 4096;

enum class LineRead
{
    Ok,
    End,
    TooLong
};

// Read one line, without its newline, stopping at max_line_length characters
// so that no line costs more memory than that. End means nothing was left to
// read; a last line without a newline still reads as Ok.
LineRead read_line(std::istream& input, std::string& line)
{
    line.clear();
    bool any = false;
    char c{};
    while (input.get(c))
    {
        any = true;
        if (c == '\n')
        {
            return LineRead::Ok;
        }
        if (line.size() == max_line_length)
        {
            return LineRead::TooLong;
        }
        line.push_back(c);
    }
    return any ? LineRead::Ok : LineRead::End;
}

// How much of a token from the file an error message repeats. Every key and
// type name is under 16 characters, so a typo always shows whole, while a
// 4096-byte token cannot fill the screen (issue #28).
constexpr std::size_t max_quoted_bytes = 64;

// A token from the file, made safe to print inside an error message, quotes
// included. The app prints errors to a terminal, and a token is any run of
// non-whitespace bytes, so it can carry ESC and BEL sequences that retitle the
// window, clear the screen, or erase the error and print a fake success line
// (issue #28). So every byte outside printable ASCII is written as \xNN, two
// lowercase hex digits, and so are the backslash and the quote, which leaves the
// text between the quotes unambiguous: a backslash there always starts an
// escape, and the closing quote cannot be forged.
//
// Bytes from 0x80 up are escaped too. No key or type is non-ASCII, so this never
// disguises a correctly typed token, and the non-ASCII characters are exactly
// the ones a raw echo hides or forges: a C1 control (0x9B is CSI to an 8-bit
// terminal), a bidi override that reorders the line as displayed, and the
// invisible or look-alike ones (a byte-order mark, issue #32; a no-break space,
// which makes `dt 0.001` one token; a Cyrillic letter in `type`).
//
// Only the first max_quoted_bytes bytes of the token are shown, counted before
// escaping so an escape is never cut, followed by the token's full length.
std::string quoted(const std::string& token)
{
    static constexpr char hex[] = "0123456789abcdef";
    const std::size_t shown = std::min(token.size(), max_quoted_bytes);
    std::string out = "'";
    for (std::size_t i = 0; i < shown; ++i)
    {
        const unsigned char c = static_cast<unsigned char>(token[i]);
        if (c < 0x20 || c >= 0x7F || c == '\\' || c == '\'')
        {
            out += "\\x";
            out += hex[c >> 4];
            out += hex[c & 0x0F];
        }
        else
        {
            out += static_cast<char>(c);
        }
    }
    out += '\'';
    if (token.size() > shown)
    {
        out += "... (" + std::to_string(token.size()) + " bytes)";
    }
    return out;
}

// The work one step of this scenario does, in pair tests. Every pairwise domain
// tests each pair of bodies once; the rigid domains also test each body against
// each ground plane; springs is linear, a body and a spring at a time.
double step_work(const Scenario& scenario, std::size_t bodies, std::size_t springs)
{
    const double n = static_cast<double>(bodies);
    if (scenario.type == ScenarioType::Springs || scenario.type == ScenarioType::Springs3D)
    {
        return n + static_cast<double>(springs);
    }
    const double planes = static_cast<double>(scenario.rigid_settings.ground.size() +
                                              scenario.rigid3d_settings.ground.size());
    return n * (n - 1.0) / 2.0 + n * planes + n;
}

// A `lattice` or `lattice3` line (M38), held until the whole file has been read
// and only then expanded into particles, so the result does not depend on
// whether `bounds` comes before or after it. `nz` stays 1 for the 2D key.
struct LatticeRequest
{
    int nx{0};
    int ny{0};
    int nz{1};
    math::Real mass{};
    math::Real radius{};
    math::Real speed{};
    std::size_t line{0};
};
} // namespace

ScenarioParseResult parse_scenario(std::istream& input)
{
    Scenario scenario;
    std::string line;
    std::size_t line_number = 0;
    // `type` selects which domain keys are legal, so it has to be settled
    // before any of them appear.
    bool saw_domain_key = false;
    std::optional<LatticeRequest> lattice;

    for (;;)
    {
        const LineRead read = read_line(input, line);
        if (read == LineRead::End)
        {
            break;
        }
        ++line_number;

        // Issue #32: Windows editors, and PowerShell 5.1's `Set-Content -Encoding
        // UTF8`, start a UTF-8 file with a byte-order mark, which used to become
        // part of the first key and fail a valid file with an invisible one. A
        // UTF-16 file (PowerShell 5.1's `>`) failed the same way, NUL bytes and
        // all, so it is now refused by name. Only a mark at the very start of
        // the input is a byte-order mark: anywhere else its bytes are text, and
        // an error shows them escaped. This runs before the length check, so a
        // UTF-16 file with no 0x0A byte in its first 4096 still gets the
        // encoding error; the mark's three bytes count toward line 1's length.
        if (line_number == 1)
        {
            if (line.starts_with("\xEF\xBB\xBF"))
            {
                line.erase(0, 3);
            }
            else if (line.starts_with("\xFF\xFE") || line.starts_with("\xFE\xFF"))
            {
                return make_error(line_number,
                                  "scenario files must be UTF-8 or ASCII (this one "
                                  "starts with a UTF-16 byte-order mark)");
            }
        }

        if (read == LineRead::TooLong)
        {
            return make_error(line_number, "line is longer than 4096 characters");
        }

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
                                  "charges, nbody3d, rigid3d, charges3d, "
                                  "springs3d or particles3d");
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
            else if (value == "springs3d")
            {
                scenario.type = ScenarioType::Springs3D;
            }
            else if (value == "particles3d")
            {
                scenario.type = ScenarioType::Particles3D;
            }
            else
            {
                return make_error(line_number, "unknown type " + quoted(value));
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
            if (scenario.steps > max_steps)
            {
                return make_error(line_number, "steps exceeds the maximum (100000000)");
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
                scenario.type != ScenarioType::Particles3D &&
                scenario.type != ScenarioType::Rigid &&
                scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number,
                                  "restitution belongs to type particles, "
                                  "particles3d, rigid or rigid3d");
            }
            saw_domain_key = true;
            math::Real value{};
            if (!(tokens >> value))
            {
                return make_error(line_number, "restitution requires a number");
            }
            // Write only the field the declared domain actually reads.
            if (scenario.type == ScenarioType::Rigid3D)
            {
                scenario.rigid3d_settings.restitution = value;
            }
            else if (scenario.type == ScenarioType::Particles)
            {
                scenario.particle_settings.restitution = value;
            }
            else if (scenario.type == ScenarioType::Particles3D)
            {
                scenario.particle3d_settings.restitution = value;
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
            // describe the one box. The moments are taken in the box's OWN axis
            // order, the order `half_extents` is in: the principal moments from
            // mass_properties_3d come back sorted, which swapped x and z for any
            // box with hx < hz (issue #21). A degenerate box gives zero moments,
            // whose invalid body the domain's validate() then rejects.
            const math::Vec3 half_extents{hx, hy, hz};
            rigid::RigidBody3D body;
            body.mass = mass;
            body.inertia = rigid::box_inertia(half_extents, mass);
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
            if (scenario.type != ScenarioType::Rigid3D &&
                scenario.type != ScenarioType::Particles3D)
            {
                return make_error(line_number,
                                  "gravity3 belongs to type rigid3d or particles3d");
            }
            saw_domain_key = true;
            math::Real gx{};
            math::Real gy{};
            math::Real gz{};
            if (!(tokens >> gx >> gy >> gz))
            {
                return make_error(line_number, "gravity3 requires: gx gy gz");
            }
            if (scenario.type == ScenarioType::Particles3D)
            {
                scenario.particle3d_settings.gravity = math::Vec3{gx, gy, gz};
            }
            else
            {
                scenario.rigid3d_settings.gravity = math::Vec3{gx, gy, gz};
            }
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
            if (scenario.rigid3d_settings.ground.size() > max_planes)
            {
                return make_error(line_number, "scenario exceeds the maximum plane count (1000)");
            }
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
            if (scenario.type != ScenarioType::Particles &&
                scenario.type != ScenarioType::Particles3D &&
                scenario.type != ScenarioType::Rigid &&
                scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number,
                                  "friction belongs to type particles, "
                                  "particles3d, rigid or rigid3d");
            }
            saw_domain_key = true;
            math::Real value{};
            if (!(tokens >> value))
            {
                return make_error(line_number, "friction requires a value");
            }
            // Write only the field the declared domain actually reads.
            if (scenario.type == ScenarioType::Rigid3D)
            {
                scenario.rigid3d_settings.friction = value;
            }
            else if (scenario.type == ScenarioType::Particles)
            {
                scenario.particle_settings.friction = value;
            }
            else if (scenario.type == ScenarioType::Particles3D)
            {
                scenario.particle3d_settings.friction = value;
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
            if (scenario.rigid_settings.ground.size() > max_planes)
            {
                return make_error(line_number, "scenario exceeds the maximum plane count (1000)");
            }
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
        else if (key == "bounds3")
        {
            if (scenario.type != ScenarioType::Particles3D)
            {
                return make_error(line_number, "bounds3 belongs to type particles3d");
            }
            saw_domain_key = true;
            math::Real min_x{};
            math::Real min_y{};
            math::Real min_z{};
            math::Real max_x{};
            math::Real max_y{};
            math::Real max_z{};
            if (!(tokens >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z))
            {
                return make_error(line_number,
                                  "bounds3 requires: minx miny minz maxx maxy maxz");
            }
            scenario.particle3d_settings.bounds = collide::Aabb3{
                math::Vec3{min_x, min_y, min_z}, math::Vec3{max_x, max_y, max_z}};
        }
        else if (key == "particle3")
        {
            if (scenario.type != ScenarioType::Particles3D)
            {
                return make_error(line_number, "particle3 belongs to type particles3d");
            }
            saw_domain_key = true;
            math::Real mass{};
            math::Real radius{};
            math::Real px{};
            math::Real py{};
            math::Real pz{};
            math::Real vx{};
            math::Real vy{};
            math::Real vz{};
            if (!(tokens >> mass >> radius >> px >> py >> pz >> vx >> vy >> vz))
            {
                return make_error(line_number,
                                  "particle3 requires: mass radius px py pz vx vy vz");
            }
            scenario.particle_list3d.push_back(particles::Particle3D{
                math::Vec3{px, py, pz}, math::Vec3{vx, vy, vz}, mass, radius});
        }
        else if (key == "lattice" || key == "lattice3")
        {
            const bool three = key == "lattice3";
            if (three ? scenario.type != ScenarioType::Particles3D
                      : scenario.type != ScenarioType::Particles)
            {
                return make_error(line_number, three
                                                   ? "lattice3 belongs to type particles3d"
                                                   : "lattice belongs to type particles");
            }
            saw_domain_key = true;
            if (lattice)
            {
                return make_error(line_number, "only one lattice line is allowed");
            }
            LatticeRequest request;
            request.line = line_number;
            const bool counts_read = three ? static_cast<bool>(tokens >> request.nx >> request.ny >> request.nz)
                                           : static_cast<bool>(tokens >> request.nx >> request.ny);
            if (!counts_read || !(tokens >> request.mass >> request.radius >> request.speed))
            {
                return make_error(line_number,
                                  three ? "lattice3 requires: nx ny nz mass radius speed"
                                        : "lattice requires: nx ny mass radius speed");
            }
            if (request.nx < 1 || request.ny < 1 || request.nz < 1)
            {
                return make_error(line_number, "lattice counts must be at least 1");
            }
            // The cap is checked HERE, before anything is allocated: this is the
            // one key that turns a few bytes of input into many entities. Each
            // count is bounded first so the product cannot overflow.
            const std::size_t nx = static_cast<std::size_t>(request.nx);
            const std::size_t ny = static_cast<std::size_t>(request.ny);
            const std::size_t nz = static_cast<std::size_t>(request.nz);
            if (nx > max_entities || ny > max_entities || nz > max_entities ||
                nx * ny * nz > max_entities)
            {
                return make_error(line_number,
                                  "lattice exceeds the maximum entity count (100000)");
            }
            lattice = request;
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
        else if (key == "rigid_box")
        {
            if (scenario.type != ScenarioType::Rigid)
            {
                return make_error(line_number, "rigid_box belongs to type rigid");
            }
            saw_domain_key = true;
            rigid::RigidBody2D body;
            math::Real hx{};
            math::Real hy{};
            math::Real px{};
            math::Real py{};
            math::Real vx{};
            math::Real vy{};
            if (!(tokens >> body.mass >> hx >> hy >> px >> py >> body.angle >> vx >>
                  vy >> body.angular_velocity))
            {
                return make_error(line_number,
                                  "rigid_box requires: mass hx hy px py angle vx vy "
                                  "omega");
            }
            body.half_extents = math::Vec2{hx, hy};
            // Inertia of a solid box about its centre, m (w^2 + h^2) / 12 with
            // w = 2 hx and h = 2 hy, computed here as the loader does for the 3D
            // box so the file need not carry it.
            body.inertia = body.mass *
                           (math::Real{4} * hx * hx + math::Real{4} * hy * hy) /
                           math::Real{12};
            body.position = math::Vec2{px, py};
            body.velocity = math::Vec2{vx, vy};
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
        else if (key == "spring_body3")
        {
            if (scenario.type != ScenarioType::Springs3D)
            {
                return make_error(line_number, "spring_body3 belongs to type springs3d");
            }
            saw_domain_key = true;
            springs::SpringBody3D body;
            math::Real px{};
            math::Real py{};
            math::Real pz{};
            math::Real vx{};
            math::Real vy{};
            math::Real vz{};
            if (!(tokens >> body.mass >> px >> py >> pz >> vx >> vy >> vz))
            {
                return make_error(line_number,
                                  "spring_body3 requires: mass px py pz vx vy vz");
            }
            body.position = math::Vec3{px, py, pz};
            body.velocity = math::Vec3{vx, vy, vz};
            scenario.spring_bodies3d.push_back(body);
        }
        else if (key == "spring")
        {
            if (scenario.type != ScenarioType::Springs &&
                scenario.type != ScenarioType::Springs3D)
            {
                return make_error(line_number,
                                  "spring belongs to type springs or springs3d");
            }
            saw_domain_key = true;
            springs::Spring spring;
            if (!(tokens >> spring.a >> spring.b >> spring.rest_length >>
                  spring.stiffness >> spring.damping))
            {
                return make_error(
                    line_number, "spring requires: a b rest_length stiffness damping");
            }
            // The Spring type is dimension-agnostic, so it goes to whichever
            // network the scenario's dimension selected.
            if (scenario.type == ScenarioType::Springs3D)
            {
                scenario.spring_network3d.add(spring);
            }
            else
            {
                scenario.spring_network.add(spring);
            }
        }
        else
        {
            return make_error(line_number, "unknown key " + quoted(key));
        }

        if (!at_end_of_directive(tokens))
        {
            std::string extra;
            tokens >> extra;
            return make_error(line_number, "unexpected extra token " + quoted(extra));
        }

        // The entity cap is checked as each line is read, so an oversized file
        // fails on the line that crossed it, before allocating any further
        // (issue #29). The lattice is the exception: its own line checks its
        // count, and its particles are added, and counted, after the file.
        if (entity_count(scenario) > max_entities)
        {
            return make_error(line_number, "scenario exceeds the maximum entity count (100000)");
        }
    }

    // --- #7: read_line stops on a read error as well as on end of input, so a
    //     truncated read would otherwise be reported as a successful parse.
    if (input.bad())
    {
        return make_error(line_number, "read error");
    }

    // M38: expand the lattice now that the bounds are final. Its particles are
    // appended after any listed one by one. The spacing check needs the bounds,
    // which is why it lives here rather than on the lattice line; an invalid
    // box, or an invalid radius, is left for validate() to reject, so each rule
    // lives in one place.
    if (lattice)
    {
        const math::Real diameter = math::Real{2} * lattice->radius;
        if (scenario.type == ScenarioType::Particles3D)
        {
            const collide::Aabb3& bounds = scenario.particle3d_settings.bounds;
            if (bounds.is_valid() &&
                ((bounds.max.x - bounds.min.x) / lattice->nx < diameter ||
                 (bounds.max.y - bounds.min.y) / lattice->ny < diameter ||
                 (bounds.max.z - bounds.min.z) / lattice->nz < diameter))
            {
                return make_error(lattice->line,
                                  "lattice spacing is smaller than a particle diameter");
            }
            const std::vector<particles::Particle3D> grid = particles::lattice_gas3d(
                bounds, lattice->nx, lattice->ny, lattice->nz, lattice->mass,
                lattice->radius, lattice->speed);
            scenario.particle_list3d.insert(scenario.particle_list3d.end(), grid.begin(),
                                            grid.end());
        }
        else
        {
            const collide::Aabb& bounds = scenario.particle_settings.bounds;
            if (bounds.is_valid() &&
                ((bounds.max.x - bounds.min.x) / lattice->nx < diameter ||
                 (bounds.max.y - bounds.min.y) / lattice->ny < diameter))
            {
                return make_error(lattice->line,
                                  "lattice spacing is smaller than a particle diameter");
            }
            const std::vector<particles::Particle2D> grid =
                particles::lattice_gas(bounds, lattice->nx, lattice->ny, lattice->mass,
                                       lattice->radius, lattice->speed);
            scenario.particle_list.insert(scenario.particle_list.end(), grid.begin(),
                                          grid.end());
        }
    }

    // Defence-in-depth cap on total scenario size (see max_entities above): a
    // hostile or fat-fingered file cannot make the app allocate or step over an
    // unbounded number of bodies and springs. Every line was checked as it was
    // read; this catches a lattice that, with the listed entities, goes over,
    // and names the lattice's line.
    const std::size_t total_entities = entity_count(scenario);
    if (total_entities > max_entities)
    {
        return make_error(lattice ? lattice->line : line_number,
                          "scenario exceeds the maximum entity count (100000)");
    }

    // The work budget and the report cap (issue #25). Both depend on values
    // that may sit on different lines, so neither names one. A report costs
    // about a step too, since its diagnostics include an O(N^2) potential in
    // the nbody and charges domains, so the work counts the reports as well as
    // the steps: otherwise a `steps 0` file of 1e5 bodies counted as no work,
    // while its one report ran 5e9 pair tests.
    const std::size_t springs =
        scenario.spring_network.springs().size() + scenario.spring_network3d.springs().size();
    const long long reports =
        scenario.output_every > 0
            ? static_cast<long long>(scenario.steps) / scenario.output_every + 1
            : 1;
    const double work = (static_cast<double>(scenario.steps) + static_cast<double>(reports)) *
                        step_work(scenario, total_entities - springs, springs);
    if (work > max_work)
    {
        return ScenarioParseResult{
            false,
            "scenario exceeds the work budget: steps plus reports, times pair tests "
            "per step, must be at most 1e11 (fewer steps or fewer bodies)",
            {}};
    }
    if (reports > max_reports)
    {
        return ScenarioParseResult{
            false,
            "scenario would print more than 100000 reports (raise output_every)",
            {}};
    }

    return ScenarioParseResult{true, "", std::move(scenario)};
}

std::size_t entity_count(const Scenario& scenario)
{
    // Only one domain's lists are ever populated, but summing them all keeps
    // this one centralised count.
    return scenario.bodies.size() + scenario.bodies3d.size() +
           scenario.particle_list.size() + scenario.particle_list3d.size() +
           scenario.rigid_bodies.size() + scenario.rigid_bodies3d.size() +
           scenario.charge_list.size() + scenario.charge_list3d.size() +
           scenario.spring_bodies.size() + scenario.spring_bodies3d.size() +
           scenario.spring_network.springs().size() +
           scenario.spring_network3d.springs().size();
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
