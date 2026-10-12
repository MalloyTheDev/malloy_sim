#include <malloy/scenario/scenario.hpp>

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <istream>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>

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
// their product, the work a run will do, counted in WORK UNITS of about 10 ns of
// Release time each: one cheap pair test (an nbody or particle pair), plus one
// per body per step. Since issue #46 the costlier tests are weighted by what
// they cost (see the weights below), so 1e11 units is at most about 20 minutes
// in Release whatever the file holds (caged piles at the limit ran 5 to 17
// minutes under the #46 weights; the M40 rigid weights are worst-case bounds,
// of which the slowest piles measured reached about half), and up to about
// three hours in Debug, where a unit costs up to about 100 ns (the same caged
// piles: 38 minutes to 2 hours). The largest shipped template uses 3.8e8
// (dropped_bodies.scn, a rigid one), the largest other 3e7 (a lattice gas).
constexpr double max_work = 1e11;

// What each test costs, in work units, for the rigid domains, whose tests
// differ in kind. Issue #46 first weighted them by what detection cost; since
// M40 a contact also feeds the per-point solver, which dominates. A step
// detects, then gives each contact point a fixed cost (its effective masses,
// the warm start, the prepass, the position projection and the warm-start
// cache) and, per sweep, a residual check and a sweep. Each phase sweeps until
// its tolerance or until kMaxSweeps, 100 (src/rigid/rigid.cpp, rigid3d.cpp), and
// a step has two phases when restitution is above 0 (impact, then rest) and one
// at 0. So the worst case a file can hold, every pair overlapping and every
// phase stopped by the cap, costs per test
//
//   detection + (most points it reports) x (fixed + phases x 100 x sweep)
//
// with the points the contact point cap counts below: 8 for a 3D box against a
// box or a plane, 2 for a 2D box pair, 4 for a 2D box on a ground, and 1 for any
// other pair or plane. Measured in Release on fresh piles in which every pair
// overlaps deeply (bodies within 0.1 of one another, approaching at up to 1,
// friction 0.5, restitution 0.5 and 0), with every phase forced to the cap, and
// split by running the same piles with the cap at 0 and with detection alone.
// Each figure is the slowest measured, rounded up into units of 10 ns:
//
//   detection, per overlapping test: 3D box-box 4.53 us (460), box-plane 1.99
//     (200), box-sphere at most 0.81 (80, the average of a mixed pile, box
//     pairs included), sphere pair 0.43 (45), sphere-plane 0.42 (45); 2D
//     box-box 0.80 (85), box-ground 0.72 (75), disc pair 0.31 (35),
//     disc-ground 0.19 (20)
//   fixed, per point: 3D, a contact with a box 0.72 us (75), spheres only 0.46
//     (50); 2D 0.73 (75)
//   sweep, per point and sweep, in the rest phase, the dearer of the two. It
//     grows with the contact list as it outgrows the caches; measured from a
//     few hundred points to 2.6e5 (the mixed-phase cost grew 0.12 to 0.17 us
//     for 3D boxes, 0.046 to 0.071 for spheres, 0.020 to 0.035 for 2D), scaled
//     by that growth: 3D, a contact with a box 0.20 to 0.28 us (28), spheres
//     only 0.084 to 0.127 (13); 2D 0.036 to 0.062 (6.5)
//
// With two phases that gives 45860 for a 3D box pair, 45600 for a box-plane
// test, 5755 for a box-sphere pair, 2695 for a sphere pair or sphere-plane test,
// and in 2D 2835 for a box pair, 5575 for a box-ground test, 1410 for a disc
// pair and 1395 for a disc-ground test; with one phase 23460, 23200, 2955,
// 1395, 1535, 2975, 760 and 745. A test that cannot make contact costs 1 (1.4 ns
// in 3D, 2.2 ns in 2D): a pair or plane test where a body has no collider, two
// immovable 2D bodies, an immovable 2D body against a ground, or a 2D box with
// a disc (issue #48), all of which the worlds skip before detecting. A rigid
// body costs 15 units a step in 3D and 10 in 2D: its integration, validation,
// and entries in the phases' per-body arrays (0.08 us measured, 0.15 us in the
// solver's own figures). A body is classified as its world classifies it: a box
// when every half-extent is positive, otherwise round when its radius is above
// 0, and in 2D immovable when its mass and inertia are both infinite.
//
// The bound is far above a typical step. A box resting on a plane took 2.5 us
// a step in 3D and 0.9 us in 2D, a three-box stack 11 us and 3 us (under one
// sweep a step), and a five-box 3D stack twisted 0.3 rad a level 91 us (20
// sweeps), where its bound is 6.9 ms. Debug ran these piles 3 to 8 times slower
// than Release (a forced-cap 3D box pile 240 us a test against 70 us, a 2D disc
// pile 18 us against 2.5 us), so a unit stays under the 100 ns allowed above
// wherever the work is large; only a body's own cost runs about 30 times
// slower, which matters in small scenes, far inside the budget.
constexpr double sweeps_per_phase = 100.0; // kMaxSweeps in src/rigid
constexpr double rigid3_body_work = 15.0;
constexpr double rigid2_body_work = 10.0;
constexpr double no_contact_work = 1.0;
constexpr double box3_box3_detect = 460.0;
constexpr double box3_plane_detect = 200.0;
constexpr double box3_sphere3_detect = 80.0;
constexpr double sphere3_sphere3_detect = 45.0;
constexpr double sphere3_plane_detect = 45.0;
constexpr double box3_point_fixed = 75.0;
constexpr double sphere3_point_fixed = 50.0;
constexpr double box3_point_sweep = 28.0;
constexpr double sphere3_point_sweep = 13.0;
constexpr double box2_box2_detect = 85.0;
constexpr double box2_ground_detect = 75.0;
constexpr double disc2_pair_detect = 35.0;
constexpr double disc2_ground_detect = 20.0;
constexpr double point2_fixed = 75.0;
constexpr double point2_sweep = 6.5;

// How many diagnostic reports a run may print. Output volume is bounded only by
// this, and a report with its ASCII frame is about 2 KB, so 1e5 reports is about
// 200 MB (issue #25). Raise `output_every` to report less often.
constexpr long long max_reports = 100000;

// Ground planes (`ground`, `plane3`) are iterated by every body every step, and
// were counted by no cap at all (issue #26). A container needs a handful.
constexpr std::size_t max_planes = 1000;

// The most contact points a rigid world could have to hold at once (M40). The
// per-point contact solver keeps every point of a step in memory together, and
// keeps the last step's as warm-start state, so memory grows with the pairs that
// touch, which neither the entity cap nor the work budget bounds: a file of a
// few thousand coincident bodies, inside both, could ask for tens of GB. This
// caps the POTENTIAL count, every pair and every body-plane test that can make
// contact, at the most points that contact can report (see
// potential_contact_points below), so it holds however the bodies end up
// arranged.
//
// The value is a memory ceiling of 2 GiB over what a point costs. At the peak
// of a step a point is held in the step's list (a 3D ContactPoint is 328 bytes,
// a 2D one 200), in at most one manifold record (32 bytes), and in the world's
// warm-start cache and that cache's rollback snapshot (a WarmContact is 80
// bytes in 3D, 56 in 2D): 520 bytes in 3D, 344 in 2D (MSVC x64). A std::vector
// can hold up to twice its size while it grows (libstdc++ doubles; MSVC adds
// half), so a point costs at most 1040 bytes, 1 KiB. 2 GiB / 1 KiB is 2097152,
// rounded down to 2000000. A body's own per-step arrays (a 3D solver body is
// 160 bytes) are bounded by the entity cap instead. The rigid templates hold at
// most 54 (dropped_bodies.scn). A pile of 1994 spheres in a box of six planes,
// or of 1997 discs in three grounds, is the largest of its kind that loads.
constexpr std::uint64_t max_contact_points = 2000000;

// The most points one contact can report, by kind, as the collide layer builds
// them, margin-aware queries included. In 3D a box face clipped against another
// box is a polygon of at most 8 vertices (a quad cut by four side planes, each
// cut adding at most one), and a box near a plane has at most its 8 corners; an
// edge-edge contact, and every other 3D pair or plane, is one point. In 2D a
// clipped box edge has 2 ends and a box has 4 corners for a ground line; a disc
// is one point. A pair the world skips before detecting (two immovable 2D
// bodies, an immovable 2D body against a ground, a 2D box with a disc, issue
// #48) makes none.
constexpr std::uint64_t box3_box3_points = 8;
constexpr std::uint64_t box3_plane_points = 8;
constexpr std::uint64_t box2_box2_points = 2;
constexpr std::uint64_t box2_ground_points = 4;
constexpr std::uint64_t single_points = 1;

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

// How many rigid bodies of each collider kind a scenario holds, as its world
// classifies them: a box when every half-extent is positive, otherwise a round
// body (a sphere, or a disc in 2D) when its radius is above 0. A body that is
// neither collides with nothing and is not counted. In 2D a body is also
// immovable when its mass and inertia are both infinite (rigid_static): two
// immovable bodies, or one and a ground, are never detected. `seen` is how many
// bodies have been classified so far, so the parser counts as it reads, in time
// proportional to the file rather than to its square.
struct ColliderCount
{
    std::uint64_t boxes{0};          // in 2D, the movable ones
    std::uint64_t rounds{0};         // in 2D, the movable ones
    std::uint64_t immovable_boxes{0};  // 2D only
    std::uint64_t immovable_rounds{0}; // 2D only
    std::size_t seen{0};
};

// Classifies the bodies added since the last call. Only the list of the
// scenario's own type is read; `type` is settled before any body line.
void count_colliders(const Scenario& scenario, ColliderCount& count)
{
    if (scenario.type == ScenarioType::Rigid3D)
    {
        for (; count.seen < scenario.rigid_bodies3d.size(); ++count.seen)
        {
            const rigid::RigidBody3D& body = scenario.rigid_bodies3d[count.seen];
            if (body.half_extents.x > math::Real{0} && body.half_extents.y > math::Real{0} &&
                body.half_extents.z > math::Real{0})
            {
                ++count.boxes;
            }
            else if (body.radius > math::Real{0})
            {
                ++count.rounds;
            }
        }
    }
    else if (scenario.type == ScenarioType::Rigid)
    {
        for (; count.seen < scenario.rigid_bodies.size(); ++count.seen)
        {
            const rigid::RigidBody2D& body = scenario.rigid_bodies[count.seen];
            const bool immovable = body.is_static();
            if (body.half_extents.x > math::Real{0} && body.half_extents.y > math::Real{0})
            {
                ++(immovable ? count.immovable_boxes : count.boxes);
            }
            else if (body.radius > math::Real{0})
            {
                ++(immovable ? count.immovable_rounds : count.rounds);
            }
        }
    }
}

// Pairs among `n`, and pairs with at least one of `n` beside `others`.
std::uint64_t pairs_of(std::uint64_t n) { return n > 0 ? n * (n - 1) / 2 : 0; }
std::uint64_t pairs_with(std::uint64_t n, std::uint64_t others) { return pairs_of(n) + n * others; }

// The potential contact points of a scenario with these colliders (see
// max_contact_points): each pair that can make contact, and each collider that
// can against each ground plane, at the most points that kind of contact
// reports. Exact in unsigned 64-bit arithmetic: inside the parser the entity
// cap keeps a count at or below 1e5, so the largest term is 8 x 1e5 x 1e5 / 2 =
// 4e10, and overflow would need about 2e9 bodies.
std::uint64_t contact_points(const Scenario& scenario, const ColliderCount& count)
{
    const std::uint64_t boxes = count.boxes;
    const std::uint64_t rounds = count.rounds;
    if (scenario.type == ScenarioType::Rigid3D)
    {
        const std::uint64_t planes = scenario.rigid3d_settings.ground.size();
        return box3_box3_points * pairs_of(boxes) +
               single_points * (boxes * rounds + pairs_of(rounds)) +
               planes * (box3_plane_points * boxes + single_points * rounds);
    }
    if (scenario.type == ScenarioType::Rigid)
    {
        // An immovable body counts only beside a movable one, and never
        // against a ground; a box and a disc are skipped (issue #48).
        const std::uint64_t planes = scenario.rigid_settings.ground.size();
        return box2_box2_points * pairs_with(boxes, count.immovable_boxes) +
               single_points * pairs_with(rounds, count.immovable_rounds) +
               planes * (box2_ground_points * boxes + single_points * rounds);
    }
    return 0;
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
    // M39: which walls a `wall_temperature` line has already named, one bit per
    // wall, so naming one twice (most likely a typo for its opposite) is an
    // error rather than a silent overwrite.
    unsigned named_walls = 0;
    // Issue #48: the first `rigid_box` line, and the first disc that collides
    // (a `rigid_body` or `rigid_static` with a radius above 0), with the key
    // that declared it. A 2D box and disc pass through each other, because box
    // against disc is not built yet, so a file that holds both is refused on
    // the line where the second kind first appears.
    std::size_t first_box_line = 0;
    std::size_t first_disc_line = 0;
    const char* first_disc_key = "";
    // M40: the rigid colliders read so far, for the contact point cap.
    ColliderCount colliders;

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
        else if (key == "force")
        {
            // M28's constant applied force, the translational half of the
            // wrench whose rotational half is `torque`, and read exactly as
            // torque is (issue #47): one world-frame vector for every body, a
            // later line replacing an earlier one, and a value whose square
            // overflows left for Rigid3DSettings::is_valid to refuse.
            if (scenario.type != ScenarioType::Rigid3D)
            {
                return make_error(line_number, "force belongs to type rigid3d");
            }
            saw_domain_key = true;
            math::Real fx{};
            math::Real fy{};
            math::Real fz{};
            if (!(tokens >> fx >> fy >> fz))
            {
                return make_error(line_number, "force requires: fx fy fz");
            }
            scenario.rigid3d_settings.force = math::Vec3{fx, fy, fz};
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
        else if (key == "seed")
        {
            if (scenario.type != ScenarioType::Particles &&
                scenario.type != ScenarioType::Particles3D)
            {
                return make_error(line_number, "seed belongs to type particles or particles3d");
            }
            saw_domain_key = true;
            std::string text;
            if (!(tokens >> text))
            {
                return make_error(line_number, "seed requires a whole number");
            }
            // from_chars, not >>: stream extraction into an unsigned type reads
            // "-1" as 18446744073709551615. from_chars takes plain decimal digits
            // only, refusing a sign, a fraction, an exponent, a hex prefix and
            // anything out of range, and the whole token must be consumed.
            std::uint64_t value{};
            const char* const first = text.data();
            const char* const last = first + text.size();
            const std::from_chars_result digits = std::from_chars(first, last, value);
            if (digits.ec != std::errc{} || digits.ptr != last)
            {
                return make_error(line_number,
                                  "seed must be a whole number from 0 to "
                                  "18446744073709551615, not " + quoted(text));
            }
            if (scenario.type == ScenarioType::Particles3D)
            {
                scenario.particle3d_settings.seed = value;
            }
            else
            {
                scenario.particle_settings.seed = value;
            }
        }
        else if (key == "wall_temperature")
        {
            if (scenario.type != ScenarioType::Particles &&
                scenario.type != ScenarioType::Particles3D)
            {
                return make_error(line_number,
                                  "wall_temperature belongs to type particles or particles3d");
            }
            saw_domain_key = true;
            std::string wall;
            if (!(tokens >> wall))
            {
                return make_error(line_number, "wall_temperature requires: wall temperature");
            }
            // The wall names are the impulse tally's. The name is resolved
            // before the temperature is read, so a bad name is reported as such
            // whatever follows it.
            const bool three = scenario.type == ScenarioType::Particles3D;
            particles::WallTemperature2D& w2 = scenario.particle_settings.wall_temperature;
            particles::WallTemperature3D& w3 = scenario.particle3d_settings.wall_temperature;
            const std::pair<const char*, math::Real*> walls[] = {
                {"min_x", three ? &w3.min_x : &w2.min_x},
                {"max_x", three ? &w3.max_x : &w2.max_x},
                {"min_y", three ? &w3.min_y : &w2.min_y},
                {"max_y", three ? &w3.max_y : &w2.max_y},
                {"min_z", three ? &w3.min_z : nullptr},
                {"max_z", three ? &w3.max_z : nullptr},
            };
            math::Real* slot = nullptr;
            unsigned bit = 0;
            for (std::size_t i = 0; i < std::size(walls); ++i)
            {
                if (walls[i].second != nullptr && wall == walls[i].first)
                {
                    slot = walls[i].second;
                    bit = 1u << i;
                }
            }
            if (slot == nullptr)
            {
                return make_error(line_number,
                                  "unknown wall " + quoted(wall) +
                                      (three ? " (walls are min_x, max_x, min_y, max_y, "
                                               "min_z and max_z)"
                                             : " (walls are min_x, max_x, min_y and max_y)"));
            }
            if ((named_walls & bit) != 0u)
            {
                // `wall` is one of the known names here, so it is safe to echo.
                return make_error(line_number, "wall_temperature names " + wall + " twice");
            }
            // A negative temperature parses and is left for validate() to
            // refuse, so that rule lives in one place. A value that is not a
            // finite number (nan, inf, or one too large for a double) does not
            // parse at all.
            math::Real temperature{};
            if (!(tokens >> temperature))
            {
                return make_error(line_number, "wall_temperature requires: wall temperature");
            }
            named_walls |= bit;
            *slot = temperature;
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
            if (body.radius > math::Real{0} && first_disc_line == 0)
            {
                first_disc_line = line_number;
                first_disc_key = "rigid_body";
            }
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
            if (body.radius > math::Real{0} && first_disc_line == 0)
            {
                first_disc_line = line_number;
                first_disc_key = "rigid_static";
            }
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
            if (first_box_line == 0)
            {
                first_box_line = line_number;
            }
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

        // Issue #48: both kinds are now present, and this line brought in the
        // second, so it is the line where the mix first occurs. An interim
        // refusal until box against disc has its own milestone. A body of
        // radius 0 collides with nothing, so it does not count. Only `type
        // rigid` reads these keys, so only it can get here.
        if (first_box_line != 0 && first_disc_line != 0)
        {
            const bool disc_here = first_disc_line == line_number;
            return make_error(
                line_number,
                "box against disc is not supported yet: this " +
                    std::string(disc_here ? first_disc_key : "rigid_box") +
                    " would pass through the " +
                    std::string(disc_here ? "rigid_box" : first_disc_key) + " on line " +
                    std::to_string(disc_here ? first_box_line : first_disc_line) +
                    " (a disc with a radius above 0 cannot share a scenario with a box)");
        }

        // The entity cap is checked as each line is read, so an oversized file
        // fails on the line that crossed it, before allocating any further
        // (issue #29). The lattice is the exception: its own line checks its
        // count, and its particles are added, and counted, after the file.
        if (entity_count(scenario) > max_entities)
        {
            return make_error(line_number, "scenario exceeds the maximum entity count (100000)");
        }

        // M40: so is the contact point cap, which grows with each colliding
        // body and each ground plane, so the error names the line that crossed
        // it. Only the bodies this line added are classified.
        count_colliders(scenario, colliders);
        if (contact_points(scenario, colliders) > max_contact_points)
        {
            return make_error(line_number,
                              "scenario exceeds the maximum contact point count (2000000)");
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
    const long long reports =
        scenario.output_every > 0
            ? static_cast<long long>(scenario.steps) / scenario.output_every + 1
            : 1;
    const double work = (static_cast<double>(scenario.steps) + static_cast<double>(reports)) *
                        step_work(scenario);
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

// The work one step does, in work units (see max_work): a unit per body, plus
// its tests. Springs is linear, a body and a spring at a time. The 2D nbody
// world evaluates every ORDERED pair, N(N-1) force terms (src/nbody/nbody.cpp),
// and was counted as N(N-1)/2, half its cost (issue #46); every other pairwise
// domain tests each pair once. The rigid domains weight each pair, and each
// body's test against each ground plane, by its kind. Public since issue #84,
// so the fuzzers bound a run with this formula and not with a copy of it.
double step_work(const Scenario& scenario)
{
    const std::size_t springs =
        scenario.spring_network.springs().size() + scenario.spring_network3d.springs().size();
    const double n = static_cast<double>(entity_count(scenario) - springs);
    if (scenario.type == ScenarioType::Springs || scenario.type == ScenarioType::Springs3D)
    {
        return n + static_cast<double>(springs);
    }
    if (scenario.type == ScenarioType::NBody)
    {
        return n * (n - 1.0) + n;
    }
    if (scenario.type == ScenarioType::Rigid3D || scenario.type == ScenarioType::Rigid)
    {
        // The weights above: detection plus, for every point the contact can
        // report, its fixed cost and the sweeps of each phase at the cap.
        ColliderCount count;
        count_colliders(scenario, count);
        const bool three = scenario.type == ScenarioType::Rigid3D;
        const math::Real restitution =
            three ? scenario.rigid3d_settings.restitution : scenario.rigid_settings.restitution;
        const double phases = restitution > math::Real{0} ? 2.0 : 1.0;
        const double sweeps = phases * sweeps_per_phase;
        const double boxes = static_cast<double>(count.boxes);
        const double rounds = static_cast<double>(count.rounds);
        const double all_pairs = n * (n - 1.0) / 2.0;
        if (three)
        {
            const double planes = static_cast<double>(scenario.rigid3d_settings.ground.size());
            const double box_point = box3_point_fixed + sweeps * box3_point_sweep;
            const double round_point = sphere3_point_fixed + sweeps * sphere3_point_sweep;
            const double box_pairs = static_cast<double>(pairs_of(count.boxes));
            const double mixed_pairs = boxes * rounds;
            const double round_pairs = static_cast<double>(pairs_of(count.rounds));
            const double box_points = static_cast<double>(box3_box3_points);
            const double plane_points = static_cast<double>(box3_plane_points);
            return rigid3_body_work * n +
                   no_contact_work * (all_pairs - box_pairs - mixed_pairs - round_pairs) +
                   box_pairs * (box3_box3_detect + box_points * box_point) +
                   mixed_pairs * (box3_sphere3_detect + box_point) +
                   round_pairs * (sphere3_sphere3_detect + round_point) +
                   planes * (boxes * (box3_plane_detect + plane_points * box_point) +
                             rounds * (sphere3_plane_detect + round_point) +
                             (n - boxes - rounds) * no_contact_work);
        }
        const double planes = static_cast<double>(scenario.rigid_settings.ground.size());
        const double point = point2_fixed + sweeps * point2_sweep;
        const double box_pairs = static_cast<double>(pairs_with(count.boxes, count.immovable_boxes));
        const double round_pairs =
            static_cast<double>(pairs_with(count.rounds, count.immovable_rounds));
        const double box_points = static_cast<double>(box2_box2_points);
        const double ground_points = static_cast<double>(box2_ground_points);
        return rigid2_body_work * n + no_contact_work * (all_pairs - box_pairs - round_pairs) +
               box_pairs * (box2_box2_detect + box_points * point) +
               round_pairs * (disc2_pair_detect + point) +
               planes * (boxes * (box2_ground_detect + ground_points * point) +
                         rounds * (disc2_ground_detect + point) +
                         (n - boxes - rounds) * no_contact_work);
    }
    return n * (n - 1.0) / 2.0 + n;
}

std::uint64_t potential_contact_points(const Scenario& scenario)
{
    ColliderCount count;
    count_colliders(scenario, count);
    return contact_points(scenario, count);
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
