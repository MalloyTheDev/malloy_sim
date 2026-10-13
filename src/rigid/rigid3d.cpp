#include <malloy/rigid/rigid.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include <malloy/collide/collide.hpp>
#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>

namespace malloy::rigid
{
namespace
{
// ============================================================================
// M40: the per-point iterative contact solver (issue #73).
//
// Every constant below is defined once, here, with its derivation; the step
// contract in rigid3d_world.hpp states each of them. The measurements quoted
// are from the M40 research prototype and the tests that pin them.
// ============================================================================

// Speculative margin, as a fraction of the smaller body's size: a pair still
// apart by up to this is kept as a contact with a negative penetration (the
// gap), and the rest phase lets it close the gap but not pass it. It must
// exceed any motion of a body at rest (the rest tolerance below leaves at most
// about 1e-8 g dt of velocity, under 1e-12 of a size per step) and the slop,
// so a resting point inside the slop band is always inside the margin; and it
// must be invisible, since a bounce can start from a gap this wide. Without it
// a resting corner rounds clear of the floor, drops out of the manifold, and
// falls for a step: a three-box stack then drifts 3.9e-5 of a size in 10 s
// and a box on a slope creeps 2.2e-5 (0.39 at dt 0.01). Any value from 1e-6 to
// 1e-2 measured the same; 1e-4 sits two decades inside both ends.
constexpr math::Real kMarginFraction = 1e-4;

// Allowed penetration, as a fraction of the smaller body's size, below which
// the position projection does not push. At least 1e10 times position rounding
// (about 1e-16 of a size), so the projection never chases rounding (a slop of
// 0 cost 18 sweeps a step instead of 2); and no larger than the precision a
// template states for a resting height (bouncing_sphere's position_z, 1e-6 at
// radius 0.5). Dynamics measured the same from 1e-9 to 1e-2.
constexpr math::Real kSlopFraction = 1e-6;

// Fraction of the remaining penetration (beyond the slop) one projection push
// removes, Box2D's position factor. Below 1 so several points that share a
// body do not overshoot each other. Measured the same from 0.1 to 1.
constexpr math::Real kProjectionFactor = 0.2;

// Projection passes per step, and the largest single push as a fraction of the
// size (Box2D's 3 and 0.2 m for metre-sized bodies). Measured the same from 1
// to 10 passes; the cap only matters for a deep overlap, which it spreads over
// a few steps instead of launching a body.
constexpr int kProjectionPasses = 3;
constexpr math::Real kMaxCorrectionFraction = 0.2;

// Sweep termination: stop when no point's constraint is violated by more than
// this times the velocity scale (the largest speed in contact, or one step of
// the external acceleration, whichever is larger). Rest: 1e-8. At 1e-6 a
// symmetric drop kept a 1.4e-9 rad/s spin (box_drop.scn checks 1e-9); 1e-8
// gives 2.7e-11 for about 2 sweeps a step at rest. Impacts, which are rare and
// violent and decide whether a symmetric landing stays symmetric: 1e-12
// (about 60 sweeps for a frictional landing). Both stop at 100 sweeps, a cap
// short stacks never reach at rest; a stack of eight twisted boxes needs about
// 300 and is beyond this milestone.
constexpr math::Real kRestTolerance = 1e-8;
constexpr math::Real kImpactTolerance = 1e-12;
constexpr int kMaxSweeps = 100;

// Restitution threshold, as a multiple of dt (|g| + |F| / m_min): a point
// approaching slower than this at mid-step is not an impact, and the rest phase
// stops it instead. A free body's resting contact approaches at |a| dt / 2 at
// mid-step, and one that fell clear for a step at 1.5 |a| dt. But a point of a
// body PIVOTING on another contact falls faster than the body: about a pivot P
// the far point Q accelerates at |a| d_G d_Q / (I_P / m) <= |a| d_Q / d_G
// (d the horizontal distances from P, I_P >= m d_G^2), and a box's centre of
// mass lies over the middle of its face, so the corner opposite the pivot has
// d_Q = 2 d_G and falls at up to 2 |a|. Fell clear for a step, it approaches at
// 1.5 x 2 = 3 |a| dt; four clears that with the 4/3 margin that two had over
// 1.5. At two, a box nudged at e 1 and dt 0.002 kept rattling on alternate
// corners, its strikes resolved a step late and so never slow enough to rest:
// 6 of 48 nudges with friction 0.5 and 32 of 48 frictionless (up to 1.5e-2 J
// after 15 s). At four, all 48 rest with friction and the frictionless ones are
// below 2e-10 J. Issue #88 changed the rotation this was tuned on, so it was
// re-measured on a reconstruction (that setup is not recorded here: a cube and
// a 0.5 x 0.3 x 0.2 box, 48 spin nudges each, e 1, dt 0.002, 15 s). Above
// 2e-10 J at the end: cube 12 and 24 of 48 (friction 0.5 and 0), box 0 and 16
// on the old rotation; cube 12 and 24, box 0 and 8 on the new. The
// reconstruction does not reproduce the counts above even on the old rotation,
// so it differs from the original setup; old against new is like for like.
constexpr math::Real kThresholdFactor = 4.0;

// Warm-start match radius for two boxes, as a fraction of body a's size: the
// nearest last-step point within it, in a's own frame, lends its impulses.
// Clipped manifold points of a resting pair repeat to rounding; a new point
// this far from every old one is a new contact and starts from 0.
constexpr math::Real kWarmMatchFraction = 0.05;

// A manifold is LEVEL when its depths agree to within this fraction of the
// boxes' size, the collide layer's touch-tolerance convention: its centroid
// weights are then equal rather than depth-weighted (issue #41).
constexpr math::Real kLevelTolerance = 1e-9;

// The friction disc projection (project_to_disc) stops when the projected
// impulse is within this fraction of the disc's radius, a few roundings, or
// after this many iterations. Newton from the bracket's low end converges in
// about 3 to 5 for any tangent effective mass; the cap is the bisection fallback
// that halves a bracket at most k_max / k_min wide to rounding.
constexpr math::Real kDiscTolerance = 1e-14;
constexpr int kDiscIterations = 60;

// The world-frame inverse inertia applied to a vector: R I^-1 R^T x. The stored
// inertia is the diagonal of principal moments in the body frame, so this
// rotates x into the body frame, divides component-wise, and rotates back. For
// an isotropic (spherical) body it reduces to x / I.
math::Vec3 inverse_inertia_world(const RigidBody3D& body, const math::Vec3& x)
{
    const math::Vec3 body_frame = math::rotate(math::conjugate(body.orientation), x);
    const math::Vec3 scaled{body_frame.x / body.inertia.x,
                            body_frame.y / body.inertia.y,
                            body_frame.z / body.inertia.z};
    return math::rotate(body.orientation, scaled);
}

// A body collides as an oriented box when all three half-extents are positive;
// otherwise it uses the sphere radius. Zero half-extents (the default) is not a
// box, so a pre-M30 body keeps its sphere behaviour.
bool body_is_box(const RigidBody3D& body)
{
    return body.half_extents.x > math::Real{0} &&
           body.half_extents.y > math::Real{0} && body.half_extents.z > math::Real{0};
}

// The size every contact tolerance scales with: a box's largest half-width, or
// a sphere's radius (0 for a body that does not collide).
math::Real collider_size(const RigidBody3D& body)
{
    if (body_is_box(body))
    {
        return std::fmax(std::fmax(body.half_extents.x, body.half_extents.y),
                         body.half_extents.z);
    }
    return body.radius;
}

// One point of a contact manifold, for one step. `b` is the other body, or the
// plane index when `ground` (an immovable participant with zero inverse mass).
struct ContactPoint
{
    std::size_t a{0};
    std::size_t b{0};
    bool ground{false};
    int key{-1};                // warm-start feature (see WarmContact)
    math::Vec3 normal{};        // unit, from a toward b
    math::Vec3 arm_a{};         // a's centre of mass to the point (world)
    math::Vec3 arm_b{};
    math::Vec3 local{};         // the point in a's body frame
    math::Vec3 anchor_a{};      // projection anchors: a's surface point in its
    math::Vec3 anchor_b{};      // frame (a sphere's centre), b's likewise or,
    math::Real radius_a{0.0};   // for the ground, the world point on the plane;
    math::Real radius_b{0.0};   // a sphere adds its radius along the normal
    math::Real penetration{0.0}; // signed: negative is a speculative gap
    math::Real weight{0.0};     // share of the manifold's centroid
    math::Real normal_mass{0.0}; // k along the normal
    math::Vec3 tangent1{};
    math::Vec3 tangent2{};
    math::Real k11{0.0};        // the 2x2 tangent block of k
    math::Real k12{0.0};
    math::Real k22{0.0};
    math::Real target{0.0};     // the normal velocity this phase drives toward
    math::Real normal_impulse{0.0};
    math::Vec3 tangent_impulse{};
    bool active{true};
};

// The points [first, first + count) of one body pair's manifold.
struct Manifold
{
    std::size_t first{0};
    std::size_t count{0};
    math::Real scale{0.0};      // the smaller body's size: margin, slop, push cap
    math::Real level_size{0.0}; // the larger: the level test
};

// A body's velocities while a phase runs: the world-frame angular velocity at
// the start of the phase plus the change the phase has made, kept apart so a
// body's own body-frame spin is only ever adjusted by that change.
struct SolverBody
{
    math::Vec3 velocity{};
    math::Vec3 spin{};          // world frame, at the phase start
    math::Vec3 spin_change{};   // world frame, accumulated this phase
    math::Mat3 inverse_inertia{};
    math::Real inverse_mass{0.0};
    bool touched{false};
};

struct Contacts
{
    std::vector<ContactPoint> points;
    std::vector<Manifold> manifolds;
};

void add_point(const std::vector<RigidBody3D>& bodies, Contacts& contacts, std::size_t a,
               std::size_t b, bool ground, int key, const math::Vec3& normal,
               const math::Vec3& point, math::Real penetration, const math::Vec3& arm_a,
               const math::Vec3& arm_b)
{
    ContactPoint c;
    c.a = a;
    c.b = b;
    c.ground = ground;
    c.key = key;
    c.normal = normal;
    c.penetration = penetration;
    c.arm_a = arm_a;
    c.arm_b = arm_b;
    const RigidBody3D& body_a = bodies[a];
    c.local = math::rotate(math::conjugate(body_a.orientation), point - body_a.position);
    // Each body's own surface point: the contact point is midway, so a's lies
    // half the penetration along the normal, b's half against it.
    const math::Vec3 surface_a = point + normal * (penetration / math::Real{2});
    const math::Vec3 surface_b = point - normal * (penetration / math::Real{2});
    if (body_is_box(body_a))
    {
        c.anchor_a = math::rotate(math::conjugate(body_a.orientation),
                                  surface_a - body_a.position);
    }
    else
    {
        c.radius_a = body_a.radius; // the centre, plus the radius along the normal
    }
    if (ground)
    {
        c.anchor_b = surface_b;
    }
    else
    {
        const RigidBody3D& body_b = bodies[b];
        if (body_is_box(body_b))
        {
            c.anchor_b = math::rotate(math::conjugate(body_b.orientation),
                                      surface_b - body_b.position);
        }
        else
        {
            c.radius_b = body_b.radius;
        }
    }
    contacts.points.push_back(c);
}

// Closes the manifold opened at `first`: its centroid weights are each point's
// penetration share (negatives count as 0), or equal when the manifold is level
// or nothing penetrates, as reduce_manifold weighted them before M40.
void close_manifold(Contacts& contacts, std::size_t first, math::Real scale,
                    math::Real level_size)
{
    const std::size_t count = contacts.points.size() - first;
    if (count == 0)
    {
        return;
    }
    math::Real deepest = -std::numeric_limits<math::Real>::infinity();
    math::Real shallowest = std::numeric_limits<math::Real>::infinity();
    math::Real total = 0.0;
    for (std::size_t k = first; k < contacts.points.size(); ++k)
    {
        const math::Real depth = contacts.points[k].penetration;
        deepest = std::fmax(deepest, depth);
        shallowest = std::fmin(shallowest, depth);
        total += std::fmax(math::Real{0}, depth);
    }
    const bool level = !(deepest - shallowest > kLevelTolerance * level_size);
    for (std::size_t k = first; k < contacts.points.size(); ++k)
    {
        ContactPoint& c = contacts.points[k];
        c.weight = (level || !(total > math::Real{0}))
                       ? math::Real{1} / static_cast<math::Real>(count)
                       : std::fmax(math::Real{0}, c.penetration) / total;
    }
    contacts.manifolds.push_back(Manifold{first, count, scale, level_size});
}

// Every contact of the step, at the current poses, in the fixed order of the
// step contract.
void detect(const std::vector<RigidBody3D>& bodies, const Rigid3DSettings& settings,
            Contacts& contacts)
{
    const std::size_t count = bodies.size();
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            const RigidBody3D& a = bodies[i];
            const RigidBody3D& b = bodies[j];
            const math::Real size_a = collider_size(a);
            const math::Real size_b = collider_size(b);
            if (!(size_a > math::Real{0}) || !(size_b > math::Real{0}))
            {
                continue; // a body without a collider takes part in nothing
            }
            const math::Real scale = std::fmin(size_a, size_b);
            const math::Real level_size = std::fmax(size_a, size_b);
            const math::Real margin = kMarginFraction * scale;
            const std::size_t first = contacts.points.size();
            if (body_is_box(a) && body_is_box(b))
            {
                const std::vector<collide::Contact3> hits = collide::contacts(
                    collide::Box3{a.position, a.half_extents, a.orientation},
                    collide::Box3{b.position, b.half_extents, b.orientation}, margin);
                for (const collide::Contact3& hit : hits)
                {
                    add_point(bodies, contacts, i, j, false, -1, hit.normal, hit.point,
                              hit.penetration, hit.point - a.position,
                              hit.point - b.position);
                }
            }
            else if (body_is_box(a) || body_is_box(b))
            {
                const RigidBody3D& box = body_is_box(a) ? a : b;
                const RigidBody3D& sphere = body_is_box(a) ? b : a;
                const std::optional<collide::Contact3> hit = collide::contact(
                    collide::Box3{box.position, box.half_extents, box.orientation},
                    collide::Sphere{sphere.position, sphere.radius}, margin);
                if (hit)
                {
                    // The query's normal runs from the box toward the sphere;
                    // flip it when the box is b, so it runs from a toward b.
                    const math::Vec3 normal = body_is_box(a) ? hit->normal : -hit->normal;
                    add_point(bodies, contacts, i, j, false, 0, normal, hit->point,
                              hit->penetration, hit->point - a.position,
                              hit->point - b.position);
                }
            }
            else
            {
                const std::optional<collide::Contact3> hit =
                    collide::contact(collide::Sphere{a.position, a.radius},
                                     collide::Sphere{b.position, b.radius}, margin);
                if (hit)
                {
                    // Both arms reach the ONE contact point, so the equal and
                    // opposite impulses conserve angular momentum (issue #35).
                    add_point(bodies, contacts, i, j, false, 0, hit->normal, hit->point,
                              hit->penetration, hit->point - a.position,
                              hit->point - b.position);
                }
            }
            close_manifold(contacts, first, scale, level_size);
        }
    }
    for (std::size_t i = 0; i < count; ++i)
    {
        const RigidBody3D& body = bodies[i];
        const math::Real size = collider_size(body);
        if (!(size > math::Real{0}))
        {
            continue;
        }
        const math::Real margin = kMarginFraction * size;
        for (std::size_t k = 0; k < settings.ground.size(); ++k)
        {
            const collide::Plane3& plane = settings.ground[k];
            const std::size_t first = contacts.points.size();
            if (body_is_box(body))
            {
                const std::vector<collide::Contact3> hits = collide::contacts(
                    collide::Box3{body.position, body.half_extents, body.orientation},
                    plane, margin);
                for (const collide::Contact3& hit : hits)
                {
                    // The corner this point came from, by the signs of its body
                    // coordinates: the warm-start key a resting corner keeps.
                    const math::Vec3 corner =
                        hit.point + hit.normal * (hit.penetration / math::Real{2});
                    const math::Vec3 local = math::rotate(
                        math::conjugate(body.orientation), corner - body.position);
                    const int key = (local.x > math::Real{0} ? 1 : 0) |
                                    (local.y > math::Real{0} ? 2 : 0) |
                                    (local.z > math::Real{0} ? 4 : 0);
                    add_point(bodies, contacts, i, k, true, key, hit.normal, hit.point,
                              hit.penetration, hit.point - body.position, math::Vec3{});
                }
            }
            else
            {
                const std::optional<collide::Contact3> hit = collide::contact(
                    collide::Sphere{body.position, body.radius}, plane, margin);
                if (hit)
                {
                    // The arm reaches the sphere's own surface point, exactly
                    // radius along the normal (M22/M23), which keeps rolling
                    // exact whatever the penetration.
                    add_point(bodies, contacts, i, k, true, 0, hit->normal, hit->point,
                              hit->penetration, hit->normal * body.radius, math::Vec3{});
                }
            }
            close_manifold(contacts, first, size, size);
        }
    }
}

// The world inverse inertia R diag(1/I) R^T, computed once per phase.
math::Mat3 world_inverse_inertia(const RigidBody3D& body)
{
    const math::Mat3 r = math::to_mat3(body.orientation);
    const math::Vec3 inverse{math::Real{1} / body.inertia.x, math::Real{1} / body.inertia.y,
                             math::Real{1} / body.inertia.z};
    return r * math::diagonal3(inverse) * math::transpose(r);
}

std::vector<SolverBody> solver_bodies(const std::vector<RigidBody3D>& bodies,
                                      const Contacts& contacts)
{
    std::vector<SolverBody> solver(bodies.size());
    std::vector<bool> ready(bodies.size(), false);
    const auto prepare = [&](std::size_t i)
    {
        if (ready[i])
        {
            return;
        }
        ready[i] = true;
        const RigidBody3D& body = bodies[i];
        SolverBody& s = solver[i];
        s.velocity = body.velocity;
        s.spin = math::rotate(body.orientation, body.angular_velocity);
        s.inverse_mass = math::Real{1} / body.mass;
        s.inverse_inertia = world_inverse_inertia(body);
    };
    for (const ContactPoint& c : contacts.points)
    {
        prepare(c.a);
        if (!c.ground)
        {
            prepare(c.b);
        }
    }
    return solver;
}

// Velocity of a participant's material point at `arm`; zero for the ground.
math::Vec3 point_velocity(const std::vector<SolverBody>& solver, std::size_t i, bool ground,
                          const math::Vec3& arm)
{
    if (ground)
    {
        return math::Vec3{};
    }
    const SolverBody& s = solver[i];
    return s.velocity + math::cross(s.spin + s.spin_change, arm);
}

math::Vec3 relative_velocity(const std::vector<SolverBody>& solver, const ContactPoint& c,
                             const math::Vec3& arm_a, const math::Vec3& arm_b)
{
    return point_velocity(solver, c.b, c.ground, arm_b) -
           point_velocity(solver, c.a, false, arm_a);
}

// Apply `impulse` at `arm`: v += J/m, w += I^-1 (r x J). A zero impulse leaves
// the body untouched, so it is not written back.
void apply(std::vector<SolverBody>& solver, std::size_t i, bool ground,
           const math::Vec3& impulse, const math::Vec3& arm)
{
    if (ground || (impulse.x == math::Real{0} && impulse.y == math::Real{0} &&
                   impulse.z == math::Real{0}))
    {
        return;
    }
    SolverBody& s = solver[i];
    s.velocity += impulse * s.inverse_mass;
    s.spin_change += s.inverse_inertia * math::cross(arm, impulse);
    s.touched = true;
}

// Equal and opposite: -impulse on a, +impulse on b, at their own arms.
void apply_pair(std::vector<SolverBody>& solver, const ContactPoint& c,
                const math::Vec3& impulse, const math::Vec3& arm_a, const math::Vec3& arm_b)
{
    apply(solver, c.a, false, -impulse, arm_a);
    apply(solver, c.b, c.ground, impulse, arm_b);
}

// The effective-mass form k(d1, d2) = (1/m_a + 1/m_b) d1.d2
// + d1.[(I_a^-1 (r_a x d2)) x r_a] + the same for b. k(n, n) is the normal's
// effective mass; it is symmetric, so the tangent block is too.
math::Real effective_mass(const std::vector<SolverBody>& solver, const ContactPoint& c,
                          const math::Vec3& arm_a, const math::Vec3& arm_b,
                          const math::Vec3& d1, const math::Vec3& d2)
{
    const SolverBody& sa = solver[c.a];
    math::Real k = sa.inverse_mass * math::dot(d1, d2) +
                   math::dot(d1, math::cross(sa.inverse_inertia * math::cross(arm_a, d2), arm_a));
    if (!c.ground)
    {
        const SolverBody& sb = solver[c.b];
        k += sb.inverse_mass * math::dot(d1, d2) +
             math::dot(d1, math::cross(sb.inverse_inertia * math::cross(arm_b, d2), arm_b));
    }
    return k;
}

// A fixed tangent basis from the normal alone, so a run repeats; the friction
// solve below is a 2x2 block with a circular clamp, so the result does not
// depend on the basis beyond rounding.
void tangent_basis(const math::Vec3& n, math::Vec3& t1, math::Vec3& t2)
{
    if (std::abs(n.x) >= math::Real{0.57735026918962576}) // 1/sqrt(3)
    {
        t1 = math::normalize(math::Vec3{n.y, -n.x, math::Real{0}});
    }
    else
    {
        t1 = math::normalize(math::Vec3{math::Real{0}, n.z, -n.y});
    }
    t2 = math::cross(n, t1);
}

void prepare_points(const std::vector<SolverBody>& solver, Contacts& contacts)
{
    for (ContactPoint& c : contacts.points)
    {
        c.normal_mass = effective_mass(solver, c, c.arm_a, c.arm_b, c.normal, c.normal);
        tangent_basis(c.normal, c.tangent1, c.tangent2);
        c.k11 = effective_mass(solver, c, c.arm_a, c.arm_b, c.tangent1, c.tangent1);
        c.k12 = effective_mass(solver, c, c.arm_a, c.arm_b, c.tangent1, c.tangent2);
        c.k22 = effective_mass(solver, c, c.arm_a, c.arm_b, c.tangent2, c.tangent2);
        c.normal_impulse = math::Real{0};
        c.tangent_impulse = math::Vec3{};
    }
}

// The point of the friction disc |p| <= limit closest to `p` (tangent
// coordinates) in the metric of the tangent effective mass K = [k11 k12; k12
// k22]: p(g) = (K + g 1)^-1 K p with g >= 0 the root of |p(g)| = limit.
//
// This is Coulomb's law in its maximum-dissipation form. Of every impulse the
// disc allows it leaves the least kinetic energy, and the slide it leaves,
// K (p(g) - p), equals -g p(g): exactly opposite the friction. Scaling p
// radially onto the circle (the Euclidean nearest point) agrees only when K is
// isotropic; a box's K is not, and radial scaling pushed a box sliding
// diagonally sideways (0.15 m/s after 0.3 s) and spun it.
//
// The root is bracketed in closed form: k/(k + g) grows with k, so with q = |p|
// and s = q / limit - 1, |p(g)| >= limit at g = k_min s and <= limit at
// g = k_max s (k_min, k_max the eigenvalues of K). An isotropic K closes the
// bracket to one point, so the result is then exactly the radial one. Newton on
// 1/|p(g)| - 1/limit, nearly linear in g (the trust-region secular equation),
// kept inside the bracket by bisection, then scaled onto the circle.
void project_to_disc(math::Real k11, math::Real k12, math::Real k22, math::Real limit,
                     math::Real& p1, math::Real& p2)
{
    const math::Real q = std::sqrt(p1 * p1 + p2 * p2);
    if (!(q > limit))
    {
        return;
    }
    if (!(limit > math::Real{0}))
    {
        p1 = math::Real{0};
        p2 = math::Real{0};
        return;
    }
    const math::Real b1 = k11 * p1 + k12 * p2;
    const math::Real b2 = k12 * p1 + k22 * p2;
    const math::Real half_sum = (k11 + k22) / math::Real{2};
    const math::Real radius = std::sqrt((k11 - k22) * (k11 - k22) / math::Real{4} + k12 * k12);
    const math::Real k_max = half_sum + radius;
    const math::Real k_min = (k11 * k22 - k12 * k12) / k_max; // no cancellation
    const math::Real s = q / limit - math::Real{1};
    math::Real lo = k_min * s;
    math::Real hi = k_max * s;
    math::Real g = lo;
    math::Real o1 = p1;
    math::Real o2 = p2;
    for (int iteration = 0; iteration < kDiscIterations; ++iteration)
    {
        const math::Real a = k11 + g;
        const math::Real d = k22 + g;
        const math::Real det = a * d - k12 * k12;
        o1 = (d * b1 - k12 * b2) / det;
        o2 = (-k12 * b1 + a * b2) / det;
        const math::Real length = std::sqrt(o1 * o1 + o2 * o2);
        const math::Real error = length - limit;
        if (std::abs(error) <= kDiscTolerance * limit || !(hi > lo))
        {
            break;
        }
        if (error > math::Real{0})
        {
            lo = g;
        }
        else
        {
            hi = g;
        }
        // d|p|/dg = -p.(K + g)^-1 p / |p|, so the Newton step on
        // 1/|p| - 1/limit is (1/|p| - 1/limit) |p|^3 / (p.(K + g)^-1 p).
        const math::Real m1 = (d * o1 - k12 * o2) / det;
        const math::Real m2 = (-k12 * o1 + a * o2) / det;
        const math::Real curvature = o1 * m1 + o2 * m2;
        math::Real next = curvature > math::Real{0}
                              ? g + (math::Real{1} / length - math::Real{1} / limit) * length *
                                        length * length / curvature
                              : (lo + hi) / math::Real{2};
        if (!(next > lo && next < hi))
        {
            next = (lo + hi) / math::Real{2};
        }
        g = next;
    }
    const math::Real length = std::sqrt(o1 * o1 + o2 * o2);
    p1 = length > math::Real{0} ? o1 * (limit / length) : math::Real{0};
    p2 = length > math::Real{0} ? o2 * (limit / length) : math::Real{0};
}

// The friction step of one point, or of a manifold's centroid: the 2x2 tangent
// solve for the impulse that stops the slide, then the ACCUMULATED impulse
// brought into the disc of radius `limit` by project_to_disc. Returns the change
// to apply.
math::Vec3 friction_change(const math::Vec3& relative, const math::Vec3& t1,
                           const math::Vec3& t2, math::Real k11, math::Real k12,
                           math::Real k22, math::Vec3& accumulated, math::Real limit)
{
    const math::Real determinant = k11 * k22 - k12 * k12;
    if (!(determinant > math::Real{0}))
    {
        return math::Vec3{};
    }
    const math::Real u1 = math::dot(relative, t1);
    const math::Real u2 = math::dot(relative, t2);
    const math::Real d1 = -(k22 * u1 - k12 * u2) / determinant;
    const math::Real d2 = -(-k12 * u1 + k11 * u2) / determinant;
    math::Vec3 next = accumulated + t1 * d1 + t2 * d2;
    const math::Real length = math::length(next);
    if (length > limit)
    {
        if (!(limit > math::Real{0}))
        {
            next = math::Vec3{};
        }
        else if (k12 == math::Real{0} && k11 == k22)
        {
            next = next * (limit / length); // isotropic: the radial point is the closest
        }
        else
        {
            math::Real p1 = math::dot(next, t1);
            math::Real p2 = math::dot(next, t2);
            project_to_disc(k11, k12, k22, limit, p1, p2);
            next = t1 * p1 + t2 * p2;
        }
    }
    const math::Vec3 change = next - accumulated;
    accumulated = next;
    return change;
}

// One impulse at the manifold's weighted centroid before the per-point sweeps.
// For a level face landing flat it is the whole answer, applied in one piece,
// so the landing stays exactly symmetric (the four corners resolved one at a
// time would spin the box); otherwise it is a starting guess the sweeps refine.
// Its share is booked to each point's accumulated impulse by weight.
void prepass(std::vector<SolverBody>& solver, Contacts& contacts, math::Real friction)
{
    for (const Manifold& m : contacts.manifolds)
    {
        if (m.count < 2)
        {
            continue;
        }
        math::Real weight_sum = 0.0;
        std::size_t active = 0;
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            if (contacts.points[k].active)
            {
                weight_sum += contacts.points[k].weight;
                ++active;
            }
        }
        if (active < 2 || !(weight_sum > math::Real{0}))
        {
            continue;
        }
        math::Vec3 arm_a{};
        math::Vec3 arm_b{};
        math::Real target = 0.0;
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            const ContactPoint& c = contacts.points[k];
            if (!c.active)
            {
                continue;
            }
            const math::Real share = c.weight / weight_sum;
            arm_a += c.arm_a * share;
            arm_b += c.arm_b * share;
            target += c.target * share;
        }
        ContactPoint& head = contacts.points[m.first];
        const math::Vec3 n = head.normal;
        const math::Real mass = effective_mass(solver, head, arm_a, arm_b, n, n);
        if (!(mass > math::Real{0}))
        {
            continue;
        }
        // The centroid impulse is clamped as an accumulated impulse, like every
        // other row: it may take back what the warm start applied, but no
        // point's booked impulse may go below 0. So a warm start larger than
        // this step needs (a landing's impulse remembered into the first resting
        // step) is reduced here, in one symmetric piece, rather than left to the
        // per-point sweeps, whose order would leave a spin at the tolerance.
        const math::Real approach = math::dot(relative_velocity(solver, head, arm_a, arm_b), n);
        math::Real lambda = (target - approach) / mass;
        math::Vec3 centroid_friction{};
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            const ContactPoint& c = contacts.points[k];
            if (c.active && c.weight > math::Real{0})
            {
                lambda = std::fmax(lambda, -c.normal_impulse * weight_sum / c.weight);
            }
            if (c.active)
            {
                centroid_friction += c.tangent_impulse;
            }
        }
        apply_pair(solver, head, n * lambda, arm_a, arm_b);
        math::Real total = 0.0;
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            ContactPoint& c = contacts.points[k];
            if (c.active)
            {
                c.normal_impulse =
                    std::fmax(math::Real{0}, c.normal_impulse + c.weight / weight_sum * lambda);
                total += c.normal_impulse;
            }
        }
        if (!(friction > math::Real{0}))
        {
            continue;
        }
        // Likewise accumulated: the manifold's total friction (what the warm
        // start applied plus this) stays within friction times its total push.
        const math::Vec3 change = friction_change(
            relative_velocity(solver, head, arm_a, arm_b), head.tangent1, head.tangent2,
            effective_mass(solver, head, arm_a, arm_b, head.tangent1, head.tangent1),
            effective_mass(solver, head, arm_a, arm_b, head.tangent1, head.tangent2),
            effective_mass(solver, head, arm_a, arm_b, head.tangent2, head.tangent2),
            centroid_friction, friction * total);
        apply_pair(solver, head, change, arm_a, arm_b);
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            ContactPoint& c = contacts.points[k];
            if (c.active)
            {
                c.tangent_impulse += change * (c.weight / weight_sum);
            }
        }
    }
}

// One Gauss-Seidel sweep: per manifold, every normal row, then every friction
// row. Accumulated impulses are clamped, never the change: the normal impulse
// never below 0, the friction impulse within friction times the point's own
// accumulated normal impulse.
void sweep(std::vector<SolverBody>& solver, Contacts& contacts, math::Real friction)
{
    for (const Manifold& m : contacts.manifolds)
    {
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            ContactPoint& c = contacts.points[k];
            if (!c.active || !(c.normal_mass > math::Real{0}))
            {
                continue;
            }
            const math::Real approach =
                math::dot(relative_velocity(solver, c, c.arm_a, c.arm_b), c.normal);
            const math::Real next =
                std::fmax(math::Real{0}, c.normal_impulse + (c.target - approach) / c.normal_mass);
            const math::Real change = next - c.normal_impulse;
            c.normal_impulse = next;
            apply_pair(solver, c, c.normal * change, c.arm_a, c.arm_b);
        }
        if (!(friction > math::Real{0}))
        {
            continue;
        }
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            ContactPoint& c = contacts.points[k];
            if (!c.active)
            {
                continue;
            }
            const math::Vec3 change = friction_change(
                relative_velocity(solver, c, c.arm_a, c.arm_b), c.tangent1, c.tangent2, c.k11,
                c.k12, c.k22, c.tangent_impulse, friction * c.normal_impulse);
            apply_pair(solver, c, change, c.arm_a, c.arm_b);
        }
    }
}

// The largest constraint violation: a pushing point off its target, a slack
// point approaching faster than its target allows, a sticking point still
// sliding, or a sliding point whose friction does not oppose the slide.
math::Real residual(const std::vector<SolverBody>& solver, const Contacts& contacts,
                    math::Real friction)
{
    math::Real worst = 0.0;
    for (const ContactPoint& c : contacts.points)
    {
        if (!c.active)
        {
            continue;
        }
        const math::Vec3 relative = relative_velocity(solver, c, c.arm_a, c.arm_b);
        const math::Real approach = math::dot(relative, c.normal);
        worst = std::fmax(worst, c.normal_impulse > math::Real{0}
                                     ? std::abs(approach - c.target)
                                     : std::fmax(math::Real{0}, c.target - approach));
        if (friction > math::Real{0})
        {
            const math::Vec3 slide = relative - c.normal * approach;
            const math::Real limit = friction * c.normal_impulse;
            const math::Real held = math::length(c.tangent_impulse);
            if (held < limit * (math::Real{1} - math::Real{1e-12}))
            {
                worst = std::fmax(worst, math::length(slide));
            }
            else if (held > math::Real{0})
            {
                worst = std::fmax(worst, std::fmax(math::Real{0},
                                                   math::dot(slide, c.tangent_impulse) / held));
            }
        }
    }
    return worst;
}

// The velocity scale the tolerance is relative to: the fastest material point
// in contact, or one step of the external acceleration, whichever is larger.
math::Real velocity_scale(const std::vector<RigidBody3D>& bodies,
                          const std::vector<SolverBody>& solver, const Contacts& contacts,
                          math::Real forcing_speed)
{
    math::Real scale = forcing_speed;
    for (const ContactPoint& c : contacts.points)
    {
        for (int side = 0; side < 2; ++side)
        {
            if (side == 1 && c.ground)
            {
                continue;
            }
            const std::size_t i = side == 0 ? c.a : c.b;
            const SolverBody& s = solver[i];
            scale = std::fmax(scale, math::length(s.velocity) +
                                         math::length(s.spin) * collider_size(bodies[i]));
        }
    }
    return scale;
}

// Runs sweeps until the residual is within `tolerance`, or kMaxSweeps. The
// residual is checked BEFORE the first sweep too: when the warm start and the
// centroid prepass already solve the manifold (a symmetric landing, a stack at
// rest), no per-point sweep runs, so nothing order-dependent touches it and it
// stays exactly symmetric.
void iterate(std::vector<SolverBody>& solver, Contacts& contacts, math::Real friction,
             math::Real tolerance)
{
    for (int pass = 0; pass < kMaxSweeps; ++pass)
    {
        if (residual(solver, contacts, friction) <= tolerance)
        {
            break;
        }
        sweep(solver, contacts, friction);
    }
}

// Writes back only the bodies a nonzero impulse reached. The body-frame spin is
// adjusted by the change alone, so a centred impulse (no spin change) leaves it
// bit for bit as it was.
void write_back(std::vector<RigidBody3D>& bodies, const std::vector<SolverBody>& solver)
{
    for (std::size_t i = 0; i < bodies.size(); ++i)
    {
        const SolverBody& s = solver[i];
        if (!s.touched)
        {
            continue;
        }
        RigidBody3D& body = bodies[i];
        body.velocity = s.velocity;
        if (s.spin_change.x != math::Real{0} || s.spin_change.y != math::Real{0} ||
            s.spin_change.z != math::Real{0})
        {
            body.angular_velocity +=
                math::rotate(math::conjugate(body.orientation), s.spin_change);
        }
    }
}

// The external acceleration scale: |g| + |F| / m_min, the most one step of
// forcing can change a resting contact's normal velocity, per unit dt.
math::Real forcing_acceleration(const std::vector<RigidBody3D>& bodies,
                                const Rigid3DSettings& settings)
{
    math::Real lightest = std::numeric_limits<math::Real>::infinity();
    for (const RigidBody3D& body : bodies)
    {
        lightest = std::fmin(lightest, body.mass);
    }
    const math::Real force = math::length(settings.force);
    return math::length(settings.gravity) +
           (force > math::Real{0} ? force / lightest : math::Real{0});
}

// Phase B: the impacts. On the start-of-step velocities, before any forcing.
void impact_phase(std::vector<RigidBody3D>& bodies, const Rigid3DSettings& settings,
                  math::Real dt, Contacts& contacts)
{
    const math::Real restitution = settings.restitution;
    std::vector<SolverBody> solver = solver_bodies(bodies, contacts);
    prepare_points(solver, contacts);
    const math::Real acceleration = forcing_acceleration(bodies, settings);
    const math::Real threshold = kThresholdFactor * dt * acceleration;
    bool any = false;
    for (ContactPoint& c : contacts.points)
    {
        // delta: what this step's forcing (steps 0 and 1) will add to the
        // point's normal velocity before the rest phase, from the same start
        // state the forcing reads. Since issue #88 step (1) is the torque kick
        // alone: the gyroscopic change of the world spin happens in the free
        // rotation of step (2), after the contact phases, and reaches the next
        // step's start velocity instead. Body b's change minus body a's, as in
        // relative_velocity; the ground (b) receives no forcing.
        const auto forcing_change = [&](std::size_t i, const math::Vec3& arm)
        {
            const RigidBody3D& body = bodies[i];
            const math::Vec3& I = body.inertia;
            const math::Vec3 torque_body =
                math::rotate(math::conjugate(body.orientation), settings.torque);
            const math::Vec3 rate{torque_body.x / I.x, torque_body.y / I.y, torque_body.z / I.z};
            const math::Vec3 spin_change = math::rotate(body.orientation, rate * dt);
            return settings.gravity * dt + (settings.force / body.mass) * dt +
                   math::cross(spin_change, arm);
        };
        math::Vec3 forcing = -forcing_change(c.a, c.arm_a);
        if (!c.ground)
        {
            forcing += forcing_change(c.b, c.arm_b);
        }
        const math::Real delta = math::dot(forcing, c.normal);
        const math::Real start =
            math::dot(relative_velocity(solver, c, c.arm_a, c.arm_b), c.normal);
        const math::Real middle = start + delta / math::Real{2};
        const math::Real gap = -c.penetration;
        c.active = middle < -threshold && (!(gap > math::Real{0}) || -middle * dt >= gap);
        // After the forcing adds delta, the point leaves at -e times its
        // mid-step approach: Newton's law at the middle of the step.
        c.target = -restitution * middle - delta / math::Real{2};
        any = any || c.active;
    }
    if (!any)
    {
        return;
    }
    prepass(solver, contacts, settings.friction);
    iterate(solver, contacts, settings.friction,
            kImpactTolerance *
                velocity_scale(bodies, solver, contacts, acceleration * dt));
    write_back(bodies, solver);
}

// Matches this step's points to last step's rest-phase impulses and applies
// them. Both lists are in detection order, sorted by (ground, a, b), so one
// merge finds each manifold's old points.
void warm_start(std::vector<SolverBody>& solver, Contacts& contacts,
                const std::vector<RigidBody3D>& bodies,
                const std::vector<Rigid3DWorld::WarmContact>& cache)
{
    const auto before = [](bool ground_x, std::size_t a_x, std::size_t b_x, bool ground_y,
                           std::size_t a_y, std::size_t b_y)
    {
        if (ground_x != ground_y)
        {
            return !ground_x;
        }
        if (a_x != a_y)
        {
            return a_x < a_y;
        }
        return b_x < b_y;
    };
    std::size_t cursor = 0;
    std::vector<bool> used;
    for (const Manifold& m : contacts.manifolds)
    {
        const ContactPoint& head = contacts.points[m.first];
        while (cursor < cache.size() && before(cache[cursor].ground, cache[cursor].a,
                                               cache[cursor].b, head.ground, head.a, head.b))
        {
            ++cursor;
        }
        std::size_t end = cursor;
        while (end < cache.size() && cache[end].ground == head.ground &&
               cache[end].a == head.a && cache[end].b == head.b)
        {
            ++end;
        }
        used.assign(end - cursor, false);
        const math::Real radius = kWarmMatchFraction * collider_size(bodies[head.a]);
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            ContactPoint& c = contacts.points[k];
            std::size_t best = end;
            math::Real best_distance = radius;
            for (std::size_t o = cursor; o < end; ++o)
            {
                if (used[o - cursor])
                {
                    continue;
                }
                if (c.key >= 0)
                {
                    if (cache[o].key == c.key)
                    {
                        best = o;
                        break;
                    }
                    continue;
                }
                const math::Real distance = math::length(cache[o].local - c.local);
                if (distance <= best_distance)
                {
                    best_distance = distance;
                    best = o;
                }
            }
            if (best == end)
            {
                continue;
            }
            used[best - cursor] = true;
            c.normal_impulse = cache[best].normal_impulse;
            // Last step's friction impulse, in this step's tangent plane.
            c.tangent_impulse = cache[best].tangent_impulse -
                                c.normal * math::dot(cache[best].tangent_impulse, c.normal);
            apply_pair(solver, c, c.normal * c.normal_impulse + c.tangent_impulse, c.arm_a,
                       c.arm_b);
        }
        cursor = end;
    }
}

// Phase C: every point, after the forcing, warm-started.
void rest_phase(std::vector<RigidBody3D>& bodies, const Rigid3DSettings& settings,
                math::Real dt, Contacts& contacts,
                const std::vector<Rigid3DWorld::WarmContact>& cache)
{
    std::vector<SolverBody> solver = solver_bodies(bodies, contacts);
    prepare_points(solver, contacts);
    for (ContactPoint& c : contacts.points)
    {
        c.active = true;
        // Hold a touching point at 0; let a point still apart close its gap
        // exactly within the step, and no further.
        c.target = c.penetration < math::Real{0} ? c.penetration / dt : math::Real{0};
    }
    warm_start(solver, contacts, bodies, cache);
    prepass(solver, contacts, settings.friction);
    iterate(solver, contacts, settings.friction,
            kRestTolerance * velocity_scale(bodies, solver, contacts,
                                            forcing_acceleration(bodies, settings) * dt));
    write_back(bodies, solver);
}

// A small world-frame turn of a body by `angle` (a rotation vector), the
// quaternion update the integrator uses, renormalised. The projection moves the
// pose, not the motion, so the body's WORLD-frame angular momentum is kept and
// re-expressed in the turned body frame (omega = I^-1 R'^T L). Keeping the
// body-frame angular velocity instead turns L with the body: on a frictionless
// floor a tilted cube picked up |L_z| = 3e-3 and never came to rest. Keeping
// the world-frame angular velocity is the same thing for an isotropic body (the
// cube) but not for a box with two equal moments: nudged at e 1 on a
// frictionless floor it kept |L_z| = 2.5e-8, a yaw of 1.5e-7 rad/s that nothing
// there can remove; keeping L leaves 2e-13. Contact impulses conserve L, so the
// projection is the step that must not change it.
void turn(RigidBody3D& body, const math::Vec3& angle)
{
    if (angle.x == math::Real{0} && angle.y == math::Real{0} && angle.z == math::Real{0})
    {
        return;
    }
    const math::Vec3& w = body.angular_velocity;
    const math::Vec3& I = body.inertia;
    const math::Vec3 momentum =
        math::rotate(body.orientation, math::Vec3{I.x * w.x, I.y * w.y, I.z * w.z});
    const math::Quat d{math::Real{0}, angle};
    body.orientation =
        math::normalize(body.orientation + (d * body.orientation) * math::Real{0.5});
    const math::Vec3 local = math::rotate(math::conjugate(body.orientation), momentum);
    body.angular_velocity = math::Vec3{local.x / I.x, local.y / I.y, local.z / I.z};
}

math::Vec3 anchor_a(const std::vector<RigidBody3D>& bodies, const ContactPoint& c)
{
    const RigidBody3D& body = bodies[c.a];
    return body.position + math::rotate(body.orientation, c.anchor_a) + c.normal * c.radius_a;
}

math::Vec3 anchor_b(const std::vector<RigidBody3D>& bodies, const ContactPoint& c)
{
    if (c.ground)
    {
        return c.anchor_b;
    }
    const RigidBody3D& body = bodies[c.b];
    return body.position + math::rotate(body.orientation, c.anchor_b) - c.normal * c.radius_b;
}

// One projection push of `correction` (negative: the depth to remove) along the
// point's normal, at the midpoint of the two anchors, moving and turning both
// bodies by inverse mass and inertia at their CURRENT pose. Velocities are not
// touched. A sphere is only moved, never turned: its contact normal passes
// through its centre, so the push has no lever arm, and leaving the rotation
// alone keeps a sphere's spin and orientation exactly those of a sphere with no
// contact (M22's rotationally inert contact).
void push(std::vector<RigidBody3D>& bodies, const ContactPoint& c, const math::Vec3& point_a,
          const math::Vec3& point_b, math::Real correction)
{
    const math::Vec3 middle = (point_a + point_b) * math::Real{0.5};
    const math::Vec3 n = c.normal;
    RigidBody3D& a = bodies[c.a];
    const bool turn_a = body_is_box(a);
    const math::Vec3 arm_a = middle - a.position;
    math::Real mass = math::Real{1} / a.mass;
    if (turn_a)
    {
        mass += math::dot(n, math::cross(inverse_inertia_world(a, math::cross(arm_a, n)), arm_a));
    }
    bool turn_b = false;
    math::Vec3 arm_b{};
    if (!c.ground)
    {
        const RigidBody3D& b = bodies[c.b];
        turn_b = body_is_box(b);
        arm_b = middle - b.position;
        mass += math::Real{1} / b.mass;
        if (turn_b)
        {
            mass += math::dot(n, math::cross(inverse_inertia_world(b, math::cross(arm_b, n)),
                                             arm_b));
        }
    }
    if (!(mass > math::Real{0}))
    {
        return;
    }
    const math::Vec3 impulse = n * (-correction / mass);
    a.position -= impulse / a.mass;
    if (turn_a)
    {
        turn(a, -inverse_inertia_world(a, math::cross(arm_a, impulse)));
    }
    if (!c.ground)
    {
        RigidBody3D& b = bodies[c.b];
        b.position += impulse / b.mass;
        if (turn_b)
        {
            turn(b, inverse_inertia_world(b, math::cross(arm_b, impulse)));
        }
    }
}

// Phase D: nonlinear Gauss-Seidel position projection.
void project(std::vector<RigidBody3D>& bodies, const Contacts& contacts)
{
    for (int pass = 0; pass < kProjectionPasses; ++pass)
    {
        // A LEVEL manifold is lifted whole at its centroid first: one push
        // through the middle of a flat face, so the corners are not pushed one
        // at a time and a flat landing gains no spin.
        for (const Manifold& m : contacts.manifolds)
        {
            if (m.count < 2)
            {
                continue;
            }
            const math::Real slop = kSlopFraction * m.scale;
            math::Real deepest = -std::numeric_limits<math::Real>::infinity();
            math::Real shallowest = std::numeric_limits<math::Real>::infinity();
            math::Vec3 sum_a{};
            math::Vec3 sum_b{};
            for (std::size_t k = m.first; k < m.first + m.count; ++k)
            {
                const ContactPoint& c = contacts.points[k];
                const math::Vec3 pa = anchor_a(bodies, c);
                const math::Vec3 pb = anchor_b(bodies, c);
                const math::Real depth = -math::dot(pb - pa, c.normal);
                deepest = std::fmax(deepest, depth);
                shallowest = std::fmin(shallowest, depth);
                sum_a += pa;
                sum_b += pb;
            }
            if (!(deepest > slop) || deepest - shallowest > kLevelTolerance * m.level_size)
            {
                continue;
            }
            const math::Real share = math::Real{1} / static_cast<math::Real>(m.count);
            const math::Vec3 pa = sum_a * share;
            const math::Vec3 pb = sum_b * share;
            const ContactPoint& head = contacts.points[m.first];
            const math::Real separation = math::dot(pb - pa, head.normal);
            const math::Real correction = std::clamp(
                separation + slop, -kMaxCorrectionFraction * m.scale, math::Real{0});
            if (correction < math::Real{0})
            {
                push(bodies, head, pa, pb, correction);
            }
        }
        for (const Manifold& m : contacts.manifolds)
        {
            const math::Real slop = kSlopFraction * m.scale;
            for (std::size_t k = m.first; k < m.first + m.count; ++k)
            {
                const ContactPoint& c = contacts.points[k];
                const math::Vec3 pa = anchor_a(bodies, c);
                const math::Vec3 pb = anchor_b(bodies, c);
                const math::Real separation = math::dot(pb - pa, c.normal);
                const math::Real correction =
                    std::clamp(kProjectionFactor * (separation + slop),
                               -kMaxCorrectionFraction * m.scale, math::Real{0});
                if (correction < math::Real{0})
                {
                    push(bodies, c, pa, pb, correction);
                }
            }
        }
    }
}

// The rest-phase impulses, kept for the next step's warm start.
void remember(const Contacts& contacts, std::vector<Rigid3DWorld::WarmContact>& cache)
{
    cache.clear();
    for (const ContactPoint& c : contacts.points)
    {
        cache.push_back(Rigid3DWorld::WarmContact{c.a, c.b, c.ground, c.key, c.local,
                                                  c.normal_impulse, c.tangent_impulse});
    }
}

// Step (2): one body's torque-free rotation over dt (issue #88).
//
// The free Hamiltonian H = sum L_a^2 / (2 I_a), with the axes ordered by
// moment (j the smallest, i the median, k the largest, ties by index), is
// regrouped as
//
//     H = |L|^2 / (2 I_i) + c_j L_j^2 / 2 + c_k L_k^2 / 2,
//     c_j = 1/I_j - 1/I_i,  c_k = 1/I_k - 1/I_i,
//
// and stepped by the symmetric composition A_j(dt/2), A_k(dt), C(dt),
// A_j(dt/2) of the three terms' EXACT flows:
//
//   - A_a(s), the flow of c_a L_a^2 / 2: L_a is constant along it, so the
//     body turns by phi = s c_a L_a about its axis e_a (q <- q r, with
//     r = (cos phi/2, sin phi/2 e_a)) while L_b turns by -phi about e_a;
//   - C(s), the flow of |L|^2 / (2 I_i): L_b is unchanged, and the body turns
//     by s |L| / I_i about L_b.
//
// Every piece depends on the body-frame momentum alone, so every piece keeps
// the WORLD angular momentum R L_b, and |L| is unchanged; the composition is
// time-reversible and second order. The median base makes it exact for every
// sphere and every axisymmetric body, since then c_j or c_k is 0 and the two
// remaining flows commute. Energy is not conserved for three distinct
// moments; its error stays bounded in the measured cases, which is not a
// general guarantee (rigid3d_world.hpp).
//
// Three rules keep the stored angular velocity bit for bit where it was exact
// before (spheres, eigenspace states, principal-axis spins, a body at rest):
//
//   R1. a turn whose angle is exactly 0 is skipped, which is bit-identical to
//       applying it;
//   R2. A_a is applied as a plane rotation of the two components across e_a,
//       so L_a is never touched;
//   R3. only the components of L_b the flow changed are converted back to
//       angular velocity.
//
// The orientation is not bit-identical to the forward-Euler update it
// replaces: it is turned by exact rotations rather than by q + q (0, w) dt / 2.
void free_rotation(RigidBody3D& body, math::Real dt)
{
    const std::array<math::Real, 3> moment{body.inertia.x, body.inertia.y, body.inertia.z};
    const math::Vec3& w = body.angular_velocity;
    const std::array<math::Real, 3> start{moment[0] * w.x, moment[1] * w.y, moment[2] * w.z};
    std::array<math::Real, 3> momentum = start;

    std::array<std::size_t, 3> order{0, 1, 2};
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return moment[a] < moment[b]; });
    const std::size_t j = order[0];
    const std::size_t i = order[1];
    const std::size_t k = order[2];
    const math::Real base = math::Real{1} / moment[i];
    const math::Real c_j = math::Real{1} / moment[j] - base;
    const math::Real c_k = math::Real{1} / moment[k] - base;

    math::Quat q = body.orientation;
    const auto axis_turn = [&](std::size_t a, math::Real c, math::Real s)
    {
        const math::Real phi = s * c * momentum[a];
        if (phi == math::Real{0})
        {
            return; // R1
        }
        const math::Real cos_half = std::cos(phi * math::Real{0.5});
        const math::Real sin_half = std::sin(phi * math::Real{0.5});
        std::array<math::Real, 3> axis{};
        axis[a] = sin_half;
        q = q * math::Quat{cos_half, math::Vec3{axis[0], axis[1], axis[2]}};
        // R2: L_b turns by -phi about e_a, in the plane of the other two axes,
        // by the matrix of the same quaternion (cos phi = c^2 - s^2 and
        // sin phi = 2 s c), so the momentum and the body turn together and the
        // turn costs two trigonometric calls rather than four.
        const std::size_t b = (a + 1) % 3;
        const std::size_t c2 = (a + 2) % 3;
        const math::Real cos_phi = cos_half * cos_half - sin_half * sin_half;
        const math::Real sin_phi = math::Real{2} * sin_half * cos_half;
        const math::Real lb = momentum[b];
        const math::Real lc = momentum[c2];
        momentum[b] = lb * cos_phi + lc * sin_phi;
        momentum[c2] = lc * cos_phi - lb * sin_phi;
    };
    const auto casimir_turn = [&](math::Real s)
    {
        const math::Vec3 axis{momentum[0], momentum[1], momentum[2]};
        const math::Real size = math::length(axis);
        const math::Real angle = s * size * base;
        if (angle == math::Real{0})
        {
            return; // R1
        }
        // from_axis_angle(axis, angle), with the length already in hand.
        const math::Real half_angle = angle * math::Real{0.5};
        q = q * math::Quat{std::cos(half_angle), axis * (std::sin(half_angle) / size)};
    };

    const math::Real half = dt * math::Real{0.5};
    axis_turn(j, c_j, half);
    axis_turn(k, c_k, dt);
    casimir_turn(dt);
    axis_turn(j, c_j, half);
    body.orientation = q;

    // R3: a component the flow left exactly as it was keeps its bits.
    math::Vec3& omega = body.angular_velocity;
    if (momentum[0] != start[0])
    {
        omega.x = momentum[0] / moment[0];
    }
    if (momentum[1] != start[1])
    {
        omega.y = momentum[1] / moment[1];
    }
    if (momentum[2] != start[2])
    {
        omega.z = momentum[2] / moment[2];
    }
}
} // namespace

bool RigidBody3D::is_valid() const
{
    if (!(mass > math::Real{0} && math::is_finite(mass) &&
          inertia.x > math::Real{0} && inertia.y > math::Real{0} &&
          inertia.z > math::Real{0} && math::is_finite(inertia) &&
          radius >= math::Real{0} && math::is_finite(radius) &&
          half_extents.x >= math::Real{0} && half_extents.y >= math::Real{0} &&
          half_extents.z >= math::Real{0} && math::is_finite(half_extents) &&
          math::is_unit(orientation) && math::is_squarable(position) &&
          math::is_squarable(velocity) && math::is_squarable(angular_velocity)))
    {
        return false;
    }

    // The body's own linear momentum m v and spin angular momentum (each moment
    // times its w) must be SQUARABLE too (issue #51). A squarable velocity is
    // not enough once the mass or a moment is huge: mass and moments 1e300 at
    // speed and spin 1e5 validated clean while E, |L| and |p| all printed inf.
    // Squarable rather than finite because the world reports the MAGNITUDES
    // |p| and |L|, which square them: mass 1e300 at speed 1e-140 has a finite
    // energy and a finite momentum, and still a |p| of inf.
    //
    // That also keeps the kinetic energy finite, with no check of its own:
    // (1/2) m v.v = (1/2) (m v).v <= (1/2) |m v| |v|, the product of two
    // squarable magnitudes, so below half the largest double, and the same for
    // (1/2) sum I w^2 <= (1/2) |I w| |w|. Every 3D mass and moment is finite,
    // so there is no infinity to guard. (The orbital angular momentum and the
    // gravity and force potentials depend on the origin and on the settings,
    // so they are not the body's own and are not bounded here.)
    const math::Vec3 spin{inertia.x * angular_velocity.x, inertia.y * angular_velocity.y,
                          inertia.z * angular_velocity.z};
    return math::is_squarable(velocity * mass) && math::is_squarable(spin);
}

bool Rigid3DSettings::is_valid() const
{
    if (!math::is_squarable(torque) || !math::is_squarable(force) ||
        !math::is_squarable(gravity))
    {
        return false;
    }
    if (!(restitution >= math::Real{0}) || !(restitution <= math::Real{1}) ||
        !math::is_finite(restitution))
    {
        return false;
    }
    if (!(friction >= math::Real{0}) || !math::is_finite(friction))
    {
        return false;
    }
    for (const collide::Plane3& plane : ground)
    {
        if (!plane.is_valid())
        {
            return false;
        }
    }
    return true;
}

math::Vec3 spin_angular_momentum(const RigidBody3D& body)
{
    // I omega in the body frame, where the inertia is diagonal, then rotated
    // into the world frame. The body-frame vector tumbles; this one does not.
    const math::Vec3 body_frame{body.inertia.x * body.angular_velocity.x,
                                body.inertia.y * body.angular_velocity.y,
                                body.inertia.z * body.angular_velocity.z};
    return math::rotate(body.orientation, body_frame);
}

math::Real rotational_energy(const RigidBody3D& body)
{
    const math::Vec3& w = body.angular_velocity;
    // (1/2) omega . (I omega), in the body frame. A scalar, so it does not
    // matter which frame it is computed in, and the body frame is where the
    // inertia is diagonal.
    return math::Real{0.5} * (body.inertia.x * w.x * w.x + body.inertia.y * w.y * w.y +
                              body.inertia.z * w.z * w.z);
}

Rigid3DWorld::Rigid3DWorld(sim_core::SimulationSettings simulation_settings,
                           std::vector<RigidBody3D> bodies,
                           Rigid3DSettings settings)
    : simulation_settings_{simulation_settings},
      settings_{settings},
      bodies_{std::move(bodies)},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus Rigid3DWorld::validate() const
{
    if (!simulation_settings_.is_valid() || // dt > 0, finite
        !settings_.is_valid())              // torque, force, gravity squarable
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const RigidBody3D& body : bodies_)
    {
        if (!body.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

sim_core::StepResult Rigid3DWorld::step()
{
    const sim_core::StepStatus status = validate();
    if (status != sim_core::StepStatus::Ok)
    {
        return sim_core::StepResult{status}; // leave state unchanged
    }

    // The state is the bodies AND the warm-start impulses; a failed step
    // restores both.
    previous_ = bodies_;
    previous_warm_ = warm_;

    const math::Real dt = simulation_settings_.dt;

    // (A) detect at the start-of-step poses, speculatively (step contract).
    Contacts contacts;
    detect(bodies_, settings_, contacts);

    // (B) the impacts, on the start-of-step velocities. Restitution 0 has none:
    //     every contact is then resolved by the rest phase alone.
    if (!contacts.points.empty() && settings_.restitution > math::Real{0})
    {
        impact_phase(bodies_, settings_, dt, contacts);
    }

    for (RigidBody3D& body : bodies_)
    {
        // (0) gravity, before the position update, which keeps this
        //     semi-implicit Euler. An acceleration, so it does not scale with
        //     mass, and every RigidBody3D is movable (the ground is geometry,
        //     not a body). Zero by default, so a pre-M22 world is unchanged.
        body.velocity += settings_.gravity * dt;

        //     Then the applied FORCE (M28), as the acceleration F/m it
        //     produces. Unlike gravity this scales with 1/mass, which is the
        //     whole difference between a force and an acceleration and is why it
        //     needs the body's mass (computed from geometry, M25 to M27). It
        //     acts through the centre of mass, so it makes no torque; the
        //     rotational half of a wrench is the torque in step (1). Zero by
        //     default, and F = 0 leaves this line adding a zero vector, so a
        //     pre-M28 world is unchanged bit for bit.
        body.velocity += (settings_.force / body.mass) * dt;

        // (1) the applied torque (M21) as ONE kick: the world angular momentum
        //     gains T dt exactly, at the start-of-step pose (issue #88). The
        //     torque is a WORLD-frame setting and the body stores its angular
        //     velocity in the body frame, so it is rotated in by the CURRENT
        //     orientation, T_body = R^-1 T_world, `rotate` by the conjugate.
        //
        //     There is no gyroscopic term here any more: the rotation itself,
        //     the term that makes a body with three distinct moments tumble,
        //     is the exact free flow of step (2), after the contact phase.
        //     A zero torque adds a zero vector, so a torque-free world is
        //     unchanged by this line.
        const math::Vec3& I = body.inertia;
        const math::Vec3 torque_body =
            math::rotate(math::conjugate(body.orientation), settings_.torque);
        const math::Vec3 rate{torque_body.x / I.x, torque_body.y / I.y, torque_body.z / I.z};
        body.angular_velocity += rate * dt;
    }

    // (C) every contact, after the forcing, warm-started from the last step.
    if (!contacts.points.empty())
    {
        rest_phase(bodies_, settings_, dt, contacts, warm_);
    }

    for (RigidBody3D& body : bodies_)
    {
        // (2) the torque-free rotation over dt, from the angular velocity the
        //     kick and the contact phases left: the exact sub-flows of
        //     `free_rotation`, which turn the orientation and carry the
        //     body-frame momentum with it, keeping the world angular momentum.
        free_rotation(body, dt);

        // (3) renormalize. Each turn is a unit quaternion, so the product
        //     stays on the unit sphere up to rounding, which this removes.
        //
        //     This is NOT the positional correction's kind of projection and
        //     costs nothing physically. A nonzero quaternion and any scaling of
        //     it represent the same rotation, since q v q^-1 is invariant under
        //     q -> k q, so renormalizing moves no body and perturbs no
        //     conserved quantity. It is required because `rotate` uses the
        //     optimized unit-only form, which is NOT scale invariant: fed a
        //     quaternion of norm k it returns p + k^2 (R p - p), scaling the
        //     displacement rather than performing the rotation.
        body.orientation = math::normalize(body.orientation);

        // (4) translation, with the solved velocity.
        body.position += body.velocity * dt;
    }

    // (D) project the remaining penetration out, positions only.
    if (!contacts.points.empty())
    {
        project(bodies_, contacts);
    }
    remember(contacts, warm_);

    for (const RigidBody3D& body : bodies_)
    {
        if (!body.is_valid())
        {
            bodies_.swap(previous_);
            warm_.swap(previous_warm_);
            return sim_core::StepResult{sim_core::StepStatus::InvalidState};
        }
    }

    step_->advance();
    return sim_core::StepResult{sim_core::StepStatus::Ok};
}

math::Vec3 total_linear_momentum3d(const std::vector<RigidBody3D>& bodies)
{
    math::Vec3 momentum{};
    for (const RigidBody3D& body : bodies)
    {
        momentum += body.velocity * body.mass;
    }
    return momentum;
}

math::Vec3 total_angular_momentum3d(const std::vector<RigidBody3D>& bodies)
{
    math::Vec3 angular{};
    for (const RigidBody3D& body : bodies)
    {
        // Spin about the centre of mass, plus the orbital term about the world
        // origin. Dropping either makes a whole class of motion look conserved
        // when it is not, which is why the 2D version sums both as well.
        //
        // No infinity guard, unlike the 2D version: there an infinite mass is a
        // legitimate state meaning an immovable body, while a RigidBody3D with
        // an infinite mass or moment is not valid, so a world that steps cannot
        // hold one. `position` is the centre of mass; there is no local offset
        // in three dimensions.
        angular += spin_angular_momentum(body);
        angular += body.mass * math::cross(body.position, body.velocity);
    }
    return angular;
}

math::Real total_kinetic_energy3d(const std::vector<RigidBody3D>& bodies)
{
    math::Real kinetic = 0.0;
    for (const RigidBody3D& body : bodies)
    {
        kinetic += math::Real{0.5} * body.mass * math::dot(body.velocity, body.velocity);
        kinetic += rotational_energy(body);
    }
    return kinetic;
}

math::Real total_potential_energy3d(const std::vector<RigidBody3D>& bodies,
                                    const math::Vec3& gravity)
{
    // Sum of -m (g . r), r the centre of mass. Zero when gravity is zero.
    math::Real potential = 0.0;
    for (const RigidBody3D& body : bodies)
    {
        potential -= body.mass * math::dot(gravity, body.position);
    }
    return potential;
}

math::Real total_force_potential3d(const std::vector<RigidBody3D>& bodies,
                                   const math::Vec3& force)
{
    // Sum of -(F . r), r the centre of mass. The applied force is uniform, not
    // mass-scaled like gravity's m g, so there is no mass factor here. Zero when
    // the force is zero.
    math::Real potential = 0.0;
    for (const RigidBody3D& body : bodies)
    {
        potential -= math::dot(force, body.position);
    }
    return potential;
}
} // namespace malloy::rigid
