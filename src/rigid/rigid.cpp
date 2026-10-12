#include <malloy/rigid/rigid.hpp>

#include <algorithm>
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
// Rotate a vector by an angle. The one place trigonometry enters this module,
// so a sign error here has exactly one home.
math::Vec2 rotate(const math::Vec2& v, math::Real angle)
{
    const math::Real c = std::cos(angle);
    const math::Real s = std::sin(angle);
    return math::Vec2{v.x * c - v.y * s, v.x * s + v.y * c};
}

// The 2D cross product, a scalar: the z component of a x b.
math::Real cross(const math::Vec2& a, const math::Vec2& b)
{
    return a.x * b.y - a.y * b.x;
}

// ============================================================================
// M40: the per-point iterative contact solver (issue #73), the 2D sibling of
// the one in rigid3d.cpp, with the same constants for the same reasons. Each is
// defined once, here, with its derivation; the step contract in rigid_world.hpp
// states them.
// ============================================================================

// Speculative margin, as a fraction of the smaller body's size: a pair still
// apart by up to this is kept as a contact with a negative penetration (the
// gap), and the rest phase lets it close the gap but not pass it. It must
// exceed any motion of a body at rest (the rest tolerance leaves at most about
// 1e-8 g dt of velocity, under 1e-12 of a size per step) and the slop, so a
// resting point inside the slop band is always inside the margin; and it must
// be invisible, since a bounce can start from a gap this wide. Without it a
// resting corner rounds clear of the floor, drops out of the manifold and falls
// for a step (a box nudged to rest flickered off its own corner after about
// 8500 steps in the research prototype). Measured the same from 1e-6 to 1e-2.
constexpr math::Real kMarginFraction = 1e-4;

// Allowed penetration, as a fraction of the smaller body's size, below which
// the position projection does not push: at least 1e10 times position rounding,
// so the projection never chases rounding, and no larger than the precision a
// template states for a resting height. Measured the same from 1e-9 to 1e-2.
constexpr math::Real kSlopFraction = 1e-6;

// Fraction of the remaining penetration one projection push removes (Box2D's
// position factor, below 1 so points sharing a body do not overshoot), the
// passes per step, and the largest push as a fraction of the size. Measured the
// same from 0.1 to 1 and from 1 to 10 passes.
constexpr math::Real kProjectionFactor = 0.2;
constexpr int kProjectionPasses = 3;
constexpr math::Real kMaxCorrectionFraction = 0.2;

// Sweep termination: stop when no point's constraint is violated by more than
// this times the velocity scale (the largest speed in contact, or one step of
// gravity, whichever is larger); impacts to the tighter tolerance; either way
// at most 100 sweeps. The 3D world measured these; see rigid3d.cpp.
constexpr math::Real kRestTolerance = 1e-8;
constexpr math::Real kImpactTolerance = 1e-12;
constexpr int kMaxSweeps = 100;

// Restitution threshold, as a multiple of dt |g| (see rigid3d.cpp, where it is
// derived and measured): a point of a box pivoting on another contact can fall
// at twice the box's acceleration, so one that fell clear for a step approaches
// at up to 1.5 x 2 = 3 |g| dt at mid-step; four times |g| dt clears that, so a
// resting or rocking body never bounces on its own weight.
constexpr math::Real kThresholdFactor = 4.0;

// Warm-start match radius for two boxes, as a fraction of body a's size.
constexpr math::Real kWarmMatchFraction = 0.05;

// A manifold is LEVEL when its depths agree to within this fraction of the
// boxes' size, the collide layer's touch-tolerance convention: its centroid
// weights are then equal rather than depth-weighted (issue #41).
constexpr math::Real kLevelTolerance = 1e-9;

// A body collides as an oriented box when both half-extents are positive;
// otherwise it uses the disc radius. Zero half-extents (the default) is not a
// box, so a pre-M36 body keeps its disc behaviour, mirroring `body_is_box` in
// the 3D world.
bool body_is_box2(const RigidBody2D& body)
{
    return body.half_extents.x > math::Real{0} && body.half_extents.y > math::Real{0};
}

// The size every contact tolerance scales with: a box's larger half-width, or a
// disc's radius (0 for a body that does not collide).
math::Real collider_size(const RigidBody2D& body)
{
    if (body_is_box2(body))
    {
        return std::fmax(body.half_extents.x, body.half_extents.y);
    }
    return body.radius;
}

// One point of a contact manifold, for one step. `b` is the other body, or the
// plane index when `ground` (an immovable participant).
struct ContactPoint
{
    std::size_t a{0};
    std::size_t b{0};
    bool ground{false};
    int key{-1};                // warm-start feature (see WarmContact)
    math::Vec2 normal{};        // unit, from a toward b
    math::Vec2 arm_a{};         // a's centre of mass to the point (world)
    math::Vec2 arm_b{};
    math::Vec2 local{};         // the point in a's body frame
    math::Vec2 anchor_a{};      // projection anchors: a's surface point in its
    math::Vec2 anchor_b{};      // frame (a disc's centre, its origin), b's
    math::Real radius_a{0.0};   // likewise or, for the ground, the world point;
    math::Real radius_b{0.0};   // a disc adds its radius along the normal
    math::Real penetration{0.0}; // signed: negative is a speculative gap
    math::Real weight{0.0};
    math::Real normal_mass{0.0};
    math::Real tangent_mass{0.0};
    math::Real target{0.0};
    math::Real normal_impulse{0.0};
    math::Real tangent_impulse{0.0};
    bool active{true};
};

struct Manifold
{
    std::size_t first{0};
    std::size_t count{0};
    math::Real scale{0.0};      // the smaller body's size: margin, slop, push cap
    math::Real level_size{0.0}; // the larger: the level test
};

// A body's velocities while a phase runs. Inverse mass and inertia are exactly 0
// for an immovable quantity, so nothing moves it.
struct SolverBody
{
    math::Vec2 velocity{};
    math::Real spin{0.0};
    math::Real inverse_mass{0.0};
    math::Real inverse_inertia{0.0};
    bool touched{false};
};

struct Contacts
{
    std::vector<ContactPoint> points;
    std::vector<Manifold> manifolds;
};

void add_point(const std::vector<RigidBody2D>& bodies, Contacts& contacts, std::size_t a,
               std::size_t b, bool ground, int key, const collide::Contact& hit)
{
    ContactPoint c;
    c.a = a;
    c.b = b;
    c.ground = ground;
    c.key = key;
    c.normal = hit.normal;
    c.penetration = hit.penetration;
    const RigidBody2D& body_a = bodies[a];
    // Arms from each centre of mass to the contact point, not from the origins
    // (M14): a disc is centred on its ORIGIN, so a contact normal through the
    // disc centre misses an offset centre of mass, and that offset is torque.
    c.arm_a = hit.point - center_of_mass(body_a);
    c.local = to_local(body_a, hit.point);
    // Each body's own surface point: the contact point is midway, so a's lies
    // half the penetration along the normal, b's half against it.
    const math::Vec2 surface_a = hit.point + hit.normal * (hit.penetration / math::Real{2});
    const math::Vec2 surface_b = hit.point - hit.normal * (hit.penetration / math::Real{2});
    if (body_is_box2(body_a))
    {
        c.anchor_a = to_local(body_a, surface_a);
    }
    else
    {
        c.radius_a = body_a.radius; // the origin, plus the radius along the normal
    }
    if (ground)
    {
        c.anchor_b = surface_b;
    }
    else
    {
        const RigidBody2D& body_b = bodies[b];
        c.arm_b = hit.point - center_of_mass(body_b);
        if (body_is_box2(body_b))
        {
            c.anchor_b = to_local(body_b, surface_b);
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

void detect(const std::vector<RigidBody2D>& bodies, const RigidSettings& settings,
            Contacts& contacts)
{
    const std::size_t count = bodies.size();
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            const RigidBody2D& a = bodies[i];
            const RigidBody2D& b = bodies[j];
            if (a.is_static() && b.is_static())
            {
                continue; // nothing can move either
            }
            const bool a_box = body_is_box2(a);
            const bool b_box = body_is_box2(b);
            if (a_box != b_box)
            {
                continue; // a box against a disc is deferred (issue #48)
            }
            const math::Real size_a = collider_size(a);
            const math::Real size_b = collider_size(b);
            if (!(size_a > math::Real{0}) || !(size_b > math::Real{0}))
            {
                continue; // a zero radius takes part in no contact
            }
            const math::Real scale = std::fmin(size_a, size_b);
            const math::Real margin = kMarginFraction * scale;
            const std::size_t first = contacts.points.size();
            if (a_box)
            {
                const std::vector<collide::Contact> hits =
                    collide::contacts(collide::Obb2{a.position, a.half_extents, a.angle},
                                      collide::Obb2{b.position, b.half_extents, b.angle}, margin);
                for (const collide::Contact& hit : hits)
                {
                    add_point(bodies, contacts, i, j, false, -1, hit);
                }
            }
            else
            {
                const std::optional<collide::Contact> hit =
                    collide::contact(collide::Circle{a.position, a.radius},
                                     collide::Circle{b.position, b.radius}, margin);
                if (hit)
                {
                    add_point(bodies, contacts, i, j, false, 0, *hit);
                }
            }
            close_manifold(contacts, first, scale, std::fmax(size_a, size_b));
        }
    }
    for (std::size_t i = 0; i < count; ++i)
    {
        const RigidBody2D& body = bodies[i];
        const math::Real size = collider_size(body);
        if (body.is_static() || !(size > math::Real{0}))
        {
            continue; // the ground cannot push what cannot move
        }
        const math::Real margin = kMarginFraction * size;
        for (std::size_t k = 0; k < settings.ground.size(); ++k)
        {
            const collide::Halfplane& plane = settings.ground[k];
            const std::size_t first = contacts.points.size();
            if (body_is_box2(body))
            {
                const std::vector<collide::Contact> hits = collide::contacts(
                    collide::Obb2{body.position, body.half_extents, body.angle}, plane, margin);
                for (const collide::Contact& hit : hits)
                {
                    // The corner this point came from, by the signs of its body
                    // coordinates: the warm-start key a resting corner keeps.
                    const math::Vec2 corner =
                        hit.point + hit.normal * (hit.penetration / math::Real{2});
                    const math::Vec2 local = to_local(body, corner);
                    const int key =
                        (local.x > math::Real{0} ? 1 : 0) | (local.y > math::Real{0} ? 2 : 0);
                    add_point(bodies, contacts, i, k, true, key, hit);
                }
            }
            else
            {
                const std::optional<collide::Contact> hit = collide::contact(
                    collide::Circle{body.position, body.radius}, plane, margin);
                if (hit)
                {
                    add_point(bodies, contacts, i, k, true, 0, *hit);
                }
            }
            close_manifold(contacts, first, size, size);
        }
    }
}

std::vector<SolverBody> solver_bodies(const std::vector<RigidBody2D>& bodies,
                                      const Contacts& contacts)
{
    std::vector<SolverBody> solver(bodies.size());
    const auto prepare = [&](std::size_t i)
    {
        const RigidBody2D& body = bodies[i];
        SolverBody& s = solver[i];
        s.velocity = body.velocity;
        s.spin = body.angular_velocity;
        s.inverse_mass = body.inverse_mass();
        s.inverse_inertia = body.inverse_inertia();
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
math::Vec2 point_velocity(const std::vector<SolverBody>& solver, std::size_t i, bool ground,
                          const math::Vec2& arm)
{
    if (ground)
    {
        return math::Vec2{};
    }
    const SolverBody& s = solver[i];
    return s.velocity + math::Vec2{-s.spin * arm.y, s.spin * arm.x};
}

math::Vec2 relative_velocity(const std::vector<SolverBody>& solver, const ContactPoint& c,
                             const math::Vec2& arm_a, const math::Vec2& arm_b)
{
    return point_velocity(solver, c.b, c.ground, arm_b) -
           point_velocity(solver, c.a, false, arm_a);
}

// Apply `impulse` at `arm`: v += J/m, w += (r x J)/I. An immovable quantity has
// an inverse of exactly 0 and is skipped rather than multiplied, so an immovable
// body is never written back; a zero impulse touches nothing.
void apply(std::vector<SolverBody>& solver, std::size_t i, bool ground,
           const math::Vec2& impulse, const math::Vec2& arm)
{
    if (ground || (impulse.x == math::Real{0} && impulse.y == math::Real{0}))
    {
        return;
    }
    SolverBody& s = solver[i];
    if (s.inverse_mass > math::Real{0})
    {
        s.velocity += impulse * s.inverse_mass;
        s.touched = true;
    }
    if (s.inverse_inertia > math::Real{0})
    {
        const math::Real turn = cross(arm, impulse) * s.inverse_inertia;
        if (turn != math::Real{0})
        {
            s.spin += turn;
            s.touched = true;
        }
    }
}

void apply_pair(std::vector<SolverBody>& solver, const ContactPoint& c,
                const math::Vec2& impulse, const math::Vec2& arm_a, const math::Vec2& arm_b)
{
    apply(solver, c.a, false, -impulse, arm_a);
    apply(solver, c.b, c.ground, impulse, arm_b);
}

// k = 1/m_a + 1/m_b + (r_a x d)^2 / I_a + (r_b x d)^2 / I_b along d (M14).
math::Real effective_mass(const std::vector<SolverBody>& solver, const ContactPoint& c,
                          const math::Vec2& arm_a, const math::Vec2& arm_b,
                          const math::Vec2& d)
{
    const SolverBody& sa = solver[c.a];
    const math::Real ca = cross(arm_a, d);
    math::Real k = sa.inverse_mass + ca * ca * sa.inverse_inertia;
    if (!c.ground)
    {
        const SolverBody& sb = solver[c.b];
        const math::Real cb = cross(arm_b, d);
        k += sb.inverse_mass + cb * cb * sb.inverse_inertia;
    }
    return k;
}

// The tangent: the normal turned a quarter counter-clockwise.
math::Vec2 tangent_of(const math::Vec2& n)
{
    return math::Vec2{-n.y, n.x};
}

void prepare_points(const std::vector<SolverBody>& solver, Contacts& contacts)
{
    for (ContactPoint& c : contacts.points)
    {
        c.normal_mass = effective_mass(solver, c, c.arm_a, c.arm_b, c.normal);
        c.tangent_mass = effective_mass(solver, c, c.arm_a, c.arm_b, tangent_of(c.normal));
        c.normal_impulse = math::Real{0};
        c.tangent_impulse = math::Real{0};
    }
}

// One impulse at the manifold's weighted centroid before the per-point sweeps:
// the whole answer for a level landing, applied in one piece so it stays exactly
// symmetric, and otherwise a starting guess the sweeps refine. Its share is
// booked to each point's accumulated impulse by weight.
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
        math::Vec2 arm_a{};
        math::Vec2 arm_b{};
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
        const ContactPoint& head = contacts.points[m.first];
        const math::Vec2 n = head.normal;
        const math::Real mass = effective_mass(solver, head, arm_a, arm_b, n);
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
        math::Real tangent_held = 0.0;
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            const ContactPoint& c = contacts.points[k];
            if (c.active && c.weight > math::Real{0})
            {
                lambda = std::fmax(lambda, -c.normal_impulse * weight_sum / c.weight);
            }
            if (c.active)
            {
                tangent_held += c.tangent_impulse;
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
        const math::Vec2 t = tangent_of(n);
        const math::Real tangent_mass = effective_mass(solver, head, arm_a, arm_b, t);
        if (!(tangent_mass > math::Real{0}))
        {
            continue;
        }
        // Likewise accumulated: the manifold's total friction (what the warm
        // start applied plus this) stays within friction times its total push.
        const math::Real limit = friction * total;
        const math::Real slide = math::dot(relative_velocity(solver, head, arm_a, arm_b), t);
        const math::Real change =
            std::clamp(tangent_held - slide / tangent_mass, -limit, limit) - tangent_held;
        apply_pair(solver, head, t * change, arm_a, arm_b);
        for (std::size_t k = m.first; k < m.first + m.count; ++k)
        {
            ContactPoint& c = contacts.points[k];
            if (c.active)
            {
                c.tangent_impulse += c.weight / weight_sum * change;
            }
        }
    }
}

// One Gauss-Seidel sweep: per manifold, every normal row, then every friction
// row. Accumulated impulses are clamped, never the change.
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
            if (!c.active || !(c.tangent_mass > math::Real{0}))
            {
                continue;
            }
            const math::Vec2 t = tangent_of(c.normal);
            const math::Real slide = math::dot(relative_velocity(solver, c, c.arm_a, c.arm_b), t);
            const math::Real limit = friction * c.normal_impulse;
            const math::Real next =
                std::clamp(c.tangent_impulse - slide / c.tangent_mass, -limit, limit);
            const math::Real change = next - c.tangent_impulse;
            c.tangent_impulse = next;
            apply_pair(solver, c, t * change, c.arm_a, c.arm_b);
        }
    }
}

// The largest constraint violation (see rigid3d.cpp).
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
        const math::Vec2 relative = relative_velocity(solver, c, c.arm_a, c.arm_b);
        const math::Real approach = math::dot(relative, c.normal);
        worst = std::fmax(worst, c.normal_impulse > math::Real{0}
                                     ? std::abs(approach - c.target)
                                     : std::fmax(math::Real{0}, c.target - approach));
        if (friction > math::Real{0} && c.tangent_mass > math::Real{0})
        {
            const math::Real slide = math::dot(relative, tangent_of(c.normal));
            const math::Real limit = friction * c.normal_impulse;
            const math::Real held = std::abs(c.tangent_impulse);
            if (held < limit * (math::Real{1} - math::Real{1e-12}))
            {
                worst = std::fmax(worst, std::abs(slide));
            }
            else if (held > math::Real{0})
            {
                worst = std::fmax(worst, std::fmax(math::Real{0},
                                                   c.tangent_impulse > math::Real{0} ? slide
                                                                                     : -slide));
            }
        }
    }
    return worst;
}

math::Real velocity_scale(const std::vector<RigidBody2D>& bodies,
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
                                         std::abs(s.spin) * collider_size(bodies[i]));
        }
    }
    return scale;
}

// Sweeps until the residual is within `tolerance`, checked before the first
// sweep too, so an already-solved symmetric manifold is left exactly as it is;
// at most kMaxSweeps.
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

void write_back(std::vector<RigidBody2D>& bodies, const std::vector<SolverBody>& solver)
{
    for (std::size_t i = 0; i < bodies.size(); ++i)
    {
        const SolverBody& s = solver[i];
        if (s.touched)
        {
            bodies[i].velocity = s.velocity;
            bodies[i].angular_velocity = s.spin;
        }
    }
}

// Phase B: the impacts, on the start-of-step velocities.
void impact_phase(std::vector<RigidBody2D>& bodies, const RigidSettings& settings,
                  math::Real dt, Contacts& contacts)
{
    const math::Real restitution = settings.restitution;
    std::vector<SolverBody> solver = solver_bodies(bodies, contacts);
    prepare_points(solver, contacts);
    const math::Real acceleration = math::length(settings.gravity);
    const math::Real threshold = kThresholdFactor * dt * acceleration;
    bool any = false;
    for (ContactPoint& c : contacts.points)
    {
        // delta: what this step's gravity will add to the point's normal
        // velocity; an infinite mass does not fall (step 1).
        const auto fall = [&](std::size_t i)
        {
            return bodies[i].has_infinite_mass() ? math::Vec2{} : settings.gravity * dt;
        };
        math::Vec2 forcing = -fall(c.a);
        if (!c.ground)
        {
            forcing += fall(c.b);
        }
        const math::Real delta = math::dot(forcing, c.normal);
        const math::Real start =
            math::dot(relative_velocity(solver, c, c.arm_a, c.arm_b), c.normal);
        const math::Real middle = start + delta / math::Real{2};
        const math::Real gap = -c.penetration;
        c.active = middle < -threshold && (!(gap > math::Real{0}) || -middle * dt >= gap);
        c.target = -restitution * middle - delta / math::Real{2};
        any = any || c.active;
    }
    if (!any)
    {
        return;
    }
    prepass(solver, contacts, settings.friction);
    iterate(solver, contacts, settings.friction,
            kImpactTolerance * velocity_scale(bodies, solver, contacts, acceleration * dt));
    write_back(bodies, solver);
}

// Matches this step's points to last step's rest-phase impulses and applies
// them (see rigid3d.cpp).
void warm_start(std::vector<SolverBody>& solver, Contacts& contacts,
                const std::vector<RigidBody2D>& bodies,
                const std::vector<RigidWorld::WarmContact>& cache)
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
            c.tangent_impulse = cache[best].tangent_impulse;
            apply_pair(solver, c,
                       c.normal * c.normal_impulse + tangent_of(c.normal) * c.tangent_impulse,
                       c.arm_a, c.arm_b);
        }
        cursor = end;
    }
}

// Phase C: every point, after the gravity, warm-started.
void rest_phase(std::vector<RigidBody2D>& bodies, const RigidSettings& settings, math::Real dt,
                Contacts& contacts, const std::vector<RigidWorld::WarmContact>& cache)
{
    std::vector<SolverBody> solver = solver_bodies(bodies, contacts);
    prepare_points(solver, contacts);
    for (ContactPoint& c : contacts.points)
    {
        c.active = true;
        c.target = c.penetration < math::Real{0} ? c.penetration / dt : math::Real{0};
    }
    warm_start(solver, contacts, bodies, cache);
    prepass(solver, contacts, settings.friction);
    iterate(solver, contacts, settings.friction,
            kRestTolerance *
                velocity_scale(bodies, solver, contacts, math::length(settings.gravity) * dt));
    write_back(bodies, solver);
}

// Moves a body's centre of mass by `shift` and turns it by `angle` about it,
// for whichever of the two it can do; the origin follows from the centre of
// mass, as in the integration.
void displace(RigidBody2D& body, const math::Vec2& shift, math::Real angle)
{
    const bool moves = !body.has_infinite_mass() && (shift.x != math::Real{0} ||
                                                     shift.y != math::Real{0});
    const bool turns = !body.has_infinite_inertia() && angle != math::Real{0};
    if (!moves && !turns)
    {
        return;
    }
    const math::Vec2 com = center_of_mass(body) + (moves ? shift : math::Vec2{});
    if (turns)
    {
        body.angle += angle;
    }
    body.position = com - rotate(body.local_center_of_mass, body.angle);
}

math::Vec2 anchor_a(const std::vector<RigidBody2D>& bodies, const ContactPoint& c)
{
    return to_world(bodies[c.a], c.anchor_a) + c.normal * c.radius_a;
}

math::Vec2 anchor_b(const std::vector<RigidBody2D>& bodies, const ContactPoint& c)
{
    if (c.ground)
    {
        return c.anchor_b;
    }
    return to_world(bodies[c.b], c.anchor_b) - c.normal * c.radius_b;
}

// The lever of a push along the normal about a body's centre of mass. For a
// box it is the anchor midpoint's arm; for a disc the push acts along the line
// through the disc's centre (its origin), so the lever is that centre's offset
// from the centre of mass, exactly 0 for a centred disc.
math::Vec2 push_arm(const RigidBody2D& body, bool box, const math::Vec2& middle)
{
    return box ? middle - center_of_mass(body) : body.position - center_of_mass(body);
}

// One projection push of `correction` (negative: the depth to remove), moving
// and turning both bodies by inverse mass and inertia at their CURRENT pose.
// Velocities are not touched.
void push(std::vector<RigidBody2D>& bodies, const ContactPoint& c, const math::Vec2& point_a,
          const math::Vec2& point_b, math::Real correction)
{
    const math::Vec2 middle = (point_a + point_b) * math::Real{0.5};
    const math::Vec2 n = c.normal;
    RigidBody2D& a = bodies[c.a];
    const math::Vec2 arm_a = push_arm(a, body_is_box2(a), middle);
    const math::Real cross_a = cross(arm_a, n);
    math::Real mass = a.inverse_mass() + cross_a * cross_a * a.inverse_inertia();
    math::Vec2 arm_b{};
    if (!c.ground)
    {
        const RigidBody2D& b = bodies[c.b];
        arm_b = push_arm(b, body_is_box2(b), middle);
        const math::Real cross_b = cross(arm_b, n);
        mass += b.inverse_mass() + cross_b * cross_b * b.inverse_inertia();
    }
    if (!(mass > math::Real{0}))
    {
        return;
    }
    const math::Vec2 impulse = n * (-correction / mass);
    displace(a, -impulse * a.inverse_mass(), -cross(arm_a, impulse) * a.inverse_inertia());
    if (!c.ground)
    {
        RigidBody2D& b = bodies[c.b];
        displace(b, impulse * b.inverse_mass(), cross(arm_b, impulse) * b.inverse_inertia());
    }
}

// Phase D: nonlinear Gauss-Seidel position projection (see rigid3d.cpp).
void project(std::vector<RigidBody2D>& bodies, const Contacts& contacts)
{
    for (int pass = 0; pass < kProjectionPasses; ++pass)
    {
        for (const Manifold& m : contacts.manifolds)
        {
            if (m.count < 2)
            {
                continue;
            }
            const math::Real slop = kSlopFraction * m.scale;
            math::Real deepest = -std::numeric_limits<math::Real>::infinity();
            math::Real shallowest = std::numeric_limits<math::Real>::infinity();
            math::Vec2 sum_a{};
            math::Vec2 sum_b{};
            for (std::size_t k = m.first; k < m.first + m.count; ++k)
            {
                const ContactPoint& c = contacts.points[k];
                const math::Vec2 pa = anchor_a(bodies, c);
                const math::Vec2 pb = anchor_b(bodies, c);
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
            const math::Vec2 pa = sum_a * share;
            const math::Vec2 pb = sum_b * share;
            const ContactPoint& head = contacts.points[m.first];
            const math::Real correction =
                std::clamp(math::dot(pb - pa, head.normal) + slop,
                           -kMaxCorrectionFraction * m.scale, math::Real{0});
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
                const math::Vec2 pa = anchor_a(bodies, c);
                const math::Vec2 pb = anchor_b(bodies, c);
                const math::Real correction =
                    std::clamp(kProjectionFactor * (math::dot(pb - pa, c.normal) + slop),
                               -kMaxCorrectionFraction * m.scale, math::Real{0});
                if (correction < math::Real{0})
                {
                    push(bodies, c, pa, pb, correction);
                }
            }
        }
    }
}

void remember(const Contacts& contacts, std::vector<RigidWorld::WarmContact>& cache)
{
    cache.clear();
    for (const ContactPoint& c : contacts.points)
    {
        cache.push_back(RigidWorld::WarmContact{c.a, c.b, c.ground, c.key, c.local,
                                                c.normal_impulse, c.tangent_impulse});
    }
}
} // namespace

bool RigidBody2D::is_valid() const
{
    // mass and inertia may be infinite, meaning immovable. NaN still fails,
    // because NaN > 0 is false.
    if (!(mass > math::Real{0} && inertia > math::Real{0} &&
          radius >= math::Real{0} && math::is_finite(radius) &&
          half_extents.x >= math::Real{0} && half_extents.y >= math::Real{0} &&
          math::is_finite(half_extents) &&
          math::is_finite(angle) && math::is_finite(angular_velocity) &&
          math::is_squarable(position) && math::is_squarable(velocity) &&
          math::is_squarable(local_center_of_mass)))
    {
        return false;
    }

    // The body's own kinetic energy must be finite too (issue #51). A squarable
    // velocity is not enough once the mass or inertia is huge: mass 1e300 at
    // speed 1e5 validated clean while every energy the world reported was inf.
    // A finite energy also bounds the momenta, m |v| and I w, since a momentum
    // overflows only if its energy does. Each half is skipped when its mass or
    // inertia is infinite, exactly as the diagnostics skip it: an immovable
    // body has no energy to report, and inf * 0 would be NaN. (The orbital
    // angular momentum and the gravity potential depend on the origin and on g,
    // so they are not the body's own and are not bounded here.)
    math::Real kinetic = 0.0;
    if (!has_infinite_mass())
    {
        kinetic += math::Real{0.5} * mass * math::dot(velocity, velocity);
    }
    if (!has_infinite_inertia())
    {
        kinetic += math::Real{0.5} * inertia * angular_velocity * angular_velocity;
    }
    return math::is_finite(kinetic);
}

bool RigidBody2D::has_infinite_mass() const
{
    return !math::is_finite(mass);
}

bool RigidBody2D::has_infinite_inertia() const
{
    return !math::is_finite(inertia);
}

bool RigidBody2D::is_static() const
{
    // AND, not OR. See the header: infinity is per quantity, and a body that
    // is immovable in only one sense is still movable in the other.
    return has_infinite_mass() && has_infinite_inertia();
}

math::Real RigidBody2D::inverse_mass() const
{
    return math::is_finite(mass) ? math::Real{1} / mass : math::Real{0};
}

math::Real RigidBody2D::inverse_inertia() const
{
    return math::is_finite(inertia) ? math::Real{1} / inertia : math::Real{0};
}

math::Vec2 center_of_mass(const RigidBody2D& body)
{
    return body.position + rotate(body.local_center_of_mass, body.angle);
}

math::Vec2 to_world(const RigidBody2D& body, const math::Vec2& local_point)
{
    return body.position + rotate(local_point, body.angle);
}

math::Vec2 to_local(const RigidBody2D& body, const math::Vec2& world_point)
{
    return rotate(world_point - body.position, -body.angle);
}

math::Vec2 velocity_at(const RigidBody2D& body, const math::Vec2& world_point)
{
    // r from the centre of mass, not from the body origin. Measuring from the
    // origin is the classic mistake, and it is invisible whenever the two
    // coincide (ADR 0007).
    const math::Vec2 r = world_point - center_of_mass(body);
    const math::Real w = body.angular_velocity;
    return body.velocity + math::Vec2{-w * r.y, w * r.x};
}

// ----------------------------------------------------------------------------
// Mass properties
// ----------------------------------------------------------------------------

namespace
{
MassProperties from_area(math::Real area, math::Real second_moment,
                         const math::Vec2& centroid, math::Real density)
{
    if (!(density > math::Real{0}) || !math::is_finite(density) ||
        !(area > math::Real{0}))
    {
        return MassProperties{};
    }
    return MassProperties{density * area, density * second_moment, centroid};
}
} // namespace

MassProperties mass_properties(const collide::Circle& shape, math::Real density)
{
    return from_area(collide::area(shape), collide::second_moment_of_area(shape),
                     collide::centroid(shape), density);
}

MassProperties mass_properties(const collide::Aabb& shape, math::Real density)
{
    return from_area(collide::area(shape), collide::second_moment_of_area(shape),
                     collide::centroid(shape), density);
}

math::Real shift_inertia(math::Real inertia_about_com, math::Real mass,
                         math::Real distance)
{
    // A negative or non-finite distance is a caller error, and squaring it
    // would quietly turn it into a plausible-looking answer.
    if (!(distance >= math::Real{0}) || !math::is_finite(distance) ||
        !math::is_finite(mass) || !math::is_finite(inertia_about_com))
    {
        return inertia_about_com;
    }
    return inertia_about_com + mass * distance * distance;
}

// ----------------------------------------------------------------------------
// World
// ----------------------------------------------------------------------------

bool RigidSettings::is_valid() const
{
    if (!(restitution >= math::Real{0}) || !(restitution <= math::Real{1}) ||
        !math::is_finite(restitution) || !math::is_finite(gravity))
    {
        return false;
    }
    if (!(friction >= math::Real{0}) || !math::is_finite(friction))
    {
        return false;
    }
    for (const collide::Halfplane& plane : ground)
    {
        if (!plane.is_valid())
        {
            return false;
        }
    }
    return true;
}

RigidWorld::RigidWorld(sim_core::SimulationSettings simulation_settings,
                       std::vector<RigidBody2D> bodies, RigidSettings settings)
    : simulation_settings_{simulation_settings},
      settings_{settings},
      bodies_{std::move(bodies)},
      step_{time::FixedStep::create(simulation_settings.dt)}
{
}

sim_core::StepStatus RigidWorld::validate() const
{
    if (!simulation_settings_.is_valid()) // dt > 0, finite
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    if (!settings_.is_valid())
    {
        return sim_core::StepStatus::InvalidSettings;
    }
    for (const RigidBody2D& body : bodies_)
    {
        if (!body.is_valid())
        {
            return sim_core::StepStatus::InvalidState;
        }
    }
    return sim_core::StepStatus::Ok;
}

bool RigidWorld::apply_impulse_at(std::size_t index, const math::Vec2& impulse,
                                  const math::Vec2& world_point)
{
    if (index >= bodies_.size() || !math::is_finite(impulse) ||
        !math::is_finite(world_point))
    {
        return false;
    }
    RigidBody2D& body = bodies_[index];
    if (!body.is_valid())
    {
        return false;
    }

    // r from the centre of mass. An impulse through the centre of mass gives
    // cross(r, J) == 0 and so leaves the spin untouched.
    const math::Vec2 r = world_point - center_of_mass(body);
    // Inverse form rather than division, so an immovable body simply does not
    // move instead of needing a special case.
    body.velocity += impulse * body.inverse_mass();
    body.angular_velocity += cross(r, impulse) * body.inverse_inertia();
    return true;
}

sim_core::StepResult RigidWorld::step()
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

    // (A) detect at the start-of-step poses, speculatively (step contract), in
    //     a fixed order so a run repeats: pairs ascending, then ground. Fixed
    //     because floating-point addition is not associative, so the order is
    //     part of the observable behaviour (docs/04).
    Contacts contacts;
    detect(bodies_, settings_, contacts);

    // (B) the impacts, on the start-of-step velocities. Restitution 0 has none.
    if (!contacts.points.empty() && settings_.restitution > math::Real{0})
    {
        impact_phase(bodies_, settings_, dt, contacts);
    }

    // (1) gravity, before the position update, which keeps this semi-implicit
    //     (symplectic) Euler and matches every other world here. It is an
    //     acceleration, so it does NOT scale with mass; a body of infinite
    //     mass is skipped explicitly, because otherwise an immovable wall would
    //     start falling. The test is on the MASS, not on is_static(): a body
    //     that merely cannot spin still falls.
    const math::Vec2 gravity = settings_.gravity;
    for (RigidBody2D& body : bodies_)
    {
        if (!body.has_infinite_mass())
        {
            body.velocity += gravity * dt;
        }
    }

    // (C) every contact, after the gravity, warm-started from the last step.
    if (!contacts.points.empty())
    {
        rest_phase(bodies_, settings_, dt, contacts, warm_);
    }

    // (2) integrate the pose with the solved velocities.
    for (RigidBody2D& body : bodies_)
    {
        // The centre of mass is what translates; the body origin follows from
        // it. Advancing the origin directly would make an offset body rotate
        // about the wrong point.
        const math::Vec2 com = center_of_mass(body) + body.velocity * dt;
        body.angle += body.angular_velocity * dt;
        body.position = com - rotate(body.local_center_of_mass, body.angle);
    }

    // (D) project the remaining penetration out, positions only.
    if (!contacts.points.empty())
    {
        project(bodies_, contacts);
    }
    remember(contacts, warm_);

    for (const RigidBody2D& body : bodies_)
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

// ----------------------------------------------------------------------------
// Diagnostics
// ----------------------------------------------------------------------------

math::Vec2 total_linear_momentum(const std::vector<RigidBody2D>& bodies)
{
    math::Vec2 momentum{};
    for (const RigidBody2D& body : bodies)
    {
        // Infinite mass with zero velocity is inf * 0 = NaN, and such a body
        // never translates, so it carries no linear momentum to report. This
        // asks about the MASS specifically: a body with infinite inertia but
        // finite mass does translate, and its momentum is real.
        if (body.has_infinite_mass())
        {
            continue;
        }
        momentum += body.mass * body.velocity;
    }
    return momentum;
}

math::Real total_angular_momentum(const std::vector<RigidBody2D>& bodies)
{
    math::Real angular = 0.0;
    for (const RigidBody2D& body : bodies)
    {
        // Spin about the centre of mass, plus the orbital term about the world
        // origin. Dropping either makes a whole class of motion look conserved
        // when it is not.
        //
        // The two terms are guarded SEPARATELY, because they depend on
        // different quantities. A body with infinite inertia and finite mass
        // has no spin term (inf * 0 is NaN, and it cannot spin) but a perfectly
        // real orbital one.
        if (!body.has_infinite_inertia())
        {
            angular += body.inertia * body.angular_velocity;
        }
        if (!body.has_infinite_mass())
        {
            angular += body.mass * cross(center_of_mass(body), body.velocity);
        }
    }
    return angular;
}

math::Real total_kinetic_energy(const std::vector<RigidBody2D>& bodies)
{
    math::Real kinetic = 0.0;
    for (const RigidBody2D& body : bodies)
    {
        // Translational and rotational energy are guarded separately, for the
        // same reason as angular momentum: an infinite quantity contributes
        // nothing and would produce NaN, but the other half is still real.
        if (!body.has_infinite_mass())
        {
            kinetic +=
                math::Real{0.5} * body.mass * math::dot(body.velocity, body.velocity);
        }
        if (!body.has_infinite_inertia())
        {
            kinetic += math::Real{0.5} * body.inertia * body.angular_velocity *
                       body.angular_velocity;
        }
    }
    return kinetic;
}

math::Real total_potential_energy(const std::vector<RigidBody2D>& bodies,
                                  const math::Vec2& gravity)
{
    math::Real potential = 0.0;
    for (const RigidBody2D& body : bodies)
    {
        if (body.has_infinite_mass())
        {
            continue; // infinite mass times a position is not a number
        }
        // U = -m (g . r) with r the centre of mass, so with g pointing down a
        // higher body has more.
        potential -= body.mass * math::dot(gravity, center_of_mass(body));
    }
    return potential;
}

math::Real total_energy(const std::vector<RigidBody2D>& bodies,
                        const math::Vec2& gravity)
{
    return total_kinetic_energy(bodies) + total_potential_energy(bodies, gravity);
}
} // namespace malloy::rigid
