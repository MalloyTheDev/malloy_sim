#include <malloy/view3d/mesh.hpp>
#include <malloy/view3d/view3d.hpp>

#include <malloy/scenario/scenario.hpp>

#include <test_check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace
{
using namespace malloy;
using namespace malloy::view3d;
using malloy::math::Quat;
using malloy::math::Real;
using malloy::math::Vec3;

const sim_core::SimulationSettings sim{0.001};
constexpr Real nan_value = std::numeric_limits<Real>::quiet_NaN();
constexpr Real inf_value = std::numeric_limits<Real>::infinity();
constexpr Real pi = 3.14159265358979323846;

// Bit-for-bit equality, for values the code computes exactly.
bool same(const Vec3& a, const Vec3& b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

bool same(const Quat& a, const Quat& b)
{
    return a.w == b.w && same(a.v, b.v);
}

rigid::RigidBody3D body_at(const Vec3& p)
{
    rigid::RigidBody3D body;
    body.position = p;
    return body;
}

// The rotation matrix of a unit quaternion, written out independently of
// math::rotate, so the corner tests do not check the library against itself.
Vec3 rotate_by_matrix(const Quat& q, const Vec3& p)
{
    const Real w = q.w;
    const Real x = q.v.x;
    const Real y = q.v.y;
    const Real z = q.v.z;
    const Vec3 c0{1 - 2 * (y * y + z * z), 2 * (x * y + w * z), 2 * (x * z - w * y)};
    const Vec3 c1{2 * (x * y - w * z), 1 - 2 * (x * x + z * z), 2 * (y * z + w * x)};
    const Vec3 c2{2 * (x * z + w * y), 2 * (y * z - w * x), 1 - 2 * (x * x + y * y)};
    return c0 * p.x + c1 * p.y + c2 * p.z;
}

Vec3 sign_of_corner(std::size_t i)
{
    return Vec3{(i & 1U) ? 1.0 : -1.0, (i & 2U) ? 1.0 : -1.0, (i & 4U) ? 1.0 : -1.0};
}

// --- 1. Exact box corners -----------------------------------------------------

int test_box_corners_exact()
{
    // 120 degrees about (1, 1, 1): x -> y, y -> z, z -> x. Every component is a
    // half, so the rotation is EXACT in binary and the corners compare with ==.
    // Off-origin centre and three unequal half-widths (docs/05, asymmetry), so
    // a transposed quaternion, a wrong operation order, a doubled half-width or
    // a missing corner each lands somewhere else.
    const Quat q{0.5, 0.5, 0.5, 0.5};
    const auto corners = box_corners(Vec3{10.0, 20.0, 30.0}, q, Vec3{1.0, 2.0, 3.0});
    for (std::size_t i = 0; i < 8; ++i)
    {
        const Vec3 s = sign_of_corner(i);
        // local (sx, 2 sy, 3 sz) turns into world (3 sz, sx, 2 sy)
        const Vec3 expected{10.0 + 3.0 * s.z, 20.0 + s.x, 30.0 + 2.0 * s.y};
        MALLOY_CHECK_TRUE(same(corners[i], expected));
    }
    // The identity leaves the box axis-aligned, again exactly.
    const auto plain = box_corners(Vec3{-1.0, 2.0, 0.5}, Quat{}, Vec3{0.25, 0.5, 2.0});
    MALLOY_CHECK_TRUE(same(plain[0], Vec3{-1.25, 1.5, -1.5}));
    MALLOY_CHECK_TRUE(same(plain[7], Vec3{-0.75, 2.5, 2.5}));
    MALLOY_CHECK_TRUE(same(plain[6], Vec3{-1.25, 2.5, 2.5}));
    return 0;
}

// --- 2. A general rotation ------------------------------------------------------

int test_rotation()
{
    // 0.7 rad about (1, -2, 0.5): not exact in binary, so compared with the
    // independent matrix to 1e-12. The rounding bound of the ~10 operations on
    // values up to ~10 is about 2e-14; the smallest defect in the mutation plan
    // moves a corner by a half-width, 0.3.
    const Quat r = math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.7);
    const Vec3 c{-4.0, 7.5, 2.25};
    const Vec3 h{0.3, 1.1, 0.6};
    const auto k = box_corners(c, r, h);
    for (std::size_t i = 0; i < 8; ++i)
    {
        const Vec3 s = sign_of_corner(i);
        const Vec3 local{s.x * h.x, s.y * h.y, s.z * h.z};
        MALLOY_CHECK_TRUE(math::approx_equal(k[i], c + rotate_by_matrix(r, local), 1e-12));
    }
    // Each edge is 2h long, along the turned body axis it belongs to.
    const Vec3 axes[3] = {rotate_by_matrix(r, Vec3{1, 0, 0}),
                          rotate_by_matrix(r, Vec3{0, 1, 0}),
                          rotate_by_matrix(r, Vec3{0, 0, 1})};
    const Real lengths[3] = {2 * h.x, 2 * h.y, 2 * h.z};
    for (std::size_t e = 0; e < box_edges.size(); ++e)
    {
        const Vec3 d = k[static_cast<std::size_t>(box_edges[e][1])] -
                       k[static_cast<std::size_t>(box_edges[e][0])];
        MALLOY_CHECK_NEAR(math::length(d), lengths[e / 4], 1e-12);
        MALLOY_CHECK_NEAR(math::dot(d, axes[e / 4]), lengths[e / 4], 1e-12);
    }
    Vec3 mean{};
    for (const Vec3& p : k)
    {
        mean += p;
    }
    MALLOY_CHECK_TRUE(math::approx_equal(mean / 8.0, c, 1e-12));
    return 0;
}

// --- 3. Faces and winding -------------------------------------------------------

int test_faces_and_winding()
{
    const Quat r = math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.7);
    const Vec3 c{-4.0, 7.5, 2.25};
    const auto k = box_corners(c, r, Vec3{0.3, 1.1, 0.6});
    for (const auto& face : box_faces)
    {
        const Vec3 p0 = k[static_cast<std::size_t>(face[0])];
        const Vec3 p1 = k[static_cast<std::size_t>(face[1])];
        const Vec3 p2 = k[static_cast<std::size_t>(face[2])];
        const Vec3 p3 = k[static_cast<std::size_t>(face[3])];
        const Vec3 n = math::cross(p1 - p0, p2 - p0);
        // Counter-clockwise seen from outside: the winding normal points away
        // from the centre. And the face is flat.
        MALLOY_CHECK_TRUE(math::dot(n, (p0 + p1 + p2 + p3) * 0.25 - c) > 0.0);
        MALLOY_CHECK_NEAR(math::dot(math::normalize(n), p3 - p0), 0.0, 1e-12);
    }
    // Every corner is in exactly three faces, so no corner is missing from the
    // table and none is used twice.
    for (int i = 0; i < 8; ++i)
    {
        int faces = 0;
        for (const auto& face : box_faces)
        {
            faces += static_cast<int>(std::count(face.begin(), face.end(), i));
        }
        MALLOY_CHECK_EQ(faces, 3);
    }
    // Every edge joins corners that differ in exactly one bit, and each of the
    // three bits names four edges.
    int per_bit[3] = {0, 0, 0};
    for (const auto& e : box_edges)
    {
        const int bits = e[0] ^ e[1];
        MALLOY_CHECK_TRUE(bits == 1 || bits == 2 || bits == 4);
        ++per_bit[bits == 1 ? 0 : (bits == 2 ? 1 : 2)];
    }
    MALLOY_CHECK_TRUE(per_bit[0] == 4 && per_bit[1] == 4 && per_bit[2] == 4);
    return 0;
}

// --- 4. Builders ----------------------------------------------------------------

int test_rigid_builder()
{
    rigid::RigidBody3D box = body_at(Vec3{1.0, -2.0, 3.0});
    box.half_extents = Vec3{0.25, 0.5, 0.75};
    box.orientation = Quat{0.5, 0.5, 0.5, 0.5};
    rigid::RigidBody3D ball = body_at(Vec3{-5.0, 4.0, 0.5});
    ball.radius = 0.4;
    ball.orientation = math::from_axis_angle(Vec3{0, 1, 0}, 0.3);
    rigid::RigidBody3D bare = body_at(Vec3{0.0, 9.0, 1.0});
    bare.inertia = Vec3{1.0, 2.0, 3.0};
    // One zero half-width: not a box collider, so the sphere radius rules.
    rigid::RigidBody3D slab = body_at(Vec3{2.0, 2.0, 2.0});
    slab.half_extents = Vec3{1.0, 0.0, 1.0};
    slab.radius = 0.2;

    rigid::Rigid3DSettings settings;
    settings.ground.push_back(collide::Plane3{Vec3{0.0, 0.0, 1.0}, -2.0});
    const rigid::Rigid3DWorld world{sim, {box, ball, bare, slab}, settings};
    const Bounds view{Vec3{1.0, 2.0, 5.0}, 3.0};
    const Frame f = build_frame(world, view);

    MALLOY_CHECK_EQ(f.dropped, std::size_t{0});
    MALLOY_CHECK_EQ(f.boxes.size(), std::size_t{2});
    MALLOY_CHECK_EQ(f.spheres.size(), std::size_t{2});
    MALLOY_CHECK_EQ(f.segments.size(), std::size_t{12}); // a triad per body
    MALLOY_CHECK_TRUE(f.boxes[0].role == Role::Body && f.boxes[0].source == 0);
    MALLOY_CHECK_TRUE(f.boxes[1].role == Role::Proxy && f.boxes[1].source == 2);
    MALLOY_CHECK_TRUE(f.spheres[0].source == 1 && f.spheres[1].source == 3);
    MALLOY_CHECK_TRUE(same(f.spheres[0].center, ball.position));
    MALLOY_CHECK_EQ(f.spheres[0].radius, 0.4);
    MALLOY_CHECK_EQ(f.spheres[1].radius, 0.2);
    MALLOY_CHECK_TRUE(same(f.boxes[0].half_extents, box.half_extents));
    MALLOY_CHECK_TRUE(same(f.boxes[0].orientation, box.orientation));
    // Corner 7 is local (0.25, 0.5, 0.75), turned to (0.75, 0.25, 0.5): exact.
    MALLOY_CHECK_TRUE(same(f.boxes[0].corners[7], Vec3{1.75, -1.75, 3.5}));
    MALLOY_CHECK_TRUE(same(f.boxes[0].corners[0], Vec3{0.25, -2.25, 2.5}));
    // The proxy is the inertia box of (1, 2, 3), turned by the identity.
    MALLOY_CHECK_TRUE(same(f.boxes[1].half_extents, inertia_box(bare.inertia, bare.mass)));

    // Triads: body x, y, z, three equal arms of 1.25 times the largest
    // half-width (0.9375 here), turned by the body's quaternion.
    const Segment& bx = f.segments[0];
    const Segment& by = f.segments[1];
    const Segment& bz = f.segments[2];
    MALLOY_CHECK_TRUE(bx.role == Role::AxisX && by.role == Role::AxisY &&
                      bz.role == Role::AxisZ);
    MALLOY_CHECK_TRUE(same(bx.a, box.position));
    MALLOY_CHECK_TRUE(same(bx.b, Vec3{1.0, -2.0 + 0.9375, 3.0})); // x -> y
    MALLOY_CHECK_TRUE(same(by.b, Vec3{1.0, -2.0, 3.0 + 0.9375})); // y -> z
    MALLOY_CHECK_TRUE(same(bz.b, Vec3{1.0 + 0.9375, -2.0, 3.0})); // z -> x
    MALLOY_CHECK_TRUE(math::dot(math::cross(bx.b - bx.a, by.b - by.a), bz.b - bz.a) > 0.0);
    // The sphere's arms are 1.25 radius = 0.5 long.
    MALLOY_CHECK_TRUE(math::approx_equal(
        f.segments[3].b, ball.position + rotate_by_matrix(ball.orientation, Vec3{0.5, 0, 0}),
        1e-12));

    // The ground patch: centred straight below the view centre, on the plane,
    // half-size twice the view radius, wound counter-clockwise from above.
    MALLOY_CHECK_EQ(f.patches.size(), std::size_t{1});
    const Patch& g = f.patches[0];
    MALLOY_CHECK_TRUE(same(g.center, Vec3{1.0, 2.0, -2.0}));
    MALLOY_CHECK_EQ(g.half_size, 6.0);
    MALLOY_CHECK_EQ(g.offset, -2.0);
    for (const Vec3& c : g.corners)
    {
        MALLOY_CHECK_EQ(c.z, -2.0);
        MALLOY_CHECK_NEAR(math::distance(c, g.center), 6.0 * std::sqrt(2.0), 1e-12);
    }
    const Vec3 wn = math::cross(g.corners[1] - g.corners[0], g.corners[2] - g.corners[0]);
    MALLOY_CHECK_TRUE(math::dot(wn, g.normal) > 0.0);
    MALLOY_CHECK_EQ(f.tick, std::uint64_t{0});
    return 0;
}

int test_inertia_box()
{
    // A uniform box's own moments give back its half-widths.
    const Vec3 h{0.5, 1.0, 1.5};
    const Real m = 2.0;
    const Vec3 inertia{m / 3 * (h.y * h.y + h.z * h.z), m / 3 * (h.x * h.x + h.z * h.z),
                       m / 3 * (h.x * h.x + h.y * h.y)};
    MALLOY_CHECK_TRUE(math::approx_equal(inertia_box(inertia, m), h, 1e-12));
    // gyroscope.scn: I_z = 2.5 > I_x + I_y = 2, which no solid has; clamped to
    // a plate a twentieth as thick as it is wide.
    const Vec3 g = inertia_box(Vec3{1.0, 1.0, 2.5}, 1.0);
    MALLOY_CHECK_TRUE(math::is_finite(g) && g.x > 0 && g.y > 0 && g.z > 0);
    MALLOY_CHECK_NEAR(g.x, std::sqrt(3.75), 1e-12);
    MALLOY_CHECK_NEAR(g.z, 0.05 * g.x, 1e-15);
    // intermediate_axis.scn: (1, 2, 3) is exactly a flat plate.
    const Vec3 p = inertia_box(Vec3{1.0, 2.0, 3.0}, 1.0, 0.1);
    MALLOY_CHECK_NEAR(p.x, std::sqrt(6.0), 1e-12);
    MALLOY_CHECK_NEAR(p.y, std::sqrt(3.0), 1e-12);
    MALLOY_CHECK_NEAR(p.z, 0.1 * std::sqrt(6.0), 1e-12);
    return 0;
}

int test_tilted_plane()
{
    const collide::Plane3 plane{math::normalize(Vec3{0.6, -0.3, 0.8}), 1.5};
    MALLOY_CHECK_TRUE(plane.is_valid());
    rigid::Rigid3DSettings settings;
    settings.ground.push_back(plane);
    const rigid::Rigid3DWorld world{sim, {body_at(Vec3{})}, settings};
    const Bounds view{Vec3{-1.0, 4.0, 2.0}, 0.7};
    const Patch g = build_frame(world, view).patches.at(0);
    MALLOY_CHECK_NEAR(math::dot(plane.normal, g.center), 1.5, 1e-12);
    // The centre is the foot of the perpendicular from the view centre.
    MALLOY_CHECK_NEAR(math::length(math::cross(view.center - g.center, plane.normal)), 0.0,
                      1e-12);
    for (const Vec3& c : g.corners)
    {
        MALLOY_CHECK_NEAR(math::dot(plane.normal, c), 1.5, 1e-12);
    }
    MALLOY_CHECK_NEAR(math::length(g.u), 1.0, 1e-15);
    MALLOY_CHECK_NEAR(math::length(g.v), 1.0, 1e-15);
    MALLOY_CHECK_NEAR(math::dot(g.u, plane.normal), 0.0, 1e-15);
    MALLOY_CHECK_TRUE(math::approx_equal(math::cross(g.u, g.v), plane.normal, 1e-15));
    const Vec3 wn = math::cross(g.corners[1] - g.corners[0], g.corners[2] - g.corners[0]);
    MALLOY_CHECK_TRUE(math::dot(wn, plane.normal) > 0.0);
    return 0;
}

int test_other_builders()
{
    // particles3d: spheres, and the twelve walls exactly on the bounds.
    particles::ParticleSettings3D ps;
    ps.bounds = collide::Aabb3{Vec3{-1.0, -2.0, -3.0}, Vec3{4.0, 5.0, 6.0}};
    particles::Particle3D a;
    a.position = Vec3{0.5, -1.5, 2.0};
    a.radius = 0.25;
    const particles::ParticleWorld3D pw{sim, ps, {a}};
    const Frame pf = build_frame(pw, Bounds{});
    MALLOY_CHECK_EQ(pf.spheres.size(), std::size_t{1});
    MALLOY_CHECK_EQ(pf.spheres[0].radius, 0.25);
    MALLOY_CHECK_EQ(pf.segments.size(), std::size_t{12});
    Real wall_length = 0;
    for (const Segment& s : pf.segments)
    {
        MALLOY_CHECK_TRUE(s.role == Role::Wall);
        for (const Vec3& p : {s.a, s.b})
        {
            MALLOY_CHECK_TRUE((p.x == -1.0 || p.x == 4.0) && (p.y == -2.0 || p.y == 5.0) &&
                              (p.z == -3.0 || p.z == 6.0));
        }
        wall_length += math::distance(s.a, s.b);
    }
    MALLOY_CHECK_EQ(wall_length, 4.0 * (5.0 + 7.0 + 9.0));

    // nbody3d: points weighted by mass.
    nbody::Body3D m1;
    m1.position = Vec3{1, 2, 3};
    m1.mass = 5.0;
    nbody::Body3D m2;
    m2.position = Vec3{-1, 0, 0.5};
    m2.mass = 0.25;
    const nbody::NBody3DWorld nw{sim, nbody::NBodySettings{1.0, 0.0}, {m1, m2}};
    const Frame nf = build_frame(nw, Bounds{});
    MALLOY_CHECK_EQ(nf.points.size(), std::size_t{2});
    MALLOY_CHECK_EQ(nf.points[0].weight, 5.0);
    MALLOY_CHECK_EQ(nf.points[1].weight, 0.25);
    MALLOY_CHECK_TRUE(same(nf.points[1].position, m2.position) && nf.points[1].source == 1);

    // springs3d: extension = length - rest, positive when stretched; 3-4-5 and
    // 5-12-13 triangles keep it exact.
    springs::SpringBody3D s0;
    s0.position = Vec3{0.0, 0.0, 0.0};
    springs::SpringBody3D s1;
    s1.position = Vec3{3.0, 4.0, 0.0};
    springs::SpringBody3D s2;
    s2.position = Vec3{3.0, 4.0, 12.0};
    springs::SpringNetwork net;
    net.add(springs::Spring{0, 1, 4.0, 1.0, 0.0});  // stretched by 1
    net.add(springs::Spring{1, 2, 13.0, 1.0, 0.0}); // compressed by 1
    net.add(springs::Spring{0, 7, 1.0, 1.0, 0.0});  // no body 7: dropped
    const springs::SpringWorld3D sw{sim, net, {s0, s1, s2}};
    const Frame sf = build_frame(sw, Bounds{});
    MALLOY_CHECK_EQ(sf.points.size(), std::size_t{3});
    MALLOY_CHECK_TRUE(sf.points[2].role == Role::Node);
    MALLOY_CHECK_EQ(sf.segments.size(), std::size_t{2});
    MALLOY_CHECK_EQ(sf.segments[0].weight, 1.0);
    MALLOY_CHECK_EQ(sf.segments[1].weight, -1.0);
    MALLOY_CHECK_TRUE(same(sf.segments[1].a, s1.position) && same(sf.segments[1].b, s2.position));
    MALLOY_CHECK_EQ(sf.dropped, std::size_t{1});

    // charges3d: points weighted by charge, and an arrow per nonzero field,
    // along it, half the view radius long. The fields differ, so a swap shows.
    charges::Charge3DSettings cs;
    cs.electric = Vec3{0.0, 0.0, 2.0};
    cs.magnetic = Vec3{3.0, 0.0, 0.0};
    charges::ChargedParticle3D q;
    q.position = Vec3{1, 1, 1};
    q.charge = -2.0;
    const charges::Charge3DWorld cw{sim, cs, {q}};
    const Bounds view{Vec3{0.0, 1.0, 2.0}, 4.0};
    const Frame cf = build_frame(cw, view);
    MALLOY_CHECK_EQ(cf.points.at(0).weight, -2.0);
    MALLOY_CHECK_EQ(cf.arrows.size(), std::size_t{2});
    MALLOY_CHECK_TRUE(cf.arrows[0].role == Role::Electric);
    MALLOY_CHECK_TRUE(math::approx_equal(cf.arrows[0].tip - cf.arrows[0].origin,
                                         Vec3{0.0, 0.0, 2.0}, 1e-12));
    MALLOY_CHECK_TRUE(cf.arrows[1].role == Role::Magnetic);
    MALLOY_CHECK_TRUE(math::approx_equal(cf.arrows[1].tip - cf.arrows[1].origin,
                                         Vec3{2.0, 0.0, 0.0}, 1e-12));
    // Beside the view centre, three quarters of its radius out.
    MALLOY_CHECK_NEAR(math::distance(cf.arrows[0].origin, view.center), 3.0, 1e-12);
    // No field, no arrow.
    const charges::Charge3DWorld none{sim, charges::Charge3DSettings{}, {q}};
    MALLOY_CHECK_TRUE(build_frame(none, view).arrows.empty());
    return 0;
}

int test_builders()
{
    if (test_rigid_builder() != 0 || test_inertia_box() != 0 || test_tilted_plane() != 0 ||
        test_other_builders() != 0)
    {
        return 1;
    }
    return 0;
}

// --- 5. Non-finite input is contained -------------------------------------------

int test_non_finite_is_contained()
{
    rigid::RigidBody3D good = body_at(Vec3{1.0, 2.0, 3.0});
    good.radius = 0.5;
    rigid::RigidBody3D bad = body_at(Vec3{nan_value, 0.0, 0.0});
    bad.radius = 0.5;
    rigid::RigidBody3D worse = body_at(Vec3{0.0, inf_value, 0.0});
    worse.half_extents = Vec3{1, 1, 1};
    // Finite fields that make a non-finite shape: no mass, no proxy box.
    rigid::RigidBody3D massless = body_at(Vec3{0.0, 0.0, 0.0});
    massless.mass = 0.0;
    rigid::Rigid3DSettings settings;
    settings.ground.push_back(collide::Plane3{Vec3{0.0, 0.0, nan_value}, 0.0});
    const rigid::Rigid3DWorld world{sim, {bad, good, worse, massless}, settings};
    MALLOY_CHECK_TRUE(world.validate() != sim_core::StepStatus::Ok);
    const Frame f = build_frame(world, Bounds{});
    MALLOY_CHECK_EQ(f.dropped, std::size_t{4});
    MALLOY_CHECK_EQ(f.spheres.size(), std::size_t{1});
    MALLOY_CHECK_EQ(f.spheres[0].source, std::size_t{1});
    MALLOY_CHECK_TRUE(f.boxes.empty());
    MALLOY_CHECK_EQ(f.segments.size(), std::size_t{3}); // the good body's triad only
    MALLOY_CHECK_TRUE(is_finite(f));
    const Bounds b = fit_bounds(f);
    MALLOY_CHECK_TRUE(math::is_finite(b.center) && math::is_finite(b.radius));

    // A non-finite view drops the patches it would place, not the bodies.
    rigid::Rigid3DSettings floor;
    floor.ground.push_back(collide::Plane3{Vec3{0.0, 0.0, 1.0}, 0.0});
    const rigid::Rigid3DWorld floored{sim, {good}, floor};
    const Frame fv = build_frame(floored, Bounds{Vec3{nan_value, 0.0, 0.0}, 1.0});
    MALLOY_CHECK_EQ(fv.dropped, std::size_t{1});
    MALLOY_CHECK_EQ(fv.spheres.size(), std::size_t{1});
    MALLOY_CHECK_TRUE(is_finite(fv));

    // The other domains drop a non-finite entry and keep the rest.
    nbody::Body3D n1;
    n1.position = Vec3{nan_value, 0, 0};
    nbody::Body3D n2;
    n2.position = Vec3{1, 1, 1};
    const Frame nf = build_frame(nbody::NBody3DWorld{sim, {1.0, 0.0}, {n1, n2}}, Bounds{});
    MALLOY_CHECK_EQ(nf.dropped, std::size_t{1});
    MALLOY_CHECK_EQ(nf.points.size(), std::size_t{1});
    MALLOY_CHECK_TRUE(is_finite(nf));
    particles::Particle3D p1;
    p1.position = Vec3{0, 0, inf_value};
    particles::Particle3D p2;
    const Frame pf = build_frame(
        particles::ParticleWorld3D{sim, particles::ParticleSettings3D{}, {p1, p2}}, Bounds{});
    MALLOY_CHECK_EQ(pf.dropped, std::size_t{1});
    MALLOY_CHECK_TRUE(is_finite(pf));
    particles::ParticleSettings3D nan_box;
    nan_box.bounds.max.x = nan_value;
    const Frame pb = build_frame(particles::ParticleWorld3D{sim, nan_box, {p2}}, Bounds{});
    MALLOY_CHECK_EQ(pb.dropped, std::size_t{1});
    MALLOY_CHECK_TRUE(pb.segments.empty() && is_finite(pb));
    charges::ChargedParticle3D c1;
    c1.position = Vec3{0, nan_value, 0};
    charges::Charge3DSettings bad_field;
    bad_field.magnetic = Vec3{inf_value, 0, 0};
    const Frame cf = build_frame(charges::Charge3DWorld{sim, bad_field, {c1}}, Bounds{});
    MALLOY_CHECK_EQ(cf.dropped, std::size_t{2});
    MALLOY_CHECK_TRUE(cf.arrows.empty() && is_finite(cf));
    springs::SpringBody3D sb1;
    sb1.position = Vec3{nan_value, 0, 0};
    springs::SpringBody3D sb2;
    springs::SpringNetwork net;
    net.add(springs::Spring{0, 1, 1.0, 1.0, 0.0});
    const Frame sf = build_frame(springs::SpringWorld3D{sim, net, {sb1, sb2}}, Bounds{});
    MALLOY_CHECK_EQ(sf.dropped, std::size_t{2}); // the body, and the spring that needs it
    MALLOY_CHECK_TRUE(is_finite(sf));

    // is_finite itself sees a NaN anywhere.
    Frame planted;
    planted.points.push_back(Point{Vec3{0, 0, 0}, nan_value, 0, Role::Mass});
    MALLOY_CHECK_FALSE(is_finite(planted));
    return 0;
}

// --- 6. fit_bounds --------------------------------------------------------------

Frame frame_of(std::vector<Sphere> spheres)
{
    Frame f;
    f.spheres = std::move(spheres);
    return f;
}

int test_fit_bounds()
{
    // The outlier is FIRST and the radii are unequal, so skipping a body or
    // ignoring a radius moves the result.
    Frame f = frame_of({Sphere{Vec3{-3.0, 0.0, 0.0}, 1.0, 0, Role::Body},
                        Sphere{Vec3{5.0, 2.0, 1.0}, 0.5, 1, Role::Body}});
    f.boxes.push_back(Box{Vec3{}, Quat{}, Vec3{},
                          box_corners(Vec3{1.0, -0.5, 0.0}, Quat{}, Vec3{0.5, 0.5, 0.5}), 2,
                          Role::Body});
    // x: -4 .. 5.5, y: -1 .. 2.5, z: -1 .. 1.5
    const Bounds b = fit_bounds(f);
    MALLOY_CHECK_TRUE(math::approx_equal(b.center, Vec3{0.75, 0.75, 0.25}, 1e-15));
    MALLOY_CHECK_NEAR(b.radius, 0.5 * std::sqrt(9.5 * 9.5 + 3.5 * 3.5 + 2.5 * 2.5), 1e-12);
    for (const Sphere& s : f.spheres)
    {
        MALLOY_CHECK_TRUE(math::distance(s.center, b.center) + s.radius <=
                          b.radius * (1 + 1e-12));
    }
    for (const Vec3& c : f.boxes[0].corners)
    {
        MALLOY_CHECK_TRUE(math::distance(c, b.center) <= b.radius * (1 + 1e-12));
    }

    // Context never enlarges the fit: a wall a kilometre away, a ground patch,
    // a field arrow, a context box.
    Frame g = f;
    g.segments.push_back(Segment{Vec3{-1000, 0, 0}, Vec3{1000, 0, 0}, 0, 0, Role::Wall});
    g.patches.push_back(Patch{Vec3{0, 0, -500}, Vec3{0, 0, 1}, -500.0, Vec3{1, 0, 0},
                              Vec3{0, 1, 0}, 900.0, {}, 0, Role::Ground});
    g.arrows.push_back(Arrow{Vec3{0, 900, 0}, Vec3{0, 901, 0}, Role::Electric});
    const Bounds bg = fit_bounds(g);
    MALLOY_CHECK_TRUE(same(bg.center, b.center) && bg.radius == b.radius);
    // And every subject kind does count: a triad arm, a spring, a point.
    for (const Role role : {Role::AxisX, Role::Spring})
    {
        Frame h = f;
        h.segments.push_back(Segment{Vec3{0, 0, 0}, Vec3{0, 0, 40}, 0, 0, role});
        MALLOY_CHECK_TRUE(fit_bounds(h).radius > 20.0);
    }
    Frame pt = f;
    pt.points.push_back(Point{Vec3{0, -60, 0}, 1.0, 0, Role::Charge});
    MALLOY_CHECK_TRUE(fit_bounds(pt).radius > 30.0);

    // The degenerate cases, as fit_viewport's.
    const Bounds empty = fit_bounds(Frame{});
    MALLOY_CHECK_TRUE(same(empty.center, Vec3{}) && empty.radius == 1.0);
    Frame one;
    one.points.push_back(Point{Vec3{7, 8, 9}, 1.0, 0, Role::Mass});
    const Bounds single = fit_bounds(one);
    MALLOY_CHECK_TRUE(same(single.center, Vec3{7, 8, 9}) && single.radius == 1.0);
    one.points.push_back(Point{Vec3{nan_value, 0, 0}, 1.0, 1, Role::Mass});
    const Bounds with_nan = fit_bounds(one);
    MALLOY_CHECK_TRUE(same(with_nan.center, Vec3{7, 8, 9}) && with_nan.radius == 1.0);
    Frame all_bad;
    all_bad.points.push_back(Point{Vec3{inf_value, 0, 0}, 1.0, 0, Role::Mass});
    const Bounds none = fit_bounds(all_bad);
    MALLOY_CHECK_TRUE(same(none.center, Vec3{}) && none.radius == 1.0);
    // Finite points whose extent overflows: the documented default, which
    // then does not contain them (as fit_viewport, issue #5).
    Frame huge;
    huge.points.push_back(Point{Vec3{-1e308, 0, 0}, 1.0, 0, Role::Mass});
    huge.points.push_back(Point{Vec3{1e308, 0, 0}, 1.0, 1, Role::Mass});
    const Bounds over = fit_bounds(huge);
    MALLOY_CHECK_TRUE(same(over.center, Vec3{}) && over.radius == 1.0);
    return 0;
}

// --- 7. merge -------------------------------------------------------------------

int test_merge()
{
    const Bounds a{Vec3{0, 0, 0}, 1.0};
    const Bounds b{Vec3{5, 1, -2}, 0.5};
    const Bounds m = merge(a, b);
    MALLOY_CHECK_TRUE(math::distance(m.center, a.center) + a.radius <= m.radius * (1 + 1e-12));
    MALLOY_CHECK_TRUE(math::distance(m.center, b.center) + b.radius <= m.radius * (1 + 1e-12));
    MALLOY_CHECK_NEAR(m.radius, (std::sqrt(30.0) + 1.5) * 0.5, 1e-12);
    const Bounds inner{Vec3{0.2, 0, 0}, 0.3};
    MALLOY_CHECK_TRUE(same(merge(a, inner).center, a.center) && merge(a, inner).radius == 1.0);
    MALLOY_CHECK_TRUE(same(merge(inner, a).center, a.center) && merge(inner, a).radius == 1.0);
    // A bad side is ignored rather than spreading a NaN.
    const Bounds broken{Vec3{nan_value, 0, 0}, 1.0};
    MALLOY_CHECK_TRUE(same(merge(a, broken).center, a.center));
    MALLOY_CHECK_TRUE(same(merge(broken, b).center, b.center));
    return 0;
}

// --- 8. project -----------------------------------------------------------------

int test_project()
{
    // Camera on -y looking at the origin, z up, 90 degree vertical field.
    Camera cam;
    cam.target = Vec3{0, 0, 0};
    cam.position = Vec3{0, -10, 0};
    cam.up = Vec3{0, 0, 1};
    cam.fovy = pi / 2;
    const Real t = std::tan(cam.fovy / 2);
    const Vec3 ndc = project(cam, 2.0, Vec3{2.0, 0.0, 1.0});
    // +x is on the RIGHT of the screen for a right-handed, z-up world.
    MALLOY_CHECK_NEAR(ndc.x, 2.0 / (10.0 * 2.0 * t), 1e-15);
    MALLOY_CHECK_NEAR(ndc.y, 1.0 / (10.0 * t), 1e-15);
    MALLOY_CHECK_NEAR(ndc.z, 10.0, 1e-15);
    MALLOY_CHECK_TRUE(ndc.x > 0.0 && ndc.y > 0.0);
    // How big a pixel is at a depth: 5 px of a 500 px viewport at depth 10.
    MALLOY_CHECK_NEAR(world_size_of_pixels(cam, 500.0, 10.0, 5.0), 0.2 * t, 1e-15);
    return 0;
}

// --- 9. fit_camera --------------------------------------------------------------

// Every subject extreme point of a frame projects within the bound, and
// between the clip planes.
int check_framed(const Frame& f, const Camera& cam, Real aspect, Real bound)
{
    std::vector<Vec3> extremes;
    const Vec3 fwd = math::normalize(cam.target - cam.position);
    const Vec3 side = math::normalize(math::cross(fwd, cam.up));
    const Vec3 up = math::cross(side, fwd);
    for (const Sphere& s : f.spheres)
    {
        for (const Vec3& d : {side, up, fwd, -side, -up, -fwd})
        {
            extremes.push_back(s.center + d * s.radius);
        }
    }
    for (const Box& b : f.boxes)
    {
        extremes.insert(extremes.end(), b.corners.begin(), b.corners.end());
    }
    for (const Segment& s : f.segments)
    {
        if (is_subject(s.role))
        {
            extremes.push_back(s.a);
            extremes.push_back(s.b);
        }
    }
    for (const Point& p : f.points)
    {
        extremes.push_back(p.position);
    }
    for (const Vec3& p : extremes)
    {
        const Vec3 ndc = project(cam, aspect, p);
        MALLOY_CHECK_TRUE(std::abs(ndc.x) <= bound + 1e-12);
        MALLOY_CHECK_TRUE(std::abs(ndc.y) <= bound + 1e-12);
        MALLOY_CHECK_TRUE(ndc.z > cam.near_plane && ndc.z < cam.far_plane);
    }
    return 0;
}

// The derived bound: subjects lie within the fitted sphere of radius R, which
// is framed as R (1 + m) filling the narrower half-angle theta, so none sits
// further off axis than asin(sin(theta) / (1 + m)).
Real ndc_bound(const ViewAngles& a)
{
    const Real half_y = a.fovy / 2;
    const Real theta = std::min(half_y, std::atan(a.aspect * std::tan(half_y)));
    return std::tan(std::asin(std::sin(theta) / (1 + a.margin))) / std::tan(theta);
}

Frame sample_frame()
{
    Frame f = frame_of({Sphere{Vec3{-3.0, 0.5, 0.0}, 1.0, 0, Role::Body},
                        Sphere{Vec3{5.0, 2.0, 1.0}, 0.5, 1, Role::Body}});
    f.boxes.push_back(Box{Vec3{}, Quat{}, Vec3{},
                          box_corners(Vec3{1.0, -0.5, 4.0},
                                      math::from_axis_angle(Vec3{1, 1, 0}, 0.4),
                                      Vec3{0.5, 1.5, 0.25}),
                          2, Role::Body});
    f.points.push_back(Point{Vec3{0.0, -6.0, -2.0}, 1.0, 3, Role::Mass});
    return f;
}

bool finite_camera(const Camera& c)
{
    return math::is_finite(c.position) && math::is_finite(c.target) &&
           math::is_finite(c.near_plane) && math::is_finite(c.far_plane) &&
           c.near_plane > 0 && c.far_plane > c.near_plane;
}

int test_fit_camera()
{
    const Frame f = sample_frame();
    const Bounds b = fit_bounds(f);
    for (const Real aspect : {16.0 / 9.0, 1.0, 0.5})
    {
        ViewAngles a;
        a.aspect = aspect;
        const Camera cam = fit_camera(b, a);
        // The framed sphere exactly fills the narrower field: its distance to
        // the nearer side planes is R (1 + margin), to rounding.
        const Real dist = math::distance(cam.position, cam.target);
        const Real half_y = a.fovy / 2;
        const Real half_x = std::atan(aspect * std::tan(half_y));
        MALLOY_CHECK_NEAR(dist * std::sin(std::min(half_x, half_y)) / (b.radius * (1 + a.margin)),
                          1.0, 1e-12);
        MALLOY_CHECK_TRUE(same(cam.target, b.center));
        MALLOY_CHECK_TRUE(same(cam.up, Vec3{0, 0, 1}));
        if (check_framed(f, cam, aspect, ndc_bound(a)) != 0)
        {
            return 1;
        }
        // The framed sphere's own rim touches the narrower edge.
        const Vec3 fwd = math::normalize(cam.target - cam.position);
        const Vec3 side = math::normalize(math::cross(fwd, cam.up));
        const Vec3 up = math::cross(side, fwd);
        const Vec3 limiting = (half_x < half_y) ? side : up;
        Real most = 0;
        for (int k = 0; k <= 2000; ++k)
        {
            const Real phi = pi * k / 2000.0;
            const Vec3 rim = b.center + (limiting * std::cos(phi) - fwd * std::sin(phi)) *
                                            (b.radius * (1 + a.margin));
            const Vec3 ndc = project(cam, aspect, rim);
            most = std::max(most, std::abs(half_x < half_y ? ndc.x : ndc.y));
        }
        MALLOY_CHECK_NEAR(most, 1.0, 1e-6);
    }
    // The default view's derived bound, in numbers: 0.8958 at any aspect of at
    // least 1, with the default 45 degree field and 0.1 margin.
    MALLOY_CHECK_NEAR(ndc_bound(ViewAngles{}), 0.8958508244213581, 1e-12);

    // Guards: nothing non-finite ever comes out.
    ViewAngles bad;
    bad.fovy = nan_value;
    bad.aspect = -1.0;
    bad.margin = inf_value;
    bad.yaw = nan_value;
    bad.pitch = inf_value;
    MALLOY_CHECK_TRUE(finite_camera(fit_camera(b, bad)));
    MALLOY_CHECK_TRUE(finite_camera(fit_camera(Bounds{Vec3{nan_value, 0, 0}, 1.0}, ViewAngles{})));
    MALLOY_CHECK_TRUE(finite_camera(fit_camera(Bounds{Vec3{}, 0.0}, ViewAngles{})));
    MALLOY_CHECK_TRUE(finite_camera(fit_camera(Bounds{Vec3{}, inf_value}, ViewAngles{})));
    const Camera too_big = fit_camera(Bounds{Vec3{}, 1e308}, ViewAngles{});
    MALLOY_CHECK_TRUE(finite_camera(too_big));
    MALLOY_CHECK_TRUE(same(too_big.target, Vec3{}) &&
                      math::distance(too_big.position, too_big.target) < 10.0);
    // Looking straight down is clamped, not a NaN basis.
    ViewAngles down;
    down.pitch = pi / 2;
    const Camera top = fit_camera(b, down);
    MALLOY_CHECK_TRUE(finite_camera(top));
    MALLOY_CHECK_TRUE(math::is_finite(project(top, down.aspect, b.center + Vec3{1, 0, 0})));
    return 0;
}

// --- 10. Orbit and resize -------------------------------------------------------

int test_orbit_and_resize()
{
    // A sphere fit is the same in every direction, so orbiting to any yaw and
    // pitch, or resizing the window to any aspect, never clips a subject.
    const Frame f = sample_frame();
    const Bounds b = fit_bounds(f);
    for (int i = 0; i < 24; ++i)
    {
        for (int j = -3; j <= 3; ++j)
        {
            for (const Real aspect : {0.3, 0.75, 1.0, 16.0 / 9.0, 3.0})
            {
                ViewAngles a;
                a.yaw = 2 * pi * i / 24;
                a.pitch = 0.5 * j;
                a.aspect = aspect;
                if (check_framed(f, fit_camera(b, a), aspect, ndc_bound(a)) != 0)
                {
                    std::cerr << "clipped at yaw " << a.yaw << " pitch " << a.pitch
                              << " aspect " << aspect << "\n";
                    return 1;
                }
            }
        }
    }
    // The default view looks along no world axis and in no world plane, so
    // motion along any axis shows on screen.
    const Camera cam = fit_camera(b, ViewAngles{});
    const Vec3 d = math::normalize(cam.position - cam.target);
    MALLOY_CHECK_TRUE(std::abs(d.x) > 0.3 && std::abs(d.y) > 0.3 && std::abs(d.z) > 0.3);

    // Orbit turns the eye about the target and changes nothing else.
    OrbitView view = fit_view(b, ViewAngles{});
    const OrbitView start = view;
    orbit(view, 10.0, -5.0, 0.01);
    MALLOY_CHECK_NEAR(view.yaw, start.yaw - 0.1, 1e-15);
    MALLOY_CHECK_NEAR(view.pitch, start.pitch - 0.05, 1e-15);
    MALLOY_CHECK_TRUE(same(view.target, start.target) && view.distance == start.distance);
    orbit(view, 0.0, 1e6, 0.01);
    MALLOY_CHECK_EQ(view.pitch, pitch_limit);
    orbit(view, 0.0, -1e6, 0.01);
    MALLOY_CHECK_EQ(view.pitch, -pitch_limit);
    orbit(view, nan_value, 0.0, 0.01);
    MALLOY_CHECK_TRUE(math::is_finite(view.yaw));

    // Pan: the point at the target's depth follows the cursor exactly. In a
    // 800 x 600 window, dragging by (+37, -21) pixels (screen y grows down)
    // moves the old target from the centre to (437, 279).
    {
        const Real w = 800.0;
        const Real h = 600.0;
        OrbitView p = start;
        const Vec3 old_target = p.target;
        pan(p, 37.0, -21.0, h);
        MALLOY_CHECK_NEAR(p.distance, start.distance, 0.0);
        const Vec3 ndc = project(to_camera(p, b), w / h, old_target);
        MALLOY_CHECK_NEAR((ndc.x + 1) * 0.5 * w, 437.0, 1e-9);
        MALLOY_CHECK_NEAR((1 - ndc.y) * 0.5 * h, 279.0, 1e-9);
        pan(p, 1.0, 1.0, 0.0); // no viewport, no pan
        MALLOY_CHECK_TRUE(math::is_finite(p.target));
    }

    // Zoom scales the distance by 0.85 a step, within the limits.
    {
        OrbitView z = start;
        zoom(z, 2.0, 0.1, 1000.0);
        MALLOY_CHECK_NEAR(z.distance, start.distance * 0.85 * 0.85, 1e-12);
        MALLOY_CHECK_TRUE(same(z.target, start.target));
        zoom(z, 1000.0, 0.1, 1000.0);
        MALLOY_CHECK_EQ(z.distance, 0.1);
        zoom(z, -1000.0, 0.1, 1000.0);
        MALLOY_CHECK_EQ(z.distance, 1000.0);
    }

    // The clip planes follow the scene through any zoom and pan, inside it or
    // outside it: positive, ordered, at most a thousand to one, and holding
    // the whole sphere whenever it is wholly in front.
    for (const Real distance : {0.01, 0.5, 3.0, b.radius, 2.0 * b.radius, 50.0, 1e5})
    {
        for (const Real shift : {0.0, 3.0, -40.0})
        {
            OrbitView z = start;
            z.distance = distance;
            z.target = z.target + Vec3{shift, 0.5 * shift, 0.0};
            const Camera c = to_camera(z, b);
            MALLOY_CHECK_TRUE(finite_camera(c));
            MALLOY_CHECK_TRUE(c.far_plane / c.near_plane <= 1000.0 * (1 + 1e-12));
            const Vec3 fwd = math::normalize(c.target - c.position);
            const Real middle = math::dot(fwd, b.center - c.position);
            // Wholly in front, and not grazing the near floor of far / 1000:
            // then the near plane is in front of it.
            if (middle - b.radius > 0.01 * (middle + b.radius))
            {
                MALLOY_CHECK_TRUE(c.near_plane < middle - b.radius);
            }
            MALLOY_CHECK_TRUE(c.far_plane > middle + b.radius);
        }
    }
    return 0;
}

// --- 11. Pacer ------------------------------------------------------------------

int test_pacer()
{
    const Real dt = 0.001;
    // The same wall time, cut into frames four ways, owes the same steps to
    // within one, the floor at the boundary; never a drift.
    const Real wall = 2.0;
    std::vector<std::vector<Real>> schedules;
    for (const Real hz : {30.0, 60.0, 144.0})
    {
        schedules.emplace_back(static_cast<std::size_t>(std::llround(wall * hz)), 1.0 / hz);
    }
    std::vector<Real> jitter;
    std::uint64_t s = 0x9E3779B97F4A7C15ull;
    Real sum = 0;
    while (sum < wall)
    {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        const Real step = 0.002 + 0.038 * static_cast<Real>(s >> 11) / 9007199254740992.0;
        const Real take = std::min(step, wall - sum);
        jitter.push_back(take);
        sum += take;
    }
    schedules.push_back(jitter);
    for (const auto& schedule : schedules)
    {
        Pacer pacer{dt, 100};
        pacer.rebase(0.0, 0);
        std::uint64_t tick = 0;
        Real now = 0;
        for (const Real frame : schedule)
        {
            now += frame;
            const std::uint64_t n = pacer.due(now, tick);
            MALLOY_CHECK_TRUE(n <= 100);
            tick += n;
        }
        MALLOY_CHECK_TRUE(tick + 1 >= 2000 && tick <= 2001);
    }

    // A stall is capped and not repaid: the run falls behind instead.
    {
        Pacer pacer{dt, 100};
        pacer.rebase(0.0, 0);
        MALLOY_CHECK_EQ(pacer.due(1.0, 0), std::uint64_t{100});
        MALLOY_CHECK_TRUE(pacer.behind());
        const std::uint64_t next = pacer.due(1.0 + 0.016, 100);
        MALLOY_CHECK_TRUE(next >= 15 && next <= 17);
        MALLOY_CHECK_FALSE(pacer.behind());
    }
    // Pause owes nothing, and resuming does not jump.
    {
        Pacer pacer{dt, 1000};
        pacer.rebase(0.0, 0);
        const std::uint64_t tick = pacer.due(0.5, 0);
        pacer.set_paused(true, 0.5, tick);
        MALLOY_CHECK_TRUE(pacer.paused());
        MALLOY_CHECK_EQ(pacer.due(10.0, tick), std::uint64_t{0});
        pacer.set_paused(false, 10.0, tick);
        const std::uint64_t n = pacer.due(10.1, tick);
        MALLOY_CHECK_TRUE(n >= 99 && n <= 101);
    }
    // Speed: doubling halves the wall time per step from the moment it
    // changes; nonsense speeds are ignored and extremes clamped.
    {
        Pacer pacer{dt, 100000};
        pacer.rebase(0.0, 0);
        const std::uint64_t tick = pacer.due(1.0, 0);
        pacer.set_speed(2.0, 1.0, tick);
        const std::uint64_t n = pacer.due(2.0, tick);
        MALLOY_CHECK_TRUE(n >= 1999 && n <= 2001);
        pacer.set_speed(0.0, 2.0, tick + n);
        pacer.set_speed(-1.0, 2.0, tick + n);
        pacer.set_speed(nan_value, 2.0, tick + n);
        MALLOY_CHECK_EQ(pacer.speed(), 2.0);
        pacer.set_speed(1e9, 2.0, tick + n);
        MALLOY_CHECK_EQ(pacer.speed(), 64.0);
        pacer.set_speed(1e-9, 2.0, tick + n);
        MALLOY_CHECK_EQ(pacer.speed(), 1.0 / 64.0);
    }
    // A clock that runs backwards, or a NaN one, owes nothing.
    {
        Pacer pacer{dt, 100};
        pacer.rebase(5.0, 10);
        MALLOY_CHECK_EQ(pacer.due(4.0, 10), std::uint64_t{0});
        MALLOY_CHECK_EQ(pacer.due(nan_value, 10), std::uint64_t{0});
    }
    // The run ends at end_tick, cap or no cap, and holds there.
    {
        Pacer pacer{dt, 100, 250};
        pacer.rebase(0.0, 0);
        std::uint64_t tick = 0;
        Real now = 0.0;
        for (int frame = 0; frame < 40; ++frame)
        {
            now += 0.05; // 50 steps a frame
            tick += pacer.due(now, tick);
            MALLOY_CHECK_TRUE(tick <= 250);
        }
        MALLOY_CHECK_EQ(tick, std::uint64_t{250});
        MALLOY_CHECK_EQ(pacer.due(1000.0, 250), std::uint64_t{0});
        Pacer stall{dt, 100, 30};
        stall.rebase(0.0, 0);
        MALLOY_CHECK_EQ(stall.due(5.0, 0), std::uint64_t{30}); // the end comes first
        MALLOY_CHECK_FALSE(stall.behind());
    }
    return 0;
}

// --- 14. Meshes -----------------------------------------------------------------

// Exact identity of a point, for matching shared vertices bit for bit.
std::array<unsigned char, sizeof(Real) * 3> key_of(const Vec3& p)
{
    std::array<unsigned char, sizeof(Real) * 3> k{};
    std::memcpy(k.data(), &p.x, sizeof(Real));
    std::memcpy(k.data() + sizeof(Real), &p.y, sizeof(Real));
    std::memcpy(k.data() + 2 * sizeof(Real), &p.z, sizeof(Real));
    return k;
}

// Closed and consistently oriented: every directed edge appears exactly once
// and its reverse exactly once.
template <typename Triangles>
bool closed(const Triangles& triangles)
{
    using Key = std::array<unsigned char, sizeof(Real) * 3>;
    std::vector<std::pair<Key, Key>> edges;
    for (const Triangle& t : triangles)
    {
        edges.emplace_back(key_of(t.a), key_of(t.b));
        edges.emplace_back(key_of(t.b), key_of(t.c));
        edges.emplace_back(key_of(t.c), key_of(t.a));
    }
    std::sort(edges.begin(), edges.end());
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        if (i + 1 < edges.size() && edges[i] == edges[i + 1])
        {
            return false;
        }
        if (!std::binary_search(edges.begin(), edges.end(),
                                std::make_pair(edges[i].second, edges[i].first)))
        {
            return false;
        }
    }
    return true;
}

// The signed volume enclosed, by tetrahedra from `origin`.
template <typename Triangles>
Real volume(const Triangles& triangles, const Vec3& origin)
{
    Real v = 0;
    for (const Triangle& t : triangles)
    {
        v += math::dot(t.a - origin, math::cross(t.b - origin, t.c - origin)) / 6.0;
    }
    return v;
}

int test_box_mesh()
{
    const Vec3 c{-4.0, 7.5, 2.25};
    const Vec3 h{0.3, 1.1, 0.6};
    const Quat r = math::from_axis_angle(Vec3{1.0, -2.0, 0.5}, 0.7);
    const Box box{c, r, h, box_corners(c, r, h), 0, Role::Body};
    const std::array<Triangle, 12> mesh = box_triangles(box);
    for (const Triangle& t : mesh)
    {
        const Vec3 wound = math::cross(t.b - t.a, t.c - t.a);
        MALLOY_CHECK_TRUE(math::dot(wound, t.normal) > 0.0);             // CCW
        MALLOY_CHECK_TRUE(math::dot((t.a + t.b + t.c) / 3.0 - c, t.normal) > 0.0); // outward
        MALLOY_CHECK_NEAR(math::length(t.normal), 1.0, 1e-15);
        MALLOY_CHECK_TRUE(math::approx_equal(math::normalize(wound), t.normal, 1e-12));
        // Every vertex is one of the box's own corners, bit for bit.
        for (const Vec3& p : {t.a, t.b, t.c})
        {
            MALLOY_CHECK_TRUE(std::any_of(box.corners.begin(), box.corners.end(),
                                          [&](const Vec3& k) { return same(k, p); }));
        }
    }
    MALLOY_CHECK_TRUE(closed(mesh));
    // The enclosed volume is the box's, 8 hx hy hz, to a relative 1e-12 (the
    // rounding of twelve triple products of unit-sized values is ~1e-15).
    MALLOY_CHECK_NEAR(volume(mesh, c) / (8 * h.x * h.y * h.z), 1.0, 1e-12);
    return 0;
}

int test_sphere_mesh()
{
    const Vec3 c{3.0, -2.0, 7.0};
    const Real r = 0.75;
    const SphereDetail details[] = {{16, 24}, {10, 16}, {6, 10}, {4, 6}, {2, 3}, {1, 1}};
    for (const SphereDetail& d : details)
    {
        const std::vector<Triangle> mesh = sphere_triangles(c, r, d);
        const int rings = std::max(d.rings, 2);
        const int slices = std::max(d.slices, 3);
        MALLOY_CHECK_EQ(mesh.size(), static_cast<std::size_t>(slices * (2 * rings - 2)));
        Real inner = r;
        for (const Triangle& t : mesh)
        {
            // Every vertex on the sphere. Each coordinate carries a few ulps of
            // the centre's magnitude (~7, so ~2e-15) and of sin and cos times
            // r; 1e-14 is about three times that bound, and a radius error a
            // mutation could make is 1e-1 or more.
            for (const Vec3& p : {t.a, t.b, t.c})
            {
                MALLOY_CHECK_NEAR(math::distance(p, c), r, 1e-14);
            }
            MALLOY_CHECK_TRUE(math::dot((t.a + t.b + t.c) / 3.0 - c, t.normal) > 0.0); // outward
            MALLOY_CHECK_TRUE(math::dot(math::cross(t.b - t.a, t.c - t.a), t.normal) > 0.0);
            inner = std::min(inner, math::dot(t.normal, t.a - c));
        }
        MALLOY_CHECK_TRUE(closed(mesh));
        // As fine as its detail says: each face's corners lie within half a
        // band (pi / 2 rings) and half a slice (pi / slices) of the face's
        // middle direction, so every face plane is at least
        // r cos(pi / 2 rings) cos(pi / slices) from the centre. A mesh cut
        // coarser than it claims falls below that.
        const Real sag_bound = r * std::cos(pi / (2 * rings)) * std::cos(pi / slices);
        MALLOY_CHECK_TRUE(inner >= sag_bound * (1 - 1e-12));
        // Inscribed: the mesh holds the ball of its nearest face plane and lies
        // inside the true sphere.
        const Real ball = 4.0 / 3.0 * pi * r * r * r;
        const Real enclosed = volume(mesh, c);
        MALLOY_CHECK_TRUE(enclosed < ball);
        MALLOY_CHECK_TRUE(enclosed >= 4.0 / 3.0 * pi * inner * inner * inner);
    }

    // Level of detail by count.
    const auto is = [](const SphereDetail& d, int rings, int slices) {
        return d.rings == rings && d.slices == slices;
    };
    MALLOY_CHECK_TRUE(is(sphere_detail(1), 16, 24) && is(sphere_detail(64), 16, 24));
    MALLOY_CHECK_TRUE(is(sphere_detail(65), 10, 16) && is(sphere_detail(512), 10, 16));
    MALLOY_CHECK_TRUE(is(sphere_detail(513), 6, 10) && is(sphere_detail(4096), 6, 10));
    MALLOY_CHECK_TRUE(is(sphere_detail(4097), 4, 6));
    return 0;
}

int test_patch_mesh_and_grid()
{
    // A tilted plane, and two patches of the same size centred at different
    // places on it, as a re-fitted camera would make.
    const collide::Plane3 plane{math::normalize(Vec3{0.3, -0.4, 0.85}), -1.25};
    rigid::Rigid3DSettings settings;
    settings.ground.push_back(plane);
    const rigid::Rigid3DWorld world{sim, {body_at(Vec3{})}, settings};
    const Patch first = build_frame(world, Bounds{Vec3{0.3, 0.1, 2.0}, 4.3}).patches.at(0);
    const Patch second = build_frame(world, Bounds{Vec3{7.77, -3.21, 0.5}, 4.3}).patches.at(0);

    const std::array<Triangle, 2> tris = patch_triangles(first);
    for (const Triangle& t : tris)
    {
        MALLOY_CHECK_TRUE(math::dot(math::cross(t.b - t.a, t.c - t.a), plane.normal) > 0.0);
        MALLOY_CHECK_TRUE(same(t.normal, plane.normal));
    }
    MALLOY_CHECK_NEAR(math::length(math::cross(tris[0].b - tris[0].a, tris[0].c - tris[0].a)) +
                          math::length(math::cross(tris[1].b - tris[1].a, tris[1].c - tris[1].a)),
                      2.0 * (2 * first.half_size) * (2 * first.half_size), 1e-9);

    const Real spacing = grid_spacing(2 * first.half_size, 20);
    // 2 x 8.6 / 20 = 0.86, which is 8.6 x 10^-1, and 8.6 rounds to the nice 10.
    MALLOY_CHECK_EQ(spacing, 1.0);
    const Vec3 anchor = plane.normal * plane.offset;
    for (const Patch* p : {&first, &second})
    {
        const std::vector<Segment> grid = patch_grid(*p);
        std::size_t along_u = 0;
        std::size_t along_v = 0;
        for (const Segment& line : grid)
        {
            MALLOY_CHECK_TRUE(line.role == Role::Ground);
            for (const Vec3& q : {line.a, line.b})
            {
                MALLOY_CHECK_NEAR(math::dot(plane.normal, q), plane.offset, 1e-12);
                MALLOY_CHECK_TRUE(std::abs(math::dot(q - p->center, p->u)) <= p->half_size + 1e-9);
                MALLOY_CHECK_TRUE(std::abs(math::dot(q - p->center, p->v)) <= p->half_size + 1e-9);
            }
            // A line runs along u or along v, and sits at a whole number of
            // spacings from the PLANE's anchor, wherever the patch is.
            const Vec3 d = line.b - line.a;
            const bool runs_along_v = std::abs(math::dot(d, p->u)) < 1e-9;
            const Real at = math::dot(line.a - anchor, runs_along_v ? p->u : p->v) / spacing;
            MALLOY_CHECK_NEAR(at, std::round(at), 1e-9);
            ++(runs_along_v ? along_u : along_v);
        }
        // About twenty a way: the nice spacing is between 0.57 and 1.43 times
        // the raw one, so 2h / spacing + 1 lies in [14, 36].
        MALLOY_CHECK_TRUE(along_u >= 13 && along_u <= 36);
        MALLOY_CHECK_TRUE(along_v >= 13 && along_v <= 36);
    }
    // The nice spacings, exactly.
    MALLOY_CHECK_EQ(grid_spacing(20.0, 20), 1.0);
    MALLOY_CHECK_EQ(grid_spacing(30.0, 20), 2.0);
    MALLOY_CHECK_EQ(grid_spacing(100.0, 20), 5.0);
    MALLOY_CHECK_EQ(grid_spacing(160.0, 20), 10.0);
    MALLOY_CHECK_EQ(grid_spacing(0.0, 20), 0.0);
    MALLOY_CHECK_EQ(grid_spacing(nan_value, 20), 0.0);
    MALLOY_CHECK_EQ(grid_spacing(10.0, 0), 0.0);
    // A patch too far out for its spacing to resolve draws no grid, and comes
    // back at once rather than looping.
    Patch far = first;
    far.center = far.center + far.u * 1e18;
    MALLOY_CHECK_TRUE(patch_grid(far).empty());
    Patch broken = first;
    broken.center.x = nan_value;
    MALLOY_CHECK_TRUE(patch_grid(broken).empty());
    return 0;
}

int test_arrow_and_shade()
{
    const Arrow arrow{Vec3{1.0, 2.0, 3.0}, Vec3{1.0, 2.0, 8.0}, Role::Magnetic};
    const std::vector<Segment> lines = arrow_lines(arrow);
    MALLOY_CHECK_EQ(lines.size(), std::size_t{5});
    MALLOY_CHECK_TRUE(same(lines[0].a, arrow.origin) && same(lines[0].b, arrow.tip));
    for (std::size_t i = 1; i < lines.size(); ++i)
    {
        MALLOY_CHECK_TRUE(same(lines[i].a, arrow.tip) && lines[i].role == Role::Magnetic);
        // A head a fifth of the shaft (1.0) back from the tip.
        MALLOY_CHECK_NEAR(lines[i].b.z, 7.0, 1e-12);
        MALLOY_CHECK_NEAR(math::distance(lines[i].b, Vec3{1.0, 2.0, 7.0}), 0.35, 1e-12);
    }
    MALLOY_CHECK_TRUE(arrow_lines(Arrow{Vec3{1, 1, 1}, Vec3{1, 1, 1}, Role::Electric}).empty());

    const Vec3 n{0.0, 0.0, 1.0};
    MALLOY_CHECK_EQ(shade(n, n, 0.25), 1.0);
    MALLOY_CHECK_EQ(shade(n, -n, 0.25), 0.25);
    MALLOY_CHECK_EQ(shade(n, Vec3{1.0, 0.0, 0.0}, 0.25), 0.25);
    MALLOY_CHECK_NEAR(shade(n, Vec3{std::sqrt(0.75), 0.0, 0.5}, 0.25), 0.625, 1e-15);
    return 0;
}

int test_meshes()
{
    if (test_box_mesh() != 0 || test_sphere_mesh() != 0 || test_patch_mesh_and_grid() != 0 ||
        test_arrow_and_shade() != 0)
    {
        return 1;
    }
    return 0;
}

// --- 12 and 13. Every 3D template -----------------------------------------------

bool equal_state(const rigid::Rigid3DWorld& a, const rigid::Rigid3DWorld& b)
{
    if (a.bodies().size() != b.bodies().size() || a.tick_count() != b.tick_count() ||
        a.warm_contacts().size() != b.warm_contacts().size())
    {
        return false;
    }
    for (std::size_t i = 0; i < a.bodies().size(); ++i)
    {
        const auto& x = a.bodies()[i];
        const auto& y = b.bodies()[i];
        if (!same(x.position, y.position) || !same(x.velocity, y.velocity) ||
            !same(x.orientation, y.orientation) ||
            !same(x.angular_velocity, y.angular_velocity))
        {
            return false;
        }
    }
    for (std::size_t i = 0; i < a.warm_contacts().size(); ++i)
    {
        const auto& x = a.warm_contacts()[i];
        const auto& y = b.warm_contacts()[i];
        if (x.a != y.a || x.b != y.b || x.ground != y.ground ||
            x.normal_impulse != y.normal_impulse || !same(x.tangent_impulse, y.tangent_impulse) ||
            !same(x.local, y.local))
        {
            return false;
        }
    }
    return true;
}

bool equal_state(const particles::ParticleWorld3D& a, const particles::ParticleWorld3D& b)
{
    if (a.particles().size() != b.particles().size() || a.tick_count() != b.tick_count() ||
        a.random_state() != b.random_state())
    {
        return false;
    }
    for (std::size_t i = 0; i < a.particles().size(); ++i)
    {
        if (!same(a.particles()[i].position, b.particles()[i].position) ||
            !same(a.particles()[i].velocity, b.particles()[i].velocity))
        {
            return false;
        }
    }
    const auto& wa = a.wall_impulse();
    const auto& wb = b.wall_impulse();
    const auto& ha = a.wall_heat();
    const auto& hb = b.wall_heat();
    return wa.min_x == wb.min_x && wa.max_x == wb.max_x && wa.min_y == wb.min_y &&
           wa.max_y == wb.max_y && wa.min_z == wb.min_z && wa.max_z == wb.max_z &&
           ha.min_x == hb.min_x && ha.max_x == hb.max_x && ha.min_y == hb.min_y &&
           ha.max_y == hb.max_y && ha.min_z == hb.min_z && ha.max_z == hb.max_z;
}

bool equal_state(const nbody::NBody3DWorld& a, const nbody::NBody3DWorld& b)
{
    if (a.bodies().size() != b.bodies().size() || a.tick_count() != b.tick_count())
    {
        return false;
    }
    for (std::size_t i = 0; i < a.bodies().size(); ++i)
    {
        if (!same(a.bodies()[i].position, b.bodies()[i].position) ||
            !same(a.bodies()[i].velocity, b.bodies()[i].velocity))
        {
            return false;
        }
    }
    return true;
}

bool equal_state(const springs::SpringWorld3D& a, const springs::SpringWorld3D& b)
{
    if (a.bodies().size() != b.bodies().size() || a.tick_count() != b.tick_count() ||
        a.forces().size() != b.forces().size())
    {
        return false;
    }
    for (std::size_t i = 0; i < a.bodies().size(); ++i)
    {
        if (!same(a.bodies()[i].position, b.bodies()[i].position) ||
            !same(a.bodies()[i].velocity, b.bodies()[i].velocity))
        {
            return false;
        }
    }
    for (std::size_t i = 0; i < a.forces().size(); ++i)
    {
        if (!same(a.forces()[i], b.forces()[i]))
        {
            return false;
        }
    }
    return true;
}

bool equal_state(const charges::Charge3DWorld& a, const charges::Charge3DWorld& b)
{
    if (a.particles().size() != b.particles().size() || a.tick_count() != b.tick_count())
    {
        return false;
    }
    for (std::size_t i = 0; i < a.particles().size(); ++i)
    {
        if (!same(a.particles()[i].position, b.particles()[i].position) ||
            !same(a.particles()[i].velocity, b.particles()[i].velocity))
        {
            return false;
        }
    }
    return true;
}

// Group 12: a finite frame at step 0 and at the last step, and with the view
// grown (never shrunk) after every step, as the viewer's default auto-fit
// does, every subject inside the view at every step.
template <typename World>
int run_framed(const std::string& name, World world, int steps)
{
    if (world.validate() != sim_core::StepStatus::Ok)
    {
        std::cerr << name << ": invalid\n";
        return 1;
    }
    const ViewAngles angles;
    Bounds view = fit_bounds(build_frame(world, Bounds{}));
    Frame frame = build_frame(world, view);
    MALLOY_CHECK_TRUE(is_finite(frame));
    MALLOY_CHECK_EQ(frame.dropped, std::size_t{0});
    MALLOY_CHECK_TRUE(!frame.spheres.empty() || !frame.boxes.empty() || !frame.points.empty());
    Camera cam = fit_camera(view, angles);
    if (check_framed(frame, cam, angles.aspect, ndc_bound(angles) + 1e-9) != 0)
    {
        std::cerr << name << ": not framed at step 0\n";
        return 1;
    }
    for (int k = 1; k <= steps; ++k)
    {
        if (!world.step().ok())
        {
            std::cerr << name << ": step " << k << " failed\n";
            return 1;
        }
        frame = build_frame(world, view);
        const Bounds grown = merge(view, fit_bounds(frame));
        if (!same(grown.center, view.center) || grown.radius != view.radius)
        {
            view = grown;
            cam = fit_camera(view, angles);
        }
        if (check_framed(frame, cam, angles.aspect, ndc_bound(angles) + 1e-9) != 0)
        {
            std::cerr << name << ": a subject left the view at step " << k << '\n';
            return 1;
        }
    }
    MALLOY_CHECK_TRUE(is_finite(frame));
    MALLOY_CHECK_EQ(frame.dropped, std::size_t{0});
    MALLOY_CHECK_EQ(frame.tick, static_cast<std::uint64_t>(steps));
    return 0;
}

// Group 13: the template stepped in lock step with a twin that never builds a
// frame. Building frames, fitting and moving the camera must leave the
// simulation bit-identical at every step.
template <typename World>
int run_untouched(const std::string& name, World viewed, World plain, int steps)
{
    Bounds view = fit_bounds(build_frame(viewed, Bounds{}));
    OrbitView orbit_view = fit_view(view, ViewAngles{});
    for (int k = 1; k <= steps; ++k)
    {
        if (!viewed.step().ok() || !plain.step().ok())
        {
            std::cerr << name << ": step " << k << " failed\n";
            return 1;
        }
        const Frame frame = build_frame(viewed, view);
        view = merge(view, fit_bounds(frame));
        orbit(orbit_view, 1.0, 0.5, 0.001);
        const Camera cam = to_camera(orbit_view, view);
        for (const Box& b : frame.boxes)
        {
            (void)box_triangles(b);
        }
        for (const Patch& p : frame.patches)
        {
            (void)patch_grid(p);
        }
        if (!frame.spheres.empty())
        {
            (void)sphere_triangles(frame.spheres.front().center, frame.spheres.front().radius,
                                   sphere_detail(frame.spheres.size()));
        }
        (void)cam;
        if (!equal_state(viewed, plain))
        {
            std::cerr << name << ": building frames changed the simulation at step " << k
                      << '\n';
            return 1;
        }
    }
    return 0;
}

std::vector<std::filesystem::path> templates()
{
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(MALLOY_SCENARIO_DIR))
    {
        if (entry.path().extension() == ".scn")
        {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

// Runs `check` on every 3D template, giving it two fresh worlds built from the
// parsed scenario, and counts them.
template <typename Check>
int every_3d_template(const Check& check, int& count)
{
    using scenario::ScenarioType;
    count = 0;
    for (const std::filesystem::path& path : templates())
    {
        const scenario::ScenarioParseResult r = scenario::parse_scenario_file(path.string());
        MALLOY_CHECK_TRUE(r.ok);
        const scenario::Scenario& s = r.scenario;
        const std::string name = path.stem().string();
        int rc = 0;
        switch (s.type)
        {
        case ScenarioType::Rigid3D:
            rc = check(name,
                       rigid::Rigid3DWorld{s.simulation, s.rigid_bodies3d, s.rigid3d_settings},
                       rigid::Rigid3DWorld{s.simulation, s.rigid_bodies3d, s.rigid3d_settings},
                       s.steps);
            break;
        case ScenarioType::Particles3D:
            rc = check(name,
                       particles::ParticleWorld3D{s.simulation, s.particle3d_settings,
                                                  s.particle_list3d},
                       particles::ParticleWorld3D{s.simulation, s.particle3d_settings,
                                                  s.particle_list3d},
                       s.steps);
            break;
        case ScenarioType::NBody3D:
            rc = check(name, nbody::NBody3DWorld{s.simulation, s.nbody_settings, s.bodies3d},
                       nbody::NBody3DWorld{s.simulation, s.nbody_settings, s.bodies3d},
                       s.steps);
            break;
        case ScenarioType::Springs3D:
            rc = check(name,
                       springs::SpringWorld3D{s.simulation, s.spring_network3d,
                                              s.spring_bodies3d},
                       springs::SpringWorld3D{s.simulation, s.spring_network3d,
                                              s.spring_bodies3d},
                       s.steps);
            break;
        case ScenarioType::Charges3D:
            rc = check(name,
                       charges::Charge3DWorld{s.simulation, s.charge3d_settings, s.charge_list3d},
                       charges::Charge3DWorld{s.simulation, s.charge3d_settings, s.charge_list3d},
                       s.steps);
            break;
        case ScenarioType::NBody:
        case ScenarioType::Particles:
        case ScenarioType::Rigid:
        case ScenarioType::Springs:
        case ScenarioType::Charges:
            continue; // 2D: not drawn by the 3D view
        }
        if (rc != 0)
        {
            return rc;
        }
        ++count;
    }
    return 0;
}

int test_templates_framed()
{
    int count = 0;
    const auto check = [](const std::string& name, auto world, auto /*twin*/, int steps) {
        return run_framed(name, std::move(world), steps);
    };
    if (every_3d_template(check, count) != 0)
    {
        return 1;
    }
    // Every 3D type is covered, and the sweep is not vacuous: 15 templates at
    // M40, and a new one only adds to the count.
    MALLOY_CHECK_TRUE(count >= 15);
    return 0;
}

int test_templates_untouched()
{
    int count = 0;
    const auto check = [](const std::string& name, auto viewed, auto plain, int steps) {
        return run_untouched(name, std::move(viewed), std::move(plain), steps);
    };
    if (every_3d_template(check, count) != 0)
    {
        return 1;
    }
    MALLOY_CHECK_TRUE(count >= 15);
    return 0;
}
} // namespace

int main(int argc, char** argv)
{
    // One group by name, for a quick run while working on it.
    const std::string only = argc > 1 ? argv[1] : "";
    const std::pair<const char*, int (*)()> groups[] = {
        {"box_corners", test_box_corners_exact},
        {"rotation", test_rotation},
        {"faces_and_winding", test_faces_and_winding},
        {"builders", test_builders},
        {"non_finite", test_non_finite_is_contained},
        {"fit_bounds", test_fit_bounds},
        {"merge", test_merge},
        {"project", test_project},
        {"fit_camera", test_fit_camera},
        {"orbit_and_resize", test_orbit_and_resize},
        {"pacer", test_pacer},
        {"templates_framed", test_templates_framed},
        {"templates_untouched", test_templates_untouched},
        {"meshes", test_meshes},
    };
    int failures = 0;
    for (const auto& [name, test] : groups)
    {
        if (!only.empty() && only != name)
        {
            continue;
        }
        if (test() != 0)
        {
            std::cerr << "FAILED: " << name << '\n';
            ++failures;
        }
    }
    if (failures != 0)
    {
        std::cerr << failures << " test group(s) failed\n";
        return 1;
    }
    std::cout << "malloy_view3d_tests passed\n";
    return 0;
}
