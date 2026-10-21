#include <malloy/collide/collide.hpp>

#include <malloy/math/math.hpp>
#include <test_check.hpp>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>

using malloy::collide::Aabb;
using malloy::collide::Circle;
using malloy::collide::Halfplane;
using malloy::collide::Contact;
using malloy::collide::area;
using malloy::collide::centroid;
using malloy::collide::second_moment_of_area;
using malloy::collide::contact;
using malloy::collide::overlaps;
using malloy::math::Real;
using malloy::math::Vec2;

int main()
{
    const Real eps = 1e-12;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    // --- Shape validation ---
    {
        const Circle good{Vec2{1.0, 2.0}, 3.0};
        const Circle point_circle{Vec2{}, 0.0}; // zero radius is a point, legal
        const Circle negative{Vec2{}, -1.0};
        const Circle nan_radius{Vec2{}, nan};
        const Circle inf_radius{Vec2{}, inf};
        const Circle nan_center{Vec2{nan, 0.0}, 1.0};

        MALLOY_CHECK_TRUE(good.is_valid());
        MALLOY_CHECK_TRUE(point_circle.is_valid());
        MALLOY_CHECK_FALSE(negative.is_valid());
        MALLOY_CHECK_FALSE(nan_radius.is_valid());
        MALLOY_CHECK_FALSE(inf_radius.is_valid());
        MALLOY_CHECK_FALSE(nan_center.is_valid());
    }
    {
        const Aabb good{Vec2{0.0, 0.0}, Vec2{2.0, 3.0}};
        const Aabb degenerate{Vec2{1.0, 1.0}, Vec2{1.0, 1.0}}; // a point, legal
        const Aabb inverted_x{Vec2{2.0, 0.0}, Vec2{1.0, 3.0}};
        const Aabb inverted_y{Vec2{0.0, 5.0}, Vec2{2.0, 3.0}};
        const Aabb nan_corner{Vec2{0.0, nan}, Vec2{2.0, 3.0}};
        const Aabb inf_corner{Vec2{0.0, 0.0}, Vec2{inf, 3.0}};

        MALLOY_CHECK_TRUE(good.is_valid());
        MALLOY_CHECK_TRUE(degenerate.is_valid());
        MALLOY_CHECK_FALSE(inverted_x.is_valid());
        MALLOY_CHECK_FALSE(inverted_y.is_valid());
        MALLOY_CHECK_FALSE(nan_corner.is_valid());
        MALLOY_CHECK_FALSE(inf_corner.is_valid());
    }

    // --- Circle/circle: separated, touching, overlapping. ---
    {
        // Radii 1 and 2, centers 5 apart: a gap of 2.
        const Circle a{Vec2{0.0, 0.0}, 1.0};
        const Circle b{Vec2{5.0, 0.0}, 2.0};
        MALLOY_CHECK_FALSE(overlaps(a, b));
        MALLOY_CHECK_FALSE(contact(a, b).has_value());
    }
    {
        // Exactly touching at distance 3 == 1 + 2: a contact with zero
        // penetration, per the documented convention.
        const Circle a{Vec2{0.0, 0.0}, 1.0};
        const Circle b{Vec2{3.0, 0.0}, 2.0};
        MALLOY_CHECK_TRUE(overlaps(a, b));
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_NEAR(c->penetration, 0.0, eps);
        MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(c->point, Vec2(1.0, 0.0), eps); // the touch point
    }
    {
        // Overlapping: centers 2 apart, radii 1 and 2, so penetration is 1.
        const Circle a{Vec2{0.0, 0.0}, 1.0};
        const Circle b{Vec2{2.0, 0.0}, 2.0};
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_NEAR(c->penetration, 1.0, eps);
        MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
        // Midway between the surfaces: a's surface is at x=1, b's at x=0.
        MALLOY_CHECK_VEC2_NEAR(c->point, Vec2(0.5, 0.0), eps);
    }
    {
        // The normal is a unit vector on a diagonal, not just on an axis.
        // Centers (0,0) and (3,4): distance 5, radii 3 and 4, penetration 2.
        const Circle a{Vec2{0.0, 0.0}, 3.0};
        const Circle b{Vec2{3.0, 4.0}, 4.0};
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_NEAR(c->penetration, 2.0, eps);
        MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(0.6, 0.8), eps);
        MALLOY_CHECK_NEAR(malloy::math::length(c->normal), 1.0, eps);
    }

    // --- Circle/circle: the convention actually separates the shapes. ---
    {
        const Circle a{Vec2{0.0, 0.0}, 1.0};
        const Circle b{Vec2{1.5, 0.0}, 1.0};
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        // Moving b by +normal*penetration must leave them exactly touching.
        const Circle moved{b.center + c->normal * c->penetration, b.radius};
        const auto after = contact(a, moved);
        MALLOY_CHECK_TRUE(after.has_value()); // touching still counts
        MALLOY_CHECK_NEAR(after->penetration, 0.0, 1e-12);
    }

    // --- Circle/circle: coincident centers use the documented fallback
    //     instead of dividing by zero. ---
    {
        const Circle a{Vec2{2.0, 2.0}, 1.0};
        const Circle b{Vec2{2.0, 2.0}, 3.0};
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
        MALLOY_CHECK_NEAR(c->penetration, 4.0, eps); // the full radius sum
        MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));
    }

    // --- Circle/circle: two zero-radius points only touch when coincident. ---
    {
        const Circle p{Vec2{1.0, 1.0}, 0.0};
        const Circle q{Vec2{1.0, 1.0}, 0.0};
        const Circle r{Vec2{1.0, 1.000001}, 0.0};
        MALLOY_CHECK_TRUE(overlaps(p, q));
        MALLOY_CHECK_FALSE(overlaps(p, r));
    }

    // --- Invalid shapes never produce a contact. ---
    {
        const Circle valid{Vec2{}, 1.0};
        const Circle bad{Vec2{}, -1.0};
        const Aabb valid_box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const Aabb bad_box{Vec2{1.0, 1.0}, Vec2{-1.0, -1.0}};

        MALLOY_CHECK_FALSE(overlaps(valid, bad));
        MALLOY_CHECK_FALSE(contact(valid, bad).has_value());
        MALLOY_CHECK_FALSE(overlaps(bad, valid));
        MALLOY_CHECK_FALSE(contact(bad, valid).has_value());
        MALLOY_CHECK_FALSE(overlaps(valid_box, bad_box));
        MALLOY_CHECK_FALSE(contact(valid_box, bad_box).has_value());
        MALLOY_CHECK_FALSE(overlaps(bad, valid_box));
        MALLOY_CHECK_FALSE(contact(bad, valid_box).has_value());
    }

    // --- Coordinates so large the squared distance overflows are reported as
    //     separated rather than producing a non-finite contact. With small
    //     radii the ordinary comparison already rejects these, since
    //     inf > radius_sum^2. ---
    {
        const Circle a{Vec2{-1e300, 0.0}, 1.0};
        const Circle b{Vec2{1e300, 0.0}, 1.0};
        MALLOY_CHECK_FALSE(overlaps(a, b));
        MALLOY_CHECK_FALSE(contact(a, b).has_value());
    }
    {
        // Radii large enough that the squared radius sum overflows too. Now
        // inf > inf is false, so the comparison alone lets the pair through and
        // the explicit finite check is the only thing standing between the
        // caller and a contact reporting -infinity penetration.
        const Circle a{Vec2{-1e300, 0.0}, 1e300};
        const Circle b{Vec2{1e300, 0.0}, 1e300};
        MALLOY_CHECK_FALSE(overlaps(a, b));
        MALLOY_CHECK_FALSE(contact(a, b).has_value());
    }
    {
        // The same hazard in the circle/box path.
        const Circle huge{Vec2{-1e300, 0.0}, 1e300};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        MALLOY_CHECK_FALSE(overlaps(huge, box));
        MALLOY_CHECK_FALSE(contact(huge, box).has_value());
    }

    // --- Aabb/Aabb: separated, touching, overlapping. ---
    {
        const Aabb a{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
        const Aabb b{Vec2{2.0, 0.0}, Vec2{3.0, 1.0}};
        MALLOY_CHECK_FALSE(overlaps(a, b));
        MALLOY_CHECK_FALSE(contact(a, b).has_value());
    }
    {
        // Edge to edge at x == 1: touching, zero penetration.
        const Aabb a{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
        const Aabb b{Vec2{1.0, 0.0}, Vec2{2.0, 1.0}};
        MALLOY_CHECK_TRUE(overlaps(a, b));
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_NEAR(c->penetration, 0.0, eps);
        MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
    }
    {
        // Overlapping by 0.25 in x and 1.0 in y: x is the shallower axis, so
        // the normal is along x and b is to the right of a.
        const Aabb a{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
        const Aabb b{Vec2{0.75, 0.0}, Vec2{2.0, 1.0}};
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_NEAR(c->penetration, 0.25, eps);
        MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
    }
    {
        // Now y is the shallower axis, and b sits below a, so the normal is -y.
        const Aabb a{Vec2{0.0, 0.0}, Vec2{4.0, 4.0}};
        const Aabb b{Vec2{0.0, -3.5}, Vec2{4.0, 0.5}};
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_NEAR(c->penetration, 0.5, eps);
        MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(0.0, -1.0), eps);
    }
    {
        // Equal overlap on both axes: the documented tie-break picks x.
        const Aabb a{Vec2{0.0, 0.0}, Vec2{2.0, 2.0}};
        const Aabb b{Vec2{1.0, 1.0}, Vec2{3.0, 3.0}};
        const auto c = contact(a, b);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_NEAR(c->penetration, 1.0, eps);
        MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
    }
    {
        // A fully contained box still reports the shallowest way out.
        const Aabb outer{Vec2{0.0, 0.0}, Vec2{10.0, 10.0}};
        const Aabb inner{Vec2{1.0, 4.0}, Vec2{9.0, 6.0}};
        const auto c = contact(outer, inner);
        MALLOY_CHECK_TRUE(c.has_value());
        MALLOY_CHECK_TRUE(c->penetration > 0.0);
        MALLOY_CHECK_TRUE(malloy::math::is_finite(c->normal));
    }

    // --- Circle/Aabb: outside, touching, overlapping. ---
    {
        const Circle c{Vec2{5.0, 0.0}, 1.0};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        MALLOY_CHECK_FALSE(overlaps(c, box));
        MALLOY_CHECK_FALSE(contact(c, box).has_value());
    }
    {
        // Touching the right face exactly: centre at x=2, radius 1, face at x=1.
        const Circle c{Vec2{2.0, 0.0}, 1.0};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        MALLOY_CHECK_TRUE(overlaps(c, box));
        const auto k = contact(c, box);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_NEAR(k->penetration, 0.0, eps);
        // Normal points from the circle toward the box, so -x.
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(-1.0, 0.0), eps);
    }
    {
        // Overlapping the right face by 0.25.
        const Circle c{Vec2{1.75, 0.0}, 1.0};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const auto k = contact(c, box);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_NEAR(k->penetration, 0.25, eps);
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(-1.0, 0.0), eps);
    }
    {
        // Near a corner, so the normal is diagonal rather than axis aligned.
        // Corner at (1,1), centre at (4,5): distance 5, radius 6, depth 1.
        const Circle c{Vec2{4.0, 5.0}, 6.0};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const auto k = contact(c, box);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_NEAR(k->penetration, 1.0, eps);
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(-0.6, -0.8), eps);
        MALLOY_CHECK_NEAR(malloy::math::length(k->normal), 1.0, eps);
    }

    // --- Circle/Aabb: centre inside the box leaves through the nearest face. ---
    {
        // Centre at (0.8, 0), nearest face is the right one at x=1.
        const Circle c{Vec2{0.8, 0.0}, 0.5};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const auto k = contact(c, box);
        MALLOY_CHECK_TRUE(k.has_value());
        // Moving the box by +normal separates them, so to push the circle out
        // to the right the box normal is -x.
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(-1.0, 0.0), eps);
        MALLOY_CHECK_NEAR(k->penetration, 0.7, eps); // 0.2 to the face + 0.5 radius
    }
    {
        // Centre at (0, -0.9), nearest face is the bottom one at y=-1.
        const Circle c{Vec2{0.0, -0.9}, 0.25};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const auto k = contact(c, box);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(0.0, 1.0), eps);
        MALLOY_CHECK_NEAR(k->penetration, 0.35, eps); // 0.1 + 0.25
    }
    {
        // Exactly centred in a square box: all four faces tie, and the
        // documented tie-break picks +x.
        const Circle c{Vec2{0.0, 0.0}, 0.5};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const auto k = contact(c, box);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(1.0, 0.0), eps);
        MALLOY_CHECK_NEAR(k->penetration, 1.5, eps); // 1.0 to the face + 0.5
        MALLOY_CHECK_TRUE(malloy::math::is_finite(k->point));
    }

    // --- overlaps() and contact() must always agree. They are separate code
    //     paths, so a change to one could silently diverge from the other. ---
    {
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        for (int i = -30; i <= 30; ++i)
        {
            for (int j = -30; j <= 30; ++j)
            {
                const Vec2 centre{static_cast<Real>(i) * 0.1,
                                  static_cast<Real>(j) * 0.1};
                const Circle c{centre, 0.35};
                MALLOY_CHECK_EQ(overlaps(c, box), contact(c, box).has_value());

                const Circle other{Vec2{0.0, 0.0}, 0.5};
                MALLOY_CHECK_EQ(overlaps(c, other), contact(c, other).has_value());

                const Aabb moving{centre - Vec2{0.2, 0.2}, centre + Vec2{0.2, 0.2}};
                MALLOY_CHECK_EQ(overlaps(moving, box), contact(moving, box).has_value());
            }
        }
    }

    // --- Determinism: the same query twice gives bit-identical results,
    //     including in the degenerate branches (docs/04). ---
    {
        const Circle a{Vec2{1.0, 1.0}, 2.0};
        const Circle b{Vec2{1.0, 1.0}, 2.0}; // coincident: the fallback branch
        const auto first = contact(a, b);
        const auto second = contact(a, b);
        MALLOY_CHECK_TRUE(first.has_value() && second.has_value());
        MALLOY_CHECK_VEC2_NEAR(first->normal, second->normal, 0.0);
        MALLOY_CHECK_NEAR(first->penetration, second->penetration, 0.0);
        MALLOY_CHECK_VEC2_NEAR(first->point, second->point, 0.0);
    }

    // --- Area properties: closed-form values, no density and no mass. ---
    {
        const Real pi = 3.14159265358979323846;

        // Disc: A = pi R^2, J = pi R^4 / 2. R = 2 gives 4pi and 8pi, and the
        // two differ, so a copied formula cannot pass both.
        const Circle disc{Vec2{7.0, -3.0}, 2.0};
        MALLOY_CHECK_NEAR(area(disc), 4.0 * pi, 1e-12);
        MALLOY_CHECK_NEAR(second_moment_of_area(disc), 8.0 * pi, 1e-12);
        MALLOY_CHECK_VEC2_NEAR(centroid(disc), Vec2(7.0, -3.0), eps); // not the origin

        // Box 6 by 2, deliberately not square and not centred on the origin:
        // A = 12, J = w*h*(w^2+h^2)/12 = 12*40/12 = 40.
        const Aabb box{Vec2{1.0, 5.0}, Vec2{7.0, 7.0}};
        MALLOY_CHECK_NEAR(area(box), 12.0, eps);
        MALLOY_CHECK_NEAR(second_moment_of_area(box), 40.0, eps);
        MALLOY_CHECK_VEC2_NEAR(centroid(box), Vec2(4.0, 6.0), eps);

        // Swapping the box extents must change the centroid but not the area or
        // the polar second moment, which is symmetric in w and h.
        const Aabb rotated{Vec2{1.0, 5.0}, Vec2{3.0, 11.0}}; // 2 by 6
        MALLOY_CHECK_NEAR(area(rotated), 12.0, eps);
        MALLOY_CHECK_NEAR(second_moment_of_area(rotated), 40.0, eps);
        MALLOY_CHECK_VEC2_NEAR(centroid(rotated), Vec2(2.0, 8.0), eps);
    }

    // --- Degenerate and invalid shapes give 0 rather than a stray number. ---
    {
        const Circle point{Vec2{1.0, 2.0}, 0.0};
        MALLOY_CHECK_NEAR(area(point), 0.0, eps);
        MALLOY_CHECK_NEAR(second_moment_of_area(point), 0.0, eps);
        MALLOY_CHECK_VEC2_NEAR(centroid(point), Vec2(1.0, 2.0), eps);

        const Aabb flat{Vec2{1.0, 1.0}, Vec2{5.0, 1.0}}; // zero height
        MALLOY_CHECK_NEAR(area(flat), 0.0, eps);
        MALLOY_CHECK_NEAR(second_moment_of_area(flat), 0.0, eps);

        const Circle bad{Vec2{}, -1.0};
        const Aabb bad_box{Vec2{1.0, 1.0}, Vec2{-1.0, -1.0}};
        MALLOY_CHECK_NEAR(area(bad), 0.0, eps);
        MALLOY_CHECK_NEAR(second_moment_of_area(bad), 0.0, eps);
        MALLOY_CHECK_VEC2_NEAR(centroid(bad), Vec2(0.0, 0.0), eps);
        MALLOY_CHECK_NEAR(area(bad_box), 0.0, eps);
        MALLOY_CHECK_NEAR(second_moment_of_area(bad_box), 0.0, eps);
        MALLOY_CHECK_VEC2_NEAR(centroid(bad_box), Vec2(0.0, 0.0), eps);
    }

    // --- Two of the four circle-in-box faces were never reached with a
    //     distinguishing configuration: right and bottom were covered, left and
    //     top were not, so their normals could be anything. ---
    {
        // Centre at (-0.7, 0), nearest face is the left one at x = -1.
        const Circle c{Vec2{-0.7, 0.0}, 0.3};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const auto k = contact(c, box);
        MALLOY_CHECK_TRUE(k.has_value());
        // To push the circle out to the left, the box normal is +x.
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(1.0, 0.0), eps);
        MALLOY_CHECK_NEAR(k->penetration, 0.6, eps); // 0.3 to the face + 0.3 radius
    }
    {
        // Centre at (0, 0.85), nearest face is the top one at y = 1.
        const Circle c{Vec2{0.0, 0.85}, 0.2};
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const auto k = contact(c, box);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(0.0, -1.0), eps);
        MALLOY_CHECK_NEAR(k->penetration, 0.35, eps); // 0.15 + 0.2
    }

    // --- Box against box picks the axis of least penetration and then a SIGN
    //     from the relative centres. Every existing case put b to the +x or -y
    //     of a, so two of the four sign branches were never asserted. Each
    //     penetration below is different, so a case cannot pass on another's
    //     expected value. ---
    {
        const Aabb a{Vec2{0.0, 0.0}, Vec2{2.0, 2.0}};
        struct DirCase
        {
            Aabb b;
            Vec2 normal;
            Real penetration;
        };
        const DirCase cases[] = {
            // b to the +x: x overlaps least, b centre is to the right
            {Aabb{Vec2{1.5, 0.0}, Vec2{3.5, 2.0}}, Vec2{1.0, 0.0}, 0.5},
            // b to the -x
            {Aabb{Vec2{-1.6, 0.0}, Vec2{0.4, 2.0}}, Vec2{-1.0, 0.0}, 0.4},
            // b below: y overlaps least, b centre is lower
            {Aabb{Vec2{0.0, -1.7}, Vec2{2.0, 0.3}}, Vec2{0.0, -1.0}, 0.3},
            // b above
            {Aabb{Vec2{0.0, 1.4}, Vec2{2.0, 3.4}}, Vec2{0.0, 1.0}, 0.6},
        };
        for (const DirCase& c : cases)
        {
            const auto k = contact(a, c.b);
            MALLOY_CHECK_TRUE(k.has_value());
            MALLOY_CHECK_VEC2_NEAR(k->normal, c.normal, eps);
            MALLOY_CHECK_NEAR(k->penetration, c.penetration, eps);
        }
    }

    // --- M16: Halfplane validation. The normal must be a UNIT vector, because
    //     the signed distance the queries rely on is only a distance when it
    //     is, and a silently normalized (0, 0) would be indistinguishable from
    //     a deliberate direction. ---
    {
        MALLOY_CHECK_TRUE(Halfplane{}.is_valid()); // default is the floor (0,1) at 0
        MALLOY_CHECK_TRUE((Halfplane{Vec2{0.6, 0.8}, -4.25}.is_valid()));
        MALLOY_CHECK_TRUE((Halfplane{Vec2{0.0, -1.0}, 7.5}.is_valid()));

        MALLOY_CHECK_FALSE((Halfplane{Vec2{0.0, 0.0}, 0.0}.is_valid()));   // no direction
        MALLOY_CHECK_FALSE((Halfplane{Vec2{0.0, 2.0}, 0.0}.is_valid()));   // not unit
        MALLOY_CHECK_FALSE((Halfplane{Vec2{1.0, 1.0}, 0.0}.is_valid()));   // length sqrt(2)
        MALLOY_CHECK_FALSE((Halfplane{Vec2{0.0, 1.0}, nan}.is_valid()));
        MALLOY_CHECK_FALSE((Halfplane{Vec2{nan, 0.0}, 0.0}.is_valid()));
        MALLOY_CHECK_FALSE((Halfplane{Vec2{inf, 0.0}, 0.0}.is_valid()));
    }

    // --- Circle against halfplane. Touching counts, matching every other
    //     query, and an invalid shape never overlaps. ---
    {
        const Halfplane floor{Vec2{0.0, 1.0}, 0.0};
        MALLOY_CHECK_FALSE(overlaps(Circle{Vec2{3.0, 0.6}, 0.5}, floor));
        MALLOY_CHECK_TRUE(overlaps(Circle{Vec2{3.0, 0.5}, 0.5}, floor));  // exact touch
        MALLOY_CHECK_TRUE(overlaps(Circle{Vec2{3.0, 0.3}, 0.5}, floor));
        MALLOY_CHECK_TRUE(overlaps(Circle{Vec2{3.0, -9.0}, 0.5}, floor)); // fully inside
        MALLOY_CHECK_FALSE(overlaps(Circle{Vec2{0.0, 0.0}, -1.0}, floor));
        MALLOY_CHECK_FALSE((overlaps(Circle{Vec2{0.0, 0.0}, 1.0}, Halfplane{Vec2{0.0, 3.0}, 0.0})));
    }

    // --- Contact data, off-origin so a dropped term shows. Floor at y = 0, a
    //     disc of radius 0.5 whose centre sits 0.3 above it: 0.2 of overlap,
    //     and the representative point is midway between the two surfaces, so
    //     between y = -0.2 and y = 0. ---
    {
        const Halfplane floor{Vec2{0.0, 1.0}, 0.0};
        const auto k = contact(Circle{Vec2{3.0, 0.3}, 0.5}, floor);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(0.0, -1.0), eps);
        MALLOY_CHECK_NEAR(k->penetration, 0.2, eps);
        MALLOY_CHECK_VEC2_NEAR(k->point, Vec2(3.0, -0.1), eps);
    }
    {
        // Exactly touching: reported, with zero penetration, not dropped.
        const auto k = contact(Circle{Vec2{-2.0, 0.5}, 0.5}, Halfplane{});
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_NEAR(k->penetration, 0.0, 0.0);
        MALLOY_CHECK_VEC2_NEAR(k->point, Vec2(-2.0, 0.0), eps);
    }
    {
        // Clear of it: no contact.
        MALLOY_CHECK_FALSE(contact(Circle{Vec2{0.0, 0.6}, 0.5}, Halfplane{}).has_value());
    }
    {
        // The case that WOULD be degenerate for circle/circle: the centre lies
        // exactly on the line, so there is no centre-to-centre direction to
        // infer a normal from. The plane supplies it, so no fallback is used.
        const auto k = contact(Circle{Vec2{2.0, 0.0}, 0.5}, Halfplane{});
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(0.0, -1.0), eps);
        MALLOY_CHECK_NEAR(k->penetration, 0.5, eps);
        MALLOY_CHECK_VEC2_NEAR(k->point, Vec2(2.0, -0.25), eps);
    }
    {
        // A plane facing DOWN is a ceiling, and pushes the other way. The line
        // is y = 3, free space is below it.
        const Halfplane ceiling{Vec2{0.0, -1.0}, -3.0};
        const auto k = contact(Circle{Vec2{1.0, 2.7}, 0.5}, ceiling);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(0.0, 1.0), eps);
        MALLOY_CHECK_NEAR(k->penetration, 0.2, eps);
    }
    {
        // Off-axis and off-origin, so neither component can be dropped.
        // distance = 0.6*1 + 0.8*2 - 1 = 1.2, radius 1.5, so 0.3 of overlap.
        const Halfplane slope{Vec2{0.6, 0.8}, 1.0};
        const auto k = contact(Circle{Vec2{1.0, 2.0}, 1.5}, slope);
        MALLOY_CHECK_TRUE(k.has_value());
        MALLOY_CHECK_VEC2_NEAR(k->normal, Vec2(-0.6, -0.8), eps);
        MALLOY_CHECK_NEAR(k->penetration, 0.3, eps);
        // point = centre + normal * (radius - penetration/2)
        //       = (1, 2) + (-0.6, -0.8) * 1.35
        MALLOY_CHECK_VEC2_NEAR(k->point, Vec2(0.19, 0.92), eps);
    }
    {
        // An invalid shape on either side yields nothing rather than NaN.
        MALLOY_CHECK_FALSE(contact(Circle{Vec2{}, nan}, Halfplane{}).has_value());
        MALLOY_CHECK_FALSE((contact(Circle{Vec2{}, 1.0},
                                    Halfplane{Vec2{1.0, 1.0}, 0.0}).has_value()));
    }

    // --- The documented separation contract: moving the circle by
    //     -normal * penetration separates them EXACTLY, leaving them touching
    //     rather than overlapping or apart. Checked on a slanted plane, where
    //     both components matter. ---
    {
        const Halfplane slope{Vec2{0.6, 0.8}, 1.0};
        const Circle before{Vec2{1.0, 2.0}, 1.5};
        const auto k = contact(before, slope);
        MALLOY_CHECK_TRUE(k.has_value());

        const Circle after{before.center - k->normal * k->penetration, before.radius};
        const auto again = contact(after, slope);
        MALLOY_CHECK_TRUE(again.has_value()); // touching still reports
        MALLOY_CHECK_NEAR(again->penetration, 0.0, 1e-15);
    }

    // --- M22: Sphere and Plane3, the 3D siblings of Circle and Halfplane.
    //     Like circle against halfplane, sphere against plane has no degenerate
    //     case: the normal is the plane's own. Configurations are asymmetric on
    //     purpose (docs/05): a tilted plane whose normal has all three
    //     components nonzero, so a dropped or transposed axis shows. ---
    {
        using malloy::collide::contact;
        using malloy::collide::overlaps;
        using malloy::collide::Plane3;
        using malloy::collide::Sphere;
        using malloy::math::Vec3;
        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();

        // Validation. A unit normal is required, exactly as Halfplane requires.
        MALLOY_CHECK_TRUE(Plane3{}.is_valid()); // default floor (0,0,1) at 0
        MALLOY_CHECK_TRUE((Sphere{Vec3{1.0, -2.0, 3.0}, 0.5}.is_valid()));
        MALLOY_CHECK_TRUE((Sphere{Vec3{}, 0.0}.is_valid())); // a point
        MALLOY_CHECK_FALSE((Sphere{Vec3{}, -1.0}.is_valid()));
        MALLOY_CHECK_FALSE((Sphere{Vec3{0.0, nan3, 0.0}, 1.0}.is_valid()));
        MALLOY_CHECK_FALSE((Sphere{Vec3{}, inf3}.is_valid()));

        const Vec3 tilted = malloy::math::normalize(Vec3{1.0, 2.0, 2.0}); // |.|=3
        MALLOY_CHECK_TRUE((Plane3{tilted, -4.25}.is_valid()));
        MALLOY_CHECK_FALSE((Plane3{Vec3{0.0, 0.0, 0.0}, 0.0}.is_valid())); // no direction
        MALLOY_CHECK_FALSE((Plane3{Vec3{0.0, 0.0, 2.0}, 0.0}.is_valid())); // not unit
        MALLOY_CHECK_FALSE((Plane3{Vec3{1.0, 1.0, 1.0}, 0.0}.is_valid())); // length sqrt(3)
        MALLOY_CHECK_FALSE((Plane3{Vec3{0.0, 0.0, 1.0}, nan3}.is_valid()));
        MALLOY_CHECK_FALSE((Plane3{Vec3{nan3, 0.0, 0.0}, 0.0}.is_valid()));

        // overlaps: touching counts, an invalid shape never overlaps.
        const Plane3 floor{Vec3{0.0, 0.0, 1.0}, 0.0};
        MALLOY_CHECK_FALSE((overlaps(Sphere{Vec3{3.0, -1.0, 0.6}, 0.5}, floor)));
        MALLOY_CHECK_TRUE((overlaps(Sphere{Vec3{3.0, -1.0, 0.5}, 0.5}, floor))); // exact touch
        MALLOY_CHECK_TRUE((overlaps(Sphere{Vec3{3.0, -1.0, 0.3}, 0.5}, floor)));
        MALLOY_CHECK_TRUE((overlaps(Sphere{Vec3{3.0, -1.0, -9.0}, 0.5}, floor))); // fully inside
        MALLOY_CHECK_FALSE((overlaps(Sphere{Vec3{}, -1.0}, floor)));            // invalid sphere

        // contact against the default floor: separated gives nothing, touching
        // gives penetration 0, overlapping gives a depth and the plane's normal.
        MALLOY_CHECK_FALSE((contact(Sphere{Vec3{0.0, 0.0, 0.6}, 0.5}, floor).has_value()));
        {
            const auto touch = contact(Sphere{Vec3{-2.0, 1.0, 0.5}, 0.5}, floor);
            MALLOY_CHECK_TRUE(touch.has_value());
            MALLOY_CHECK_NEAR(touch->penetration, 0.0, 0.0);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(touch->normal, Vec3{0.0, 0.0, -1.0}, 0.0));
        }
        {
            // Centre 0.2 above the floor, radius 0.5, so depth 0.3.
            const auto k = contact(Sphere{Vec3{2.0, -3.0, 0.2}, 0.5}, floor);
            MALLOY_CHECK_TRUE(k.has_value());
            MALLOY_CHECK_NEAR(k->penetration, 0.3, eps);
            // Normal points from the sphere into the plane's solid side: -z.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(k->normal, Vec3{0.0, 0.0, -1.0}, 0.0));
            // Point midway between the surfaces: sphere bottom at z=-0.3, plane
            // at z=0, midpoint z=-0.15; x and y are the centre's.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                k->point, Vec3{2.0, -3.0, -0.15}, eps));
        }

        // The tilted plane, no fallback: the normal is exactly -plane.normal
        // even for a centre exactly on the surface.
        {
            const Plane3 slope{tilted, 1.0};
            const auto k = contact(Sphere{tilted * 1.0, 0.5}, slope); // centre on the surface
            MALLOY_CHECK_TRUE(k.has_value());
            MALLOY_CHECK_NEAR(k->penetration, 0.5, eps);       // radius - 0
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(k->normal, tilted * -1.0, 0.0));
        }

        // Invalid shapes yield no contact rather than throwing.
        MALLOY_CHECK_FALSE((contact(Sphere{Vec3{}, nan3}, floor).has_value()));
        MALLOY_CHECK_FALSE((contact(Sphere{Vec3{}, 1.0}, Plane3{Vec3{1.0, 1.0, 1.0}, 0.0}).has_value()));
    }

    // --- M24: Sphere against Sphere, the 3D sibling of Circle against Circle.
    //     Normal points from a toward b; coincident centres fall back to +x. ---
    {
        using malloy::collide::contact;
        using malloy::collide::overlaps;
        using malloy::collide::Sphere;
        using malloy::math::Vec3;

        // Separated: no overlap, no contact.
        MALLOY_CHECK_FALSE((overlaps(Sphere{Vec3{0.0, 0.0, 0.0}, 0.5},
                                     Sphere{Vec3{3.0, 0.0, 0.0}, 0.5})));
        MALLOY_CHECK_FALSE((contact(Sphere{Vec3{0.0, 0.0, 0.0}, 0.5},
                                    Sphere{Vec3{3.0, 0.0, 0.0}, 0.5}).has_value()));

        // Overlapping along a diagonal: normal from a to b, penetration is the
        // sum of radii minus the centre distance, point midway between surfaces.
        {
            const Vec3 ca{1.0, 2.0, 3.0};
            const Vec3 dir = malloy::math::normalize(Vec3{2.0, -1.0, 2.0}); // |.|=3
            const Vec3 cb = ca + dir * 2.5;   // centres 2.5 apart
            const auto k = contact(Sphere{ca, 1.5}, Sphere{cb, 2.0}); // sum 3.5
            MALLOY_CHECK_TRUE(k.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(k->normal, dir, 1e-15));
            MALLOY_CHECK_NEAR(k->penetration, 3.5 - 2.5, eps); // 1.0
            MALLOY_CHECK_TRUE(overlaps(Sphere{ca, 1.5}, Sphere{cb, 2.0}));
            // Point is a.radius - penetration/2 along the normal from a's centre.
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(
                k->point, ca + dir * (1.5 - 0.5), eps));
        }

        // Exactly touching gives penetration 0, not "no contact".
        {
            const auto k = contact(Sphere{Vec3{0.0, 0.0, 0.0}, 1.0},
                                   Sphere{Vec3{3.0, 0.0, 0.0}, 2.0}); // dist 3 == sum
            MALLOY_CHECK_TRUE(k.has_value());
            MALLOY_CHECK_NEAR(k->penetration, 0.0, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(k->normal, Vec3{1.0, 0.0, 0.0}, eps));
        }

        // Coincident centres: the fallback normal is +x, penetration the full
        // sum of radii, and no divide by zero.
        {
            const auto k = contact(Sphere{Vec3{5.0, -1.0, 2.0}, 1.0},
                                   Sphere{Vec3{5.0, -1.0, 2.0}, 1.5});
            MALLOY_CHECK_TRUE(k.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(k->normal, Vec3{1.0, 0.0, 0.0}, 0.0));
            MALLOY_CHECK_NEAR(k->penetration, 2.5, eps);
        }

        // An invalid sphere never overlaps or contacts.
        MALLOY_CHECK_FALSE((overlaps(Sphere{Vec3{}, -1.0}, Sphere{Vec3{}, 1.0})));
        MALLOY_CHECK_FALSE((contact(Sphere{Vec3{}, -1.0}, Sphere{Vec3{}, 1.0}).has_value()));
    }

    // --- M30: Box3 against Plane3. An oriented box has up to four corners
    //     against a plane at once, so this returns a manifold (a list), unlike
    //     the single-point sphere pairs. No separating-axis search is needed:
    //     the plane has one normal, and the eight corners tested against it are
    //     the whole story. ---
    {
        using malloy::collide::Box3;
        using malloy::collide::contacts;
        using malloy::collide::overlaps;
        using malloy::collide::Plane3;
        using malloy::math::from_axis_angle;
        using malloy::math::Quat;
        using malloy::math::Vec3;
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();
        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Plane3 floor{Vec3{0.0, 0.0, 1.0}, 0.0};

        // Validation: positive finite half-extents, a unit orientation, a finite
        // centre.
        MALLOY_CHECK_TRUE((Box3{Vec3{}, Vec3{0.5, 0.5, 0.5}, Quat{}}.is_valid()));
        MALLOY_CHECK_FALSE((Box3{Vec3{}, Vec3{0.0, 0.5, 0.5}, Quat{}}.is_valid())); // no width
        MALLOY_CHECK_FALSE((Box3{Vec3{}, Vec3{-1.0, 0.5, 0.5}, Quat{}}.is_valid()));
        MALLOY_CHECK_FALSE((Box3{Vec3{}, Vec3{0.5, inf3, 0.5}, Quat{}}.is_valid()));
        MALLOY_CHECK_FALSE(
            (Box3{Vec3{0.0, nan3, 0.0}, Vec3{0.5, 0.5, 0.5}, Quat{}}.is_valid()));
        MALLOY_CHECK_FALSE(
            (Box3{Vec3{}, Vec3{0.5, 0.5, 0.5}, Quat{2.0, Vec3{}}}.is_valid())); // not unit

        // Clear of the floor: no overlap, no contacts.
        MALLOY_CHECK_FALSE((overlaps(Box3{Vec3{0.0, 0.0, 1.0}, Vec3{0.5, 0.5, 0.5}, Quat{}}, floor)));
        MALLOY_CHECK_TRUE(
            contacts(Box3{Vec3{0.0, 0.0, 1.0}, Vec3{0.5, 0.5, 0.5}, Quat{}}, floor).empty());

        // An axis-aligned cube with its centre 0.3 above the floor, half-width
        // 0.5: its four bottom corners are 0.2 below the surface, its four top
        // corners well above. Four contacts, each 0.2 deep, normal -z.
        {
            const Box3 box{Vec3{0.0, 0.0, 0.3}, Vec3{0.5, 0.5, 0.5}, Quat{}};
            MALLOY_CHECK_TRUE(overlaps(box, floor));
            const auto hits = contacts(box, floor);
            MALLOY_CHECK_TRUE(hits.size() == 4);
            Real x_times_y = 0.0;
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_NEAR(hit.penetration, 0.2, eps);
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(hit.normal, Vec3{0.0, 0.0, -1.0}, 0.0));
                // Midway between the corner (z = -0.2) and the surface (z = 0).
                MALLOY_CHECK_NEAR(hit.point.z, -0.1, eps);
                x_times_y += hit.point.x * hit.point.y;
            }
            // The four bottom corners are (+-0.5, +-0.5), all four sign
            // combinations, so the products x*y cancel to zero. A corner set
            // that collapsed x and y onto the same axis would sum to 1 instead.
            MALLOY_CHECK_NEAR(x_times_y, 0.0, eps);
        }

        // Tilted 45 degrees about x, centre at z = 0.6: exactly two corners dip
        // below, each by 0.5*sqrt(2) - 0.6. A single dropped or transposed axis
        // would change the count or the depth.
        {
            const Box3 box{Vec3{0.0, 0.0, 0.6}, Vec3{0.5, 0.5, 0.5},
                           from_axis_angle(Vec3{1.0, 0.0, 0.0}, 0.78539816339744830961)};
            const auto hits = contacts(box, floor);
            MALLOY_CHECK_TRUE(hits.size() == 2);
            const Real expected = 0.5 * std::sqrt(2.0) - 0.6;
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_NEAR(hit.penetration, expected, 1e-9);
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(hit.normal, Vec3{0.0, 0.0, -1.0}, 0.0));
            }
        }

        // A tilted plane whose normal has all three components nonzero, so a
        // dropped axis in the corner test shows. The box straddles it.
        {
            const Vec3 n = malloy::math::normalize(Vec3{1.0, 2.0, 2.0});
            const Plane3 ramp{n, 0.0};
            const Box3 box{Vec3{}, Vec3{0.5, 0.5, 0.5}, Quat{}}; // centred on the surface
            MALLOY_CHECK_TRUE(overlaps(box, ramp));
            const auto hits = contacts(box, ramp);
            // The centre is on the surface, so four corners are inside and four
            // out; each reported normal is the plane's own, negated.
            MALLOY_CHECK_TRUE(hits.size() == 4);
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(hit.normal, n * -1.0, 1e-15));
                MALLOY_CHECK_TRUE(hit.penetration > 0.0);
            }
        }

        // Invalid inputs give an empty manifold.
        MALLOY_CHECK_TRUE(
            contacts(Box3{Vec3{}, Vec3{-1.0, 0.5, 0.5}, Quat{}}, floor).empty());
        MALLOY_CHECK_TRUE(contacts(Box3{Vec3{}, Vec3{0.5, 0.5, 0.5}, Quat{}},
                                   Plane3{Vec3{0.0, 0.0, 2.0}, 0.0})
                              .empty());
    }

    // --- M33: Box3 against Box3, by the separating-axis theorem. Unlike the
    //     box/plane manifold this reduces to a SINGLE contact: enough for a
    //     bounce, not a resting stack. One predicate drives both overlaps() and
    //     contact(), so the two must agree everywhere. Ground truth for the
    //     rotated cases was computed independently against all fifteen axes. ---
    {
        using malloy::collide::Box3;
        using malloy::collide::contact;
        using malloy::collide::overlaps;
        using malloy::math::from_axis_angle;
        using malloy::math::normalize;
        using malloy::math::Quat;
        using malloy::math::Vec3;

        const Vec3 half{0.5, 0.5, 0.5};
        const Quat none{}; // identity orientation

        // overlaps() and contact() run the SAME predicate, so every case below
        // asserts they agree (MALLOY_CHECK macros return from main on failure,
        // so this stays inline rather than in a helper).

        // Clear of each other along x: no overlap, no contact.
        {
            const Box3 a{Vec3{0.0, 0.0, 0.0}, half, none};
            const Box3 b{Vec3{2.0, 0.0, 0.0}, half, none};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_FALSE(c.has_value());
        }

        // Face of a (FaceA): two axis-aligned cubes overlapping 0.2 along x.
        // Normal is +x (from a toward b), penetration 0.2, and the contact x
        // lands midway in the overlap slab [0.3, 0.5].
        {
            const Box3 a{Vec3{0.0, 0.0, 0.0}, half, none};
            const Box3 b{Vec3{0.8, 0.0, 0.0}, half, none};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.2, eps);
            MALLOY_CHECK_NEAR(c->point.x, 0.4, eps);
            MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));

            // The same pair with b on the NEGATIVE side: the normal must flip to
            // -x (it always points from a toward b), which pins the sign
            // convention rather than leaving it to a lucky axis orientation.
            const Box3 b_neg{Vec3{-0.8, 0.0, 0.0}, half, none};
            const auto c_neg = contact(a, b_neg);
            MALLOY_CHECK_TRUE(overlaps(a, b_neg) == c_neg.has_value());
            MALLOY_CHECK_TRUE(c_neg.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c_neg->normal, Vec3{-1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c_neg->penetration, 0.2, eps);
            MALLOY_CHECK_NEAR(c_neg->point.x, -0.4, eps);
        }

        // Exact face touch: cubes one unit apart share a face. Touching counts,
        // so a contact exists with zero penetration; nudging b a hair further
        // apart separates them (the predicate is strict).
        {
            const Box3 a{Vec3{0.0, 0.0, 0.0}, half, none};
            const Box3 touch{Vec3{1.0, 0.0, 0.0}, half, none};
            const auto c = contact(a, touch);
            MALLOY_CHECK_TRUE(overlaps(a, touch) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_NEAR(c->penetration, 0.0, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->point.x, 0.5, eps); // the shared face
            const Box3 apart{Vec3{1.0 + 1e-6, 0.0, 0.0}, half, none};
            const auto c2 = contact(a, apart);
            MALLOY_CHECK_TRUE(overlaps(a, apart) == c2.has_value());
            MALLOY_CHECK_FALSE(c2.has_value());
        }

        // Face of b (FaceB): b is turned 30 degrees about z and pushed along its
        // OWN x-axis, so the tightest axis is b's face normal, not a's. This is
        // the mirror of the FaceA branch and would break if the two were
        // confused. The normal is b's local x; the penetration follows from the
        // projected radii along it.
        {
            const Real angle = 0.52359877559829887308; // pi / 6
            const Quat rot = from_axis_angle(Vec3{0.0, 0.0, 1.0}, angle);
            const Vec3 bx{std::cos(angle), std::sin(angle), 0.0}; // b local x in world
            const Real d = 1.0;
            const Box3 a{Vec3{0.0, 0.0, 0.0}, half, none};
            const Box3 b{bx * d, half, rot};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, bx, 1e-9));
            const Real expected = 0.5 * (std::cos(angle) + std::sin(angle)) + 0.5 - d;
            MALLOY_CHECK_NEAR(c->penetration, expected, 1e-9);
            MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));
        }

        // Edge against edge, the ninth family of axes: a is an axis-aligned cube,
        // b is tilted 45 degrees about (1,1,0) and offset along that same
        // diagonal, so the least axis is cross(a.z, b.z) = normalize(1,1,0). No
        // face of either box points that way, so a faces-only test could never
        // find it. Penetration is exactly sqrt(2) - d.
        {
            const Real d = 1.2;
            const Vec3 diag = normalize(Vec3{1.0, 1.0, 0.0});
            const Quat rot =
                from_axis_angle(Vec3{1.0, 1.0, 0.0}, 0.78539816339744830961); // pi / 4
            const Box3 a{Vec3{0.0, 0.0, 0.0}, half, none};
            const Box3 b{diag * d, half, rot};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, diag, 1e-9));
            MALLOY_CHECK_NEAR(c->penetration, std::sqrt(2.0) - d, 1e-9);
            MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));
        }

        // Coincident centres, different sizes: every direction overlaps, so no
        // axis separates and the least-combined-radius axis wins (here x, the
        // smallest radius sum). A fixed, repeatable fallback, the box analogue of
        // the sphere pair's +x, with no divide by zero.
        {
            const Vec3 c0{1.0, -1.0, 2.0};
            const Box3 a{c0, Vec3{0.5, 0.5, 0.5}, none};
            const Box3 b{c0, Vec3{0.3, 0.4, 0.6}, none};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.8, eps);
            MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));
        }

        // An invalid box never overlaps or contacts anything.
        {
            const Box3 good{Vec3{}, half, none};
            const Box3 bad{Vec3{}, Vec3{-1.0, 0.5, 0.5}, none};
            MALLOY_CHECK_FALSE(overlaps(good, bad));
            MALLOY_CHECK_FALSE(contact(good, bad).has_value());
            MALLOY_CHECK_FALSE(overlaps(bad, good));
            MALLOY_CHECK_FALSE(contact(bad, good).has_value());
        }
    }

    // --- M34: Box3 against Sphere, the last movable pair in 3D. No SAT: the
    //     nearest point of the box to the sphere centre (the centre clamped into
    //     the box frame) decides everything. The normal runs from that point to
    //     the centre; a centre inside the box is the degenerate case and pushes
    //     out the least-penetrated face. `overlaps` and `contact` agree on every
    //     case, and the reversed order negates the normal. ---
    {
        using malloy::collide::Box3;
        using malloy::collide::contact;
        using malloy::collide::overlaps;
        using malloy::collide::Sphere;
        using malloy::math::from_axis_angle;
        using malloy::math::normalize;
        using malloy::math::Quat;
        using malloy::math::Vec3;

        const Vec3 half{0.5, 0.5, 0.5};
        const Quat none{};
        const Box3 cube{Vec3{0.0, 0.0, 0.0}, half, none}; // axis-aligned unit cube

        // Clear of the box: no overlap, no contact.
        {
            const Sphere s{Vec3{2.0, 0.0, 0.0}, 0.5};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_FALSE(c.has_value());
        }

        // Face contact: a sphere on the +x side, its centre 1.0 out, radius 0.6.
        // The nearest point is the face centre (0.5, 0, 0), so the normal is +x
        // (box toward sphere), penetration 0.1, and the point sits midway in the
        // overlap along x (box face 0.5, sphere surface 0.4).
        {
            const Sphere s{Vec3{1.0, 0.0, 0.0}, 0.6};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.1, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->point, Vec3{0.45, 0.0, 0.0}, eps));

            // Reversed order negates the normal (sphere toward box) but keeps the
            // penetration and point.
            const auto r = contact(s, cube);
            MALLOY_CHECK_TRUE(overlaps(s, cube) == r.has_value());
            MALLOY_CHECK_TRUE(r.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(r->normal, Vec3{-1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(r->penetration, 0.1, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(r->point, Vec3{0.45, 0.0, 0.0}, eps));
        }

        // Exact touch: distance equals the radius, so a contact exists with zero
        // penetration; a hair further apart separates them (touching counts, the
        // boundary is inclusive).
        {
            // Centre 1.0 out, radius 0.5: the face at x = 0.5 and the sphere
            // surface meet exactly, and 1.0, 0.5 and 0.25 are all exact in
            // double, so this really is a zero-penetration touch.
            const Sphere s{Vec3{1.0, 0.0, 0.0}, 0.5};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_NEAR(c->penetration, 0.0, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            const Sphere apart{Vec3{1.0 + 1e-6, 0.0, 0.0}, 0.5};
            const auto c2 = contact(cube, apart);
            MALLOY_CHECK_TRUE(overlaps(cube, apart) == c2.has_value());
            MALLOY_CHECK_FALSE(c2.has_value());
        }

        // Corner contact: the nearest point is the +++ corner (0.5, 0.5, 0.5), so
        // the normal is normalize(1,1,1) and the penetration is radius minus the
        // corner distance sqrt(0.75).
        {
            const Sphere s{Vec3{1.0, 1.0, 1.0}, 1.0};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, normalize(Vec3{1.0, 1.0, 1.0}), 1e-9));
            MALLOY_CHECK_NEAR(c->penetration, 1.0 - std::sqrt(0.75), 1e-9);
        }

        // Edge contact: nearest point on the +x+y edge (0.5, 0.5, 0), normal
        // normalize(1,1,0), penetration radius minus sqrt(0.5).
        {
            const Sphere s{Vec3{1.0, 1.0, 0.0}, 0.8};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, normalize(Vec3{1.0, 1.0, 0.0}), 1e-9));
            MALLOY_CHECK_NEAR(c->penetration, 0.8 - std::sqrt(0.5), 1e-9);
        }

        // Sphere centre INSIDE the box: the degenerate case. Centre at
        // (0.2, 0, 0), radius 0.3, is 0.3 from the +x face (the nearest), so it
        // pushes out +x with penetration radius + 0.3 = 0.6.
        //
        // The point is midway between the two surfaces along the normal, as in
        // every other case (issue #57): the box's is the face point (0.5, 0, 0),
        // the sphere's deepest is (-0.1, 0, 0), so (0.2, 0, 0). Taking the
        // centre as the box's surface gave (0.05, 0, 0), half the face depth
        // short.
        {
            const Sphere s{Vec3{0.2, 0.0, 0.0}, 0.3};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.6, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->point, Vec3{0.2, 0.0, 0.0}, eps));

            // Reversed order: the normal is negated, the point is the same.
            const auto r = contact(s, cube);
            MALLOY_CHECK_TRUE(overlaps(s, cube) == r.has_value());
            MALLOY_CHECK_TRUE(r.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(r->normal, Vec3{-1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(r->point, Vec3{0.2, 0.0, 0.0}, eps));
        }

        // The mirror inside case: centre at (-0.2, 0, 0) is nearest the -x face,
        // so the push is -x, which pins the SIGN of the deep-contact normal (a
        // fixed +x would pass the case above but fail here). The point mirrors
        // too: (-0.2, 0, 0).
        {
            const Sphere s{Vec3{-0.2, 0.0, 0.0}, 0.3};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{-1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.6, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->point, Vec3{-0.2, 0.0, 0.0}, eps));
        }

        // A centre just inside a face, (0.45, 0, 0) with radius 0.3: face depth
        // 0.05, penetration 0.35, and the point midway between the face at 0.5
        // and the sphere's deepest point at 0.15, so 0.325 (issue #57 measured
        // 0.30 before the fix).
        {
            const Sphere s{Vec3{0.45, 0.0, 0.0}, 0.3};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.35, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->point, Vec3{0.325, 0.0, 0.0}, eps));
        }

        // Sphere centre exactly at the box centre: every face is equidistant, so
        // the tie breaks to axis 0 with a +x normal (a fixed, repeatable choice),
        // penetration radius + half-extent, no divide by zero. The point is
        // midway between the face at 0.5 and the sphere's deepest point at -0.3:
        // (0.1, 0, 0), where the centre as the box surface gave (-0.15, 0, 0).
        {
            const Sphere s{Vec3{0.0, 0.0, 0.0}, 0.3};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.8, eps);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->point, Vec3{0.1, 0.0, 0.0}, eps));
        }

        // A tilted box: rotated 45 degrees about z with the sphere placed along
        // the box's OWN x-axis, so the contact is a clean face and the normal is
        // that axis, normalize(1,1,0). This would break if the clamp used world
        // axes instead of the box frame.
        {
            const Quat rot = from_axis_angle(Vec3{0.0, 0.0, 1.0}, 0.78539816339744830961);
            const Vec3 bx = normalize(Vec3{1.0, 1.0, 0.0}); // box local x in world
            const Box3 tilted{Vec3{0.0, 0.0, 0.0}, half, rot};
            const Sphere s{bx * 1.0, 0.7};
            const auto c = contact(tilted, s);
            MALLOY_CHECK_TRUE(overlaps(tilted, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, bx, 1e-9));
            MALLOY_CHECK_NEAR(c->penetration, 0.2, 1e-9);
        }

        // Two properties over many random rotated boxes (a fixed, fully
        // specified generator, so the run repeats), both from issue #57.
        //
        // 1. A centre inside the box gives the midway point: halfway between
        //    the face `d` out along the normal and the sphere's deepest point
        //    `radius` back from the centre, with `d` the least face depth,
        //    recomputed here in the box frame, and the normal that face's axis.
        // 2. `overlaps` agrees with `contact` in both orders, and a centre well
        //    inside overlaps at ANY radius, zero included. Measuring the
        //    distance to the nearest point instead missed most zero-radius
        //    centres inside a rotated box (729 of 1000 in the issue): the point
        //    rebuilt from the box frame misses the centre by rounding, which
        //    exceeds a zero radius squared. Half the boxes sit near 1e6, where
        //    that rounding is largest, and the radii include 0, 1e-300 and 1e-12.
        //
        // "Well inside" keeps every coordinate within 0.9 of its half-extent,
        // so the face depth is at least 0.02, far above rounding even at 1e6,
        // and the old point's error (half that depth) far above the tolerance.
        {
            using malloy::math::conjugate;
            using malloy::math::rotate;
            std::uint64_t state = 0x2545F4914F6CDD1DULL;
            const auto next = [&state]() {
                state += 0x9E3779B97F4A7C15ULL; // splitmix64
                std::uint64_t z = state;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
                z ^= z >> 31;
                return static_cast<Real>(z >> 11) * (1.0 / 9007199254740992.0); // [0, 1)
            };
            const auto signed_unit = [&next]() { return 2.0 * next() - 1.0; };
            const Real radii[] = {0.0, 1e-300, 1e-12, 0.05, 0.3};

            int deep_cases = 0;
            int deep_missed = 0;     // a centre well inside that did not overlap
            int not_midway = 0;      // a deep contact off the midway point
            int disagreements = 0;   // overlaps != contact, either order
            int near_overlapping = 0;
            int near_clear = 0;
            for (int it = 0; it < 2000; ++it)
            {
                const Real scale = (it % 2 == 0) ? 1.0 : 1e6;
                const bool deep = (it % 4) < 2;
                const Vec3 centre{scale * signed_unit(), scale * signed_unit(),
                                  scale * signed_unit()};
                const Vec3 e{0.2 + next(), 0.2 + next(), 0.2 + next()};
                const Quat q = from_axis_angle(
                    Vec3{signed_unit(), signed_unit(), signed_unit()}, 3.14159 * signed_unit());
                const Box3 box{centre, e, q};
                // Deep: well inside. Otherwise near the surface, inside or out.
                const Real reach = deep ? 0.9 : 1.3;
                const Vec3 local{reach * e.x * signed_unit(), reach * e.y * signed_unit(),
                                 reach * e.z * signed_unit()};
                const Vec3 p = centre + rotate(q, local);

                for (const Real r : radii)
                {
                    const Sphere s{p, r};
                    const auto c = contact(box, s);
                    if (overlaps(box, s) != c.has_value() || overlaps(s, box) != c.has_value())
                    {
                        ++disagreements;
                    }
                    if (!deep)
                    {
                        if (c.has_value())
                        {
                            ++near_overlapping;
                        }
                        else
                        {
                            ++near_clear;
                        }
                        continue;
                    }
                    ++deep_cases;
                    if (!overlaps(box, s) || !c.has_value())
                    {
                        ++deep_missed;
                        continue;
                    }
                    // The least face depth and its axis, in the box frame.
                    const Vec3 l = rotate(conjugate(q), p - centre);
                    const Real depth[3] = {e.x - std::abs(l.x), e.y - std::abs(l.y),
                                           e.z - std::abs(l.z)};
                    const Real coord[3] = {l.x, l.y, l.z};
                    int k = 0;
                    for (int j = 1; j < 3; ++j)
                    {
                        if (depth[j] < depth[k])
                        {
                            k = j;
                        }
                    }
                    const Vec3 local_axis{k == 0 ? 1.0 : 0.0, k == 1 ? 1.0 : 0.0,
                                          k == 2 ? 1.0 : 0.0};
                    const Vec3 normal = rotate(q, local_axis) * (coord[k] >= 0.0 ? 1.0 : -1.0);
                    const Vec3 face = p + normal * depth[k];
                    const Vec3 deepest = p - normal * r;
                    const Vec3 midway = (face + deepest) * 0.5;
                    const Real tol = 1e-12 * scale;
                    if (!malloy::math::approx_equal(c->normal, normal, 1e-12) ||
                        !malloy::math::approx_equal(c->point, midway, tol) ||
                        !malloy::math::approx_equal(c->penetration, r + depth[k], tol))
                    {
                        ++not_midway;
                    }
                }
            }
            MALLOY_CHECK_EQ(disagreements, 0);
            MALLOY_CHECK_EQ(deep_missed, 0);
            MALLOY_CHECK_EQ(not_midway, 0);
            // The sweep really covers both regimes, so it cannot pass vacuously.
            MALLOY_CHECK_EQ(deep_cases, 1000 * 5);
            MALLOY_CHECK_TRUE(near_overlapping > 500);
            MALLOY_CHECK_TRUE(near_clear > 500);
        }

        // Invalid shapes never overlap or contact, in either order.
        {
            const Box3 bad_box{Vec3{}, Vec3{-1.0, 0.5, 0.5}, none};
            const Sphere bad_sphere{Vec3{}, -1.0};
            const Sphere ok_sphere{Vec3{}, 0.5};
            MALLOY_CHECK_FALSE(overlaps(bad_box, ok_sphere));
            MALLOY_CHECK_FALSE(contact(bad_box, ok_sphere).has_value());
            MALLOY_CHECK_FALSE(overlaps(cube, bad_sphere));
            MALLOY_CHECK_FALSE(contact(cube, bad_sphere).has_value());
            MALLOY_CHECK_FALSE(contact(bad_sphere, cube).has_value());
        }
    }

    // --- M36: Obb2 against Obb2, the 2D separating-axis test (the 2D sibling of
    //     the Box3 pair). Only the four face axes, and no edge-edge case, since
    //     in 2D an edge's separating direction is already a face normal. One
    //     predicate drives overlaps() and contact(); ground truth is
    //     hand-computed. ---
    {
        using malloy::collide::Obb2;
        using malloy::collide::contact;
        using malloy::collide::overlaps;
        using malloy::math::Vec2;

        const Vec2 half{0.5, 0.5};

        // Clear of each other along x: no overlap, no contact.
        {
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 b{Vec2{2.0, 0.0}, half, 0.0};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_FALSE(c.has_value());
        }

        // Face of a: two axis-aligned boxes overlapping 0.2 along x. Normal +x,
        // penetration 0.2, contact midway in the overlap slab. The mirror with b
        // on the negative side flips the normal to -x, pinning the sign.
        {
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 b{Vec2{0.8, 0.0}, half, 0.0};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
            MALLOY_CHECK_NEAR(c->penetration, 0.2, eps);
            MALLOY_CHECK_NEAR(c->point.x, 0.4, eps);
            MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));

            const Obb2 b_neg{Vec2{-0.8, 0.0}, half, 0.0};
            const auto c_neg = contact(a, b_neg);
            MALLOY_CHECK_TRUE(overlaps(a, b_neg) == c_neg.has_value());
            MALLOY_CHECK_TRUE(c_neg.has_value());
            MALLOY_CHECK_VEC2_NEAR(c_neg->normal, Vec2(-1.0, 0.0), eps);
            MALLOY_CHECK_NEAR(c_neg->point.x, -0.4, eps);
        }

        // Exact touch: boxes one unit apart share an edge. Touching counts;
        // nudging b a hair further apart separates them (the predicate is strict).
        {
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 touch{Vec2{1.0, 0.0}, half, 0.0};
            const auto c = contact(a, touch);
            MALLOY_CHECK_TRUE(overlaps(a, touch) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_NEAR(c->penetration, 0.0, eps);
            MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
            const Obb2 apart{Vec2{1.0 + 1e-6, 0.0}, half, 0.0};
            const auto c2 = contact(a, apart);
            MALLOY_CHECK_TRUE(overlaps(a, apart) == c2.has_value());
            MALLOY_CHECK_FALSE(c2.has_value());
        }

        // Face of b: b is turned 30 degrees and pushed along its OWN x-axis, so
        // the tightest axis is b's face normal, not a's. Normal is b's local x,
        // penetration follows from the projected radii along it. This is the
        // mirror of the FaceA branch and would break if the two were confused.
        {
            const Real angle = 0.52359877559829887308; // pi / 6
            const Vec2 bx{std::cos(angle), std::sin(angle)}; // b local x in world
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 b{bx * 1.0, half, angle};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_VEC2_NEAR(c->normal, bx, 1e-9);
            const Real expected = 0.5 * (std::cos(angle) + std::sin(angle)) + 0.5 - 1.0;
            MALLOY_CHECK_NEAR(c->penetration, expected, 1e-9);
            MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));
        }

        // Coincident centres, different sizes: every direction overlaps, so no
        // axis separates and the least-combined-width axis wins (here x, 0.8 < the
        // 0.9 along y). A fixed, repeatable fallback, the box analogue of the
        // circle pair's +x, with no divide by zero.
        {
            const Vec2 c0{1.0, -1.0};
            const Obb2 a{c0, Vec2{0.5, 0.5}, 0.0};
            const Obb2 b{c0, Vec2{0.3, 0.4}, 0.0};
            const auto c = contact(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_VEC2_NEAR(c->normal, Vec2(1.0, 0.0), eps);
            MALLOY_CHECK_NEAR(c->penetration, 0.8, eps);
        }

        // An invalid box never overlaps or contacts anything.
        {
            const Obb2 good{Vec2{}, half, 0.0};
            const Obb2 bad{Vec2{}, Vec2{-1.0, 0.5}, 0.0};
            MALLOY_CHECK_FALSE(overlaps(good, bad));
            MALLOY_CHECK_FALSE(contact(good, bad).has_value());
            MALLOY_CHECK_FALSE(overlaps(bad, good));
            MALLOY_CHECK_FALSE(contact(bad, good).has_value());
        }
    }

    // --- M37: Obb2 against a Halfplane, the 2D sibling of Box3 against a plane
    //     (M30). A box on a floor touches at up to two corners, so this returns a
    //     manifold (a list), not a single point. No separating-axis search: the
    //     plane has one normal, and the four corners tested against it are the
    //     whole story. ---
    {
        using malloy::collide::contacts;
        using malloy::collide::Halfplane;
        using malloy::collide::Obb2;
        using malloy::collide::overlaps;
        using malloy::math::Vec2;

        const Vec2 half{0.5, 0.5};
        // Floor at y = 0, normal pointing up out of the solid below it.
        const Halfplane floor{Vec2{0.0, 1.0}, 0.0};

        // Clear of the floor: no overlap, no contacts.
        {
            const Obb2 box{Vec2{0.0, 1.0}, half, 0.0};
            MALLOY_CHECK_FALSE(overlaps(box, floor));
            MALLOY_CHECK_TRUE(contacts(box, floor).empty());
        }

        // A flat square with its centre 0.3 above the floor: its two bottom
        // corners are 0.2 below the surface, its two top corners above. Two
        // contacts, each 0.2 deep, normal -y. The two bottom corners are at
        // x = +-0.5, so their x sums to zero: the centroid sits under the centre
        // of mass, which is what makes a flat drop induce no torque.
        {
            const Obb2 box{Vec2{0.0, 0.3}, half, 0.0};
            MALLOY_CHECK_TRUE(overlaps(box, floor));
            const auto hits = contacts(box, floor);
            MALLOY_CHECK_TRUE(hits.size() == 2);
            Real sum_x = 0.0;
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_NEAR(hit.penetration, 0.2, eps);
                MALLOY_CHECK_VEC2_NEAR(hit.normal, Vec2(0.0, -1.0), eps);
                MALLOY_CHECK_NEAR(hit.point.y, -0.1, eps); // midway: corner -0.2, surface 0
                sum_x += hit.point.x;
            }
            MALLOY_CHECK_NEAR(sum_x, 0.0, eps); // centroid under the centre of mass
        }

        // Tilted 45 degrees, centre at y = 0.5: only the single lowest corner
        // dips below (at y = 0.5 - sqrt(0.5)), the two side corners sit exactly on
        // y = 0.5 and the top well above. One contact, a corner landing.
        {
            const Obb2 box{Vec2{0.0, 0.5}, half, 0.78539816339744830961}; // pi/4
            const auto hits = contacts(box, floor);
            MALLOY_CHECK_TRUE(hits.size() == 1);
            MALLOY_CHECK_NEAR(hits[0].penetration, std::sqrt(0.5) - 0.5, 1e-9);
            MALLOY_CHECK_VEC2_NEAR(hits[0].normal, Vec2(0.0, -1.0), eps);
        }

        // A sloped floor whose normal has both components nonzero, so a dropped
        // axis in the corner test would show. The box straddles it.
        {
            const Vec2 n = malloy::math::normalize(Vec2{1.0, 2.0});
            const Halfplane ramp{n, 0.0};
            const Obb2 box{Vec2{}, half, 0.0}; // centred on the surface
            MALLOY_CHECK_TRUE(overlaps(box, ramp));
            const auto hits = contacts(box, ramp);
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_VEC2_NEAR(hit.normal, n * -1.0, 1e-15);
                MALLOY_CHECK_TRUE(hit.penetration > 0.0);
            }
        }

        // Invalid inputs give an empty manifold.
        {
            const Obb2 bad{Vec2{}, Vec2{-1.0, 0.5}, 0.0};
            MALLOY_CHECK_TRUE(contacts(bad, floor).empty());
            MALLOY_CHECK_FALSE(overlaps(bad, floor));
            const Obb2 good{Vec2{}, half, 0.0};
            MALLOY_CHECK_TRUE(contacts(good, Halfplane{Vec2{0.0, 2.0}, 0.0}).empty());
        }
    }

    // --- Box-box contact MANIFOLDS (review pass after M37). contact() reports a
    //     single vertex on a face corner; a rigid body that pushed there invented
    //     a lever arm for a head-on hit. contacts() clips the incident face to the
    //     reference face instead, so its centroid sits in the middle of the
    //     touch. Every case also asserts the manifold agrees with overlaps(). ---
    {
        using malloy::collide::Obb2;
        using malloy::collide::contact;
        using malloy::collide::contacts;
        using malloy::collide::overlaps;
        using malloy::math::Vec2;
        const Vec2 half{0.5, 0.5};

        // Apart: empty.
        {
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 b{Vec2{2.0, 0.0}, half, 0.0};
            MALLOY_CHECK_TRUE(contacts(a, b).empty());
            MALLOY_CHECK_FALSE(overlaps(a, b));
        }

        // Face on face, aligned: both ends of the incident edge, each 0.2 deep,
        // midway at x = 0.4, and symmetric in y so the centroid is on the line
        // of centres. The deepest equals the SAT penetration.
        {
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 b{Vec2{0.8, 0.0}, half, 0.0};
            const auto hits = contacts(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == !hits.empty());
            MALLOY_CHECK_TRUE(hits.size() == 2);
            Real sum_y = 0.0;
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_VEC2_NEAR(hit.normal, Vec2(1.0, 0.0), eps);
                MALLOY_CHECK_NEAR(hit.penetration, 0.2, eps);
                MALLOY_CHECK_NEAR(hit.point.x, 0.4, eps);
                MALLOY_CHECK_NEAR(std::abs(hit.point.y), 0.5, eps);
                sum_y += hit.point.y;
            }
            MALLOY_CHECK_NEAR(sum_y, 0.0, eps);
            MALLOY_CHECK_NEAR(hits[0].penetration, contact(a, b)->penetration, eps);
        }

        // Face on face, offset by 0.6 in y: the incident edge spans y in
        // [0.1, 1.1] and is clipped to the reference face's [-0.5, 0.5], so the
        // manifold spans the true overlap [0.1, 0.5] and its centroid is 0.3.
        {
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 b{Vec2{0.8, 0.6}, half, 0.0};
            const auto hits = contacts(a, b);
            MALLOY_CHECK_TRUE(hits.size() == 2);
            const Real lo = std::fmin(hits[0].point.y, hits[1].point.y);
            const Real hi = std::fmax(hits[0].point.y, hits[1].point.y);
            MALLOY_CHECK_NEAR(lo, 0.1, eps);
            MALLOY_CHECK_NEAR(hi, 0.5, eps);
            MALLOY_CHECK_NEAR(hits[0].penetration, 0.2, eps);
            MALLOY_CHECK_NEAR(hits[1].penetration, 0.2, eps);

            // The mirror, offset by -0.6: now the OTHER side plane does the
            // clipping, so both of the reference face's side planes are pinned.
            const Obb2 below{Vec2{0.8, -0.6}, half, 0.0};
            const auto mirrored = contacts(a, below);
            MALLOY_CHECK_TRUE(mirrored.size() == 2);
            MALLOY_CHECK_NEAR(std::fmin(mirrored[0].point.y, mirrored[1].point.y), -0.5, eps);
            MALLOY_CHECK_NEAR(std::fmax(mirrored[0].point.y, mirrored[1].point.y), -0.1, eps);
        }

        // A corner: b turned 45 degrees with its lowest corner 0.1 into a's top
        // face. The incident edge is clipped, then only the corner lies behind
        // the face: one point, at x = 0, normal +y.
        {
            const Real diag = std::sqrt(0.5);
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 b{Vec2{0.0, 0.5 + diag - 0.1}, half, 0.78539816339744830961};
            const auto hits = contacts(a, b);
            MALLOY_CHECK_TRUE(hits.size() == 1);
            MALLOY_CHECK_VEC2_NEAR(hits[0].normal, Vec2(0.0, 1.0), 1e-12);
            MALLOY_CHECK_NEAR(hits[0].penetration, 0.1, 1e-12);
            MALLOY_CHECK_NEAR(hits[0].point.x, 0.0, 1e-12);

            // Arguments swapped: the SECOND box owns the separating face (the
            // SAT's FaceB branch) and is the reference. Same corner, same depth,
            // normal reversed to point from the first box toward the second.
            const auto swapped = contacts(b, a);
            MALLOY_CHECK_TRUE(swapped.size() == 1);
            MALLOY_CHECK_VEC2_NEAR(swapped[0].normal, Vec2(0.0, -1.0), 1e-12);
            MALLOY_CHECK_NEAR(swapped[0].penetration, 0.1, 1e-12);
            MALLOY_CHECK_VEC2_NEAR(swapped[0].point, Vec2(0.0, 0.45), 1e-12); // corner 0.4, face 0.5
        }

        // A face of b wins (b turned 30 degrees, pushed along its own x), so b is
        // the reference box and a is incident. The manifold keeps contact()'s
        // normal and never reports more penetration than the SAT does.
        {
            const Real angle = 0.52359877559829887308; // pi / 6
            const Vec2 bx{std::cos(angle), std::sin(angle)};
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 b{bx * 1.0, half, angle};
            const auto single = contact(a, b);
            const auto hits = contacts(a, b);
            MALLOY_CHECK_TRUE(single.has_value());
            MALLOY_CHECK_FALSE(hits.empty());
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_VEC2_NEAR(hit.normal, single->normal, 1e-12);
                MALLOY_CHECK_TRUE(hit.penetration >= 0.0);
                MALLOY_CHECK_TRUE(hit.penetration <= single->penetration + 1e-12);
            }
        }

        // A bare touch still produces a manifold (depth zero), and an invalid
        // box produces none.
        {
            const Obb2 a{Vec2{0.0, 0.0}, half, 0.0};
            const Obb2 touch{Vec2{1.0, 0.0}, half, 0.0};
            const auto hits = contacts(a, touch);
            MALLOY_CHECK_FALSE(hits.empty());
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_NEAR(hit.penetration, 0.0, eps);
            }
            const Obb2 bad{Vec2{}, Vec2{-1.0, 0.5}, 0.0};
            MALLOY_CHECK_TRUE(contacts(a, bad).empty());
            MALLOY_CHECK_TRUE(contacts(bad, a).empty());
        }
    }
    {
        using malloy::collide::Box3;
        using malloy::collide::contact;
        using malloy::collide::contacts;
        using malloy::collide::overlaps;
        using malloy::math::from_axis_angle;
        using malloy::math::normalize;
        using malloy::math::Quat;
        using malloy::math::Vec3;
        const Vec3 half{0.5, 0.5, 0.5};
        const Quat none{};

        // Apart: empty.
        MALLOY_CHECK_TRUE(
            (contacts(Box3{Vec3{}, half, none}, Box3{Vec3{2.0, 0.0, 0.0}, half, none}).empty()));

        // Face on face, aligned: all four corners of the incident face, each 0.2
        // deep at x = 0.4, symmetric in y and z so the centroid is on the line of
        // centres.
        {
            const Box3 a{Vec3{}, half, none};
            const Box3 b{Vec3{0.8, 0.0, 0.0}, half, none};
            const auto hits = contacts(a, b);
            MALLOY_CHECK_TRUE(overlaps(a, b) == !hits.empty());
            MALLOY_CHECK_TRUE(hits.size() == 4);
            Vec3 sum{};
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(hit.normal, Vec3{1.0, 0.0, 0.0}, eps));
                MALLOY_CHECK_NEAR(hit.penetration, 0.2, eps);
                MALLOY_CHECK_NEAR(hit.point.x, 0.4, eps);
                sum += hit.point;
            }
            MALLOY_CHECK_NEAR(sum.y, 0.0, eps);
            MALLOY_CHECK_NEAR(sum.z, 0.0, eps);
        }

        // Face on face, offset by 0.6 in BOTH y and z: the incident face is
        // clipped against all four side planes of the reference face, so the
        // manifold is the overlap square y, z in [0.1, 0.5] (centroid 0.3, 0.3).
        // Offsetting both axes means a missing clip on either one shows.
        {
            const Box3 a{Vec3{}, half, none};
            const Box3 b{Vec3{0.8, 0.6, 0.6}, half, none};
            const auto hits = contacts(a, b);
            MALLOY_CHECK_TRUE(hits.size() == 4);
            Vec3 sum{};
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_TRUE(hit.point.y >= 0.1 - eps && hit.point.y <= 0.5 + eps);
                MALLOY_CHECK_TRUE(hit.point.z >= 0.1 - eps && hit.point.z <= 0.5 + eps);
                MALLOY_CHECK_NEAR(hit.penetration, 0.2, eps);
                sum += hit.point;
            }
            MALLOY_CHECK_NEAR(sum.y / 4.0, 0.3, eps);
            MALLOY_CHECK_NEAR(sum.z / 4.0, 0.3, eps);

            // The mirror, offset by -0.6 in y and z: the other two side planes do
            // the clipping, so all four side planes of the reference face are
            // pinned.
            const Box3 below{Vec3{0.8, -0.6, -0.6}, half, none};
            const auto mirrored = contacts(a, below);
            MALLOY_CHECK_TRUE(mirrored.size() == 4);
            Vec3 msum{};
            for (const auto& hit : mirrored)
            {
                MALLOY_CHECK_TRUE(hit.point.y >= -0.5 - eps && hit.point.y <= -0.1 + eps);
                MALLOY_CHECK_TRUE(hit.point.z >= -0.5 - eps && hit.point.z <= -0.1 + eps);
                msum += hit.point;
            }
            MALLOY_CHECK_NEAR(msum.y / 4.0, -0.3, eps);
            MALLOY_CHECK_NEAR(msum.z / 4.0, -0.3, eps);
        }

        // An edge resting on a face: b turned 45 degrees about z with its lowest
        // edge (parallel to z) 0.1 into a's top face. Two points, the two ends of
        // that edge clipped to the face, at x = 0 and z = +-0.5.
        {
            const Real diag = std::sqrt(0.5);
            const Box3 a{Vec3{}, half, none};
            const Box3 b{Vec3{0.0, 0.5 + diag - 0.1, 0.0}, half,
                         from_axis_angle(Vec3{0.0, 0.0, 1.0}, 0.78539816339744830961)};
            const auto hits = contacts(a, b);
            MALLOY_CHECK_TRUE(hits.size() == 2);
            for (const auto& hit : hits)
            {
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(hit.normal, Vec3{0.0, 1.0, 0.0}, 1e-12));
                MALLOY_CHECK_NEAR(hit.penetration, 0.1, 1e-12);
                MALLOY_CHECK_NEAR(hit.point.x, 0.0, 1e-12);
                MALLOY_CHECK_NEAR(std::abs(hit.point.z), 0.5, 1e-12);
            }

            // The same pair with the arguments swapped: now the SECOND box owns
            // the separating face (the SAT's FaceB branch), so it is the
            // reference. The same two points come back, with the normal
            // reversed to keep pointing from the first box toward the second.
            const auto swapped = contacts(b, a);
            MALLOY_CHECK_TRUE(swapped.size() == 2);
            for (const auto& hit : swapped)
            {
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(hit.normal, Vec3{0.0, -1.0, 0.0}, 1e-12));
                MALLOY_CHECK_NEAR(hit.penetration, 0.1, 1e-12);
                MALLOY_CHECK_NEAR(hit.point.x, 0.0, 1e-12);
                MALLOY_CHECK_NEAR(hit.point.y, 0.45, 1e-12); // midway: edge 0.4, face 0.5
            }
        }

        // Edge against edge (the M33 case): genuinely a point, so the manifold is
        // exactly contact()'s single contact.
        {
            const Vec3 diag = normalize(Vec3{1.0, 1.0, 0.0});
            const Box3 a{Vec3{}, half, none};
            const Box3 b{diag * 1.2, half,
                         from_axis_angle(Vec3{1.0, 1.0, 0.0}, 0.78539816339744830961)};
            const auto single = contact(a, b);
            const auto hits = contacts(a, b);
            MALLOY_CHECK_TRUE(single.has_value());
            MALLOY_CHECK_TRUE(hits.size() == 1);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(hits[0].point, single->point, 0.0));
            MALLOY_CHECK_NEAR(hits[0].penetration, single->penetration, 0.0);
        }

        // Invalid boxes give no manifold.
        MALLOY_CHECK_TRUE(
            (contacts(Box3{Vec3{}, Vec3{-1.0, 0.5, 0.5}, none}, Box3{Vec3{}, half, none}).empty()));
    }

    // --- Issues #22 and #23: box-box manifolds under GENERIC orientations. The
    //     tests above use identity or single-axis rotations, where every cross
    //     product and projection is exact, and so could not see either defect. ---
    {
        using malloy::collide::Box3;
        using malloy::collide::Obb2;
        using malloy::collide::contacts;
        using malloy::math::from_axis_angle;
        using malloy::math::normalize;
        using malloy::math::Quat;
        using malloy::math::rotate;
        using malloy::math::Vec3;

        // #22: face on face with a twist about the shared normal. An edge of
        // each box lies in the face plane, and their cross product IS the face
        // normal up to rounding, so the edge axis ties the face axis; it must
        // not win. The manifold is a face patch (at least three points) with
        // the face normal, and its centroid sits under b's centre.
        {
            const Quat r = from_axis_angle(Vec3{0.3, -0.7, 0.2}, 0.9);
            const Vec3 up = rotate(r, Vec3{0.0, 0.0, 1.0});
            for (int degrees = 1; degrees < 90; ++degrees)
            {
                const Real twist = degrees * 3.14159265358979323846 / 180.0;
                const Box3 a{Vec3{}, Vec3{1.0, 0.7, 0.5}, r};
                const Vec3 centre = rotate(r, Vec3{0.1, 0.05, 0.0}) + up * (0.8 - 0.004);
                const Box3 b{centre, Vec3{0.4, 0.3, 0.3},
                             normalize(r * from_axis_angle(Vec3{0.0, 0.0, 1.0}, twist))};
                const auto hits = contacts(a, b);
                MALLOY_CHECK_TRUE(hits.size() >= 3);
                Vec3 sum{};
                for (const auto& hit : hits)
                {
                    MALLOY_CHECK_TRUE(malloy::math::approx_equal(hit.normal, up, 1e-9));
                    MALLOY_CHECK_NEAR(hit.penetration, 0.004, 1e-9);
                    sum += hit.point;
                }
                const Vec3 arm = sum / static_cast<Real>(hits.size()) - centre;
                const Vec3 lateral = arm - up * malloy::math::dot(arm, up);
                MALLOY_CHECK_NEAR(malloy::math::length(lateral), 0.0, 1e-9);
                // The single contact agrees on the axis.
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(contact(a, b)->normal, up, 1e-9));
            }
        }

        // #23: a bare touch under a generic rotation. Two equal cubes, aligned
        // with each other and touching face to face to rounding, must give the
        // full four-point face, not one corner and not a fifth vertex that
        // clipping put partway along an edge lying in a side plane. At this
        // rotation clipping did both: a 5-point manifold, then a lopsided one.
        {
            for (const Real angle : {0.114, 0.702, 1.479})
            {
                const Quat r = from_axis_angle(Vec3{0.3, -0.7, 0.2}, angle);
                const Vec3 u = rotate(r, Vec3{1.0, 0.0, 0.0});
                for (const Real gap : {0.0, 1e-15, -1e-15, -0.004})
                {
                    const Box3 a{u * -0.5, Vec3{0.5, 0.5, 0.5}, r};
                    const Box3 b{u * (0.5 + gap), Vec3{0.5, 0.5, 0.5}, r};
                    const auto hits = contacts(a, b);
                    if (gap > 0.0 && hits.empty())
                    {
                        continue; // a hair apart is allowed to be apart
                    }
                    MALLOY_CHECK_TRUE(hits.size() == 4);
                    Vec3 sum{};
                    for (const auto& hit : hits)
                    {
                        MALLOY_CHECK_TRUE(hit.penetration >= 0.0);
                        sum += hit.point;
                    }
                    // The centroid lies on the line of centres: no lever arm.
                    const Vec3 mid = sum / 4.0 - (a.center + b.center) * 0.5;
                    MALLOY_CHECK_NEAR(malloy::math::length(mid - u * malloy::math::dot(mid, u)),
                                      0.0, 1e-9);
                }
            }
            // 2D: the same for two squares turned to a generic angle, both
            // endpoints of the touching face kept, symmetric about the line of
            // centres.
            for (int k = 0; k < 200; ++k)
            {
                const Real angle = 0.1 + 0.007 * k;
                const Vec2 u{std::cos(angle), std::sin(angle)};
                const Obb2 a{u * -0.5, Vec2{0.5, 0.5}, angle};
                const Obb2 b{u * 0.5, Vec2{0.5, 0.5}, angle};
                const auto hits = contacts(a, b);
                MALLOY_CHECK_TRUE(hits.size() == 2);
                MALLOY_CHECK_TRUE(hits[0].penetration >= 0.0 && hits[1].penetration >= 0.0);
                const Vec2 mid = (hits[0].point + hits[1].point) * 0.5;
                MALLOY_CHECK_NEAR(mid.x * -u.y + mid.y * u.x, 0.0, 1e-9);
            }
        }

        // A property over many random box pairs (a fixed, fully specified
        // generator, so the run repeats): whenever two boxes overlap, the single
        // contact lies within its penetration of BOTH boxes, and moving b along
        // the normal by the penetration leaves them no deeper than a touch. It
        // guards the face preference of #22 against growing too generous: a 5
        // per cent margin (ODE's) broke the first in three of 85000 such pairs,
        // the face winning where the contact was really edge on edge.
        {
            std::uint64_t state = 0x9E3779B97F4A7C15ULL;
            const auto next = [&state]() {
                state += 0x9E3779B97F4A7C15ULL; // splitmix64
                std::uint64_t z = state;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
                z ^= z >> 31;
                return static_cast<Real>(z >> 11) * (1.0 / 9007199254740992.0); // [0, 1)
            };
            const auto signed_unit = [&next]() { return 2.0 * next() - 1.0; };
            const auto distance_outside = [](const Box3& box, const Vec3& p) {
                const Vec3 d = rotate(malloy::math::conjugate(box.orientation), p - box.center);
                const Real dx = std::fmax(0.0, std::abs(d.x) - box.half_extents.x);
                const Real dy = std::fmax(0.0, std::abs(d.y) - box.half_extents.y);
                const Real dz = std::fmax(0.0, std::abs(d.z) - box.half_extents.z);
                return std::sqrt(dx * dx + dy * dy + dz * dz);
            };
            int overlapping = 0;
            for (int it = 0; it < 20000; ++it)
            {
                const Box3 a{Vec3{signed_unit(), signed_unit(), signed_unit()},
                             Vec3{0.1 + 0.9 * next(), 0.1 + 0.9 * next(), 0.1 + 0.9 * next()},
                             normalize(Quat{signed_unit(), Vec3{signed_unit(), signed_unit(),
                                                                signed_unit()}})};
                const Box3 b{Vec3{signed_unit(), signed_unit(), signed_unit()},
                             Vec3{0.1 + 0.9 * next(), 0.1 + 0.9 * next(), 0.1 + 0.9 * next()},
                             normalize(Quat{signed_unit(), Vec3{signed_unit(), signed_unit(),
                                                                signed_unit()}})};
                const auto hit = contact(a, b);
                if (!hit)
                {
                    continue;
                }
                ++overlapping;
                MALLOY_CHECK_TRUE(distance_outside(a, hit->point) <= hit->penetration + 1e-9);
                MALLOY_CHECK_TRUE(distance_outside(b, hit->point) <= hit->penetration + 1e-9);
                const Box3 moved{b.center + hit->normal * hit->penetration, b.half_extents,
                                 b.orientation};
                const auto after = contact(a, moved);
                MALLOY_CHECK_TRUE(!after || after->penetration <= 1e-9);
                MALLOY_CHECK_FALSE(contacts(a, b).empty());
            }
            MALLOY_CHECK_TRUE(overlapping > 5000); // the property was actually exercised
        }
    }

    // --- M40: the speculative (margin-aware) contact queries the iterative
    //     contact solver resolves. A pair apart by a gap within the margin is
    //     reported at penetration -gap, its point midway in the gap; with
    //     margin 0 the plain query's answer comes back; a bad margin is refused
    //     like a bad shape. ---
    {
        using malloy::collide::Box3;
        using malloy::collide::Contact3;
        using malloy::collide::contacts;
        using malloy::collide::Obb2;
        using malloy::collide::Plane3;
        using malloy::collide::Sphere;
        using malloy::math::from_axis_angle;
        using malloy::math::normalize;
        using malloy::math::Quat;
        using malloy::math::rotate;
        using malloy::math::Vec3;

        // A fixed, fully specified generator, so every run repeats.
        std::uint64_t state = 0x243F6A8885A308D3ULL;
        const auto next = [&state]() {
            state += 0x9E3779B97F4A7C15ULL; // splitmix64
            std::uint64_t z = state;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            z ^= z >> 31;
            return static_cast<Real>(z >> 11) * (1.0 / 9007199254740992.0); // [0, 1)
        };
        const auto signed_unit = [&next]() { return 2.0 * next() - 1.0; };
        const auto same2 = [](const Contact& x, const Contact& y) {
            return x.normal.x == y.normal.x && x.normal.y == y.normal.y &&
                   x.penetration == y.penetration && x.point.x == y.point.x &&
                   x.point.y == y.point.y;
        };
        const auto same3 = [](const Contact3& x, const Contact3& y) {
            return x.normal.x == y.normal.x && x.normal.y == y.normal.y &&
                   x.normal.z == y.normal.z && x.penetration == y.penetration &&
                   x.point.x == y.point.x && x.point.y == y.point.y && x.point.z == y.point.z;
        };

        // 2D, a circle 0.3 above a floor, margin 0.5: reported, depth -0.3,
        // point halfway between the circle's bottom (y 0.3) and the floor.
        {
            const Halfplane floor{Vec2{0.0, 1.0}, 0.0};
            const auto above = contact(Circle{Vec2{2.0, 1.3}, 1.0}, floor, 0.5);
            MALLOY_CHECK_TRUE(above.has_value());
            MALLOY_CHECK_NEAR(above->penetration, -0.3, 1e-15);
            MALLOY_CHECK_VEC2_NEAR(above->normal, (Vec2{0.0, -1.0}), 0.0);
            MALLOY_CHECK_VEC2_NEAR(above->point, (Vec2{2.0, 0.15}), 1e-15);
            MALLOY_CHECK_FALSE(contact(Circle{Vec2{2.0, 1.3}, 1.0}, floor, 0.25).has_value());
            MALLOY_CHECK_FALSE(contact(Circle{Vec2{2.0, 1.3}, 1.0}, floor).has_value());
            // Two circles 0.2 apart: reported at margin 0.25, not at 0.1.
            const auto pair = contact(Circle{Vec2{}, 1.0}, Circle{Vec2{2.7, 0.0}, 1.5}, 0.25);
            MALLOY_CHECK_TRUE(pair.has_value());
            MALLOY_CHECK_NEAR(pair->penetration, -0.2, 1e-15);
            MALLOY_CHECK_VEC2_NEAR(pair->point, (Vec2{1.1, 0.0}), 1e-15);
            MALLOY_CHECK_FALSE(
                contact(Circle{Vec2{}, 1.0}, Circle{Vec2{2.7, 0.0}, 1.5}, 0.1).has_value());
            // A square 0.01 clear of the floor: margin 0.02 keeps its two lower
            // corners at depth -0.01, margin 0.005 keeps neither.
            const Obb2 box{Vec2{0.0, 0.51}, Vec2{0.5, 0.5}, 0.0};
            const auto corners = contacts(box, floor, 0.02);
            MALLOY_CHECK_EQ(corners.size(), std::size_t{2});
            MALLOY_CHECK_NEAR(corners[0].penetration, -0.01, 1e-15);
            MALLOY_CHECK_NEAR(corners[0].point.y, 0.005, 1e-15);
            MALLOY_CHECK_TRUE(contacts(box, floor, 0.005).empty());
            // Two squares 0.004 apart face to face: both ends of the face, at
            // depth -0.004, with the plain query reporting nothing.
            const Obb2 left{Vec2{-0.5, 0.0}, Vec2{0.5, 0.5}, 0.0};
            const Obb2 right{Vec2{0.504, 0.1}, Vec2{0.5, 0.5}, 0.0};
            const auto face = contacts(left, right, 0.01);
            MALLOY_CHECK_EQ(face.size(), std::size_t{2});
            for (const Contact& hit : face)
            {
                MALLOY_CHECK_NEAR(hit.penetration, -0.004, 1e-15);
                MALLOY_CHECK_NEAR(hit.point.x, 0.002, 1e-15);
            }
            MALLOY_CHECK_TRUE(contacts(left, right).empty());
            MALLOY_CHECK_TRUE(contacts(left, right, 0.003).empty());
        }

        // 3D, the same four facts for spheres, a box on a plane, a box against a
        // sphere and two boxes.
        {
            const Plane3 floor{Vec3{0.0, 0.0, 1.0}, 0.0};
            const auto sphere_hit = contact(Sphere{Vec3{1.0, 2.0, 1.3}, 1.0}, floor, 0.5);
            MALLOY_CHECK_TRUE(sphere_hit.has_value());
            MALLOY_CHECK_NEAR(sphere_hit->penetration, -0.3, 1e-15);
            MALLOY_CHECK_NEAR(sphere_hit->point.z, 0.15, 1e-15);
            MALLOY_CHECK_FALSE(contact(Sphere{Vec3{1.0, 2.0, 1.3}, 1.0}, floor, 0.25).has_value());
            const auto spheres =
                contact(Sphere{Vec3{}, 1.0}, Sphere{Vec3{0.0, 2.7, 0.0}, 1.5}, 0.25);
            MALLOY_CHECK_TRUE(spheres.has_value());
            MALLOY_CHECK_NEAR(spheres->penetration, -0.2, 1e-15);
            MALLOY_CHECK_NEAR(spheres->point.y, 1.1, 1e-15);
            const Box3 cube{Vec3{0.0, 0.0, 0.51}, Vec3{0.5, 0.5, 0.5}, Quat{}};
            const auto corners = contacts(cube, floor, 0.02);
            MALLOY_CHECK_EQ(corners.size(), std::size_t{4});
            for (const Contact3& hit : corners)
            {
                MALLOY_CHECK_NEAR(hit.penetration, -0.01, 1e-15);
            }
            MALLOY_CHECK_TRUE(contacts(cube, floor, 0.005).empty());
            // The cube's top is at 1.01; a ball of radius 0.3 at 1.34 is 0.03
            // clear of it, the point halfway down the gap at 1.025.
            const auto ball = contact(cube, Sphere{Vec3{0.0, 0.0, 1.34}, 0.3}, 0.05);
            MALLOY_CHECK_TRUE(ball.has_value());
            MALLOY_CHECK_NEAR(ball->penetration, -0.03, 1e-15);
            MALLOY_CHECK_NEAR(ball->point.z, 1.025, 1e-15);
            MALLOY_CHECK_FALSE(
                contact(cube, Sphere{Vec3{0.0, 0.0, 1.41}, 0.3}, 0.05).has_value());
            const Box3 above{Vec3{0.1, 0.0, 1.514}, Vec3{0.5, 0.5, 0.5}, Quat{}};
            const auto face = contacts(cube, above, 0.01);
            MALLOY_CHECK_EQ(face.size(), std::size_t{4});
            for (const Contact3& hit : face)
            {
                MALLOY_CHECK_NEAR(hit.penetration, -0.004, 1e-15);
                MALLOY_CHECK_NEAR(hit.point.z, 1.012, 1e-15);
            }
            MALLOY_CHECK_TRUE(contacts(cube, above).empty());
            MALLOY_CHECK_TRUE(contacts(cube, above, 0.003).empty());
        }

        // A margin that is negative or not finite is refused like a bad shape.
        for (const Real bad : {-1e-9, nan, inf})
        {
            const Circle c{Vec2{}, 1.0};
            MALLOY_CHECK_FALSE(contact(c, Halfplane{Vec2{0.0, 1.0}, 0.0}, bad).has_value());
            MALLOY_CHECK_FALSE(contact(c, c, bad).has_value());
            const Obb2 square{Vec2{}, Vec2{0.5, 0.5}, 0.0};
            MALLOY_CHECK_TRUE(contacts(square, Halfplane{Vec2{0.0, 1.0}, 0.0}, bad).empty());
            MALLOY_CHECK_TRUE(contacts(square, square, bad).empty());
            const Sphere s{Vec3{}, 1.0};
            const Box3 cube{Vec3{}, Vec3{0.5, 0.5, 0.5}, Quat{}};
            MALLOY_CHECK_FALSE(contact(s, Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}, bad).has_value());
            MALLOY_CHECK_FALSE(contact(s, s, bad).has_value());
            MALLOY_CHECK_TRUE(contacts(cube, Plane3{Vec3{0.0, 0.0, 1.0}, 0.0}, bad).empty());
            MALLOY_CHECK_TRUE(contacts(cube, cube, bad).empty());
            MALLOY_CHECK_FALSE(contact(cube, s, bad).has_value());
        }

        // Margin 0 gives back exactly what the plain query gives, bit for bit,
        // over random configurations (the plain queries now delegate to these,
        // so this pins that the margin adds nothing at 0). The box/box manifold
        // is signed where the plain one clamps at 0, so its depths are compared
        // after the same clamp; in 3D only where a face wins outright (the face
        // preference is the speculative query's own, checked below).
        for (int it = 0; it < 4000; ++it)
        {
            const Vec2 p{signed_unit(), signed_unit()};
            const Halfplane plane{normalize(Vec2{signed_unit(), signed_unit() + 2.0}),
                                  0.3 * signed_unit()};
            const Circle c1{p, 0.2 + next()};
            const Circle c2{Vec2{signed_unit(), signed_unit()}, 0.2 + next()};
            const auto plain_cp = contact(c1, plane);
            const auto zero_cp = contact(c1, plane, 0.0);
            MALLOY_CHECK_EQ(plain_cp.has_value(), zero_cp.has_value());
            MALLOY_CHECK_TRUE(!plain_cp || same2(*plain_cp, *zero_cp));
            const auto plain_cc = contact(c1, c2);
            const auto zero_cc = contact(c1, c2, 0.0);
            MALLOY_CHECK_EQ(plain_cc.has_value(), zero_cc.has_value());
            MALLOY_CHECK_TRUE(!plain_cc || same2(*plain_cc, *zero_cc));
            const Obb2 b1{p, Vec2{0.1 + next(), 0.1 + next()}, 3.0 * signed_unit()};
            const Obb2 b2{Vec2{signed_unit(), signed_unit()}, Vec2{0.1 + next(), 0.1 + next()},
                          3.0 * signed_unit()};
            const auto plain_bp = contacts(b1, plane);
            const auto zero_bp = contacts(b1, plane, 0.0);
            MALLOY_CHECK_EQ(plain_bp.size(), zero_bp.size());
            for (std::size_t k = 0; k < plain_bp.size() && k < zero_bp.size(); ++k)
            {
                MALLOY_CHECK_TRUE(same2(plain_bp[k], zero_bp[k]));
            }
            const auto plain_bb = contacts(b1, b2);
            auto zero_bb = contacts(b1, b2, 0.0);
            MALLOY_CHECK_EQ(plain_bb.size(), zero_bb.size());
            for (std::size_t k = 0; k < plain_bb.size() && k < zero_bb.size(); ++k)
            {
                zero_bb[k].penetration = std::fmax(0.0, zero_bb[k].penetration);
                MALLOY_CHECK_TRUE(same2(plain_bb[k], zero_bb[k]));
            }

            const Vec3 q{signed_unit(), signed_unit(), signed_unit()};
            const Plane3 plane3{
                normalize(Vec3{signed_unit(), signed_unit(), signed_unit() + 2.0}), 0.3 * signed_unit()};
            const Sphere s1{q, 0.2 + next()};
            const Sphere s2{Vec3{signed_unit(), signed_unit(), signed_unit()}, 0.2 + next()};
            const auto plain_sp = contact(s1, plane3);
            const auto zero_sp = contact(s1, plane3, 0.0);
            MALLOY_CHECK_EQ(plain_sp.has_value(), zero_sp.has_value());
            MALLOY_CHECK_TRUE(!plain_sp || same3(*plain_sp, *zero_sp));
            const auto plain_ss = contact(s1, s2);
            const auto zero_ss = contact(s1, s2, 0.0);
            MALLOY_CHECK_EQ(plain_ss.has_value(), zero_ss.has_value());
            MALLOY_CHECK_TRUE(!plain_ss || same3(*plain_ss, *zero_ss));
            const Box3 x1{q, Vec3{0.1 + next(), 0.1 + next(), 0.1 + next()},
                          normalize(Quat{signed_unit(), Vec3{signed_unit(), signed_unit(),
                                                             signed_unit()}})};
            const Box3 x2{Vec3{signed_unit(), signed_unit(), signed_unit()},
                          Vec3{0.1 + next(), 0.1 + next(), 0.1 + next()},
                          normalize(Quat{signed_unit(), Vec3{signed_unit(), signed_unit(),
                                                             signed_unit()}})};
            const auto plain_xp = contacts(x1, plane3);
            const auto zero_xp = contacts(x1, plane3, 0.0);
            MALLOY_CHECK_EQ(plain_xp.size(), zero_xp.size());
            for (std::size_t k = 0; k < plain_xp.size() && k < zero_xp.size(); ++k)
            {
                MALLOY_CHECK_TRUE(same3(plain_xp[k], zero_xp[k]));
            }
            const auto plain_xs = contact(x1, s2);
            const auto zero_xs = contact(x1, s2, 0.0);
            MALLOY_CHECK_EQ(plain_xs.has_value(), zero_xs.has_value());
            MALLOY_CHECK_TRUE(!plain_xs || same3(*plain_xs, *zero_xs));
            const auto plain_xx = contacts(x1, x2);
            auto zero_xx = contacts(x1, x2, 0.0);
            if (plain_xx.size() >= 2) // a face won outright: same manifold
            {
                MALLOY_CHECK_EQ(plain_xx.size(), zero_xx.size());
                for (std::size_t k = 0; k < plain_xx.size() && k < zero_xx.size(); ++k)
                {
                    zero_xx[k].penetration = std::fmax(0.0, zero_xx[k].penetration);
                    MALLOY_CHECK_TRUE(same3(plain_xx[k], zero_xx[k]));
                }
            }
            MALLOY_CHECK_EQ(plain_xx.empty(), zero_xx.empty());
        }

        // The face preference (kFacePreference), and its proof. Issue #22's fix
        // kept a 5 per cent face margin OUT of the single-point contact, because
        // there a face winning where the contact is really edge on edge put the
        // point outside the boxes (3 of 85000 random pairs). The speculative
        // manifold applies it only where the face manifold keeps a point, and a
        // kept point lies on the incident box and inside the reference box's
        // side slabs. The property: every point of every manifold lies within
        // |penetration| (plus rounding, 1e-9 of the size) of BOTH boxes, over
        // 100000 configurations scaled across twelve orders of magnitude, half
        // fully random and half near-flat stacks (twisted, tilted up to 0.08
        // rad, overhanging, a touch to 0.02 deep), at margin 0 and at the
        // solver's 1e-4 of the size. The preference must actually fire (it
        // turns a plain one-point edge contact into a face patch) and never
        // breaks the property. Research run: 400000 configurations, 89338
        // firings, 0 violations, the worst point at 0.5 of its bound (a
        // midpoint, by construction).
        {
            const auto distance_outside = [](const Box3& box, const Vec3& p) {
                const Vec3 d = rotate(malloy::math::conjugate(box.orientation), p - box.center);
                const Real dx = std::fmax(0.0, std::abs(d.x) - box.half_extents.x);
                const Real dy = std::fmax(0.0, std::abs(d.y) - box.half_extents.y);
                const Real dz = std::fmax(0.0, std::abs(d.z) - box.half_extents.z);
                return std::sqrt(dx * dx + dy * dy + dz * dz);
            };
            int manifolds = 0;
            int fired = 0;
            for (int it = 0; it < 100000; ++it)
            {
                const Real scale = std::pow(10.0, -6.0 + 12.0 * next());
                Box3 a;
                Box3 b;
                if (it % 2 == 0)
                {
                    a = Box3{Vec3{signed_unit(), signed_unit(), signed_unit()} * scale,
                             Vec3{0.1 + 0.9 * next(), 0.1 + 0.9 * next(), 0.1 + 0.9 * next()} *
                                 scale,
                             normalize(Quat{signed_unit(), Vec3{signed_unit(), signed_unit(),
                                                                signed_unit()}})};
                    b = Box3{Vec3{signed_unit(), signed_unit(), signed_unit()} * scale,
                             Vec3{0.1 + 0.9 * next(), 0.1 + 0.9 * next(), 0.1 + 0.9 * next()} *
                                 scale,
                             normalize(Quat{signed_unit(), Vec3{signed_unit(), signed_unit(),
                                                                signed_unit()}})};
                }
                else
                {
                    const Quat r = normalize(
                        Quat{signed_unit(), Vec3{signed_unit(), signed_unit(), signed_unit()}});
                    const Vec3 ha{0.2 + 0.8 * next(), 0.2 + 0.8 * next(), 0.2 + 0.8 * next()};
                    const Vec3 hb{0.2 + 0.8 * next(), 0.2 + 0.8 * next(), 0.2 + 0.8 * next()};
                    const Real tilt = 0.08 * signed_unit();
                    const Vec3 tilt_axis{signed_unit(), signed_unit(), 0.0};
                    const Real twist = 3.14159 * signed_unit();
                    const Quat q = normalize(r * from_axis_angle(tilt_axis, tilt) *
                                             from_axis_angle(Vec3{0.0, 0.0, 1.0}, twist));
                    const Vec3 up = rotate(r, Vec3{0.0, 0.0, 1.0});
                    const Vec3 lateral =
                        rotate(r, Vec3{signed_unit() * ha.x, signed_unit() * ha.y, 0.0});
                    const Real depth = 0.02 * signed_unit();
                    a = Box3{Vec3{}, ha * scale, r};
                    b = Box3{(lateral + up * (ha.z + hb.z - depth)) * scale, hb * scale, q};
                }
                const Real size =
                    std::fmax(std::fmax(std::fmax(a.half_extents.x, a.half_extents.y),
                                        a.half_extents.z),
                              std::fmax(std::fmax(b.half_extents.x, b.half_extents.y),
                                        b.half_extents.z));
                for (const Real margin : {0.0, 1e-4 * size})
                {
                    const auto hits = contacts(a, b, margin);
                    if (hits.empty())
                    {
                        continue;
                    }
                    ++manifolds;
                    if (contacts(a, b).size() == 1 && hits.size() > 1)
                    {
                        ++fired;
                    }
                    for (const Contact3& hit : hits)
                    {
                        const Real bound = std::abs(hit.penetration) + 1e-9 * size;
                        MALLOY_CHECK_TRUE(distance_outside(a, hit.point) <= bound);
                        MALLOY_CHECK_TRUE(distance_outside(b, hit.point) <= bound);
                    }
                }
            }
            MALLOY_CHECK_TRUE(manifolds > 100000); // exercised at both margins
            MALLOY_CHECK_TRUE(fired > 10000);      // the preference really fired
        }
    }

    std::cout << "malloy_collide_tests passed\n";
    return 0;
}
