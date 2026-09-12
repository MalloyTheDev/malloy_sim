#include <malloy/collide/collide.hpp>

#include <malloy/math/math.hpp>
#include <test_check.hpp>

#include <cmath>
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
        {
            const Sphere s{Vec3{0.2, 0.0, 0.0}, 0.3};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.6, eps);
            MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));
        }

        // The mirror inside case: centre at (-0.2, 0, 0) is nearest the -x face,
        // so the push is -x, which pins the SIGN of the deep-contact normal (a
        // fixed +x would pass the case above but fail here).
        {
            const Sphere s{Vec3{-0.2, 0.0, 0.0}, 0.3};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(overlaps(cube, s) == c.has_value());
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{-1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.6, eps);
        }

        // Sphere centre exactly at the box centre: every face is equidistant, so
        // the tie breaks to axis 0 with a +x normal (a fixed, repeatable choice),
        // penetration radius + half-extent, no divide by zero.
        {
            const Sphere s{Vec3{0.0, 0.0, 0.0}, 0.3};
            const auto c = contact(cube, s);
            MALLOY_CHECK_TRUE(c.has_value());
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(c->normal, Vec3{1.0, 0.0, 0.0}, eps));
            MALLOY_CHECK_NEAR(c->penetration, 0.8, eps);
            MALLOY_CHECK_TRUE(malloy::math::is_finite(c->point));
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

    std::cout << "malloy_collide_tests passed\n";
    return 0;
}
