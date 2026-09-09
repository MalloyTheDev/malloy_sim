#include <malloy/collide/collide.hpp>

#include <malloy/math/math.hpp>
#include <test_check.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>

using malloy::collide::Aabb;
using malloy::collide::Circle;
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

    std::cout << "malloy_collide_tests passed\n";
    return 0;
}
