#include <malloy/charges/charges.hpp>

#include <test_check.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using malloy::charges::Charge3DSettings;
using malloy::charges::Charge3DWorld;
using malloy::charges::ChargedParticle2D;
using malloy::charges::ChargedParticle3D;
using malloy::charges::ChargeSettings;
using malloy::charges::ChargeWorld;
using malloy::charges::total_energy;
using malloy::charges::total_energy3d;
using malloy::charges::total_kinetic_energy;
using malloy::charges::total_kinetic_energy3d;
using malloy::charges::total_momentum;
using malloy::charges::total_momentum3d;
using malloy::charges::total_potential_energy;
using malloy::charges::total_potential_energy3d;
using malloy::math::Real;
using malloy::math::Vec2;
using malloy::math::Vec3;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;

namespace
{
// A timestep that divides a cyclotron period exactly, so a full turn is a whole
// number of steps and "did it come back" is a question with an exact answer.
constexpr int steps_per_turn = 1000;
const Real two_pi = Real{6.283185307179586476925286766559};
const Real turn_dt = two_pi / static_cast<Real>(steps_per_turn);

ChargedParticle2D charge_at(Vec2 position, Vec2 velocity, Real mass, Real charge)
{
    ChargedParticle2D p;
    p.position = position;
    p.velocity = velocity;
    p.mass = mass;
    p.charge = charge;
    return p;
}
} // namespace

int main()
{
    const Real eps = 1e-12;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    // --- Validation. Mass must be positive; charge must not. ---
    {
        MALLOY_CHECK_TRUE(charge_at(Vec2{1.0, 2.0}, Vec2{3.0, 4.0}, 2.5, -7.0).is_valid());
        MALLOY_CHECK_TRUE(charge_at(Vec2{}, Vec2{}, 1.0, 0.0).is_valid());  // neutral
        MALLOY_CHECK_TRUE(charge_at(Vec2{}, Vec2{}, 1.0, -3.0).is_valid()); // negative

        MALLOY_CHECK_FALSE(charge_at(Vec2{}, Vec2{}, 0.0, 1.0).is_valid());
        MALLOY_CHECK_FALSE(charge_at(Vec2{}, Vec2{}, -1.0, 1.0).is_valid());
        MALLOY_CHECK_FALSE(charge_at(Vec2{}, Vec2{}, nan, 1.0).is_valid());
        MALLOY_CHECK_FALSE(charge_at(Vec2{}, Vec2{}, 1.0, nan).is_valid());
        MALLOY_CHECK_FALSE(charge_at(Vec2{}, Vec2{}, 1.0, inf).is_valid());
        MALLOY_CHECK_FALSE(charge_at(Vec2{nan, 0.0}, Vec2{}, 1.0, 1.0).is_valid());
        MALLOY_CHECK_FALSE(charge_at(Vec2{}, Vec2{0.0, inf}, 1.0, 1.0).is_valid());

        MALLOY_CHECK_TRUE(ChargeSettings{}.is_valid());
        ChargeSettings bad_k;
        bad_k.k = nan;
        MALLOY_CHECK_FALSE(bad_k.is_valid());
        ChargeSettings bad_e;
        bad_e.electric = Vec2{inf, 0.0};
        MALLOY_CHECK_FALSE(bad_e.is_valid());
        ChargeSettings bad_b;
        bad_b.magnetic = nan;
        MALLOY_CHECK_FALSE(bad_b.is_valid());
        ChargeSettings bad_soft;
        bad_soft.softening = -1.0;
        MALLOY_CHECK_FALSE(bad_soft.is_valid());

        ChargeWorld invalid{SimulationSettings{0.001}, bad_b,
                            {charge_at(Vec2{}, Vec2{}, 1.0, 1.0)}};
        MALLOY_CHECK_TRUE(invalid.validate() == StepStatus::InvalidSettings);
        MALLOY_CHECK_FALSE(invalid.step().ok());
    }

    // --- THE MILESTONE'S INVARIANT. A magnetic force is always perpendicular
    //     to the velocity, so it does no work and CANNOT change a particle's
    //     speed. Here that is exact, not approximate: the rotation preserves
    //     the length of the velocity to the last bit.
    //
    //     The obvious implementation does not. Applying the force as another
    //     kick, v += (q/m)(v x B) dt, adds a vector perpendicular to v, so the
    //     two are the legs of a right triangle and the speed becomes
    //     |v| * sqrt(1 + alpha^2) with alpha = q*b*dt/m. Here alpha is 6.2832e-3,
    //     so that is a growth of 1.9739e-5 per step, and over 10000 steps a
    //     factor of exp(0.19739) = 1.2182. A cyclotron orbit integrated that
    //     way spirals outward by 22 per cent in ten turns.
    //
    //     This assertion is nine orders tighter than that error. ---
    {
        ChargeSettings s;
        s.k = 0.0; // no pairwise term: this is about the field alone
        s.magnetic = 1.0;

        ChargeWorld w{SimulationSettings{turn_dt}, s,
                      {charge_at(Vec2{}, Vec2{3.0, 0.0}, 1.0, 1.0)}};
        for (int i = 0; i < 10 * steps_per_turn; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        //     The bound is not zero, and the reason is worth being precise
        //     about. A rotation multiplies |v|^2 by cos^2 + sin^2, which is 1
        //     in exact arithmetic but rounds to within about 2 ulp of it, so
        //     |v|^2 can drift by at most N * 2 ulp = 4.4e-12 relative over
        //     10000 steps. Measured here: -5.6e-13 on the speed, comfortably
        //     inside that. The rotation is exact; double precision is not.
        const Real speed = malloy::math::length(w.particles()[0].velocity);
        MALLOY_CHECK_NEAR(speed, 3.0, 1e-11);

        // And the kinetic energy with it, since nothing else is acting. It
        // goes as |v|^2, so it carries twice the relative drift.
        MALLOY_CHECK_NEAR(total_kinetic_energy(w.particles()), 0.5 * 9.0, 1e-10);
    }

    // --- The orbit closes. One turn is exactly `steps_per_turn` steps because
    //     dt was chosen to divide the period, and the velocity has then been
    //     rotated by exactly -2*pi. The displacement over a full turn is
    //     dt * sum(R^n) for n = 1..N, and sum(R^n) is zero when R^N is the
    //     identity, so the particle returns to where it started. ---
    {
        ChargeSettings s;
        s.k = 0.0;
        s.magnetic = 1.0;

        const Vec2 start{2.0, -1.0};
        ChargeWorld w{SimulationSettings{turn_dt}, s,
                      {charge_at(start, Vec2{3.0, 0.0}, 1.0, 1.0)}};
        for (int i = 0; i < steps_per_turn; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].position, start, 1e-12);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, Vec2(3.0, 0.0), 1e-12);
    }

    // --- The cyclotron PERIOD does not depend on speed. Two particles with the
    //     same charge-to-mass ratio and very different speeds both return to
    //     their starting point after the same number of steps, tracing circles
    //     of very different size. ---
    {
        ChargeSettings s;
        s.k = 0.0; // independent particles, so they cannot influence each other
        s.magnetic = 1.0;

        ChargeWorld w{SimulationSettings{turn_dt}, s,
                      {charge_at(Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, 1.0, 1.0),
                       charge_at(Vec2{50.0, 0.0}, Vec2{9.0, 0.0}, 1.0, 1.0)}};
        for (int i = 0; i < steps_per_turn; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].position, Vec2(0.0, 0.0), 1e-12);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[1].position, Vec2(50.0, 0.0), 1e-12);
    }

    // --- The cyclotron RADIUS is m|v|/(|q|b), so unlike the period it does
    //     depend on mass and charge.
    //
    //     Measured without needing to locate the centre at all. The discrete
    //     orbit is a regular polygon: x_n = c + R^n u, so two vertices half a
    //     turn apart differ by R^n (R^(N/2) - I) u, and R^(N/2) is a rotation
    //     by pi, which is -I. The separation is therefore exactly 2|u|, the
    //     diameter, wherever the centre happens to be.
    //
    //     |u| is the polygon's circumradius, |v| dt / (2 sin(|theta|/2)), which
    //     exceeds the true cyclotron radius by (theta/2)/sin(theta/2), about
    //     1 + theta^2/24. Both are derived below rather than written down.
    //
    //     Note the continuum guiding centre, position + (m/qb)(v.y, -v.x), is
    //     NOT the polygon's centre: it is off by about |v| dt / 2, which is
    //     0.0094 here. Measuring the radius from it gives 2.99994 rather than
    //     3.0000049, and looks like a failing radius rather than a wrong
    //     centre. ---
    {
        ChargeSettings s;
        s.k = 0.0;
        s.magnetic = 1.0;

        const Real speed = 3.0;

        struct Case
        {
            Real mass;
            Real charge;
            Real radius; // m|v| / (|q| b)
        };
        const Case cases[] = {{1.0, 1.0, 3.0}, {2.0, 1.0, 6.0}, {1.0, 2.0, 1.5}};

        for (const Case& c : cases)
        {
            const Vec2 start{7.0, -2.0}; // off origin, so nothing cancels
            ChargeWorld w{SimulationSettings{turn_dt}, s,
                          {charge_at(start, Vec2(speed, 0.0), c.mass, c.charge)}};

            // Rotation per step, and therefore the steps in half a turn.
            const Real theta = c.charge * s.magnetic * turn_dt / c.mass;
            const int half_turn =
                static_cast<int>(std::lround(3.14159265358979323846 / theta));

            for (int i = 0; i < half_turn; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }

            const Real diameter = malloy::math::distance(w.particles()[0].position, start);
            const Real circumradius = speed * turn_dt / (2.0 * std::sin(theta / 2.0));
            MALLOY_CHECK_NEAR(diameter, 2.0 * circumradius, 1e-9);

            // And that circumradius really is the cyclotron radius, inflated by
            // the polygon factor and by nothing else.
            MALLOY_CHECK_NEAR(circumradius,
                              c.radius * (theta / 2.0) / std::sin(theta / 2.0), 1e-12);
        }
    }

    // --- Handedness. A positive charge in a positive out-of-plane field turns
    //     CLOCKWISE, and a negative one turns the other way. This is invisible
    //     to speed and invisible to radius, so it needs its own assertion: a
    //     sign error in the 2D cross product would pass both of the tests
    //     above (docs/05, asymmetric configurations). ---
    {
        ChargeSettings s;
        s.k = 0.0;
        s.magnetic = 1.0;

        ChargeWorld w{SimulationSettings{turn_dt}, s,
                      {charge_at(Vec2{}, Vec2{1.0, 0.0}, 1.0, 1.0),
                       charge_at(Vec2{100.0, 0.0}, Vec2{1.0, 0.0}, 1.0, -1.0)}};
        MALLOY_CHECK_TRUE(w.step().ok());

        // Moving along +x, the positive charge is pushed toward -y.
        MALLOY_CHECK_TRUE(w.particles()[0].velocity.y < 0.0);
        // The negative one is pushed the opposite way by the same field.
        MALLOY_CHECK_TRUE(w.particles()[1].velocity.y > 0.0);
        // Mirror images, to the last bit.
        MALLOY_CHECK_NEAR(w.particles()[0].velocity.y, -w.particles()[1].velocity.y, 1e-15);
    }

    // --- E cross B drift. In crossed fields the guiding centre marches at
    //     |E|/b in the direction (E.y/b, -E.x/b), and in the continuum that
    //     velocity depends on NEITHER the charge NOR the mass.
    //
    //     Discretely it depends on both, and the departure is derivable rather
    //     than a tolerance. The velocity map has a fixed point at
    //     v_d = (I - R)^-1 R a dt, which is the continuum drift ROTATED by
    //     half a step's angle and inflated by the same (theta/2)/sin(theta/2)
    //     polygon factor as the cyclotron radius. Both corrections carry
    //     theta = -q b dt / m, so the independence from q and m survives only
    //     to order dt^2.
    //
    //     Asserted exactly, per particle, at 1e-9. The two particles below
    //     gyrate at rates differing by a factor of three, and their drifts
    //     differ in the seventh digit, which is the point. ---
    {
        ChargeSettings s;
        s.k = 0.0;
        s.magnetic = 1.0;
        s.electric = Vec2{0.0, 2.0};

        struct Case
        {
            Vec2 start;
            Real mass;
            Real charge;
        };
        // Released from rest, so each traces a cycloid about a drifting centre.
        const Case cases[] = {{Vec2{0.0, 0.0}, 1.0, 1.0}, {Vec2{0.0, 40.0}, 3.0, 1.0}};

        // Three turns of the light particle, which is exactly one of the heavy
        // one, so both cycloids close and the displacement is the drift alone.
        const int steps = 3 * steps_per_turn;
        const Real elapsed = static_cast<Real>(steps) * turn_dt;

        ChargeWorld w{SimulationSettings{turn_dt}, s,
                      {charge_at(cases[0].start, Vec2{}, cases[0].mass, cases[0].charge),
                       charge_at(cases[1].start, Vec2{}, cases[1].mass, cases[1].charge)}};
        for (int i = 0; i < steps; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }

        const Vec2 continuum{s.electric.y / s.magnetic, -s.electric.x / s.magnetic};
        for (std::size_t i = 0; i < 2; ++i)
        {
            const Real theta =
                -cases[i].charge * s.magnetic * turn_dt / cases[i].mass;
            const Real half = theta / 2.0;
            const Real c = std::cos(half);
            const Real sn = std::sin(half);
            const Real polygon = half / sn;
            const Vec2 drift{(continuum.x * c - continuum.y * sn) * polygon,
                             (continuum.x * sn + continuum.y * c) * polygon};

            MALLOY_CHECK_VEC2_NEAR(w.particles()[i].position,
                                   cases[i].start + drift * elapsed, 1e-9);
        }

        // The continuum claim, correctly qualified: the two drifts agree to
        // order dt^2, which here is about 3e-6 relative.
        MALLOY_CHECK_NEAR(w.particles()[0].position.x / elapsed, continuum.x, 1e-4);
        MALLOY_CHECK_NEAR(
            (w.particles()[1].position.y - 40.0) / elapsed, continuum.y, 1e-2);
    }

    // --- Coulomb is SIGNED, which is the whole difference from gravity. Like
    //     charges repel, unlike attract, and a neutral particle does neither. ---
    {
        ChargeSettings s;
        s.k = 2.0;

        // Two positives at distance 2: force on the right one is
        // +k q1 q2 / r^2 = 2*3*5/4 = 7.5, pointing away, so its acceleration is
        // 7.5 / 1.5 = 5.0 in +x.
        ChargeWorld like{SimulationSettings{0.001}, s,
                         {charge_at(Vec2{-1.0, 0.0}, Vec2{}, 1.0, 3.0),
                          charge_at(Vec2{1.0, 0.0}, Vec2{}, 1.5, 5.0)}};
        const auto a = like.compute_electric_accelerations();
        MALLOY_CHECK_VEC2_NEAR(a[1], Vec2(5.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(a[0], Vec2(-7.5, 0.0), eps); // 7.5 / 1.0, toward -x

        // Flip one sign and the pair attracts instead, with the same magnitude.
        ChargeWorld unlike{SimulationSettings{0.001}, s,
                           {charge_at(Vec2{-1.0, 0.0}, Vec2{}, 1.0, 3.0),
                            charge_at(Vec2{1.0, 0.0}, Vec2{}, 1.5, -5.0)}};
        const auto b = unlike.compute_electric_accelerations();
        MALLOY_CHECK_VEC2_NEAR(b[1], Vec2(-5.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(b[0], Vec2(7.5, 0.0), eps);

        // A neutral particle neither pushes nor is pushed by the pair term.
        ChargeWorld neutral{SimulationSettings{0.001}, s,
                            {charge_at(Vec2{-1.0, 0.0}, Vec2{}, 1.0, 3.0),
                             charge_at(Vec2{1.0, 0.0}, Vec2{}, 1.5, 0.0)}};
        const auto c = neutral.compute_electric_accelerations();
        MALLOY_CHECK_VEC2_NEAR(c[0], Vec2(0.0, 0.0), 0.0);
        MALLOY_CHECK_VEC2_NEAR(c[1], Vec2(0.0, 0.0), 0.0);
    }

    // --- A uniform ELECTRIC field is not gravity. Gravity is an acceleration
    //     and every body falls at the same rate; an electric field is a force
    //     per unit charge, so the acceleration is (q/m)E and a heavier charge
    //     accelerates LESS. ---
    {
        ChargeSettings s;
        s.k = 0.0;
        s.electric = Vec2{4.0, 0.0};

        ChargeWorld w{SimulationSettings{0.5}, s,
                      {charge_at(Vec2{}, Vec2{}, 1.0, 2.0),    // q/m = 2
                       charge_at(Vec2{}, Vec2{}, 8.0, 2.0),    // q/m = 0.25
                       charge_at(Vec2{}, Vec2{}, 1.0, -2.0)}}; // opposite sign
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, Vec2(4.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[1].velocity, Vec2(0.5, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[2].velocity, Vec2(-4.0, 0.0), eps);
    }

    // --- The pair interaction conserves momentum exactly: both members of a
    //     pair are written in the same visit, so the two forces are equal and
    //     opposite by construction. A uniform field does NOT conserve it, and
    //     should not: the field is external. ---
    {
        ChargeSettings s;
        s.k = 1.5;

        const std::vector<ChargedParticle2D> start = {
            charge_at(Vec2{-1.5, 0.3}, Vec2{0.4, -0.2}, 1.0, 2.0),
            charge_at(Vec2{0.7, -0.4}, Vec2{-0.1, 0.5}, 2.5, -1.0),
            charge_at(Vec2{1.9, 1.1}, Vec2{0.25, 0.15}, 0.75, 3.0)};
        const Vec2 p0 = total_momentum(start);

        ChargeWorld w{SimulationSettings{0.0005}, s, start};
        for (int i = 0; i < 2000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_VEC2_NEAR(total_momentum(w.particles()), p0, 1e-12);
    }

    // --- Potential energy carries the sign too. Two like charges have POSITIVE
    //     potential energy, which gravity never does, and it is the energy that
    //     drives them apart. ---
    {
        ChargeSettings s;
        s.k = 2.0;

        // U = k q1 q2 / r = 2 * 3 * 5 / 2 = 15, positive for like charges.
        const std::vector<ChargedParticle2D> like = {
            charge_at(Vec2{-1.0, 0.0}, Vec2{}, 1.0, 3.0),
            charge_at(Vec2{1.0, 0.0}, Vec2{}, 1.0, 5.0)};
        MALLOY_CHECK_NEAR(total_potential_energy(like, s), 15.0, eps);

        const std::vector<ChargedParticle2D> unlike = {
            charge_at(Vec2{-1.0, 0.0}, Vec2{}, 1.0, 3.0),
            charge_at(Vec2{1.0, 0.0}, Vec2{}, 1.0, -5.0)};
        MALLOY_CHECK_NEAR(total_potential_energy(unlike, s), -15.0, eps);

        // The uniform-field term is -q (E . r), the same shape as gravity's.
        ChargeSettings field;
        field.k = 0.0;
        field.electric = Vec2{0.0, -10.0};
        const std::vector<ChargedParticle2D> high = {
            charge_at(Vec2{0.0, 3.0}, Vec2{}, 1.0, 2.0)};
        MALLOY_CHECK_NEAR(total_potential_energy(high, field), 60.0, eps);

        // A magnetic field has no potential at all, because it does no work.
        ChargeSettings magnetic_only;
        magnetic_only.k = 0.0;
        magnetic_only.magnetic = 7.0;
        MALLOY_CHECK_NEAR(total_potential_energy(high, magnetic_only), 0.0, 0.0);
    }

    // --- Softening: two coincident LIKE charges would otherwise repel with
    //     infinite force. This is the charge domain's version of the N-body
    //     coincident-body case, and it matters more here, because nothing stops
    //     two charges of the same sign from being placed on top of each other. ---
    {
        ChargeSettings s;
        s.k = 1.0;
        s.softening = 0.5;

        ChargeWorld w{SimulationSettings{0.001}, s,
                      {charge_at(Vec2{}, Vec2{}, 1.0, 1.0),
                       charge_at(Vec2{}, Vec2{}, 1.0, 1.0)}};
        for (int i = 0; i < 500; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_TRUE(malloy::math::is_finite(w.particles()[0].position));
        MALLOY_CHECK_TRUE(malloy::math::is_finite(w.particles()[1].position));

        // Softening enters the denominator SQUARED, exactly as in docs/04. With
        // r = 0 and softening 0.5, the separation used is 0.5, so the force is
        // k q1 q2 / 0.25 = 4, and with equal masses the pair is symmetric.
        const auto a = w.compute_electric_accelerations();
        MALLOY_CHECK_TRUE(malloy::math::is_finite(a[0]));

        // With no softening at all, a coincident pair contributes nothing
        // rather than NaN.
        ChargeSettings hard;
        hard.k = 1.0;
        ChargeWorld sharp{SimulationSettings{0.001}, hard,
                          {charge_at(Vec2{}, Vec2{}, 1.0, 1.0),
                           charge_at(Vec2{}, Vec2{}, 1.0, 1.0)}};
        MALLOY_CHECK_TRUE(sharp.step().ok());
        MALLOY_CHECK_VEC2_NEAR(sharp.particles()[0].velocity, Vec2(0.0, 0.0), 0.0);
    }

    // --- Softening enters the denominator SQUARED, the same contract as the
    //     N-body domain states in docs/04, and it needs a case where the
    //     softening is comparable to the separation. With r = 1 and softening
    //     2, r2 becomes 1 + 4 = 5, so the force is 5^-1.5 = 0.0894427. If the
    //     softening were added unsquared, r2 would be 3 and the force would be
    //     3^-1.5 = 0.1924501, which is not a subtle difference.
    //
    //     This is the same number the N-body test pins, because it is the same
    //     contract on the same shape of denominator. ---
    {
        ChargeSettings s;
        s.k = 1.0;
        s.softening = 2.0;

        ChargeWorld w{SimulationSettings{0.001}, s,
                      {charge_at(Vec2{0.0, 0.0}, Vec2{}, 1.0, 1.0),
                       charge_at(Vec2{1.0, 0.0}, Vec2{}, 1.0, 1.0)}};
        const auto a = w.compute_electric_accelerations();
        MALLOY_CHECK_VEC2_NEAR(a[1], Vec2(0.0894427190999916, 0.0), 1e-15);
        MALLOY_CHECK_VEC2_NEAR(a[0], Vec2(-0.0894427190999916, 0.0), 1e-15);

        // And the energy uses the same softening, so that total_energy is a
        // quantity the dynamics conserve rather than one they approach:
        // U = k q1 q2 / sqrt(5).
        MALLOY_CHECK_NEAR(total_potential_energy(w.particles(), s),
                          0.4472135954999579, 1e-15);

        // Unsoftened, the same pair would give 1.0 exactly.
        ChargeSettings hard;
        hard.k = 1.0;
        MALLOY_CHECK_NEAR(total_potential_energy(w.particles(), hard), 1.0, 1e-15);
    }

    // --- A world with no magnetic field behaves as a pure electrostatic one,
    //     and the rotation is skipped rather than applied as a zero angle, so
    //     the velocity is bit-identical to the unrotated value. ---
    {
        ChargeSettings s;
        s.k = 1.0;
        ChargeWorld w{SimulationSettings{0.01}, s,
                      {charge_at(Vec2{-1.0, 0.0}, Vec2{0.3, 0.7}, 1.0, 1.0),
                       charge_at(Vec2{1.0, 0.0}, Vec2{-0.3, 0.2}, 1.0, 1.0)}};
        MALLOY_CHECK_TRUE(w.step().ok());
        // 0.7 came only from the initial velocity: nothing acts along y here.
        MALLOY_CHECK_NEAR(w.particles()[0].velocity.y, 0.7, 0.0);
    }

    // --- Energy under the magnetic field alone is exactly conserved, because
    //     the rotation changes no speed and there is no potential to change. ---
    {
        ChargeSettings s;
        s.k = 0.0;
        s.magnetic = 2.5;

        const std::vector<ChargedParticle2D> start = {
            charge_at(Vec2{1.0, 2.0}, Vec2{0.3, -0.7}, 1.5, 2.0),
            charge_at(Vec2{-4.0, 0.5}, Vec2{-1.1, 0.9}, 0.5, -3.0)};
        const Real e0 = total_energy(start, s);

        ChargeWorld w{SimulationSettings{0.0005}, s, start};
        for (int i = 0; i < 4000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        // Same bound as the speed test: |v|^2 drifts by at most N * 2 ulp
        // through the rounding of cos^2 + sin^2, which over 4000 steps is
        // 1.8e-12 relative, or 1.7e-12 on this energy.
        MALLOY_CHECK_NEAR(total_energy(w.particles(), s), e0, 1e-11);
    }

    // --- Determinism: the same initial state and the same timestep give the
    //     same run (docs/04). ---
    {
        ChargeSettings s;
        s.k = 1.25;
        s.magnetic = 0.6;
        s.electric = Vec2{0.2, -0.1};

        const std::vector<ChargedParticle2D> start = {
            charge_at(Vec2{0.4, -0.3}, Vec2{0.0, 0.2}, 2.0, 1.0),
            charge_at(Vec2{-0.6, 0.5}, Vec2{0.1, 0.0}, 1.0, -2.0),
            charge_at(Vec2{0.1, 0.9}, Vec2{-0.2, 0.1}, 0.5, 3.0)};

        ChargeWorld a{SimulationSettings{0.001}, s, start};
        ChargeWorld b{SimulationSettings{0.001}, s, start};
        for (int i = 0; i < 1000; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        for (std::size_t i = 0; i < a.particles().size(); ++i)
        {
            MALLOY_CHECK_VEC2_NEAR(a.particles()[i].position, b.particles()[i].position, 0.0);
            MALLOY_CHECK_VEC2_NEAR(a.particles()[i].velocity, b.particles()[i].velocity, 0.0);
        }
        MALLOY_CHECK_EQ(a.tick_count(), std::uint64_t{1000});
        MALLOY_CHECK_NEAR(a.elapsed_time(), 1.0, 1e-12);
    }

    // --- A failing step leaves the state untouched and reports rather than
    //     throwing, like every other world here. ---
    {
        ChargeSettings s;
        ChargeWorld w{SimulationSettings{-1.0}, s,
                      {charge_at(Vec2{1.0, 2.0}, Vec2{}, 1.0, 1.0)}};
        MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidSettings);
        MALLOY_CHECK_FALSE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].position, Vec2(1.0, 2.0), 0.0);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
    }

    // --- Issue #17: a softening large enough that its SQUARE overflows is
    //     refused, in both domains that have one.
    //
    //     It has to be refused rather than tolerated, because the failure is
    //     invisible in every diagnostic. A squared softening of infinity makes
    //     every separation infinite, so every force is exactly zero and every
    //     potential is exactly zero, and a simulation in which nothing happens
    //     conserves everything perfectly. Measured before the fix: two like
    //     charges that should fly apart sat still for 2000 steps while the
    //     energy column read 0.00000000e+00, when the correct value was 0.5. ---
    {
        ChargeSettings ok;
        ok.softening = 1.0e154; // squares to 1e308, still finite
        MALLOY_CHECK_TRUE(ok.is_valid());

        ChargeSettings over;
        over.softening = 1.4e154; // squares to just over DBL_MAX
        MALLOY_CHECK_FALSE(over.is_valid());

        ChargeSettings absurd;
        absurd.softening = 1.0e200;
        MALLOY_CHECK_FALSE(absurd.is_valid());

        // Still refused for the older reasons.
        ChargeSettings negative;
        negative.softening = -1.0;
        MALLOY_CHECK_FALSE(negative.is_valid());
        ChargeSettings not_finite;
        not_finite.softening = inf;
        MALLOY_CHECK_FALSE(not_finite.is_valid());
    }
    {
        ChargeSettings over;
        over.k = 1.0;
        over.softening = 1.0e200;
        ChargeWorld w{SimulationSettings{0.01}, over,
                      {charge_at(Vec2{-1.0, 0.0}, Vec2{}, 1.0, 1.0),
                       charge_at(Vec2{1.0, 0.0}, Vec2{}, 1.0, 1.0)}};
        MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidSettings);
        MALLOY_CHECK_FALSE(w.step().ok());
    }

    // --- Issue #16: a position or velocity whose SQUARE overflows is refused.
    //     Finite is not enough. Everything that squares a vector reaches
    //     infinity above about 1.34e154, so accepting state up to 1.8e308 left
    //     a window in which a world validated clean while every energy it
    //     reported was inf. ---
    {
        const Real too_big = 1.4e154; // finite, and its square is not
        MALLOY_CHECK_TRUE(charge_at(Vec2{1.0, 2.0}, Vec2{}, 1.0, 1.0).is_valid());
        MALLOY_CHECK_FALSE(charge_at(Vec2{too_big, 0.0}, Vec2{}, 1.0, 1.0).is_valid());
        MALLOY_CHECK_FALSE(charge_at(Vec2{}, Vec2{0.0, too_big}, 1.0, 1.0).is_valid());
    }

    // --- M29: charged particles in three dimensions. The Coulomb and electric
    //     terms are the 2D ones with a third component; the magnetic field
    //     becomes a VECTOR, and its turn is a rotation about the field
    //     direction, so the velocity component ALONG the field is carried and
    //     the motion is helical. A magnetic force still does no work. ---
    {
        using malloy::math::approx_equal;

        const auto charge3 = [](Vec3 position, Vec3 velocity, Real mass, Real charge) {
            ChargedParticle3D p;
            p.position = position;
            p.velocity = velocity;
            p.mass = mass;
            p.charge = charge;
            return p;
        };
        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();
        const Real too_big3 = 1.4e154; // finite, but its square overflows

        // --- Validation. ---
        {
            MALLOY_CHECK_TRUE(charge3(Vec3{1.0, 2.0, 3.0}, Vec3{}, 1.0, -2.0).is_valid());
            MALLOY_CHECK_FALSE(charge3(Vec3{}, Vec3{}, 0.0, 1.0).is_valid()); // mass 0
            MALLOY_CHECK_FALSE(charge3(Vec3{}, Vec3{}, 1.0, nan3).is_valid());
            MALLOY_CHECK_FALSE(charge3(Vec3{0.0, too_big3, 0.0}, Vec3{}, 1.0, 1.0).is_valid());
            MALLOY_CHECK_FALSE(charge3(Vec3{}, Vec3{0.0, 0.0, too_big3}, 1.0, 1.0).is_valid());
            // Settings: k finite, fields squarable, softening non-negative and
            // squarable.
            Charge3DSettings good;
            good.electric = Vec3{1.0, -2.0, 0.5};
            good.magnetic = Vec3{0.0, 0.0, 3.0};
            MALLOY_CHECK_TRUE(good.is_valid());
            Charge3DSettings bad_field;
            bad_field.magnetic = Vec3{too_big3, 0.0, 0.0}; // square overflows |B|
            MALLOY_CHECK_FALSE(bad_field.is_valid());
            Charge3DSettings bad_efield;
            bad_efield.electric = Vec3{0.0, inf3, 0.0};
            MALLOY_CHECK_FALSE(bad_efield.is_valid());
            Charge3DSettings bad_soft;
            bad_soft.softening = -1.0;
            MALLOY_CHECK_FALSE(bad_soft.is_valid());
            // A world with bad settings refuses to step.
            Charge3DWorld world{SimulationSettings{0.01}, bad_field, {charge3(Vec3{}, Vec3{}, 1.0, 1.0)}};
            MALLOY_CHECK_TRUE(world.step().status == StepStatus::InvalidSettings);
        }

        // --- The electric field is a force per charge, divided by mass: a
        //     heavier charge accelerates less and an opposite charge the other
        //     way, the whole difference from gravity, now in 3D. ---
        {
            const Real dt = 0.001;
            Charge3DSettings settings;
            settings.k = 0.0; // field only
            settings.electric = Vec3{2.0, -1.0, 3.0};
            Charge3DWorld world{
                SimulationSettings{dt}, settings,
                {charge3(Vec3{}, Vec3{}, 1.0, 1.0),   // light, positive
                 charge3(Vec3{}, Vec3{}, 4.0, 1.0),   // heavy, positive
                 charge3(Vec3{}, Vec3{}, 1.0, -1.0)}};// light, negative
            for (int i = 0; i < 100; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            const auto& p = world.particles();
            const Real t = 100 * dt;
            // a = (q/m) E, so v = (q/m) E t.
            MALLOY_CHECK_TRUE(approx_equal(p[0].velocity, settings.electric * t, 1e-9));
            MALLOY_CHECK_TRUE(approx_equal(p[1].velocity, settings.electric * (t / 4.0), 1e-9));
            MALLOY_CHECK_TRUE(approx_equal(p[2].velocity, settings.electric * -t, 1e-9));
        }

        // --- A magnetic force does no work: even with a velocity component
        //     ALONG a TILTED field (so the rotation is genuinely 3D, not an
        //     axis-aligned special case), the speed is preserved and the energy
        //     is purely kinetic and constant. ---
        {
            Charge3DSettings settings;
            settings.k = 0.0;
            settings.magnetic = Vec3{1.0, 2.0, -0.5}; // tilted
            Charge3DWorld world{SimulationSettings{0.01}, settings,
                                {charge3(Vec3{}, Vec3{0.7, -1.3, 0.4}, 1.5, 2.0)}};
            const Real speed0 = malloy::math::length(world.particles().front().velocity);
            const Real energy0 = total_energy3d(world.particles(), settings);
            Real worst_speed = 0.0;
            Real worst_energy = 0.0;
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                const Real speed = malloy::math::length(world.particles().front().velocity);
                worst_speed = std::fmax(worst_speed, std::abs(speed - speed0));
                worst_energy = std::fmax(
                    worst_energy, std::abs(total_energy3d(world.particles(), settings) - energy0));
            }
            MALLOY_CHECK_TRUE(worst_speed < 1e-12);
            MALLOY_CHECK_TRUE(worst_energy < 1e-12);
        }

        // --- Helical motion: with the field along z, the velocity component
        //     along z is carried EXACTLY (a rotation about z leaves the z axis
        //     untouched, to the last bit), while the perpendicular part circles.
        //     A full cyclotron period returns the velocity. This is the drift
        //     plus circle that two dimensions cannot express. ---
        {
            Charge3DSettings settings;
            settings.k = 0.0;
            settings.magnetic = Vec3{0.0, 0.0, 1.0}; // b = 1, so omega = q/m
            // v_perp = 3, v_par = 4: speed 5, radius m v_perp/(q b) = 3, period
            // 2 pi m /(q b) = one turn in steps_per_turn steps at turn_dt.
            Charge3DWorld world{SimulationSettings{turn_dt}, settings,
                                {charge3(Vec3{}, Vec3{3.0, 0.0, 4.0}, 1.0, 1.0)}};
            const Vec3 v0 = world.particles().front().velocity;
            for (int i = 0; i < steps_per_turn; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                // The along-field component never moves, at any step, exactly.
                MALLOY_CHECK_NEAR(world.particles().front().velocity.z, 4.0, 0.0);
            }
            // After one full period the velocity has come back (to rounding).
            MALLOY_CHECK_TRUE(approx_equal(world.particles().front().velocity, v0, 1e-12));
            // And the guiding centre has drifted along z by v_par * period.
            const Real period = turn_dt * static_cast<Real>(steps_per_turn);
            MALLOY_CHECK_NEAR(world.particles().front().position.z, 4.0 * period, 1e-9);
        }

        // --- The turn rate and sense match the field: a positive charge moving
        //     +x in a +z field turns toward -y, the same convention as 2D, and
        //     the angle per step is -(q b / m) dt. One quarter turn lands +x on
        //     -y. ---
        {
            Charge3DSettings settings;
            settings.k = 0.0;
            settings.magnetic = Vec3{0.0, 0.0, 1.0};
            Charge3DWorld world{SimulationSettings{turn_dt}, settings,
                                {charge3(Vec3{}, Vec3{2.0, 0.0, 0.0}, 1.0, 1.0)}};
            for (int i = 0; i < steps_per_turn / 4; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            // A quarter turn clockwise: +x velocity is now along -y.
            const Vec3 v = world.particles().front().velocity;
            MALLOY_CHECK_TRUE(approx_equal(v, Vec3{0.0, -2.0, 0.0}, 1e-3));
        }

        // --- It reduces to the 2D domain: the same particles in the z = 0
        //     plane, with B along z and E in the plane, must trace the same xy
        //     path as ChargeWorld, and never leave the plane. This ties the 3D
        //     integrator to the proven 2D one. ---
        {
            const Real dt = 0.01;
            const Real b = 0.8;
            const Vec2 e2{0.5, -0.3};
            ChargeSettings s2;
            s2.k = 1.0;
            s2.electric = e2;
            s2.magnetic = b;
            s2.softening = 0.25;
            Charge3DSettings s3;
            s3.k = 1.0;
            s3.electric = Vec3{e2.x, e2.y, 0.0};
            s3.magnetic = Vec3{0.0, 0.0, b};
            s3.softening = 0.25;
            ChargeWorld w2{SimulationSettings{dt}, s2,
                           {charge_at(Vec2{-1.0, 0.0}, Vec2{0.0, 1.0}, 1.0, 1.0),
                            charge_at(Vec2{1.0, 0.0}, Vec2{0.0, -1.0}, 2.0, -1.5)}};
            Charge3DWorld w3{SimulationSettings{dt}, s3,
                             {charge3(Vec3{-1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, 1.0, 1.0),
                              charge3(Vec3{1.0, 0.0, 0.0}, Vec3{0.0, -1.0, 0.0}, 2.0, -1.5)}};
            for (int i = 0; i < 500; ++i)
            {
                MALLOY_CHECK_TRUE(w2.step().ok());
                MALLOY_CHECK_TRUE(w3.step().ok());
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto& a = w2.particles()[k];
                    const auto& c = w3.particles()[k];
                    MALLOY_CHECK_TRUE(approx_equal(c.position, Vec3{a.position.x, a.position.y, 0.0}, 1e-10));
                    MALLOY_CHECK_TRUE(approx_equal(c.velocity, Vec3{a.velocity.x, a.velocity.y, 0.0}, 1e-10));
                }
            }
        }

        // --- Pairwise Coulomb in 3D: two like charges off-axis repel along
        //     their 3D separation, and the equal-and-opposite forces keep total
        //     momentum exactly zero from a symmetric rest start. ---
        {
            Charge3DSettings settings;
            settings.k = 1.0;
            settings.magnetic = Vec3{}; // no field, pure Coulomb
            Charge3DWorld world{SimulationSettings{0.001}, settings,
                                {charge3(Vec3{-1.0, -1.0, -1.0}, Vec3{}, 1.0, 1.0),
                                 charge3(Vec3{1.0, 1.0, 1.0}, Vec3{}, 1.0, 1.0)}};
            for (int i = 0; i < 200; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            const auto& p = world.particles();
            // Started at rest, equal masses, so momentum stays zero and they
            // move apart symmetrically along the body diagonal.
            MALLOY_CHECK_TRUE(approx_equal(total_momentum3d(p), Vec3{}, 1e-12));
            MALLOY_CHECK_TRUE(malloy::math::length(p[1].position) >
                              malloy::math::length(Vec3{1.0, 1.0, 1.0})); // pushed out
            MALLOY_CHECK_TRUE(approx_equal(p[0].position, p[1].position * -1.0, 1e-12));
        }

        // --- The diagnostics pin their own factors, before any step. Kinetic
        //     energy is 0.5 m |v|^2, momentum is the mass-weighted velocity sum,
        //     and the potential is the pair term k q_i q_j / r plus the field
        //     term -q (E . r). ---
        {
            const std::vector<ChargedParticle3D> two = {
                charge3(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 2.0, 2.0}, 3.0, 1.0),
                charge3(Vec3{5.0, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0}, 4.0, 1.0)};
            // KE = 0.5*3*(1+4+4) + 0.5*4*(1) = 13.5 + 2 = 15.5.
            MALLOY_CHECK_NEAR(total_kinetic_energy3d(two), 15.5, 1e-12);
            // p = 3*(1,2,2) + 4*(-1,0,0) = (-1, 6, 6).
            MALLOY_CHECK_TRUE(approx_equal(total_momentum3d(two), Vec3{-1.0, 6.0, 6.0}, 1e-12));
            // Pair only (k=1, E=0): k q0 q1 / r = 1/5.
            Charge3DSettings pair_only; // k defaults to 1
            MALLOY_CHECK_NEAR(total_potential_energy3d(two, pair_only), 1.0 / 5.0, 1e-12);
            // Field only: U = -q0 (E.r0) - q1 (E.r1) = 0 - 1*(2*5) = -10.
            Charge3DSettings field_only;
            field_only.k = 0.0;
            field_only.electric = Vec3{2.0, 0.0, 0.0};
            MALLOY_CHECK_NEAR(total_potential_energy3d(two, field_only), -10.0, 1e-12);
        }
    }

    std::cout << "malloy_charges_tests passed\n";
    return 0;
}
