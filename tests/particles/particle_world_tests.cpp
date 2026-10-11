#include <malloy/particles/particles.hpp>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <test_check.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

using malloy::collide::Aabb;
using malloy::collide::Aabb3;
using malloy::math::Real;
using malloy::math::Vec2;
using malloy::math::Vec3;
using malloy::particles::Particle2D;
using malloy::particles::Particle3D;
using malloy::particles::ParticleSettings;
using malloy::particles::ParticleSettings3D;
using malloy::particles::ParticleWorld;
using malloy::particles::ParticleWorld3D;
using malloy::particles::total_potential_energy;
using malloy::particles::total_potential_energy3d;
using malloy::particles::total_energy;
using malloy::particles::total_energy3d;
using malloy::particles::total_kinetic_energy;
using malloy::particles::total_kinetic_energy3d;
using malloy::particles::total_momentum;
using malloy::particles::total_momentum3d;
using malloy::particles::lattice_gas;
using malloy::particles::lattice_gas3d;
using malloy::particles::pressure;
using malloy::particles::pressure3d;
using malloy::particles::slow_fraction;
using malloy::particles::slow_fraction3d;
using malloy::particles::speed_ratio;
using malloy::particles::speed_ratio3d;
using malloy::particles::temperature;
using malloy::particles::temperature3d;
using malloy::particles::WallImpulse2D;
using malloy::particles::WallImpulse3D;
using malloy::sim_core::SimulationSettings;
using malloy::sim_core::StepStatus;

namespace
{
// A roomy box, so tests that care about pair collisions are not perturbed by
// the walls.
Aabb big_box()
{
    return Aabb{Vec2{-100.0, -100.0}, Vec2{100.0, 100.0}};
}
} // namespace

int main()
{
    const Real eps = 1e-12;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    // --- Validation: settings ---
    {
        const ParticleSettings good{1.0, big_box()};
        const ParticleSettings inelastic{0.0, big_box()};
        const ParticleSettings too_bouncy{1.5, big_box()};
        const ParticleSettings negative{-0.1, big_box()};
        const ParticleSettings nan_restitution{nan, big_box()};
        const ParticleSettings bad_box{1.0, Aabb{Vec2{1.0, 1.0}, Vec2{-1.0, -1.0}}};

        MALLOY_CHECK_TRUE(good.is_valid());
        MALLOY_CHECK_TRUE(inelastic.is_valid()); // 0 is legal
        MALLOY_CHECK_FALSE(too_bouncy.is_valid());
        MALLOY_CHECK_FALSE(negative.is_valid());
        MALLOY_CHECK_FALSE(nan_restitution.is_valid());
        MALLOY_CHECK_FALSE(bad_box.is_valid());
    }

    // --- Validation: particles ---
    {
        const Particle2D good{Vec2{}, Vec2{}, 1.0, 0.5};
        const Particle2D point{Vec2{}, Vec2{}, 1.0, 0.0}; // zero radius is legal
        const Particle2D zero_mass{Vec2{}, Vec2{}, 0.0, 0.5};
        const Particle2D negative_radius{Vec2{}, Vec2{}, 1.0, -0.5};
        const Particle2D inf_mass{Vec2{}, Vec2{}, inf, 0.5};
        const Particle2D nan_position{Vec2{nan, 0.0}, Vec2{}, 1.0, 0.5};

        MALLOY_CHECK_TRUE(good.is_valid());
        MALLOY_CHECK_TRUE(point.is_valid());
        MALLOY_CHECK_FALSE(zero_mass.is_valid());
        MALLOY_CHECK_FALSE(negative_radius.is_valid());
        MALLOY_CHECK_FALSE(inf_mass.is_valid());
        MALLOY_CHECK_FALSE(nan_position.is_valid());
    }

    // --- A particle too big for the box is rejected, because clamping it
    //     against both walls at once would fight every step. ---
    {
        const std::vector<Particle2D> big = {Particle2D{Vec2{}, Vec2{}, 1.0, 5.0}};
        ParticleWorld w{SimulationSettings{0.01},
                        ParticleSettings{1.0, Aabb{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}}},
                        big};
        MALLOY_CHECK_TRUE(w.validate() == StepStatus::InvalidState);
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
    }

    // --- Failed validation leaves the state untouched. ---
    {
        const std::vector<Particle2D> before = {
            Particle2D{Vec2{0.0, 0.0}, Vec2{1.0, 2.0}, 1.0, 0.5}};
        ParticleWorld w{SimulationSettings{0.0}, ParticleSettings{1.0, big_box()},
                        before};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].position, before[0].position, 0.0);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, before[0].velocity, 0.0);
    }

    // --- Free flight: with nothing to hit, a particle moves in a straight line
    //     at exactly velocity * dt per step. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.0, 0.0}, Vec2{1.0, 2.0}, 1.0, 0.1}};
        ParticleWorld w{SimulationSettings{0.5}, ParticleSettings{1.0, big_box()}, p};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].position, Vec2(0.5, 1.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, Vec2(1.0, 2.0), eps);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{1});
    }

    // --- Head-on equal-mass elastic collision: the classic result is that the
    //     two particles exchange velocities exactly. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{-1.0, 0.0}, Vec2{1.0, 0.0}, 1.0, 0.5},
            Particle2D{Vec2{1.0, 0.0}, Vec2{-1.0, 0.0}, 1.0, 0.5}};
        ParticleWorld w{SimulationSettings{0.5}, ParticleSettings{1.0, big_box()}, p};

        // One step brings them to (-0.5, 0) and (0.5, 0): exactly touching.
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, Vec2(-1.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[1].velocity, Vec2(1.0, 0.0), eps);
    }

    // --- Restitution 0 is perfectly inelastic: the pair stops separating and
    //     ends with a common velocity along the normal. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{-0.4, 0.0}, Vec2{1.0, 0.0}, 1.0, 0.5},
            Particle2D{Vec2{0.4, 0.0}, Vec2{-1.0, 0.0}, 1.0, 0.5}};
        ParticleWorld w{SimulationSettings{0.001}, ParticleSettings{0.0, big_box()}, p};
        MALLOY_CHECK_TRUE(w.step().ok());
        const Real closing = w.particles()[1].velocity.x - w.particles()[0].velocity.x;
        MALLOY_CHECK_NEAR(closing, 0.0, eps); // no separation at all
    }

    // --- INVARIANT: particle/particle collisions conserve momentum at every
    //     restitution, because the impulse is equal and opposite. ---
    {
        for (const Real restitution : {0.0, 0.35, 1.0})
        {
            const std::vector<Particle2D> p = {
                Particle2D{Vec2{-0.4, 0.1}, Vec2{2.0, -0.5}, 1.0, 0.5},
                Particle2D{Vec2{0.4, -0.1}, Vec2{-1.0, 0.25}, 3.0, 0.5}};
            ParticleWorld w{SimulationSettings{0.001},
                            ParticleSettings{restitution, big_box()}, p};
            const Vec2 before = total_momentum(p);
            for (int i = 0; i < 200; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            MALLOY_CHECK_VEC2_NEAR(total_momentum(w.particles()), before, 1e-12);
        }
    }

    // --- INVARIANT: restitution 1 also conserves kinetic energy, and anything
    //     below it strictly removes energy. This is what makes restitution
    //     itself testable rather than just plumbed through. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{-0.4, 0.1}, Vec2{2.0, -0.5}, 1.0, 0.5},
            Particle2D{Vec2{0.4, -0.1}, Vec2{-1.0, 0.25}, 3.0, 0.5}};
        const Real before = total_kinetic_energy(p);

        ParticleWorld elastic{SimulationSettings{0.001},
                              ParticleSettings{1.0, big_box()}, p};
        for (int i = 0; i < 200; ++i)
        {
            MALLOY_CHECK_TRUE(elastic.step().ok());
        }
        MALLOY_CHECK_NEAR(total_kinetic_energy(elastic.particles()), before, 1e-12);

        ParticleWorld lossy{SimulationSettings{0.001},
                            ParticleSettings{0.5, big_box()}, p};
        for (int i = 0; i < 200; ++i)
        {
            MALLOY_CHECK_TRUE(lossy.step().ok());
        }
        MALLOY_CHECK_TRUE(total_kinetic_energy(lossy.particles()) < before);
    }

    // --- Walls: an elastic bounce reverses the normal component and preserves
    //     speed, while the tangential component is untouched. ---
    {
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.4, 0.0}, Vec2{1.0, 0.5}, 1.0, 0.5}};
        ParticleWorld w{SimulationSettings{0.5}, ParticleSettings{1.0, box}, p};
        MALLOY_CHECK_TRUE(w.step().ok());
        // Would reach x = 0.9, so its surface at 1.4 is past the wall: clamped
        // back to 0.5 and the x velocity reverses. y is untouched.
        MALLOY_CHECK_NEAR(w.particles()[0].position.x, 0.5, eps);
        MALLOY_CHECK_NEAR(w.particles()[0].velocity.x, -1.0, eps);
        MALLOY_CHECK_NEAR(w.particles()[0].velocity.y, 0.5, eps);
    }
    {
        // Restitution 0.5 halves the rebound speed.
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.4, 0.0}, Vec2{1.0, 0.0}, 1.0, 0.5}};
        ParticleWorld w{SimulationSettings{0.5}, ParticleSettings{0.5, box}, p};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_NEAR(w.particles()[0].velocity.x, -0.5, eps);
    }

    // --- An overlapping pair that is ALREADY separating must not be impulsed
    //     again. A second impulse conserves both momentum and energy (it just
    //     flips the pair back together), so only checking the velocities
    //     directly can catch it. ---
    {
        // Overlapping by 0.5 and already moving apart at relative speed 2.
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{-0.25, 0.0}, Vec2{-1.0, 0.0}, 1.0, 0.5},
            Particle2D{Vec2{0.25, 0.0}, Vec2{1.0, 0.0}, 1.0, 0.5}};
        ParticleWorld w{SimulationSettings{0.001}, ParticleSettings{1.0, big_box()}, p};
        MALLOY_CHECK_TRUE(w.step().ok());
        // Only the positional correction should act; velocities stay put.
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, Vec2(-1.0, 0.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[1].velocity, Vec2(1.0, 0.0), eps);
    }

    // --- All four walls, not just one. Each is a separate branch, so testing
    //     only the right wall leaves three untested. ---
    {
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        struct WallCase
        {
            Vec2 position;
            Vec2 velocity;
            Vec2 expected_position;
            Vec2 expected_velocity;
        };
        const WallCase cases[] = {
            // left wall: reaches x = -0.9, surface past -1, clamped to -0.5
            {Vec2{-0.4, 0.0}, Vec2{-1.0, 0.25}, Vec2{-0.5, 0.125}, Vec2{1.0, 0.25}},
            // right wall
            {Vec2{0.4, 0.0}, Vec2{1.0, 0.25}, Vec2{0.5, 0.125}, Vec2{-1.0, 0.25}},
            // bottom wall
            {Vec2{0.0, -0.4}, Vec2{0.25, -1.0}, Vec2{0.125, -0.5}, Vec2{0.25, 1.0}},
            // top wall
            {Vec2{0.0, 0.4}, Vec2{0.25, 1.0}, Vec2{0.125, 0.5}, Vec2{0.25, -1.0}},
        };
        for (const WallCase& c : cases)
        {
            const std::vector<Particle2D> p = {Particle2D{c.position, c.velocity, 1.0, 0.5}};
            ParticleWorld w{SimulationSettings{0.5}, ParticleSettings{1.0, box}, p};
            MALLOY_CHECK_TRUE(w.step().ok());
            MALLOY_CHECK_VEC2_NEAR(w.particles()[0].position, c.expected_position, eps);
            MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, c.expected_velocity, eps);
        }
    }

    // --- Containment: over a long run in a small box, nothing escapes. This is
    //     the property a user would actually notice breaking. ---
    {
        const Aabb box{Vec2{-2.0, -2.0}, Vec2{2.0, 2.0}};
        std::vector<Particle2D> p;
        for (int i = 0; i < 8; ++i)
        {
            const Real f = static_cast<Real>(i);
            p.push_back(Particle2D{Vec2{-1.5 + f * 0.4, -1.0 + f * 0.25},
                                   Vec2{1.0 - f * 0.3, -0.7 + f * 0.2}, 1.0 + f * 0.5,
                                   0.15});
        }
        ParticleWorld w{SimulationSettings{0.005}, ParticleSettings{0.9, box}, p};
        for (int i = 0; i < 4000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        for (const Particle2D& particle : w.particles())
        {
            MALLOY_CHECK_TRUE(malloy::math::is_finite(particle.position));
            MALLOY_CHECK_TRUE(particle.position.x - particle.radius >= box.min.x - 1e-9);
            MALLOY_CHECK_TRUE(particle.position.x + particle.radius <= box.max.x + 1e-9);
            MALLOY_CHECK_TRUE(particle.position.y - particle.radius >= box.min.y - 1e-9);
            MALLOY_CHECK_TRUE(particle.position.y + particle.radius <= box.max.y + 1e-9);
        }
    }

    // --- No energy created from nothing: with restitution 1 a crowded box must
    //     not gain kinetic energy over a long run. A double impulse on an
    //     already separating pair is the classic way this goes wrong. ---
    {
        const Aabb box{Vec2{-2.0, -2.0}, Vec2{2.0, 2.0}};
        std::vector<Particle2D> p;
        for (int i = 0; i < 9; ++i)
        {
            const Real f = static_cast<Real>(i);
            p.push_back(Particle2D{Vec2{-1.2 + f * 0.3, -1.2 + f * 0.28},
                                   Vec2{0.9 - f * 0.2, -0.6 + f * 0.15}, 1.0, 0.2});
        }
        const Real before = total_kinetic_energy(p);
        ParticleWorld w{SimulationSettings{0.002}, ParticleSettings{1.0, box}, p};
        for (int i = 0; i < 3000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_TRUE(total_kinetic_energy(w.particles()) <= before * 1.0000001);
    }

    // --- After an elastic collision the pair must actually end up apart and
    //     moving apart. Re-impulsing an already separating pair still conserves
    //     momentum AND energy (it just flips them back together), so only a
    //     separation check catches it: the particles would stick and jitter. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{-0.6, 0.0}, Vec2{1.0, 0.0}, 1.0, 0.5},
            Particle2D{Vec2{0.6, 0.0}, Vec2{-1.0, 0.0}, 1.0, 0.5}};
        ParticleWorld w{SimulationSettings{0.01}, ParticleSettings{1.0, big_box()}, p};
        for (int i = 0; i < 200; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const Vec2 delta = w.particles()[1].position - w.particles()[0].position;
        MALLOY_CHECK_TRUE(malloy::math::length(delta) > 1.0); // past the radius sum
        const Vec2 relative = w.particles()[1].velocity - w.particles()[0].velocity;
        MALLOY_CHECK_TRUE(malloy::math::dot(relative, delta) > 0.0); // moving apart
    }

    // --- Positional correction is split by inverse mass, so the light particle
    //     moves further out of an overlap. Equal masses make the split
    //     invisible, and it never touches velocity, so neither the momentum nor
    //     the energy invariant can see it. ---
    {
        // Overlapping by 0.5 and at rest, so only the correction acts.
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{-0.25, 0.0}, Vec2{}, 1.0, 0.5},  // light
            Particle2D{Vec2{0.25, 0.0}, Vec2{}, 9.0, 0.5}};  // nine times heavier
        ParticleWorld w{SimulationSettings{0.001}, ParticleSettings{1.0, big_box()}, p};
        MALLOY_CHECK_TRUE(w.step().ok());

        const Real light_moved = std::abs(w.particles()[0].position.x + 0.25);
        const Real heavy_moved = std::abs(w.particles()[1].position.x - 0.25);
        // Inverse masses 1 and 1/9, so the light one takes 9/10 of the push.
        MALLOY_CHECK_NEAR(light_moved, 0.45, 1e-12);
        MALLOY_CHECK_NEAR(heavy_moved, 0.05, 1e-12);
        // And the pair ends exactly touching: the total push equals the depth.
        MALLOY_CHECK_NEAR(light_moved + heavy_moved, 0.5, 1e-12);
    }

    // --- Determinism: identical worlds stepped identically stay identical. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{-0.4, 0.1}, Vec2{2.0, -0.5}, 1.0, 0.3},
            Particle2D{Vec2{0.4, -0.1}, Vec2{-1.0, 0.25}, 2.0, 0.3},
            Particle2D{Vec2{0.0, 0.6}, Vec2{0.1, -1.5}, 0.5, 0.3}};
        ParticleWorld a{SimulationSettings{0.003},
                        ParticleSettings{0.8, Aabb{Vec2{-2.0, -2.0}, Vec2{2.0, 2.0}}},
                        p};
        ParticleWorld b{SimulationSettings{0.003},
                        ParticleSettings{0.8, Aabb{Vec2{-2.0, -2.0}, Vec2{2.0, 2.0}}},
                        p};
        for (int i = 0; i < 500; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        for (std::size_t i = 0; i < p.size(); ++i)
        {
            MALLOY_CHECK_VEC2_NEAR(a.particles()[i].position, b.particles()[i].position,
                                   0.0);
            MALLOY_CHECK_VEC2_NEAR(a.particles()[i].velocity, b.particles()[i].velocity,
                                   0.0);
        }
    }

    // --- Diagnostics on an empty world. ---
    {
        const std::vector<Particle2D> none;
        MALLOY_CHECK_VEC2_NEAR(total_momentum(none), Vec2(0.0, 0.0), eps);
        MALLOY_CHECK_NEAR(total_kinetic_energy(none), 0.0, eps);
    }

    // --- Elapsed simulation time comes from malloy_time, shared with the other
    //     domain rather than reimplemented here. dt = 0.1 separates ticks * dt
    //     from an accumulated sum. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.0, 0.0}, Vec2{0.1, 0.0}, 1.0, 0.1}};
        ParticleWorld w{SimulationSettings{0.1}, ParticleSettings{1.0, big_box()}, p};
        MALLOY_CHECK_NEAR(w.elapsed_time(), 0.0, 0.0);
        for (int i = 0; i < 10; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{10});
        MALLOY_CHECK_NEAR(w.elapsed_time(), 1.0, 1e-17);
    }

    // --- M12: zero gravity must reproduce the pre-gravity behaviour EXACTLY.
    //     This is the backward-compatibility guard for every scenario written
    //     before the setting existed. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{-0.4, 0.1}, Vec2{2.0, -0.5}, 1.0, 0.5},
            Particle2D{Vec2{0.4, -0.1}, Vec2{-1.0, 0.25}, 3.0, 0.5}};

        ParticleSettings unset{0.9, big_box()};        // gravity defaults to zero
        ParticleSettings explicit_zero{0.9, big_box()};
        explicit_zero.gravity = Vec2{0.0, 0.0};

        ParticleWorld a{SimulationSettings{0.002}, unset, p};
        ParticleWorld b{SimulationSettings{0.002}, explicit_zero, p};
        for (int i = 0; i < 400; ++i)
        {
            MALLOY_CHECK_TRUE(a.step().ok());
            MALLOY_CHECK_TRUE(b.step().ok());
        }
        for (std::size_t i = 0; i < p.size(); ++i)
        {
            MALLOY_CHECK_VEC2_NEAR(a.particles()[i].position, b.particles()[i].position,
                                   0.0);
            MALLOY_CHECK_VEC2_NEAR(a.particles()[i].velocity, b.particles()[i].velocity,
                                   0.0);
        }
        // And with no gravity the potential term vanishes, so total energy is
        // exactly the kinetic energy.
        MALLOY_CHECK_NEAR(total_potential_energy(a.particles(), Vec2{}), 0.0, eps);
        MALLOY_CHECK_NEAR(total_energy(a.particles(), Vec2{}),
                          total_kinetic_energy(a.particles()), 0.0);
    }

    // --- Gravity is applied BEFORE the position update, which is what makes
    //     this semi-implicit Euler. One hand-computable step pins the order:
    //     with g = (0, -10) and dt = 0.5 the velocity becomes -5 and the
    //     position moves by -2.5, not by 0 (explicit Euler) and not by -1.25
    //     (the exact half-a-t-squared solution). ---
    {
        ParticleSettings s{1.0, big_box()};
        s.gravity = Vec2{0.0, -10.0};
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.0, 0.0}, Vec2{}, 2.0, 0.1}};
        ParticleWorld w{SimulationSettings{0.5}, s, p};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, Vec2(0.0, -5.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].position, Vec2(0.0, -2.5), eps);
    }

    // --- Gravity accelerates every particle equally regardless of mass, and
    //     both components act. ---
    {
        ParticleSettings s{1.0, big_box()};
        s.gravity = Vec2{3.0, -4.0};
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.0, 0.0}, Vec2{}, 1.0, 0.1},
            Particle2D{Vec2{20.0, 0.0}, Vec2{}, 50.0, 0.1}}; // far apart, no contact
        ParticleWorld w{SimulationSettings{0.25}, s, p};
        MALLOY_CHECK_TRUE(w.step().ok());
        MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, Vec2(0.75, -1.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particles()[1].velocity, Vec2(0.75, -1.0), eps);
    }

    // --- INVARIANT: in free flight momentum changes by exactly (sum m) * g * dt
    //     per step. Gravity is an external force, so momentum is not conserved,
    //     but its change is exact rather than approximate. ---
    {
        ParticleSettings s{1.0, big_box()};
        s.gravity = Vec2{0.6, -2.5};
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.0, 0.0}, Vec2{1.0, 2.0}, 1.5, 0.1},
            Particle2D{Vec2{30.0, 0.0}, Vec2{-0.5, 0.25}, 2.5, 0.1}};
        const Real total_mass = 4.0;
        const Real dt = 0.001;
        const int steps = 500;

        ParticleWorld w{SimulationSettings{dt}, s, p};
        const Vec2 p0 = total_momentum(p);
        for (int i = 0; i < steps; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const Vec2 expected =
            p0 + s.gravity * (total_mass * dt * static_cast<Real>(steps));
        MALLOY_CHECK_VEC2_NEAR(total_momentum(w.particles()), expected, 1e-10);
    }

    // --- INVARIANT: total energy is NOT conserved under a constant field. The
    //     documented behaviour is an exact loss of (1/2)(sum m)|g|^2 dt^2 per
    //     free-flight step, so that is what is asserted. Claiming conservation
    //     here would be claiming something false. ---
    {
        ParticleSettings s{1.0, big_box()};
        s.gravity = Vec2{0.0, -3.0};
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.0, 5.0}, Vec2{1.0, 0.0}, 2.0, 0.1}};
        const Real dt = 0.01;
        const int steps = 300;

        ParticleWorld w{SimulationSettings{dt}, s, p};
        const Real e0 = total_energy(p, s.gravity);
        for (int i = 0; i < steps; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const Real per_step = 0.5 * 2.0 * 9.0 * dt * dt;   // (1/2) m |g|^2 dt^2
        const Real expected = e0 - per_step * static_cast<Real>(steps);
        MALLOY_CHECK_NEAR(total_energy(w.particles(), s.gravity), expected, 1e-10);

        // The drift is real, not noise: it must be far larger than the tolerance.
        MALLOY_CHECK_TRUE(per_step * static_cast<Real>(steps) > 1e-6);
    }

    // --- Potential energy: higher is more, with gravity pointing down. ---
    {
        const Vec2 g{0.0, -10.0};
        const std::vector<Particle2D> high = {
            Particle2D{Vec2{0.0, 3.0}, Vec2{}, 2.0, 0.1}};
        const std::vector<Particle2D> low = {
            Particle2D{Vec2{0.0, 1.0}, Vec2{}, 2.0, 0.1}};
        MALLOY_CHECK_NEAR(total_potential_energy(high, g), 60.0, eps); // m*g*h
        MALLOY_CHECK_NEAR(total_potential_energy(low, g), 20.0, eps);
        MALLOY_CHECK_TRUE(total_potential_energy(high, g) >
                          total_potential_energy(low, g));
        // A horizontal field acts on x instead, so the sign convention is pinned
        // on both axes rather than just the one gravity usually uses.
        MALLOY_CHECK_NEAR(total_potential_energy(high, Vec2{-10.0, 0.0}), 0.0, eps);
        const std::vector<Particle2D> right = {
            Particle2D{Vec2{4.0, 0.0}, Vec2{}, 2.0, 0.1}};
        MALLOY_CHECK_NEAR(total_potential_energy(right, Vec2{-10.0, 0.0}), 80.0, eps);
    }

    // --- A ball dropped onto the floor bounces, and with restitution below 1
    //     it settles instead of bouncing forever. ---
    {
        const Aabb box{Vec2{-5.0, 0.0}, Vec2{5.0, 10.0}};
        ParticleSettings s{0.5, box};
        s.gravity = Vec2{0.0, -9.81};
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{0.0, 8.0}, Vec2{0.0, 0.0}, 1.0, 0.25}};
        ParticleWorld w{SimulationSettings{0.001}, s, p};
        for (int i = 0; i < 20000; ++i)
        {
            MALLOY_CHECK_TRUE(w.step().ok());
        }
        const Particle2D& ball = w.particles()[0];
        // Resting on the floor, within a radius of it, and barely moving.
        MALLOY_CHECK_TRUE(ball.position.y - ball.radius < box.min.y + 0.05);
        MALLOY_CHECK_TRUE(std::abs(ball.velocity.y) < 0.5);
        MALLOY_CHECK_TRUE(malloy::math::is_finite(ball.position));
    }

    // --- Non-finite gravity is rejected rather than poisoning every particle. ---
    {
        ParticleSettings bad{1.0, big_box()};
        bad.gravity = Vec2{0.0, nan};
        MALLOY_CHECK_FALSE(bad.is_valid());
        ParticleSettings bad_inf{1.0, big_box()};
        bad_inf.gravity = Vec2{inf, 0.0};
        MALLOY_CHECK_FALSE(bad_inf.is_valid());

        const std::vector<Particle2D> p = {
            Particle2D{Vec2{}, Vec2{}, 1.0, 0.1}};
        ParticleWorld w{SimulationSettings{0.01}, bad, p};
        MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
        MALLOY_CHECK_EQ(w.tick_count(), std::uint64_t{0});
    }

    // --- Settings accessors, previously untested public API. Every field is a
    //     different value so a getter returning the wrong member shows. ---
    {
        ParticleSettings s{0.375, Aabb{Vec2{-3.0, -4.0}, Vec2{5.0, 6.0}}};
        s.gravity = Vec2{0.25, -1.5};
        ParticleWorld w{SimulationSettings{0.0025}, s, {}};
        MALLOY_CHECK_NEAR(w.simulation_settings().dt, 0.0025, eps);
        MALLOY_CHECK_NEAR(w.particle_settings().restitution, 0.375, eps);
        MALLOY_CHECK_VEC2_NEAR(w.particle_settings().gravity, Vec2(0.25, -1.5), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particle_settings().bounds.min, Vec2(-3.0, -4.0), eps);
        MALLOY_CHECK_VEC2_NEAR(w.particle_settings().bounds.max, Vec2(5.0, 6.0), eps);
    }

    // --- The four-wall table above runs at restitution 1, where the
    //     multiplier is a no-op: deleting `* restitution` from the left or top
    //     branch changes nothing there. Only the right wall was ever pinned at
    //     a value other than 1. This runs all four at 0.4, with a different
    //     tangential speed in each case so a swapped component shows too. ---
    {
        const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
        struct WallCase
        {
            Vec2 position;
            Vec2 velocity;
            Vec2 expected_position;
            Vec2 expected_velocity;
        };
        const WallCase cases[] = {
            // left: -1.0 in becomes +0.4 out, and the tangential 0.25 survives
            {Vec2{-0.4, 0.0}, Vec2{-1.0, 0.25}, Vec2{-0.5, 0.125}, Vec2{0.4, 0.25}},
            // right
            {Vec2{0.4, 0.0}, Vec2{1.0, -0.75}, Vec2{0.5, -0.375}, Vec2{-0.4, -0.75}},
            // bottom
            {Vec2{0.0, -0.4}, Vec2{0.6, -1.0}, Vec2{0.3, -0.5}, Vec2{0.6, 0.4}},
            // top
            {Vec2{0.0, 0.4}, Vec2{-0.3, 1.0}, Vec2{-0.15, 0.5}, Vec2{-0.3, -0.4}},
        };
        for (const WallCase& c : cases)
        {
            const std::vector<Particle2D> p = {Particle2D{c.position, c.velocity, 1.0, 0.5}};
            ParticleWorld w{SimulationSettings{0.5}, ParticleSettings{0.4, box}, p};
            MALLOY_CHECK_TRUE(w.step().ok());
            MALLOY_CHECK_VEC2_NEAR(w.particles()[0].position, c.expected_position, eps);
            MALLOY_CHECK_VEC2_NEAR(w.particles()[0].velocity, c.expected_velocity, eps);
        }
    }

    // --- Issue #16: a position or velocity whose SQUARE overflows is refused.
    //     Finite is not enough. Everything that squares a vector reaches
    //     infinity above about 1.34e154, so accepting state up to 1.8e308 left
    //     a window in which a world validated clean while every energy it
    //     reported was inf. ---
    {
        const Real too_big = 1.4e154; // finite, and its square is not
        MALLOY_CHECK_TRUE((Particle2D{Vec2{1.0, 2.0}, Vec2{}, 1.0, 0.5}.is_valid()));
        MALLOY_CHECK_FALSE((Particle2D{Vec2{too_big, 0.0}, Vec2{}, 1.0, 0.5}.is_valid()));
        MALLOY_CHECK_FALSE((Particle2D{Vec2{}, Vec2{0.0, too_big}, 1.0, 0.5}.is_valid()));
    }

    // --- Issue #20: the wall clamp MOVES ENERGY, and this is the derivation.
    //
    //     Clamping a particle back to a wall changes its position without
    //     changing its velocity. Against gravity that adds potential energy:
    //     m |g| delta, once per clamping event, where delta is how far it had
    //     overshot. The project had characterized exactly one integrator error
    //     term, the (1/2) m |g|^2 dt^2 free-flight loss, and this one is
    //     comparable in size and opposite in sign.
    //
    //     The whole step is in closed form here, with restitution 1 so the
    //     impulse changes no kinetic energy at all:
    //
    //       gravity     v1 = v0 - g dt = -6 - 1  = -7
    //       position    y1 = y0 + v1 dt = 1 - 0.7 = 0.3
    //       clamp       surface is 0.5, so delta  = 0.2
    //
    //       E0 = 0.5*2*36 + 2*10*1.0 = 56
    //       E1 = 0.5*2*49 + 2*10*0.5 = 59
    //
    //     and the change of +3 is exactly the sum of the two terms:
    //     free flight -(1/2)*2*100*0.01 = -1, plus clamp 2*10*0.2 = +4. ---
    {
        const Real mass = 2.0;
        const Real radius = 0.5;
        const Real gravity = 10.0;
        const Real dt = 0.1;

        // Wide in x and tall above, so the floor is the only wall in play.
        const Aabb box{Vec2{-100.0, 0.0}, Vec2{100.0, 100.0}};
        ParticleSettings settings;
        settings.restitution = 1.0;
        settings.bounds = box;
        settings.gravity = Vec2{0.0, -gravity};

        const std::vector<Particle2D> start = {
            Particle2D{Vec2{0.0, 1.0}, Vec2{0.0, -6.0}, mass, radius}};

        const Real energy_before = total_energy(start, settings.gravity);
        MALLOY_CHECK_NEAR(energy_before, 56.0, eps);

        ParticleWorld w{SimulationSettings{dt}, settings, start};
        MALLOY_CHECK_TRUE(w.step().ok());

        // It really was clamped, and to the surface rather than somewhere else.
        MALLOY_CHECK_NEAR(w.particles()[0].position.y, radius, eps);
        // Restitution 1, so the speed is untouched and only the sign flipped.
        MALLOY_CHECK_NEAR(w.particles()[0].velocity.y, 7.0, eps);

        const Real energy_after = total_energy(w.particles(), settings.gravity);
        MALLOY_CHECK_NEAR(energy_after, 59.0, eps);

        const Real free_flight = -0.5 * mass * gravity * gravity * dt * dt;
        const Real clamp = mass * gravity * 0.2;
        MALLOY_CHECK_NEAR(energy_after - energy_before, free_flight + clamp, eps);

        // And the clamp term dominates by a factor of four here, which is the
        // point: it is not a rounding-level correction.
        MALLOY_CHECK_TRUE(std::abs(clamp) > 3.0 * std::abs(free_flight));
    }

    // --- Issue #20: angular momentum, which this domain had no way to report.
    //     A hand sum, and then the reason it is not conserved. ---
    {
        const std::vector<Particle2D> p = {
            Particle2D{Vec2{2.0, 0.0}, Vec2{0.0, 3.0}, 1.5, 0.1},
            Particle2D{Vec2{0.0, -4.0}, Vec2{5.0, 0.0}, 2.0, 0.1}};
        // 1.5*(2*3 - 0*0) + 2.0*(0*0 - (-4)*5) = 9 + 40 = 49
        MALLOY_CHECK_NEAR(total_angular_momentum(p), 49.0, eps);

        // Zero for a particle heading straight at the origin, whatever its
        // speed: the position and velocity are parallel.
        const std::vector<Particle2D> radial = {
            Particle2D{Vec2{3.0, 4.0}, Vec2{-6.0, -8.0}, 1.0, 0.1}};
        MALLOY_CHECK_NEAR(total_angular_momentum(radial), 0.0, eps);
    }

    // --- M32: colliding particles in three dimensions. The 3D sibling of the
    //     M10 domain: spheres collide (the M22/M24 geometry) and are confined to
    //     a 3D box. Elastic collisions conserve momentum and energy; walls carry
    //     momentum away; and a z = 0 run matches the 2D world. ---
    {
        using malloy::math::approx_equal;
        const auto particle3 = [](Vec3 position, Vec3 velocity, Real mass, Real radius) {
            Particle3D p;
            p.position = position;
            p.velocity = velocity;
            p.mass = mass;
            p.radius = radius;
            return p;
        };
        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();
        const Aabb3 room{Vec3{-100.0, -100.0, -100.0}, Vec3{100.0, 100.0, 100.0}};

        // --- Validation. ---
        {
            MALLOY_CHECK_TRUE(particle3(Vec3{}, Vec3{}, 1.0, 0.5).is_valid());
            MALLOY_CHECK_TRUE(particle3(Vec3{}, Vec3{}, 1.0, 0.0).is_valid()); // a point
            MALLOY_CHECK_FALSE(particle3(Vec3{}, Vec3{}, 0.0, 0.5).is_valid());
            MALLOY_CHECK_FALSE(particle3(Vec3{}, Vec3{}, 1.0, -0.5).is_valid());
            MALLOY_CHECK_FALSE(particle3(Vec3{0.0, 0.0, nan3}, Vec3{}, 1.0, 0.5).is_valid());
            ParticleSettings3D good;
            good.bounds = room;
            MALLOY_CHECK_TRUE(good.is_valid());
            ParticleSettings3D bad_rest;
            bad_rest.restitution = 1.5;
            MALLOY_CHECK_FALSE(bad_rest.is_valid());
            ParticleSettings3D bad_gravity;
            bad_gravity.gravity = Vec3{0.0, inf3, 0.0};
            MALLOY_CHECK_FALSE(bad_gravity.is_valid());
            ParticleSettings3D bad_bounds; // min > max on an axis
            bad_bounds.bounds = Aabb3{Vec3{1.0, 1.0, 1.0}, Vec3{-1.0, 1.0, 1.0}};
            MALLOY_CHECK_FALSE(bad_bounds.is_valid());
            // A particle wider than the box is rejected (it would fight two walls).
            ParticleSettings3D small_box;
            small_box.bounds = Aabb3{Vec3{-0.4, -0.4, -0.4}, Vec3{0.4, 0.4, 0.4}};
            ParticleWorld3D tight{SimulationSettings{0.01}, small_box,
                                  {particle3(Vec3{}, Vec3{}, 1.0, 0.5)}};
            MALLOY_CHECK_TRUE(tight.step().status == StepStatus::InvalidState);
        }

        // --- An elastic head-on collision conserves momentum AND kinetic
        //     energy, and the equal masses swap velocities. No walls, no
        //     gravity, restitution 1. ---
        {
            ParticleSettings3D settings;
            settings.restitution = 1.0;
            settings.bounds = room;
            ParticleWorld3D world{SimulationSettings{0.001}, settings,
                                  {particle3(Vec3{-2.0, 0.0, 0.0}, Vec3{3.0, 0.0, 0.0}, 1.0, 0.5),
                                   particle3(Vec3{2.0, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0}, 1.0, 0.5)}};
            const Vec3 p0 = total_momentum3d(world.particles());
            const Real e0 = total_kinetic_energy3d(world.particles());
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(approx_equal(total_momentum3d(world.particles()), p0, 1e-9));
            MALLOY_CHECK_NEAR(total_kinetic_energy3d(world.particles()), e0, 1e-9);
            // Equal masses swap along the line of centres: 3 -> -1 and -1 -> 3.
            MALLOY_CHECK_NEAR(world.particles()[0].velocity.x, -1.0, 1e-9);
            MALLOY_CHECK_NEAR(world.particles()[1].velocity.x, 3.0, 1e-9);
        }

        // --- An inelastic collision keeps momentum but sheds energy. ---
        {
            ParticleSettings3D settings;
            settings.restitution = 0.5;
            settings.bounds = room;
            ParticleWorld3D world{SimulationSettings{0.001}, settings,
                                  {particle3(Vec3{-2.0, 0.0, 0.0}, Vec3{3.0, 0.0, 0.0}, 1.0, 0.5),
                                   particle3(Vec3{2.0, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0}, 1.0, 0.5)}};
            const Vec3 p0 = total_momentum3d(world.particles());
            const Real e0 = total_kinetic_energy3d(world.particles());
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
            }
            MALLOY_CHECK_TRUE(approx_equal(total_momentum3d(world.particles()), p0, 1e-9));
            MALLOY_CHECK_TRUE(total_kinetic_energy3d(world.particles()) < e0 - 0.1);
        }

        // --- The box confines a particle, on every axis including the new z:
        //     one aimed diagonally at a corner stays inside and reverses on each
        //     wall it meets. ---
        {
            ParticleSettings3D settings;
            settings.restitution = 1.0;
            settings.bounds = Aabb3{Vec3{-1.0, -1.0, -1.0}, Vec3{1.0, 1.0, 1.0}};
            ParticleWorld3D world{SimulationSettings{0.01}, settings,
                                  {particle3(Vec3{}, Vec3{0.7, -0.5, 0.9}, 1.0, 0.2)}};
            for (int i = 0; i < 4000; ++i)
            {
                MALLOY_CHECK_TRUE(world.step().ok());
                const Vec3& x = world.particles()[0].position;
                // Never escapes the box (allowing the sphere radius).
                MALLOY_CHECK_TRUE(x.x >= -1.0 - 1e-9 && x.x <= 1.0 + 1e-9);
                MALLOY_CHECK_TRUE(x.y >= -1.0 - 1e-9 && x.y <= 1.0 + 1e-9);
                MALLOY_CHECK_TRUE(x.z >= -1.0 - 1e-9 && x.z <= 1.0 + 1e-9);
            }
            // Elastic walls preserve the speed exactly.
            MALLOY_CHECK_NEAR(malloy::math::length(world.particles()[0].velocity),
                              std::sqrt(0.49 + 0.25 + 0.81), 1e-12);
        }

        // --- It reduces to the 2D domain: the same particles in the z = 0
        //     plane, with gravity in the plane, trace the same xy path as
        //     ParticleWorld and never leave it. ---
        {
            const Real dt = 0.005;
            ParticleSettings s2;
            s2.restitution = 0.8;
            s2.gravity = Vec2{0.3, -2.0};
            s2.bounds = Aabb{Vec2{-5.0, -5.0}, Vec2{5.0, 5.0}};
            ParticleSettings3D s3;
            s3.restitution = 0.8;
            s3.gravity = Vec3{0.3, -2.0, 0.0};
            s3.bounds = Aabb3{Vec3{-5.0, -5.0, -5.0}, Vec3{5.0, 5.0, 5.0}};
            ParticleWorld w2{SimulationSettings{dt}, s2,
                             {Particle2D{Vec2{-1.0, 2.0}, Vec2{1.0, 0.0}, 1.0, 0.4},
                              Particle2D{Vec2{1.0, 2.0}, Vec2{-1.0, 0.5}, 2.0, 0.4}}};
            ParticleWorld3D w3{SimulationSettings{dt}, s3,
                               {particle3(Vec3{-1.0, 2.0, 0.0}, Vec3{1.0, 0.0, 0.0}, 1.0, 0.4),
                                particle3(Vec3{1.0, 2.0, 0.0}, Vec3{-1.0, 0.5, 0.0}, 2.0, 0.4)}};
            for (int i = 0; i < 800; ++i)
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

        // --- The diagnostics pin their factors, before any step. ---
        {
            const std::vector<Particle3D> two = {
                particle3(Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 2.0, 0.0}, 3.0, 0.5),
                particle3(Vec3{0.0, 0.0, 4.0}, Vec3{-1.0, 0.0, 0.0}, 2.0, 0.5)};
            // p = 3*(0,2,0) + 2*(-1,0,0) = (-2, 6, 0).
            MALLOY_CHECK_TRUE(approx_equal(total_momentum3d(two), Vec3{-2.0, 6.0, 0.0}, 1e-12));
            // KE = 0.5*3*4 + 0.5*2*1 = 6 + 1 = 7.
            MALLOY_CHECK_NEAR(total_kinetic_energy3d(two), 7.0, 1e-12);
            // U = sum of -m (g . r). With g = (0,0,-10) the only nonzero term is
            // the second body at z = 4: -2 * (-10 * 4) = 80.
            MALLOY_CHECK_NEAR(total_potential_energy3d(two, Vec3{0.0, 0.0, -10.0}), 80.0, 1e-12);
        }
    }

    // --- M35: Coulomb friction for particle contacts, in 2D and 3D. A
    //     tangential impulse clamped to `friction` times the normal impulse damps
    //     the slide; a particle has no orientation, so it imparts no spin. It
    //     applies to both particle/particle contacts and wall bounces. Friction
    //     defaults to zero, which is the old normal-only behaviour, so every
    //     earlier test in this file still passes unchanged. ---
    {
        using malloy::math::approx_equal;
        const auto particle3 = [](Vec3 position, Vec3 velocity, Real mass, Real radius) {
            Particle3D q;
            q.position = position;
            q.velocity = velocity;
            q.mass = mass;
            q.radius = radius;
            return q;
        };
        const Real inf3 = std::numeric_limits<Real>::infinity();
        const Real nan3 = std::numeric_limits<Real>::quiet_NaN();
        const Real dt = 1e-5; // tiny, so the pre-contact move barely tilts the normal

        // Validation: a negative or non-finite friction is invalid, both
        // dimensions.
        {
            ParticleSettings s{1.0, big_box()};
            MALLOY_CHECK_TRUE(s.is_valid()); // friction 0 by default
            s.friction = 0.5;
            MALLOY_CHECK_TRUE(s.is_valid());
            s.friction = -0.1;
            MALLOY_CHECK_FALSE(s.is_valid());
            s.friction = nan3;
            MALLOY_CHECK_FALSE(s.is_valid());
            s.friction = inf3;
            MALLOY_CHECK_FALSE(s.is_valid());

            ParticleSettings3D t;
            t.friction = -1.0;
            MALLOY_CHECK_FALSE(t.is_valid());
            t.friction = 0.3;
            MALLOY_CHECK_TRUE(t.is_valid());
        }

        // 2D particle/particle: two equal masses meet with the normal along x and
        // a tangential (y) relative velocity of -2. The normal impulse is 2 and
        // the impulse that would fully stop the slide is 1.
        {
            const std::vector<Particle2D> pp = {
                Particle2D{Vec2{-0.4, 0.0}, Vec2{1.0, 1.0}, 1.0, 0.5},
                Particle2D{Vec2{0.4, 0.0}, Vec2{-1.0, -1.0}, 1.0, 0.5}};

            // Frictionless: the tangential relative velocity is untouched.
            ParticleWorld none{SimulationSettings{dt}, ParticleSettings{1.0, big_box()}, pp};
            MALLOY_CHECK_TRUE(none.step().ok());
            MALLOY_CHECK_NEAR(none.particles()[1].velocity.y - none.particles()[0].velocity.y,
                              -2.0, 1e-3);
            MALLOY_CHECK_VEC2_NEAR(total_momentum(none.particles()), Vec2(0.0, 0.0), 1e-12);

            // Strong friction (the cone does not bind) drives the tangential
            // relative velocity to zero, and never reverses it. Momentum, being
            // internal, is still conserved.
            ParticleSettings strong{1.0, big_box()};
            strong.friction = 10.0;
            ParticleWorld stuck{SimulationSettings{dt}, strong, pp};
            MALLOY_CHECK_TRUE(stuck.step().ok());
            MALLOY_CHECK_NEAR(stuck.particles()[1].velocity.y - stuck.particles()[0].velocity.y,
                              0.0, 1e-3);
            MALLOY_CHECK_VEC2_NEAR(total_momentum(stuck.particles()), Vec2(0.0, 0.0), 1e-12);

            // Weak friction binds at the cone: friction 0.1 caps the tangential
            // impulse at 0.1 * 2 = 0.2, so the relative velocity falls from -2 to
            // -1.6, not to zero.
            ParticleSettings weak{1.0, big_box()};
            weak.friction = 0.1;
            ParticleWorld capped{SimulationSettings{dt}, weak, pp};
            MALLOY_CHECK_TRUE(capped.step().ok());
            MALLOY_CHECK_NEAR(capped.particles()[1].velocity.y - capped.particles()[0].velocity.y,
                              -1.6, 1e-3);
        }

        // 2D wall friction: a particle sliding into the floor keeps its normal
        // bounce but sheds tangential speed, clamped to the cone.
        {
            const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
            const std::vector<Particle2D> one = {
                Particle2D{Vec2{0.0, -0.95}, Vec2{2.0, -1.0}, 1.0, 0.1}};

            ParticleWorld none{SimulationSettings{dt}, ParticleSettings{1.0, box}, one};
            MALLOY_CHECK_TRUE(none.step().ok());
            MALLOY_CHECK_NEAR(none.particles()[0].velocity.x, 2.0, 1e-3); // tangent untouched
            MALLOY_CHECK_NEAR(none.particles()[0].velocity.y, 1.0, 1e-3); // reflected, e = 1

            // Strong friction: cap is friction*(1+e)*|v_n| = 10*2*1 = 20, far
            // above the tangential speed 2, so the slide is fully removed and NOT
            // reversed (the cap is a floor at zero, not a subtraction). The normal
            // is untouched.
            ParticleSettings strong{1.0, box};
            strong.friction = 10.0;
            ParticleWorld stuck{SimulationSettings{dt}, strong, one};
            MALLOY_CHECK_TRUE(stuck.step().ok());
            MALLOY_CHECK_NEAR(stuck.particles()[0].velocity.x, 0.0, 1e-3);
            MALLOY_CHECK_NEAR(stuck.particles()[0].velocity.y, 1.0, 1e-3);

            // Weak friction: cap is 0.25*2*1 = 0.5, so 2 falls to 1.5.
            ParticleSettings weak{1.0, box};
            weak.friction = 0.25;
            ParticleWorld capped{SimulationSettings{dt}, weak, one};
            MALLOY_CHECK_TRUE(capped.step().ok());
            MALLOY_CHECK_NEAR(capped.particles()[0].velocity.x, 1.5, 1e-3);
            MALLOY_CHECK_NEAR(capped.particles()[0].velocity.y, 1.0, 1e-3);
        }

        // 3D particle/particle: the same as the 2D pair, in space. Friction
        // removes the tangential relative velocity and conserves momentum.
        {
            const Aabb3 room{Vec3{-100.0, -100.0, -100.0}, Vec3{100.0, 100.0, 100.0}};
            const std::vector<Particle3D> pp = {
                particle3(Vec3{-0.4, 0.0, 0.0}, Vec3{1.0, 1.0, 0.0}, 1.0, 0.5),
                particle3(Vec3{0.4, 0.0, 0.0}, Vec3{-1.0, -1.0, 0.0}, 1.0, 0.5)};

            ParticleSettings3D none;
            none.restitution = 1.0;
            none.bounds = room;
            ParticleWorld3D free_w{SimulationSettings{dt}, none, pp};
            MALLOY_CHECK_TRUE(free_w.step().ok());
            MALLOY_CHECK_NEAR(free_w.particles()[1].velocity.y - free_w.particles()[0].velocity.y,
                              -2.0, 1e-3);

            ParticleSettings3D strong;
            strong.restitution = 1.0;
            strong.bounds = room;
            strong.friction = 10.0;
            ParticleWorld3D stuck{SimulationSettings{dt}, strong, pp};
            MALLOY_CHECK_TRUE(stuck.step().ok());
            MALLOY_CHECK_NEAR(stuck.particles()[1].velocity.y - stuck.particles()[0].velocity.y,
                              0.0, 1e-3);
            MALLOY_CHECK_TRUE(
                approx_equal(total_momentum3d(stuck.particles()), Vec3{0.0, 0.0, 0.0}, 1e-12));
        }

        // 3D wall friction: a particle sliding into the floor (y-min) loses both
        // in-plane tangential axes (x and z) together, clamped to the cone. With
        // tangential speed sqrt(5) and cap 2, the pair is scaled toward zero.
        {
            const Aabb3 box{Vec3{-1.0, -1.0, -1.0}, Vec3{1.0, 1.0, 1.0}};
            const std::vector<Particle3D> one = {
                particle3(Vec3{0.0, -0.95, 0.0}, Vec3{2.0, -1.0, 1.0}, 1.0, 0.1)};

            ParticleSettings3D none;
            none.restitution = 1.0;
            none.bounds = box;
            ParticleWorld3D free_w{SimulationSettings{dt}, none, one};
            MALLOY_CHECK_TRUE(free_w.step().ok());
            MALLOY_CHECK_NEAR(free_w.particles()[0].velocity.x, 2.0, 1e-3);
            MALLOY_CHECK_NEAR(free_w.particles()[0].velocity.z, 1.0, 1e-3);
            MALLOY_CHECK_NEAR(free_w.particles()[0].velocity.y, 1.0, 1e-3); // reflected

            // Strong friction: cap 10*2*1 = 20 exceeds the tangential speed
            // sqrt(5), so BOTH in-plane axes are driven to zero together, not
            // past it. This pins the two-axis clamp as a floor at zero.
            ParticleSettings3D strong;
            strong.restitution = 1.0;
            strong.bounds = box;
            strong.friction = 10.0;
            ParticleWorld3D stuck{SimulationSettings{dt}, strong, one};
            MALLOY_CHECK_TRUE(stuck.step().ok());
            MALLOY_CHECK_NEAR(stuck.particles()[0].velocity.x, 0.0, 1e-3);
            MALLOY_CHECK_NEAR(stuck.particles()[0].velocity.z, 0.0, 1e-3);
            MALLOY_CHECK_NEAR(stuck.particles()[0].velocity.y, 1.0, 1e-3);

            // Weak friction binds at the cone: cap 0.5*2*1 = 1 is below sqrt(5),
            // so the (x, z) tangent is scaled toward zero by (sqrt5 - 1)/sqrt5,
            // keeping its direction.
            ParticleSettings3D weak;
            weak.restitution = 1.0;
            weak.bounds = box;
            weak.friction = 0.5;
            ParticleWorld3D capped{SimulationSettings{dt}, weak, one};
            MALLOY_CHECK_TRUE(capped.step().ok());
            const Real s5 = std::sqrt(5.0);
            const Real scale = (s5 - 1.0) / s5; // cap 1 removed from a length-sqrt5 tangent
            MALLOY_CHECK_NEAR(capped.particles()[0].velocity.x, 2.0 * scale, 1e-3);
            MALLOY_CHECK_NEAR(capped.particles()[0].velocity.z, 1.0 * scale, 1e-3);
            MALLOY_CHECK_NEAR(capped.particles()[0].velocity.y, 1.0, 1e-3);
        }
    }

    // --- M38: kinetic theory. The gas quantities are pure functions of the
    //     state and of the wall tally, and the tally is exact bookkeeping that
    //     changes no trajectory. ---
    {

        // --- The statistics in closed form, on hand-built lists. The masses
        //     differ, which pins temperature as mass-weighted and the two
        //     velocity statistics as unweighted. ---
        {
            const std::vector<Particle2D> two = {
                Particle2D{Vec2{}, Vec2{3.0, 4.0}, 2.0, 0.1}, // v^2 25, KE 25
                Particle2D{Vec2{}, Vec2{1.0, 0.0}, 4.0, 0.1}, // v^2 1, KE 2
            };
            // T = 2 KE / (d N) = 2 * 27 / (2 * 2)
            MALLOY_CHECK_NEAR(temperature(two), 13.5, 1e-12);
            // <v^4> / <v^2>^2 = ((625 + 1) / 2) / 13^2
            MALLOY_CHECK_NEAR(speed_ratio(two), 313.0 / 169.0, 1e-12);
            // the mean v^2 is 13, and one of the two is below it
            MALLOY_CHECK_NEAR(slow_fraction(two), 0.5, 1e-12);

            const std::vector<Particle3D> three = {
                Particle3D{Vec3{}, Vec3{1.0, 2.0, 2.0}, 2.0, 0.1}, // v^2 9, KE 9
                Particle3D{Vec3{}, Vec3{0.0, 0.0, 1.0}, 1.0, 0.1}, // v^2 1, KE 0.5
                Particle3D{Vec3{}, Vec3{0.0, 1.0, 1.0}, 1.0, 0.1}, // v^2 2, KE 1
            };
            // T = 2 * 10.5 / (3 * 3); <v^4> = 86/3 over <v^2>^2 = 16; 2 of 3 below 4
            MALLOY_CHECK_NEAR(temperature3d(three), 21.0 / 9.0, 1e-12);
            MALLOY_CHECK_NEAR(speed_ratio3d(three), 86.0 / 48.0, 1e-12);
            MALLOY_CHECK_NEAR(slow_fraction3d(three), 2.0 / 3.0, 1e-12);

            // Empty and motionless lists give 0 rather than dividing by zero.
            MALLOY_CHECK_EQ(temperature(std::vector<Particle2D>{}), 0.0);
            MALLOY_CHECK_EQ(speed_ratio(std::vector<Particle2D>{}), 0.0);
            MALLOY_CHECK_EQ(slow_fraction(std::vector<Particle2D>{}), 0.0);
            MALLOY_CHECK_EQ(temperature3d(std::vector<Particle3D>{}), 0.0);
            MALLOY_CHECK_EQ(speed_ratio3d(std::vector<Particle3D>{}), 0.0);
            MALLOY_CHECK_EQ(slow_fraction3d(std::vector<Particle3D>{}), 0.0);
            const std::vector<Particle2D> still = {
                Particle2D{Vec2{}, Vec2{}, 1.0, 0.1},
                Particle2D{Vec2{1.0, 0.0}, Vec2{}, 2.0, 0.1},
            };
            MALLOY_CHECK_EQ(temperature(still), 0.0);
            MALLOY_CHECK_EQ(speed_ratio(still), 0.0);
            MALLOY_CHECK_EQ(slow_fraction(still), 0.0); // nothing is below a mean of 0
        }

        // --- Pressure is total impulse over elapsed time and total wall: the
        //     perimeter in 2D, the surface area in 3D. Every wall's impulse is
        //     distinct, so a wall left out of the sum shows. Before any time has
        //     passed there is no pressure at all, not a 0. ---
        {
            const Aabb box{Vec2{-1.0, -2.0}, Vec2{1.0, 1.0}}; // 2 by 3, perimeter 10
            const WallImpulse2D impulse{1.0, 2.0, 3.0, 4.0};   // total 10
            const auto p = pressure(impulse, box, 2.0);
            MALLOY_CHECK_TRUE(p.has_value());
            MALLOY_CHECK_NEAR(*p, 10.0 / (2.0 * 10.0), 1e-15);
            MALLOY_CHECK_FALSE(pressure(impulse, box, 0.0).has_value());
            MALLOY_CHECK_FALSE(pressure(impulse, box, -1.0).has_value());
            MALLOY_CHECK_FALSE(pressure(impulse, box, nan).has_value());
            MALLOY_CHECK_FALSE(pressure(impulse, box, inf).has_value());

            const Aabb3 box3{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 2.0, 3.0}}; // area 2(2+6+3)
            const WallImpulse3D impulse3{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};  // total 21
            const auto p3 = pressure3d(impulse3, box3, 0.5);
            MALLOY_CHECK_TRUE(p3.has_value());
            MALLOY_CHECK_NEAR(*p3, 21.0 / (0.5 * 22.0), 1e-15);
            MALLOY_CHECK_FALSE(pressure3d(impulse3, box3, 0.0).has_value());
            MALLOY_CHECK_FALSE(pressure3d(impulse3, box3, -1.0).has_value());
            MALLOY_CHECK_FALSE(pressure3d(impulse3, box3, nan).has_value());
            MALLOY_CHECK_FALSE(pressure3d(impulse3, box3, inf).has_value());
        }

        // --- Each bounce adds m (1 + e) times the incoming normal speed to the
        //     wall it hit and to no other. One particle per wall, each with a
        //     distinct mass and speed, at restitution 0.5. Friction is on, and
        //     its tangential impulse is NOT tallied: it pushes along a wall,
        //     not into it, so the normal figures are unchanged by it. ---
        {
            const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
            ParticleSettings settings;
            settings.restitution = 0.5;
            settings.bounds = box;
            settings.friction = 0.3;
            const std::vector<Particle2D> four = {
                Particle2D{Vec2{-0.85, -0.5}, Vec2{-3.0, 0.5}, 2.0, 0.1}, // min_x 2*1.5*3
                Particle2D{Vec2{0.85, 0.5}, Vec2{1.0, -0.5}, 3.0, 0.1},   // max_x 3*1.5*1
                Particle2D{Vec2{0.5, -0.85}, Vec2{0.5, -2.0}, 1.0, 0.1},  // min_y 1*1.5*2
                Particle2D{Vec2{-0.5, 0.85}, Vec2{-0.5, 4.0}, 1.0, 0.1},  // max_y 1*1.5*4
            };
            ParticleWorld w{SimulationSettings{0.1}, settings, four};
            MALLOY_CHECK_TRUE(w.step().ok());
            MALLOY_CHECK_NEAR(w.wall_impulse().min_x, 9.0, 1e-12);
            MALLOY_CHECK_NEAR(w.wall_impulse().max_x, 4.5, 1e-12);
            MALLOY_CHECK_NEAR(w.wall_impulse().min_y, 3.0, 1e-12);
            MALLOY_CHECK_NEAR(w.wall_impulse().max_y, 6.0, 1e-12);

            // 3D: one particle driven into a corner hits three walls in one
            // step, each tallied on its own axis (friction 0, so no wall's
            // bounce alters the normal speed another one sees).
            ParticleSettings3D settings3;
            settings3.restitution = 0.5;
            settings3.bounds = Aabb3{Vec3{-1.0, -1.0, -1.0}, Vec3{1.0, 1.0, 1.0}};
            ParticleWorld3D w3{SimulationSettings{0.1}, settings3,
                               {Particle3D{Vec3{0.85, -0.85, 0.85}, Vec3{1.0, -2.0, 3.0}, 2.0, 0.1}}};
            MALLOY_CHECK_TRUE(w3.step().ok());
            MALLOY_CHECK_NEAR(w3.wall_impulse().max_x, 3.0, 1e-12); // 2 * 1.5 * 1
            MALLOY_CHECK_NEAR(w3.wall_impulse().min_y, 6.0, 1e-12); // 2 * 1.5 * 2
            MALLOY_CHECK_NEAR(w3.wall_impulse().max_z, 9.0, 1e-12); // 2 * 1.5 * 3
            MALLOY_CHECK_EQ(w3.wall_impulse().min_x, 0.0);
            MALLOY_CHECK_EQ(w3.wall_impulse().max_y, 0.0);
            MALLOY_CHECK_EQ(w3.wall_impulse().min_z, 0.0);
        }

        // --- One point particle between the x walls, in closed form. dt = 1/64
        //     and speed 1/2 put every position on a multiple of 1/128, exact in
        //     binary. Landing exactly ON a wall is not a bounce; the next step
        //     overshoots and is clamped back, so the particle spends one step
        //     there and a crossing of the 2-wide box takes 256 + 1 steps. The
        //     max-x wall is hit at steps 129, 643, ..., the min-x wall at 386,
        //     each hit delivering 2 m v = 3. That lag is the clamp's: the mean
        //     force is m v^2 / (L + v dt), not m v^2 / L. ---
        {
            const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
            ParticleWorld w{SimulationSettings{1.0 / 64.0}, ParticleSettings{1.0, box},
                            {Particle2D{Vec2{0.0, 0.0}, Vec2{0.5, 0.0}, 3.0, 0.0}}};
            const auto run_to = [&w](int step) {
                while (static_cast<int>(w.tick_count()) < step)
                {
                    if (!w.step().ok())
                    {
                        return false;
                    }
                }
                return true;
            };
            MALLOY_CHECK_TRUE(run_to(128));
            MALLOY_CHECK_EQ(w.wall_impulse().max_x, 0.0);
            MALLOY_CHECK_TRUE(run_to(129));
            MALLOY_CHECK_EQ(w.wall_impulse().max_x, 3.0);
            MALLOY_CHECK_TRUE(run_to(385));
            MALLOY_CHECK_EQ(w.wall_impulse().min_x, 0.0);
            MALLOY_CHECK_TRUE(run_to(386));
            MALLOY_CHECK_EQ(w.wall_impulse().min_x, 3.0);
            MALLOY_CHECK_TRUE(run_to(642));
            MALLOY_CHECK_EQ(w.wall_impulse().max_x, 3.0);
            MALLOY_CHECK_TRUE(run_to(643));
            MALLOY_CHECK_EQ(w.wall_impulse().max_x, 6.0);
            MALLOY_CHECK_EQ(w.wall_impulse().min_y + w.wall_impulse().max_y, 0.0);
            // 9 of impulse over 643/64 of time and a perimeter of 8.
            MALLOY_CHECK_NEAR(*pressure(w.wall_impulse(), box, w.elapsed_time()),
                              72.0 / 643.0, 1e-15);
            // One period is 514 steps: 2 m v per period on the max-x wall is
            // exactly m v^2 / (L + v dt).
            MALLOY_CHECK_NEAR(3.0 / (514.0 / 64.0), 3.0 * 0.25 / (2.0 + 0.5 / 64.0), 1e-15);
        }

        // --- Momentum balance, exact: walls are the only external agent here
        //     (no gravity, friction 0), and each wall's tally is the momentum it
        //     gave, so p(t) - p(0) = (J_minx - J_maxx, J_miny - J_maxy). Unequal
        //     masses and radii, a non-square box and restitution 0.9 keep it
        //     asymmetric. Run twice: the tally is as repeatable as the state. ---
        {
            const Aabb box{Vec2{-3.0, -2.0}, Vec2{3.0, 2.0}};
            std::vector<Particle2D> gas = lattice_gas(box, 6, 4, 1.0, 0.15, 1.0);
            for (std::size_t i = 0; i < gas.size(); i += 3)
            {
                gas[i].mass = 2.5;
                gas[i].radius = 0.25;
            }
            const Vec2 p0 = total_momentum(gas);
            ParticleWorld a{SimulationSettings{0.01}, ParticleSettings{0.9, box}, gas};
            ParticleWorld b{SimulationSettings{0.01}, ParticleSettings{0.9, box}, gas};
            for (int i = 0; i < 1500; ++i)
            {
                MALLOY_CHECK_TRUE(a.step().ok());
                MALLOY_CHECK_TRUE(b.step().ok());
            }
            const WallImpulse2D& j = a.wall_impulse();
            MALLOY_CHECK_TRUE(j.min_x > 0.0 && j.max_x > 0.0 && j.min_y > 0.0 && j.max_y > 0.0);
            const Vec2 dp = total_momentum(a.particles()) - p0;
            MALLOY_CHECK_NEAR(dp.x, j.min_x - j.max_x, 1e-9);
            MALLOY_CHECK_NEAR(dp.y, j.min_y - j.max_y, 1e-9);
            MALLOY_CHECK_EQ(j.min_x, b.wall_impulse().min_x);
            MALLOY_CHECK_EQ(j.max_x, b.wall_impulse().max_x);
            MALLOY_CHECK_EQ(j.min_y, b.wall_impulse().min_y);
            MALLOY_CHECK_EQ(j.max_y, b.wall_impulse().max_y);
            MALLOY_CHECK_EQ(a.particles()[5].position.x, b.particles()[5].position.x);
        }
        {
            // 3D, now WITH gravity, which adds M g t to the balance: every step
            // adds m g dt to every velocity, walls or not.
            const Aabb3 box{Vec3{-2.0, -1.5, -1.0}, Vec3{2.0, 1.5, 1.0}};
            std::vector<Particle3D> gas = lattice_gas3d(box, 4, 3, 2, 1.0, 0.15, 1.0);
            Real total_mass = 0.0;
            for (std::size_t i = 0; i < gas.size(); ++i)
            {
                if (i % 3 == 0)
                {
                    gas[i].mass = 2.5;
                    gas[i].radius = 0.25;
                }
                total_mass += gas[i].mass;
            }
            ParticleSettings3D settings;
            settings.restitution = 0.9;
            settings.bounds = box;
            settings.gravity = Vec3{0.2, -0.5, 0.3};
            const Vec3 p0 = total_momentum3d(gas);
            ParticleWorld3D w{SimulationSettings{0.01}, settings, gas};
            for (int i = 0; i < 1500; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            const WallImpulse3D& j = w.wall_impulse();
            const Vec3 pushed = settings.gravity * (total_mass * w.elapsed_time());
            const Vec3 dp = total_momentum3d(w.particles()) - p0;
            MALLOY_CHECK_NEAR(dp.x, j.min_x - j.max_x + pushed.x, 1e-9);
            MALLOY_CHECK_NEAR(dp.y, j.min_y - j.max_y + pushed.y, 1e-9);
            MALLOY_CHECK_NEAR(dp.z, j.min_z - j.max_z + pushed.z, 1e-9);
            MALLOY_CHECK_TRUE(j.min_z > 0.0 && j.max_z > 0.0);
        }

        // --- Point particles never collide (a contact needs the centres within
        //     r_a + r_b = 0), so each bounces between the walls on its own and
        //     the gas is exactly ideal: P A = N T = KE in 2D, P V = N T =
        //     (2/3) KE in 3D. The error is bounded, not fitted: a particle's
        //     hits on a pair of walls are within one of their mean count, each
        //     hit carrying 2 m |v_axis|, and the clamp lag lowers each rate by
        //     at most v dt / L. So
        //       |ratio - 1| <= sum(2 m |v_axis|) L / (4 t KE) + v dt / L. ---
        {
            const Real side = 1.0;
            const Aabb box{Vec2{0.0, 0.0}, Vec2{side, side}};
            const std::vector<Particle2D> points = lattice_gas(box, 4, 4, 1.0, 0.0, 1.0);
            ParticleWorld w{SimulationSettings{0.01}, ParticleSettings{1.0, box}, points};
            for (int i = 0; i < 20000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            Real lumps = 0.0;
            for (const Particle2D& p : points)
            {
                lumps += 2.0 * p.mass * (std::abs(p.velocity.x) + std::abs(p.velocity.y));
            }
            const Real ke = total_kinetic_energy(points);
            const Real t = w.elapsed_time();
            const Real bound = lumps * side / (4.0 * t * ke) + 1.0 * 0.01 / side;
            const Real ratio = *pressure(w.wall_impulse(), box, t) * side * side / ke;
            MALLOY_CHECK_TRUE(bound < 0.02); // the test is only meaningful if tight
            MALLOY_CHECK_NEAR(ratio, 1.0, bound);
        }
        {
            const Real side = 1.0;
            const Aabb3 box{Vec3{0.0, 0.0, 0.0}, Vec3{side, side, side}};
            const std::vector<Particle3D> points = lattice_gas3d(box, 4, 2, 2, 1.0, 0.0, 1.0);
            ParticleWorld3D w{SimulationSettings{0.01}, ParticleSettings3D{1.0, box}, points};
            for (int i = 0; i < 20000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            Real lumps = 0.0;
            for (const Particle3D& p : points)
            {
                lumps += 2.0 * p.mass *
                         (std::abs(p.velocity.x) + std::abs(p.velocity.y) + std::abs(p.velocity.z));
            }
            const Real ke = total_kinetic_energy3d(points);
            const Real t = w.elapsed_time();
            const Real bound = lumps * side / (4.0 * t * ke) + 1.0 * 0.01 / side;
            const Real ratio =
                *pressure3d(w.wall_impulse(), box, t) * side * side * side / (2.0 / 3.0 * ke);
            MALLOY_CHECK_TRUE(bound < 0.025);
            MALLOY_CHECK_NEAR(ratio, 1.0, bound);
        }

        // --- Equipartition. A gas started moving almost entirely along x (a
        //     small y component lets rows meet off-axis) shares its energy
        //     until KE_x / KE is 1/2. At equilibrium that ratio is Beta(N/2,
        //     N/2), standard deviation 1 / (2 sqrt(N + 1)); the check allows
        //     three of those. Temperature is conserved throughout (elastic). ---
        {
            const Aabb box{Vec2{-4.0, -4.0}, Vec2{4.0, 4.0}};
            std::vector<Particle2D> gas = lattice_gas(box, 8, 8, 1.0, 0.2, 1.0);
            for (std::size_t i = 0; i < gas.size(); ++i)
            {
                const Real along = i % 2 == 0 ? 1.0 : -1.0;
                const Real across = (i / 2) % 2 == 0 ? 0.2 : -0.2;
                gas[i].velocity = Vec2{along, across};
            }
            const auto x_share = [](const std::vector<Particle2D>& list) {
                Real kx = 0.0;
                for (const Particle2D& p : list)
                {
                    kx += 0.5 * p.mass * p.velocity.x * p.velocity.x;
                }
                return kx / total_kinetic_energy(list);
            };
            MALLOY_CHECK_TRUE(x_share(gas) > 0.95);
            const Real t0 = temperature(gas);
            ParticleWorld w{SimulationSettings{0.01}, ParticleSettings{1.0, box}, gas};
            for (int i = 0; i < 4000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            const Real n = static_cast<Real>(gas.size());
            MALLOY_CHECK_NEAR(x_share(w.particles()), 0.5, 3.0 / (2.0 * std::sqrt(n + 1.0)));
            MALLOY_CHECK_NEAR(temperature(w.particles()), t0, 1e-12 * t0);
            // Relaxed: the slow fraction is the 2D Maxwell-Boltzmann 1 - 1/e,
            // within three binomial standard deviations for this N.
            const Real expected = 1.0 - std::exp(-1.0);
            MALLOY_CHECK_NEAR(slow_fraction(w.particles()), expected,
                              3.0 * std::sqrt(expected * (1.0 - expected) / n));
        }

        // --- Equipartition across masses: a mixture of masses 1 and 4 started
        //     at the SAME speed holds four times the energy in its heavy half.
        //     Collisions bring the mean kinetic energy per particle of the two
        //     species to the same value, the temperature both share. A snapshot
        //     of each half (32 disks, 64 degrees of freedom) scatters by
        //     sqrt(2/64) relatively, so their ratio by about 0.25; averaging
        //     over the second half of the run, with at least 10 independent
        //     samples (one per two collision times), brings three standard
        //     deviations to 0.24. ---
        {
            const Aabb box{Vec2{-4.0, -4.0}, Vec2{4.0, 4.0}};
            std::vector<Particle2D> gas = lattice_gas(box, 8, 8, 1.0, 0.2, 1.0);
            for (std::size_t i = 0; i < gas.size(); i += 2)
            {
                gas[i].mass = 4.0;
            }
            const auto species_ratio = [](const std::vector<Particle2D>& list) {
                Real heavy = 0.0;
                Real light = 0.0;
                for (const Particle2D& p : list)
                {
                    const Real ke = 0.5 * p.mass * malloy::math::dot(p.velocity, p.velocity);
                    (p.mass > 2.0 ? heavy : light) += ke;
                }
                return std::pair<Real, Real>{heavy, light};
            };
            const auto [heavy0, light0] = species_ratio(gas);
            MALLOY_CHECK_NEAR(heavy0 / light0, 4.0, 1e-12);
            ParticleWorld w{SimulationSettings{0.01}, ParticleSettings{1.0, box}, gas};
            Real heavy_sum = 0.0;
            Real light_sum = 0.0;
            for (int i = 1; i <= 4000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                if (i > 2000 && i % 20 == 0)
                {
                    const auto [heavy, light] = species_ratio(w.particles());
                    heavy_sum += heavy;
                    light_sum += light;
                }
            }
            MALLOY_CHECK_NEAR(heavy_sum / light_sum, 1.0, 0.24);
        }

        // --- Pressure does not hinge on the timestep: halving dt (and doubling
        //     the steps) moves it by less than the run-to-run scatter allows.
        //     Different dt means different collisions, so this is a statistical
        //     comparison: three standard deviations of the difference (about
        //     0.005 each at this size) plus 0.01 for any first-order bias. ---
        {
            const Aabb box{Vec2{-4.0, -4.0}, Vec2{4.0, 4.0}};
            const std::vector<Particle2D> gas = lattice_gas(box, 8, 8, 1.0, 0.2, 1.0);
            ParticleWorld coarse{SimulationSettings{0.01}, ParticleSettings{1.0, box}, gas};
            ParticleWorld fine{SimulationSettings{0.005}, ParticleSettings{1.0, box}, gas};
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(coarse.step().ok());
            }
            for (int i = 0; i < 4000; ++i)
            {
                MALLOY_CHECK_TRUE(fine.step().ok());
            }
            const Real p_coarse = *pressure(coarse.wall_impulse(), box, coarse.elapsed_time());
            const Real p_fine = *pressure(fine.wall_impulse(), box, fine.elapsed_time());
            MALLOY_CHECK_NEAR(p_coarse, p_fine, 0.03);
            MALLOY_CHECK_TRUE(p_coarse > 0.55); // and clearly above the ideal 0.5
        }

        // --- A step rolled back for non-finite state rolls the tally back too.
        //     A huge (finite) gravity makes the velocity unsquarable in one
        //     step, AFTER the particle has bounced off the max-x wall in that
        //     same step; the bounce must not survive the rollback. ---
        {
            const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
            ParticleSettings settings;
            settings.bounds = box;
            settings.gravity = Vec2{1e300, 0.0};
            ParticleWorld w{SimulationSettings{1.0}, settings,
                            {Particle2D{Vec2{0.5, 0.0}, Vec2{1.0, 0.0}, 1.0, 0.1}}};
            MALLOY_CHECK_TRUE(w.validate() == StepStatus::Ok);
            MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
            MALLOY_CHECK_EQ(w.wall_impulse().max_x, 0.0);
            MALLOY_CHECK_EQ(w.particles()[0].position.x, 0.5);

            ParticleSettings3D settings3;
            settings3.bounds = Aabb3{Vec3{-1.0, -1.0, -1.0}, Vec3{1.0, 1.0, 1.0}};
            settings3.gravity = Vec3{0.0, 0.0, 1e300};
            ParticleWorld3D w3{SimulationSettings{1.0}, settings3,
                               {Particle3D{Vec3{0.0, 0.0, 0.5}, Vec3{0.0, 0.0, 1.0}, 1.0, 0.1}}};
            MALLOY_CHECK_TRUE(w3.validate() == StepStatus::Ok);
            MALLOY_CHECK_TRUE(w3.step().status == StepStatus::InvalidState);
            MALLOY_CHECK_EQ(w3.wall_impulse().max_z, 0.0);
            MALLOY_CHECK_EQ(w3.particles()[0].position.z, 0.5);
        }

        // --- An empty world: steps fine, delivers no impulse, and after time
        //     has passed its pressure exists and is 0. ---
        {
            const Aabb box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
            ParticleWorld w{SimulationSettings{0.1}, ParticleSettings{1.0, box}, {}};
            MALLOY_CHECK_FALSE(pressure(w.wall_impulse(), box, w.elapsed_time()).has_value());
            MALLOY_CHECK_TRUE(w.step().ok());
            const auto p = pressure(w.wall_impulse(), box, w.elapsed_time());
            MALLOY_CHECK_TRUE(p.has_value());
            MALLOY_CHECK_EQ(*p, 0.0);
            MALLOY_CHECK_EQ(temperature(w.particles()), 0.0);
        }

        // --- The lattice builder: cell centres in row-major order (x fastest),
        //     one mass, radius and speed, and directions in exactly opposite
        //     pairs, pair k at k golden angles, so the momentum is exactly zero
        //     for an even count and one particle's for an odd one. ---
        {
            const Aabb box{Vec2{-2.0, -1.0}, Vec2{2.0, 1.0}}; // 1 by 1 cells
            const std::vector<Particle2D> grid = lattice_gas(box, 4, 2, 1.5, 0.2, 2.0);
            MALLOY_CHECK_EQ(grid.size(), std::size_t{8});
            MALLOY_CHECK_VEC2_NEAR(grid[0].position, (Vec2{-1.5, -0.5}), 1e-15);
            MALLOY_CHECK_VEC2_NEAR(grid[1].position, (Vec2{-0.5, -0.5}), 1e-15);
            MALLOY_CHECK_VEC2_NEAR(grid[5].position, (Vec2{-0.5, 0.5}), 1e-15);
            MALLOY_CHECK_VEC2_NEAR(grid[7].position, (Vec2{1.5, 0.5}), 1e-15);
            for (std::size_t i = 0; i < grid.size(); ++i)
            {
                MALLOY_CHECK_EQ(grid[i].mass, 1.5);
                MALLOY_CHECK_EQ(grid[i].radius, 0.2);
                MALLOY_CHECK_NEAR(malloy::math::length(grid[i].velocity), 2.0, 1e-12);
                if (i % 2 == 1)
                {
                    MALLOY_CHECK_EQ(grid[i].velocity.x, -grid[i - 1].velocity.x);
                    MALLOY_CHECK_EQ(grid[i].velocity.y, -grid[i - 1].velocity.y);
                }
            }
            MALLOY_CHECK_EQ(total_momentum(grid).x, 0.0);
            MALLOY_CHECK_EQ(total_momentum(grid).y, 0.0);
            // pair 0 points along +x; pair 1 a golden angle, about 137.5 degrees, on
            const Real golden = 3.14159265358979323846 * (3.0 - std::sqrt(5.0));
            MALLOY_CHECK_VEC2_NEAR(grid[0].velocity, (Vec2{2.0, 0.0}), 1e-15);
            MALLOY_CHECK_VEC2_NEAR(grid[2].velocity,
                                   (Vec2{2.0 * std::cos(golden), 2.0 * std::sin(golden)}), 1e-12);
            // An odd count leaves the last particle unpaired.
            const std::vector<Particle2D> odd = lattice_gas(box, 3, 1, 1.0, 0.2, 1.0);
            MALLOY_CHECK_EQ(odd.size(), std::size_t{3});
            MALLOY_CHECK_VEC2_NEAR(total_momentum(odd), odd[2].velocity, 1e-15);
            // Degenerate counts build nothing.
            MALLOY_CHECK_TRUE(lattice_gas(box, 0, 3, 1.0, 0.2, 1.0).empty());
            MALLOY_CHECK_TRUE(lattice_gas(box, 3, -1, 1.0, 0.2, 1.0).empty());
            // Directions cover the circle evenly: on a 20 by 20 lattice the x
            // share of the kinetic energy is within 0.01 of a half.
            const Aabb wide{Vec2{-10.0, -10.0}, Vec2{10.0, 10.0}};
            const std::vector<Particle2D> many = lattice_gas(wide, 20, 20, 1.0, 0.1, 1.0);
            Real kx = 0.0;
            for (const Particle2D& p : many)
            {
                kx += 0.5 * p.velocity.x * p.velocity.x;
            }
            MALLOY_CHECK_NEAR(kx / total_kinetic_energy(many), 0.5, 0.01);
        }
        {
            const Aabb3 box{Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 3.0, 4.0}}; // 1 by 1 by 1 cells
            const std::vector<Particle3D> grid = lattice_gas3d(box, 2, 3, 4, 1.5, 0.2, 2.0);
            MALLOY_CHECK_EQ(grid.size(), std::size_t{24});
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(grid[0].position, Vec3{0.5, 0.5, 0.5}, 1e-15));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(grid[1].position, Vec3{1.5, 0.5, 0.5}, 1e-15));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(grid[2].position, Vec3{0.5, 1.5, 0.5}, 1e-15));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(grid[6].position, Vec3{0.5, 0.5, 1.5}, 1e-15));
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(grid[23].position, Vec3{1.5, 2.5, 3.5}, 1e-15));
            for (std::size_t i = 0; i < grid.size(); ++i)
            {
                MALLOY_CHECK_EQ(grid[i].mass, 1.5);
                MALLOY_CHECK_EQ(grid[i].radius, 0.2);
                MALLOY_CHECK_NEAR(malloy::math::length(grid[i].velocity), 2.0, 1e-12);
                if (i % 2 == 1)
                {
                    MALLOY_CHECK_EQ(grid[i].velocity.x, -grid[i - 1].velocity.x);
                    MALLOY_CHECK_EQ(grid[i].velocity.y, -grid[i - 1].velocity.y);
                    MALLOY_CHECK_EQ(grid[i].velocity.z, -grid[i - 1].velocity.z);
                }
            }
            const Vec3 p = total_momentum3d(grid);
            MALLOY_CHECK_TRUE(p.x == 0.0 && p.y == 0.0 && p.z == 0.0);
            // The documented directions: pair k at azimuth k golden angles and
            // height z = 1 - 2 frac((k + 1/2) sqrt 2). Pinned for pairs 0 and 1,
            // because a skewed z sequence would still look isotropic once each
            // pair is mirrored, so the energy shares below cannot catch it.
            const Real golden3 = 3.14159265358979323846 * (3.0 - std::sqrt(5.0));
            for (int k = 0; k < 2; ++k)
            {
                const Real s = (k + 0.5) * std::sqrt(2.0);
                const Real z = 1.0 - 2.0 * (s - std::floor(s));
                const Real ring = std::sqrt(1.0 - z * z);
                const Vec3 expected{2.0 * ring * std::cos(k * golden3),
                                    2.0 * ring * std::sin(k * golden3), 2.0 * z};
                MALLOY_CHECK_TRUE(malloy::math::approx_equal(grid[2 * k].velocity, expected, 1e-12));
            }
            const std::vector<Particle3D> odd = lattice_gas3d(box, 1, 1, 3, 1.0, 0.2, 1.0);
            MALLOY_CHECK_TRUE(malloy::math::approx_equal(total_momentum3d(odd), odd[2].velocity, 1e-15));
            MALLOY_CHECK_TRUE(lattice_gas3d(box, 2, 2, 0, 1.0, 0.2, 1.0).empty());
            // Directions cover the sphere evenly: on a 10 by 10 by 10 lattice
            // each axis holds a third of the kinetic energy to within 0.005.
            const Aabb3 cube{Vec3{0.0, 0.0, 0.0}, Vec3{10.0, 10.0, 10.0}};
            const std::vector<Particle3D> many = lattice_gas3d(cube, 10, 10, 10, 1.0, 0.1, 1.0);
            Vec3 share{};
            for (const Particle3D& q : many)
            {
                share += Vec3{q.velocity.x * q.velocity.x, q.velocity.y * q.velocity.y,
                              q.velocity.z * q.velocity.z};
            }
            const Real all = share.x + share.y + share.z;
            MALLOY_CHECK_NEAR(share.x / all, 1.0 / 3.0, 0.005);
            MALLOY_CHECK_NEAR(share.y / all, 1.0 / 3.0, 0.005);
            MALLOY_CHECK_NEAR(share.z / all, 1.0 / 3.0, 0.005);
        }
    }

    // --- M39: thermal walls. A wall held at a temperature re-emits each particle
    //     it meets from that temperature's distribution, drawing from the
    //     world's seeded generator, so heat flows between the walls and the
    //     gas. The generator and samplers first, then the bookkeeping, which is
    //     exact, then the physics, which is statistical with derived tolerances
    //     (docs/05). The physics centres and tolerances were measured across
    //     40 to 60 seeds by an independent prototype before this was built. ---
    {
        using malloy::particles::emission_normal_speed;
        using malloy::particles::emission_tangent;
        using malloy::particles::SplitMix64;
        using malloy::particles::WallTemperature2D;
        using malloy::particles::WallTemperature3D;
        const Real pi = 3.14159265358979323846;
        const std::uint64_t gamma = SplitMix64::gamma;

        // The generator: SplitMix64 is pinned by known answers (they match an
        // independent implementation, Java's SplittableRandom), at compile
        // time and at run time, and its uniform deviate stays strictly inside
        // (0, 1). After n draws its state is exactly seed + n gamma.
        {
            static_assert([] {
                SplitMix64 g{0};
                return g.next();
            }() == 0xe220a8397b1dcdafULL);
            static_assert(SplitMix64::to_unit_open(0) == 0x1p-53);
            static_assert(SplitMix64::to_unit_open(~std::uint64_t{0}) == 1.0 - 0x1p-53);
            SplitMix64 zero{0};
            MALLOY_CHECK_TRUE(zero.next() == 0xe220a8397b1dcdafULL);
            MALLOY_CHECK_TRUE(zero.next() == 0x6e789e6aa1b965f4ULL);
            MALLOY_CHECK_TRUE(zero.next() == 0x06c45d188009454fULL);
            SplitMix64 other{1234567};
            MALLOY_CHECK_TRUE(other.next() == 0x599ed017fb08fc85ULL);
            MALLOY_CHECK_TRUE(other.next() == 0x2c73f08458540fa5ULL);
            MALLOY_CHECK_TRUE(other.next() == 0x883ebce5a3f27c77ULL);
            SplitMix64 counted{42};
            for (int i = 0; i < 1000; ++i)
            {
                const Real u = counted.uniform();
                MALLOY_CHECK_TRUE(u > 0.0 && u < 1.0);
            }
            MALLOY_CHECK_TRUE(counted.state() == 42 + 1000 * gamma);
        }

        // The samplers' distributions, a million emissions at T = 2, m = 0.5
        // (sigma^2 = T/m = 4; m != 1 so a missing division by m shows). Each
        // tolerance is three standard errors of the mean in question. The
        // normal speed is Rayleigh: mean sigma sqrt(pi/2) = 2.5066 (a
        // half-Gaussian would give sigma sqrt(2/pi) = 1.596), mean square
        // 2 sigma^2, and v^2 below 2T/m with probability 1 - 1/e (0.843 for a
        // half-Gaussian). The tangential pair is two independent N(0, sigma^2).
        {
            SplitMix64 g{7};
            const Real t = 2.0;
            const Real m = 0.5;
            const int count = 1000000;
            Real normal = 0.0;
            Real normal_sq = 0.0;
            Real tangent = 0.0;
            Real tangent_sq = 0.0;
            Real product = 0.0;
            int slow = 0;
            bool all_positive = true;
            for (int i = 0; i < count; ++i)
            {
                const Real vn = emission_normal_speed(g, t, m);
                const Vec2 vt = emission_tangent(g, t, m);
                all_positive = all_positive && vn > 0.0 && std::isfinite(vn);
                normal += vn;
                normal_sq += vn * vn;
                slow += vn * vn < 2.0 * t / m ? 1 : 0;
                tangent += vt.x + vt.y;
                tangent_sq += vt.x * vt.x + vt.y * vt.y;
                product += vt.x * vt.y;
            }
            const Real k = static_cast<Real>(count);
            MALLOY_CHECK_TRUE(all_positive);
            MALLOY_CHECK_NEAR(normal / k, 2.0 * std::sqrt(pi / 2.0), 0.0039); // sd 1.310
            MALLOY_CHECK_NEAR(normal_sq / k, 8.0, 0.024);                      // sd 8
            MALLOY_CHECK_NEAR(slow / k, 1.0 - std::exp(-1.0), 0.00145);        // binomial
            MALLOY_CHECK_NEAR(tangent / (2.0 * k), 0.0, 0.0043);               // sd 2
            MALLOY_CHECK_NEAR(tangent_sq / (2.0 * k), 4.0, 0.012);             // sd 5.66
            MALLOY_CHECK_NEAR(product / k, 0.0, 0.012);                        // sd 4
        }

        // The samplers' exact values, not just their distributions: a rewrite
        // that keeps the law (log(1 - u) for log(u), cos and sin swapped, u1
        // and u2 swapped) would still change every thermal run, and ADR 0010
        // makes that a deliberate act. The tolerance is relative 1e-12, which
        // absorbs a last-bit difference in log, cos or sin between C runtimes.
        {
            SplitMix64 g{0};
            const Real golden[2][3] = {
                {0.99630577286538602, 2.5572672056599703, 0.42868116209283502},
                {0.48621283223457967, -1.9772082492486569, 3.7442027606630823}};
            for (const auto& row : golden)
            {
                const Real vn = emission_normal_speed(g, 2.0, 0.5);
                const Vec2 vt = emission_tangent(g, 2.0, 0.5);
                MALLOY_CHECK_NEAR(vn, row[0], 1e-12 * std::abs(row[0]));
                MALLOY_CHECK_NEAR(vt.x, row[1], 1e-12 * std::abs(row[1]));
                MALLOY_CHECK_NEAR(vt.y, row[2], 1e-12 * std::abs(row[2]));
            }
        }

        // Validation: every wall's temperature is checked on its own, so a
        // validator that looked only at one wall fails. NaN, an infinity and a
        // negative are refused; -0.0 is allowed and, not above 0, specular.
        {
            Real WallTemperature2D::* const walls2[] = {
                &WallTemperature2D::min_x, &WallTemperature2D::max_x,
                &WallTemperature2D::min_y, &WallTemperature2D::max_y};
            for (Real WallTemperature2D::* const wall : walls2)
            {
                for (const Real bad : {nan, inf, -1.0})
                {
                    ParticleSettings s;
                    s.wall_temperature.*wall = bad;
                    MALLOY_CHECK_FALSE(s.is_valid());
                }
                ParticleSettings fine;
                fine.wall_temperature.*wall = -0.0;
                MALLOY_CHECK_TRUE(fine.is_valid());
                fine.wall_temperature.*wall = 3.5;
                MALLOY_CHECK_TRUE(fine.is_valid());
            }
            Real WallTemperature3D::* const walls3[] = {
                &WallTemperature3D::min_x, &WallTemperature3D::max_x,
                &WallTemperature3D::min_y, &WallTemperature3D::max_y,
                &WallTemperature3D::min_z, &WallTemperature3D::max_z};
            for (Real WallTemperature3D::* const wall : walls3)
            {
                for (const Real bad : {nan, inf, -1.0})
                {
                    ParticleSettings3D s;
                    s.wall_temperature.*wall = bad;
                    MALLOY_CHECK_FALSE(s.is_valid());
                }
                ParticleSettings3D fine;
                fine.wall_temperature.*wall = -0.0;
                MALLOY_CHECK_TRUE(fine.is_valid());
            }
            // And a world with an invalid temperature refuses to step.
            ParticleSettings hot_nan;
            hot_nan.wall_temperature.max_y = nan;
            ParticleWorld w{SimulationSettings{0.01}, hot_nan,
                            {Particle2D{Vec2{}, Vec2{1.0, 0.0}, 1.0, 0.1}}};
            MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidSettings);
        }

        // No thermal wall, no draw: a box full of bounces leaves the stream at
        // its seed, and two different seeds give bit-identical runs.
        // Determinism: with thermal walls the same seed repeats the run exactly
        // and a different seed does not, and the stream has moved by a whole
        // number of thermal bounces, three draws each.
        {
            const Aabb box{Vec2{-2.0, -2.0}, Vec2{2.0, 2.0}};
            const std::vector<Particle2D> gas = lattice_gas(box, 4, 4, 1.0, 0.1, 1.0);
            ParticleSettings plain{1.0, box};
            plain.seed = 1;
            ParticleSettings plain_other = plain;
            plain_other.seed = 0xDEADBEEFULL;
            ParticleWorld a{SimulationSettings{0.01}, plain, gas};
            ParticleWorld b{SimulationSettings{0.01}, plain_other, gas};
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(a.step().ok());
                MALLOY_CHECK_TRUE(b.step().ok());
            }
            MALLOY_CHECK_TRUE(a.random_state() == 1);
            MALLOY_CHECK_TRUE(b.random_state() == 0xDEADBEEFULL);
            MALLOY_CHECK_TRUE(a.wall_impulse().min_x > 0.0); // it did bounce
            for (std::size_t i = 0; i < gas.size(); ++i)
            {
                MALLOY_CHECK_TRUE(a.particles()[i].position.x == b.particles()[i].position.x);
                MALLOY_CHECK_TRUE(a.particles()[i].velocity.y == b.particles()[i].velocity.y);
            }
            MALLOY_CHECK_TRUE(a.wall_heat().min_x == 0.0 && a.wall_heat().max_y == 0.0);

            ParticleSettings warm = plain;
            warm.wall_temperature = WallTemperature2D{2.0, 0.5, 0.0, 1.0};
            warm.seed = 99;
            ParticleSettings warm_other = warm;
            warm_other.seed = 100;
            ParticleWorld c{SimulationSettings{0.01}, warm, gas};
            ParticleWorld d{SimulationSettings{0.01}, warm, gas};
            ParticleWorld e{SimulationSettings{0.01}, warm_other, gas};
            for (int i = 0; i < 2000; ++i)
            {
                MALLOY_CHECK_TRUE(c.step().ok());
                MALLOY_CHECK_TRUE(d.step().ok());
                MALLOY_CHECK_TRUE(e.step().ok());
            }
            MALLOY_CHECK_TRUE(c.random_state() == d.random_state());
            bool same = true;
            bool differs = false;
            for (std::size_t i = 0; i < gas.size(); ++i)
            {
                same = same && c.particles()[i].position.x == d.particles()[i].position.x &&
                       c.particles()[i].velocity.y == d.particles()[i].velocity.y;
                differs = differs || c.particles()[i].position.x != e.particles()[i].position.x;
            }
            MALLOY_CHECK_TRUE(same);
            MALLOY_CHECK_TRUE(differs);
            MALLOY_CHECK_TRUE(c.wall_heat().min_x == d.wall_heat().min_x);
            MALLOY_CHECK_TRUE(c.wall_impulse().max_y == d.wall_impulse().max_y);
            MALLOY_CHECK_TRUE(c.wall_heat().min_y == 0.0); // the specular wall
            // Draws so far, n = (state - seed) / gamma mod 2^64, via gamma's
            // inverse (Newton's iteration doubles the correct low bits each time).
            std::uint64_t inverse = gamma;
            for (int i = 0; i < 6; ++i)
            {
                inverse *= 2 - gamma * inverse;
            }
            MALLOY_CHECK_TRUE(gamma * inverse == 1);
            const std::uint64_t draws = (c.random_state() - 99) * inverse;
            MALLOY_CHECK_TRUE(draws > 0 && draws % 3 == 0 && draws < 3000000);
        }

        // One thermal bounce, in closed form against the public samplers: the
        // world draws the normal speed, then the tangential pair, from its
        // stream (three draws), points the normal speed into the box, and puts
        // the particle on the wall. The tallies get exactly m (|v_in| + v_out)
        // and the change in kinetic energy, computed the same way. Restitution
        // and friction do not touch a thermal bounce.
        {
            const Aabb box{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
            for (const bool damped : {false, true})
            {
                ParticleSettings s{1.0, box};
                s.wall_temperature.min_x = 2.0;
                s.seed = 2024;
                if (damped)
                {
                    s.restitution = 0.3;
                    s.friction = 0.7;
                }
                ParticleWorld w{SimulationSettings{0.1}, s,
                                {Particle2D{Vec2{0.05, 0.5}, Vec2{-1.0, 0.3}, 0.5, 0.0}}};
                MALLOY_CHECK_TRUE(w.step().ok());
                SplitMix64 copy{2024};
                const Real vn = emission_normal_speed(copy, 2.0, 0.5);
                const Vec2 vt = emission_tangent(copy, 2.0, 0.5);
                const Particle2D& p = w.particles()[0];
                MALLOY_CHECK_TRUE(p.velocity.x == vn); // into the box, +x
                MALLOY_CHECK_TRUE(p.velocity.y == vt.x);
                MALLOY_CHECK_TRUE(p.position.x == 0.0); // on the wall
                MALLOY_CHECK_TRUE(w.random_state() == copy.state());
                MALLOY_CHECK_TRUE(copy.state() == 2024 + 3 * gamma);
                MALLOY_CHECK_TRUE(w.wall_impulse().min_x == 0.5 * (1.0 + vn));
                const Real before = 0.5 * 0.5 * (1.0 * 1.0 + 0.3 * 0.3);
                const Real after = 0.5 * 0.5 * (vn * vn + vt.x * vt.x);
                MALLOY_CHECK_NEAR(w.wall_heat().min_x, after - before, 1e-15);
                MALLOY_CHECK_TRUE(w.wall_heat().max_x == 0.0 && w.wall_heat().min_y == 0.0);
            }
            // A max wall emits toward -x, and a y wall sets x as the tangent.
            ParticleSettings s{1.0, box};
            s.wall_temperature.max_x = 1.5;
            s.wall_temperature.min_y = 0.25;
            s.seed = 77;
            ParticleWorld w{SimulationSettings{0.1}, s,
                            {Particle2D{Vec2{0.97, 0.5}, Vec2{1.0, 0.0}, 2.0, 0.0},
                             Particle2D{Vec2{0.5, 0.04}, Vec2{0.2, -1.0}, 3.0, 0.0}}};
            MALLOY_CHECK_TRUE(w.step().ok());
            SplitMix64 copy{77};
            const Real vn0 = emission_normal_speed(copy, 1.5, 2.0);
            const Vec2 vt0 = emission_tangent(copy, 1.5, 2.0);
            const Real vn1 = emission_normal_speed(copy, 0.25, 3.0);
            const Vec2 vt1 = emission_tangent(copy, 0.25, 3.0);
            MALLOY_CHECK_TRUE(w.particles()[0].velocity.x == -vn0);
            MALLOY_CHECK_TRUE(w.particles()[0].velocity.y == vt0.x);
            MALLOY_CHECK_TRUE(w.particles()[1].velocity.y == vn1);
            MALLOY_CHECK_TRUE(w.particles()[1].velocity.x == vt1.x);
            MALLOY_CHECK_TRUE(w.random_state() == copy.state());
            MALLOY_CHECK_TRUE(w.wall_impulse().max_x == 2.0 * (1.0 + vn0));
            MALLOY_CHECK_TRUE(w.wall_impulse().min_y == 3.0 * (1.0 + vn1));
        }
        {
            // 3D: an x wall sets the tangents y then z, a z wall x then y.
            const Aabb3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};
            ParticleSettings3D s;
            s.bounds = box;
            s.wall_temperature.min_x = 2.0;
            s.wall_temperature.max_z = 0.5;
            s.seed = 31337;
            ParticleWorld3D w{SimulationSettings{0.1}, s,
                              {Particle3D{Vec3{0.05, 0.5, 0.5}, Vec3{-1.0, 0.2, 0.1}, 0.5, 0.0},
                               Particle3D{Vec3{0.5, 0.5, 0.97}, Vec3{0.1, 0.0, 1.0}, 2.0, 0.0}}};
            MALLOY_CHECK_TRUE(w.step().ok());
            SplitMix64 copy{31337};
            const Real vn0 = emission_normal_speed(copy, 2.0, 0.5);
            const Vec2 vt0 = emission_tangent(copy, 2.0, 0.5);
            const Real vn1 = emission_normal_speed(copy, 0.5, 2.0);
            const Vec2 vt1 = emission_tangent(copy, 0.5, 2.0);
            const Vec3 v0 = w.particles()[0].velocity;
            const Vec3 v1 = w.particles()[1].velocity;
            MALLOY_CHECK_TRUE(v0.x == vn0 && v0.y == vt0.x && v0.z == vt0.y);
            MALLOY_CHECK_TRUE(v1.x == vt1.x && v1.y == vt1.y && v1.z == -vn1);
            MALLOY_CHECK_TRUE(w.random_state() == copy.state());
            MALLOY_CHECK_TRUE(copy.state() == 31337 + 6 * gamma);
            MALLOY_CHECK_TRUE(w.wall_impulse().min_x == 0.5 * (1.0 + vn0));
            MALLOY_CHECK_TRUE(w.wall_impulse().max_z == 2.0 * (1.0 + vn1));
            const Real before = 0.5 * 0.5 * (1.0 * 1.0 + 0.2 * 0.2 + 0.1 * 0.1);
            const Real after = 0.5 * 0.5 * (vn0 * vn0 + vt0.x * vt0.x + vt0.y * vt0.y);
            MALLOY_CHECK_NEAR(w.wall_heat().min_x, after - before, 1e-15);
            MALLOY_CHECK_TRUE(w.wall_heat().max_x == 0.0 && w.wall_heat().min_z == 0.0);
        }

        // Every wall at once: one particle into each wall in the same step,
        // each wall at its own temperature and each particle of its own mass,
        // checked component by component, and every impulse and heat slot,
        // against the samplers drawn in particle order. Each wall has its own
        // call site, so this pins the normal sign, the tangent axes and the
        // tally slot at every one of them, not just the few above.
        {
            const Aabb box{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
            ParticleSettings s{1.0, box};
            s.wall_temperature = WallTemperature2D{1.1, 1.3, 1.7, 1.9};
            s.seed = 404;
            const std::vector<Particle2D> movers = {
                Particle2D{Vec2{0.05, 0.4}, Vec2{-1.0, 0.25}, 1.5, 0.0},  // into min_x
                Particle2D{Vec2{0.95, 0.6}, Vec2{1.0, -0.5}, 2.5, 0.0},   // into max_x
                Particle2D{Vec2{0.3, 0.05}, Vec2{0.75, -1.0}, 3.5, 0.0},  // into min_y
                Particle2D{Vec2{0.7, 0.95}, Vec2{-0.125, 1.0}, 4.5, 0.0}, // into max_y
            };
            ParticleWorld w{SimulationSettings{0.1}, s, movers};
            MALLOY_CHECK_TRUE(w.step().ok());
            SplitMix64 copy{404};
            const Real temps[4] = {1.1, 1.3, 1.7, 1.9};
            Real impulse[4]{};
            Real heat[4]{};
            for (int i = 0; i < 4; ++i)
            {
                const Particle2D& in = movers[i];
                const Real vn = emission_normal_speed(copy, temps[i], in.mass);
                const Vec2 vt = emission_tangent(copy, temps[i], in.mass);
                const Real sign = i % 2 == 0 ? 1.0 : -1.0; // min walls emit +, max walls -
                const Vec2 out = i < 2 ? Vec2{sign * vn, vt.x} : Vec2{vt.x, sign * vn};
                const Vec2 got = w.particles()[static_cast<std::size_t>(i)].velocity;
                MALLOY_CHECK_TRUE(got.x == out.x && got.y == out.y);
                const Real speed_in = std::abs(i < 2 ? in.velocity.x : in.velocity.y);
                impulse[i] = in.mass * (speed_in + vn);
                heat[i] = 0.5 * in.mass * malloy::math::dot(out, out) -
                          0.5 * in.mass * malloy::math::dot(in.velocity, in.velocity);
            }
            MALLOY_CHECK_TRUE(w.random_state() == copy.state());
            const auto& j = w.wall_impulse();
            const auto& q = w.wall_heat();
            MALLOY_CHECK_TRUE(j.min_x == impulse[0] && j.max_x == impulse[1] &&
                              j.min_y == impulse[2] && j.max_y == impulse[3]);
            MALLOY_CHECK_NEAR(q.min_x, heat[0], 1e-14);
            MALLOY_CHECK_NEAR(q.max_x, heat[1], 1e-14);
            MALLOY_CHECK_NEAR(q.min_y, heat[2], 1e-14);
            MALLOY_CHECK_NEAR(q.max_y, heat[3], 1e-14);
        }
        {
            // 3D: each wall's tangents in ascending axis order (an x wall sets
            // y then z, a y wall x then z, a z wall x then y).
            const Aabb3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};
            ParticleSettings3D s;
            s.bounds = box;
            s.wall_temperature = WallTemperature3D{1.1, 1.3, 1.7, 1.9, 2.3, 2.9};
            s.seed = 505;
            const std::vector<Particle3D> movers = {
                Particle3D{Vec3{0.05, 0.4, 0.3}, Vec3{-1.0, 0.25, 0.5}, 1.5, 0.0},
                Particle3D{Vec3{0.95, 0.6, 0.7}, Vec3{1.0, -0.5, 0.25}, 2.5, 0.0},
                Particle3D{Vec3{0.3, 0.05, 0.4}, Vec3{0.75, -1.0, -0.25}, 3.5, 0.0},
                Particle3D{Vec3{0.7, 0.95, 0.6}, Vec3{-0.125, 1.0, 0.375}, 4.5, 0.0},
                Particle3D{Vec3{0.4, 0.3, 0.05}, Vec3{0.5, 0.25, -1.0}, 5.5, 0.0},
                Particle3D{Vec3{0.6, 0.7, 0.95}, Vec3{-0.25, 0.125, 1.0}, 6.5, 0.0},
            };
            ParticleWorld3D w{SimulationSettings{0.1}, s, movers};
            MALLOY_CHECK_TRUE(w.step().ok());
            SplitMix64 copy{505};
            const Real temps[6] = {1.1, 1.3, 1.7, 1.9, 2.3, 2.9};
            Real impulse[6]{};
            Real heat[6]{};
            for (int i = 0; i < 6; ++i)
            {
                const Particle3D& in = movers[static_cast<std::size_t>(i)];
                const Real vn = emission_normal_speed(copy, temps[i], in.mass);
                const Vec2 vt = emission_tangent(copy, temps[i], in.mass);
                const Real n = (i % 2 == 0 ? 1.0 : -1.0) * vn;
                const Vec3 out = i < 2   ? Vec3{n, vt.x, vt.y}
                                 : i < 4 ? Vec3{vt.x, n, vt.y}
                                         : Vec3{vt.x, vt.y, n};
                const Vec3 got = w.particles()[static_cast<std::size_t>(i)].velocity;
                MALLOY_CHECK_TRUE(got.x == out.x && got.y == out.y && got.z == out.z);
                const Real speed_in = std::abs(i < 2   ? in.velocity.x
                                               : i < 4 ? in.velocity.y
                                                       : in.velocity.z);
                impulse[i] = in.mass * (speed_in + vn);
                heat[i] = 0.5 * in.mass * malloy::math::dot(out, out) -
                          0.5 * in.mass * malloy::math::dot(in.velocity, in.velocity);
            }
            MALLOY_CHECK_TRUE(w.random_state() == copy.state());
            const auto& j = w.wall_impulse();
            const auto& q = w.wall_heat();
            MALLOY_CHECK_TRUE(j.min_x == impulse[0] && j.max_x == impulse[1] &&
                              j.min_y == impulse[2] && j.max_y == impulse[3] &&
                              j.min_z == impulse[4] && j.max_z == impulse[5]);
            const Real got_heat[6] = {q.min_x, q.max_x, q.min_y, q.max_y, q.min_z, q.max_z};
            for (int i = 0; i < 6; ++i)
            {
                MALLOY_CHECK_NEAR(got_heat[i], heat[i], 1e-14);
            }
        }

        // A corner: a particle crossing a thermal x wall and a thermal y wall
        // in the same step meets the x wall first (the documented order). The x
        // emission sets its velocity; if the new y component still points into
        // the y wall, the y wall emits too, from the next three draws;
        // otherwise the particle is only clamped there.
        {
            const Aabb box{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
            int emitted_twice = 0;
            for (std::uint64_t seed = 1; seed <= 8; ++seed)
            {
                ParticleSettings s{1.0, box};
                s.wall_temperature = WallTemperature2D{1.2, 0.0, 2.7, 0.0};
                s.seed = seed;
                ParticleWorld w{SimulationSettings{0.1}, s,
                                {Particle2D{Vec2{0.05, 0.03}, Vec2{-1.0, -0.9}, 1.25, 0.0}}};
                MALLOY_CHECK_TRUE(w.step().ok());
                SplitMix64 copy{seed};
                const Real vn = emission_normal_speed(copy, 1.2, 1.25);
                const Vec2 vt = emission_tangent(copy, 1.2, 1.25);
                Vec2 expected{vn, vt.x};
                if (expected.y < 0.0) // still heading into min_y
                {
                    const Real vn2 = emission_normal_speed(copy, 2.7, 1.25);
                    const Vec2 vt2 = emission_tangent(copy, 2.7, 1.25);
                    expected = Vec2{vt2.x, vn2};
                    ++emitted_twice;
                }
                const Particle2D& p = w.particles()[0];
                MALLOY_CHECK_TRUE(p.velocity.x == expected.x && p.velocity.y == expected.y);
                MALLOY_CHECK_TRUE(p.position.x == 0.0 && p.position.y == 0.0);
                MALLOY_CHECK_TRUE(w.random_state() == copy.state());
            }
            MALLOY_CHECK_TRUE(emitted_twice > 0 && emitted_twice < 8); // both branches ran
        }

        // A rolled-back step rolls the stream and the heat tally back too.
        // Particle A bounces off a thermal max_x wall in step 1, a good step
        // that draws and tallies. In step 2 particle B reaches min_x, whose
        // T/m overflows (1e308 / 1e-10), so its emission is infinite and the
        // step is refused: every tally, the stream, the positions and the tick
        // are exactly as step 1 left them, and a retry fails the same way.
        {
            const Aabb box{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
            ParticleSettings s{1.0, box};
            s.wall_temperature.min_x = 1e308;
            s.wall_temperature.max_x = 1.0;
            s.seed = 5;
            ParticleWorld w{SimulationSettings{0.1}, s,
                            {Particle2D{Vec2{0.95, 0.5}, Vec2{1.0, 0.0}, 1.0, 0.0},
                             Particle2D{Vec2{0.15, 0.3}, Vec2{-1.0, 0.0}, 1e-10, 0.0}}};
            MALLOY_CHECK_TRUE(w.validate() == StepStatus::Ok);
            MALLOY_CHECK_TRUE(w.step().ok());
            const std::uint64_t state = w.random_state();
            const Real heat = w.wall_heat().max_x;
            const Real impulse = w.wall_impulse().max_x;
            const Vec2 a = w.particles()[0].position;
            const Vec2 b = w.particles()[1].position;
            MALLOY_CHECK_TRUE(state == 5 + 3 * gamma);
            MALLOY_CHECK_TRUE(heat != 0.0);
            for (int attempt = 0; attempt < 2; ++attempt)
            {
                MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
                MALLOY_CHECK_TRUE(w.random_state() == state);
                MALLOY_CHECK_TRUE(w.wall_heat().max_x == heat);
                MALLOY_CHECK_TRUE(w.wall_heat().min_x == 0.0);
                MALLOY_CHECK_TRUE(w.wall_impulse().max_x == impulse);
                MALLOY_CHECK_TRUE(w.wall_impulse().min_x == 0.0);
                MALLOY_CHECK_TRUE(w.particles()[0].position.x == a.x);
                MALLOY_CHECK_TRUE(w.particles()[1].position.x == b.x);
                MALLOY_CHECK_TRUE(w.tick_count() == 1);
            }
        }
        {
            // The same in 3D: A bounces off a thermal max_z wall in step 1; B,
            // of mass 1e-10, reaches a min_y wall at 1e308 in step 2.
            const Aabb3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};
            ParticleSettings3D s;
            s.bounds = box;
            s.wall_temperature.min_y = 1e308;
            s.wall_temperature.max_z = 1.0;
            s.seed = 6;
            ParticleWorld3D w{SimulationSettings{0.1}, s,
                              {Particle3D{Vec3{0.5, 0.5, 0.95}, Vec3{0.0, 0.0, 1.0}, 1.0, 0.0},
                               Particle3D{Vec3{0.3, 0.15, 0.4}, Vec3{0.0, -1.0, 0.0}, 1e-10, 0.0}}};
            MALLOY_CHECK_TRUE(w.step().ok());
            const std::uint64_t state = w.random_state();
            const Real heat = w.wall_heat().max_z;
            const Real impulse = w.wall_impulse().max_z;
            const Vec3 a = w.particles()[0].position;
            const Vec3 b = w.particles()[1].position;
            MALLOY_CHECK_TRUE(state == 6 + 3 * gamma);
            MALLOY_CHECK_TRUE(heat != 0.0);
            for (int attempt = 0; attempt < 2; ++attempt)
            {
                MALLOY_CHECK_TRUE(w.step().status == StepStatus::InvalidState);
                MALLOY_CHECK_TRUE(w.random_state() == state);
                MALLOY_CHECK_TRUE(w.wall_heat().max_z == heat);
                MALLOY_CHECK_TRUE(w.wall_heat().min_y == 0.0);
                MALLOY_CHECK_TRUE(w.wall_impulse().max_z == impulse);
                MALLOY_CHECK_TRUE(w.wall_impulse().min_y == 0.0);
                MALLOY_CHECK_TRUE(w.particles()[0].position.z == a.z);
                MALLOY_CHECK_TRUE(w.particles()[1].position.y == b.y);
                MALLOY_CHECK_TRUE(w.tick_count() == 1);
            }
        }

        // The energy books balance exactly: with restitution 1, no friction and
        // no gravity, only a thermal wall changes the gas's kinetic energy, so
        // KE(t) - KE(0) is the sum of the wall heats to rounding. A specular
        // wall's heat is exactly 0. Disks of unequal mass and radius, a hot and
        // a cold wall and a third thermal wall, 1500 steps. And the x momentum
        // balances against the x walls' impulse tallies when only x walls are
        // thermal and the others specular and frictionless.
        {
            const Aabb box{Vec2{-3.0, -2.0}, Vec2{3.0, 2.0}};
            std::vector<Particle2D> gas = lattice_gas(box, 6, 4, 1.0, 0.15, 1.0);
            for (std::size_t i = 0; i < gas.size(); i += 3)
            {
                gas[i].mass = 2.5;
                gas[i].radius = 0.25;
            }
            ParticleSettings s{1.0, box};
            s.wall_temperature = WallTemperature2D{3.0, 0.5, 0.0, 1.5};
            s.seed = 11;
            ParticleWorld w{SimulationSettings{0.01}, s, gas};
            const Real ke0 = total_kinetic_energy(gas);
            for (int i = 0; i < 1500; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            const auto& q = w.wall_heat();
            MALLOY_CHECK_NEAR(total_kinetic_energy(w.particles()) - ke0,
                              q.min_x + q.max_x + q.min_y + q.max_y, 1e-9);
            MALLOY_CHECK_TRUE(q.min_y == 0.0);
            MALLOY_CHECK_TRUE(q.min_x > 1.0); // the hot wall heats a cooler gas

            ParticleSettings x_only{1.0, box};
            x_only.wall_temperature = WallTemperature2D{3.0, 0.5, 0.0, 0.0};
            x_only.seed = 12;
            ParticleWorld m{SimulationSettings{0.01}, x_only, gas};
            const Vec2 p0 = total_momentum(gas);
            for (int i = 0; i < 1500; ++i)
            {
                MALLOY_CHECK_TRUE(m.step().ok());
            }
            const auto& j = m.wall_impulse();
            MALLOY_CHECK_NEAR(total_momentum(m.particles()).x - p0.x, j.min_x - j.max_x, 1e-9);
        }
        {
            const Aabb3 box{Vec3{-2.0, -1.5, -1.0}, Vec3{2.0, 1.5, 1.0}};
            std::vector<Particle3D> gas = lattice_gas3d(box, 4, 3, 2, 1.0, 0.15, 1.0);
            for (std::size_t i = 0; i < gas.size(); i += 3)
            {
                gas[i].mass = 2.5;
                gas[i].radius = 0.25;
            }
            ParticleSettings3D s;
            s.bounds = box;
            s.wall_temperature = WallTemperature3D{3.0, 0.5, 0.0, 1.5, 2.0, 0.0};
            s.seed = 13;
            ParticleWorld3D w{SimulationSettings{0.01}, s, gas};
            const Real ke0 = total_kinetic_energy3d(gas);
            for (int i = 0; i < 1500; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
            }
            const auto& q = w.wall_heat();
            MALLOY_CHECK_NEAR(total_kinetic_energy3d(w.particles()) - ke0,
                              q.min_x + q.max_x + q.min_y + q.max_y + q.min_z + q.max_z,
                              1e-9);
            MALLOY_CHECK_TRUE(q.min_y == 0.0 && q.max_z == 0.0);

            ParticleSettings3D x_only;
            x_only.bounds = box;
            x_only.wall_temperature = WallTemperature3D{3.0, 0.5, 0.0, 0.0, 0.0, 0.0};
            x_only.seed = 14;
            ParticleWorld3D m{SimulationSettings{0.01}, x_only, gas};
            const Vec3 p0 = total_momentum3d(gas);
            for (int i = 0; i < 1500; ++i)
            {
                MALLOY_CHECK_TRUE(m.step().ok());
            }
            const auto& j = m.wall_impulse();
            MALLOY_CHECK_NEAR(total_momentum3d(m.particles()).x - p0.x, j.min_x - j.max_x, 1e-9);
        }

        // Equilibrium: with every wall thermal at T_w the gas relaxes to T_w,
        // for any masses (the diffuse wall's stationary state is the canonical
        // Maxwell-Boltzmann distribution at T_w). Point particles, half of
        // mass 4, started cold (T 0.3125) on a lattice in a unit box, time-
        // averaged after a burn-in. Tolerances: at least three times the
        // seed-to-seed spread measured across 1000 seeds, plus the fixed step's
        // bias (it raises T by about dt / (2 d E[tau]), 0.3 per cent here, and
        // moves the heavy-to-light ratio to about 0.998). They exclude the
        // unrelaxed start (0.31), a half-Gaussian emitter (about 0.50 T_w for
        // these collisionless points in 2D; in a colliding gas it would be
        // d/(d+1) T_w), and a sampler that forgot to divide by m (the
        // heavy-to-light energy ratio would be 4).
        const auto relax2 = [&](std::uint64_t seed, Real& temperature, Real& ratio) {
            const Aabb box{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
            std::vector<Particle2D> gas = lattice_gas(box, 4, 4, 1.0, 0.0, 0.5);
            for (std::size_t i = 0; i < gas.size(); i += 2)
            {
                gas[i].mass = 4.0;
            }
            ParticleSettings s{1.0, box};
            s.wall_temperature = WallTemperature2D{1.0, 1.0, 1.0, 1.0};
            s.seed = seed;
            ParticleWorld w{SimulationSettings{0.01}, s, gas};
            Real t_sum = 0.0;
            Real heavy = 0.0;
            Real light = 0.0;
            int samples = 0;
            for (int i = 1; i <= 40000; ++i)
            {
                if (!w.step().ok())
                {
                    return false;
                }
                if (i > 2000)
                {
                    t_sum += malloy::particles::temperature(w.particles());
                    for (const Particle2D& p : w.particles())
                    {
                        const Real ke = 0.5 * p.mass * malloy::math::dot(p.velocity, p.velocity);
                        (p.mass > 2.0 ? heavy : light) += ke;
                    }
                    ++samples;
                }
            }
            temperature = t_sum / samples;
            ratio = heavy / light;
            return true;
        };
        {
            Real temperature = 0.0;
            Real ratio = 0.0;
            MALLOY_CHECK_TRUE(relax2(3, temperature, ratio));
            MALLOY_CHECK_NEAR(temperature, 1.0, 0.056);
            MALLOY_CHECK_NEAR(ratio, 1.0, 0.095);
        }
        {
            const Aabb3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};
            std::vector<Particle3D> gas = lattice_gas3d(box, 4, 2, 2, 1.0, 0.0, 0.5);
            for (std::size_t i = 0; i < gas.size(); i += 2)
            {
                gas[i].mass = 4.0;
            }
            ParticleSettings3D s;
            s.bounds = box;
            s.wall_temperature = WallTemperature3D{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
            s.seed = 4;
            ParticleWorld3D w{SimulationSettings{0.01}, s, gas};
            Real t_sum = 0.0;
            Real heavy = 0.0;
            Real light = 0.0;
            int samples = 0;
            for (int i = 1; i <= 40000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                if (i > 2000)
                {
                    t_sum += malloy::particles::temperature3d(w.particles());
                    for (const Particle3D& p : w.particles())
                    {
                        const Real ke = 0.5 * p.mass * malloy::math::dot(p.velocity, p.velocity);
                        (p.mass > 2.0 ? heavy : light) += ke;
                    }
                    ++samples;
                }
            }
            MALLOY_CHECK_NEAR(t_sum / samples, 1.0, 0.04);
            MALLOY_CHECK_NEAR(heavy / light, 1.0, 0.065);
        }

        // Free-molecular heat flow, in closed form: point particles never
        // collide, so each one bounces between a hot wall (min_x, T_h = 4) and
        // a cold one (max_x, T_c = 1) on its own, the side walls specular. By
        // renewal-reward the long-run heat rate per particle is the mean energy
        // carried per round trip over the mean trip time,
        //   q = ((d+1)/2) (T_h - T_c) / (L sqrt(pi m / 2) (T_h^-1/2 + T_c^-1/2) + dt),
        // the dt being the fixed step's lag (each trip wastes half a step at
        // each end). 2.380990 in 2D and 3.174654 in 3D for L = 1, m = 1,
        // dt = 0.01. Estimated as the mean of the two walls' magnitudes,
        // (Q_hot - Q_cold) / (2 N t). Tolerance: three times the derived sigma
        // (the trip time is heavy-tailed, so it is truncated at the run length),
        // 0.072 and 0.090, plus 0.5 per cent. It excludes a wall that kept the
        // tangential velocity (0.67 q in 2D, 0.50 q in 3D), a half-Gaussian
        // normal (about 0.2 q), and no flow. The time-averaged temperature,
        // about sqrt(T_h T_c) = 2 (2.0093 and 2.0071 with the step), is checked
        // too. It is heavy-tailed (a slow, hot emission lingers), so three
        // standard deviations pass only 99.3 per cent of seeds; its tolerances,
        // 0.19 and 0.18, are the 99.73 per cent quantiles of |T - centre|
        // measured over 5000 seeds. They exclude the arithmetic mean 2.5 and
        // the harmonic 1.6.
        {
            const Aabb box{Vec2{0.0, 0.0}, Vec2{1.0, 1.0}};
            ParticleSettings s{1.0, box};
            s.wall_temperature = WallTemperature2D{4.0, 1.0, 0.0, 0.0};
            s.seed = 21;
            ParticleWorld w{SimulationSettings{0.01}, s, lattice_gas(box, 4, 4, 1.0, 0.0, 1.0)};
            Real t_sum = 0.0;
            int samples = 0;
            for (int i = 1; i <= 40000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                if (i > 2000)
                {
                    t_sum += malloy::particles::temperature(w.particles());
                    ++samples;
                }
            }
            const auto& q = w.wall_heat();
            const Real rate = (q.min_x - q.max_x) / (2.0 * 16.0 * w.elapsed_time());
            MALLOY_CHECK_NEAR(rate, 2.380990, 0.227);
            MALLOY_CHECK_NEAR(t_sum / samples, 2.0093, 0.19);
            MALLOY_CHECK_TRUE(q.min_y == 0.0 && q.max_y == 0.0);
        }
        {
            const Aabb3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};
            ParticleSettings3D s;
            s.bounds = box;
            s.wall_temperature = WallTemperature3D{4.0, 1.0, 0.0, 0.0, 0.0, 0.0};
            s.seed = 22;
            ParticleWorld3D w{SimulationSettings{0.01}, s,
                              lattice_gas3d(box, 4, 2, 2, 1.0, 0.0, 1.0)};
            Real t_sum = 0.0;
            int samples = 0;
            for (int i = 1; i <= 40000; ++i)
            {
                MALLOY_CHECK_TRUE(w.step().ok());
                if (i > 2000)
                {
                    t_sum += malloy::particles::temperature3d(w.particles());
                    ++samples;
                }
            }
            const auto& q = w.wall_heat();
            const Real rate = (q.min_x - q.max_x) / (2.0 * 16.0 * w.elapsed_time());
            MALLOY_CHECK_NEAR(rate, 3.174654, 0.286);
            MALLOY_CHECK_NEAR(t_sum / samples, 2.0071, 0.18);
        }
    }

    std::cout << "malloy_particles_tests passed\n";
    return 0;
}
