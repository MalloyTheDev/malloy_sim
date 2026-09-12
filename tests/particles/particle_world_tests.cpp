#include <malloy/particles/particles.hpp>

#include <malloy/math/math.hpp>
#include <malloy/sim_core/sim_core.hpp>
#include <test_check.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
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

    std::cout << "malloy_particles_tests passed\n";
    return 0;
}
