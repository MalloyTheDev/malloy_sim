#include <malloy/particles/thermo.hpp>

#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

#include <malloy/math/math.hpp>

namespace malloy::particles
{
namespace
{
const math::Real pi = math::Real{3.14159265358979323846};

// The golden angle, pi (3 - sqrt 5), about 137.5 degrees. Successive multiples
// of it never line up, so consecutive pairs of a lattice point unrelated ways
// while the set still covers the circle evenly.
const math::Real golden_angle = pi * (math::Real{3} - std::sqrt(math::Real{5}));

// Sum of v^2 and of v^4 over a list, the two moments every statistic here uses.
struct SpeedMoments
{
    math::Real v2{0.0};
    math::Real v4{0.0};
};

template <typename Particle>
SpeedMoments speed_moments(const std::vector<Particle>& particles)
{
    SpeedMoments moments;
    for (const Particle& p : particles)
    {
        const math::Real v2 = math::dot(p.velocity, p.velocity);
        moments.v2 += v2;
        moments.v4 += v2 * v2;
    }
    return moments;
}

template <typename Particle>
math::Real kinetic(const std::vector<Particle>& particles)
{
    math::Real energy = 0.0;
    for (const Particle& p : particles)
    {
        energy += math::Real{0.5} * p.mass * math::dot(p.velocity, p.velocity);
    }
    return energy;
}

template <typename Particle>
math::Real temperature_in(const std::vector<Particle>& particles, math::Real dimension)
{
    if (particles.empty())
    {
        return math::Real{0};
    }
    return math::Real{2} * kinetic(particles) /
           (dimension * static_cast<math::Real>(particles.size()));
}

template <typename Particle>
math::Real speed_ratio_of(const std::vector<Particle>& particles)
{
    const SpeedMoments moments = speed_moments(particles);
    if (!(moments.v2 > math::Real{0}))
    {
        return math::Real{0}; // empty, or nothing moves: no distribution to shape
    }
    // <v^4> / <v^2>^2 = N sum(v^4) / sum(v^2)^2
    return static_cast<math::Real>(particles.size()) * moments.v4 /
           (moments.v2 * moments.v2);
}

template <typename Particle>
math::Real slow_fraction_of(const std::vector<Particle>& particles)
{
    if (particles.empty())
    {
        return math::Real{0};
    }
    const math::Real mean =
        speed_moments(particles).v2 / static_cast<math::Real>(particles.size());
    std::size_t slow = 0;
    for (const Particle& p : particles)
    {
        if (math::dot(p.velocity, p.velocity) < mean)
        {
            ++slow;
        }
    }
    return static_cast<math::Real>(slow) / static_cast<math::Real>(particles.size());
}

// The fractional part, in [0, 1).
math::Real fraction(math::Real x)
{
    return x - std::floor(x);
}

// Pair k of a lattice: the even particle takes `direction`, the odd one its
// exact negative, so each pair carries zero momentum.
math::Real pair_sign(std::size_t index)
{
    return index % 2 == 0 ? math::Real{1} : math::Real{-1};
}
} // namespace

math::Real temperature(const std::vector<Particle2D>& particles)
{
    return temperature_in(particles, math::Real{2});
}

math::Real temperature3d(const std::vector<Particle3D>& particles)
{
    return temperature_in(particles, math::Real{3});
}

math::Real speed_ratio(const std::vector<Particle2D>& particles)
{
    return speed_ratio_of(particles);
}

math::Real speed_ratio3d(const std::vector<Particle3D>& particles)
{
    return speed_ratio_of(particles);
}

math::Real slow_fraction(const std::vector<Particle2D>& particles)
{
    return slow_fraction_of(particles);
}

math::Real slow_fraction3d(const std::vector<Particle3D>& particles)
{
    return slow_fraction_of(particles);
}

std::optional<math::Real> pressure(const WallImpulse2D& impulse,
                                   const collide::Aabb& bounds,
                                   math::Real elapsed)
{
    if (!(elapsed > math::Real{0}) || !math::is_finite(elapsed))
    {
        return std::nullopt; // no time has passed, so nothing has been measured
    }
    const math::Real width = bounds.max.x - bounds.min.x;
    const math::Real height = bounds.max.y - bounds.min.y;
    const math::Real perimeter = math::Real{2} * (width + height);
    const math::Real total = impulse.min_x + impulse.max_x + impulse.min_y + impulse.max_y;
    return total / (elapsed * perimeter);
}

std::optional<math::Real> pressure3d(const WallImpulse3D& impulse,
                                     const collide::Aabb3& bounds,
                                     math::Real elapsed)
{
    if (!(elapsed > math::Real{0}) || !math::is_finite(elapsed))
    {
        return std::nullopt;
    }
    const math::Real width = bounds.max.x - bounds.min.x;
    const math::Real height = bounds.max.y - bounds.min.y;
    const math::Real depth = bounds.max.z - bounds.min.z;
    const math::Real area =
        math::Real{2} * (width * height + height * depth + depth * width);
    const math::Real total = impulse.min_x + impulse.max_x + impulse.min_y +
                             impulse.max_y + impulse.min_z + impulse.max_z;
    return total / (elapsed * area);
}

std::vector<Particle2D> lattice_gas(const collide::Aabb& bounds, int nx, int ny,
                                    math::Real mass, math::Real radius,
                                    math::Real speed)
{
    std::vector<Particle2D> particles;
    if (nx < 1 || ny < 1)
    {
        return particles;
    }
    particles.reserve(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));
    const math::Real cell_x = (bounds.max.x - bounds.min.x) / static_cast<math::Real>(nx);
    const math::Real cell_y = (bounds.max.y - bounds.min.y) / static_cast<math::Real>(ny);
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const std::size_t index = particles.size();
            const math::Real angle = static_cast<math::Real>(index / 2) * golden_angle;
            const math::Vec2 direction{std::cos(angle), std::sin(angle)};
            Particle2D p;
            p.position = math::Vec2{bounds.min.x + (static_cast<math::Real>(i) + 0.5) * cell_x,
                                    bounds.min.y + (static_cast<math::Real>(j) + 0.5) * cell_y};
            p.velocity = direction * (pair_sign(index) * speed);
            p.mass = mass;
            p.radius = radius;
            particles.push_back(p);
        }
    }
    return particles;
}

math::Real emission_normal_speed(SplitMix64& random, math::Real temperature,
                                 math::Real mass)
{
    const math::Real u = random.uniform();
    return std::sqrt(math::Real{-2} * (temperature / mass) * std::log(u));
}

math::Vec2 emission_tangent(SplitMix64& random, math::Real temperature,
                            math::Real mass)
{
    // Two draws, in this order, as separate statements so the order is fixed.
    const math::Real u1 = random.uniform();
    const math::Real u2 = random.uniform();
    const math::Real radius =
        std::sqrt(math::Real{-2} * (temperature / mass) * std::log(u1));
    const math::Real angle = math::Real{2} * pi * u2;
    return math::Vec2{radius * std::cos(angle), radius * std::sin(angle)};
}

std::vector<Particle3D> lattice_gas3d(const collide::Aabb3& bounds, int nx, int ny,
                                      int nz, math::Real mass, math::Real radius,
                                      math::Real speed)
{
    std::vector<Particle3D> particles;
    if (nx < 1 || ny < 1 || nz < 1)
    {
        return particles;
    }
    particles.reserve(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                      static_cast<std::size_t>(nz));
    const math::Real cell_x = (bounds.max.x - bounds.min.x) / static_cast<math::Real>(nx);
    const math::Real cell_y = (bounds.max.y - bounds.min.y) / static_cast<math::Real>(ny);
    const math::Real cell_z = (bounds.max.z - bounds.min.z) / static_cast<math::Real>(nz);
    const math::Real sqrt2 = std::sqrt(math::Real{2});
    for (int l = 0; l < nz; ++l)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                const std::size_t index = particles.size();
                const math::Real pair = static_cast<math::Real>(index / 2);
                const math::Real azimuth = pair * golden_angle;
                const math::Real z = math::Real{1} - math::Real{2} * fraction((pair + 0.5) * sqrt2);
                const math::Real ring = std::sqrt(std::fmax(math::Real{0}, math::Real{1} - z * z));
                const math::Vec3 direction{ring * std::cos(azimuth), ring * std::sin(azimuth), z};
                Particle3D p;
                p.position =
                    math::Vec3{bounds.min.x + (static_cast<math::Real>(i) + 0.5) * cell_x,
                               bounds.min.y + (static_cast<math::Real>(j) + 0.5) * cell_y,
                               bounds.min.z + (static_cast<math::Real>(l) + 0.5) * cell_z};
                p.velocity = direction * (pair_sign(index) * speed);
                p.mass = mass;
                p.radius = radius;
                particles.push_back(p);
            }
        }
    }
    return particles;
}
} // namespace malloy::particles
