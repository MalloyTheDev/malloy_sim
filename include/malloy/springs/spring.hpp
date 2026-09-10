#pragma once

#include <cstddef>
#include <vector>

#include <malloy/math/real.hpp>

namespace malloy::springs
{
// Index of a body in the world's body list. Springs reference bodies by id
// rather than owning body state, so the topology and the state stay separable
// (ADR 0008).
using BodyId = std::size_t;

// Position of a spring within its network. Positional, so it is invalidated by
// SpringNetwork::remove, exactly like an index into a vector. Documented rather
// than made stable with a slot map, which nothing yet needs.
using SpringId = std::size_t;

// A linear spring with a damper, connecting two bodies at their centres.
//
// Force magnitude along the axis is `stiffness * extension + damping * v_axial`
// (see accumulate_spring_forces for the full derivation and sign convention).
struct Spring
{
    BodyId a{0};
    BodyId b{1};

    math::Real rest_length{1.0};
    math::Real stiffness{1.0};
    math::Real damping{0.0};

    // Valid when the two endpoints are different bodies, rest_length, stiffness
    // and damping are all non-negative, and every value is finite.
    //
    // a == b is rejected rather than tolerated: a self-spring has zero length
    // forever, so its axis is permanently undefined and no useful force exists.
    bool is_valid() const;
};

// An ordered collection of springs.
//
// Order is part of the contract, not an implementation detail. Forces are
// accumulated by traversing this list, and floating-point addition is not
// associative, so a different traversal order gives different results
// (docs/04 determinism). add appends; remove preserves the relative order of
// everything else.
class SpringNetwork
{
public:
    // Appends and returns the new spring's id, which is its position.
    SpringId add(const Spring& spring);

    // Removes one spring, shifting later ids down by one. Returns false when
    // the id is out of range, changing nothing.
    bool remove(SpringId id);

    const std::vector<Spring>& springs() const { return springs_; }
    std::size_t size() const { return springs_.size(); }
    void clear() { springs_.clear(); }

    // True when every spring is valid and every endpoint id is below
    // body_count. Endpoint validity depends on the world, so it is checked
    // here rather than on the spring itself.
    bool is_valid_for(std::size_t body_count) const;

private:
    std::vector<Spring> springs_;
};
} // namespace malloy::springs
