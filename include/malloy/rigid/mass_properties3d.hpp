#pragma once

#include <cstddef>
#include <vector>

#include <malloy/math/quat.hpp>
#include <malloy/math/vec3.hpp>
#include <malloy/rigid/rigid_body3d.hpp>

namespace malloy::rigid
{
// One solid uniform sphere, a part of a compound body. Its own density, so a
// body can be made of parts of different materials, a dense core inside a
// lighter shell.
struct SolidSphere
{
    math::Vec3 center{};
    math::Real radius{0.0};
    math::Real density{1.0};
};

// One solid uniform box, the second mass primitive (M26). Centred at `center`,
// with half-widths `half_extents` along its OWN axes, turned by `orientation`.
// Unlike a sphere, a box has three distinct principal moments and a definite
// orientation, so a single tilted box already gives an inertia tensor that is
// not diagonal in the lab frame.
struct SolidBox
{
    math::Vec3 center{};
    math::Vec3 half_extents{}; // half the width along each of the box's own axes
    math::Quat orientation{};  // the box frame in the lab frame
    math::Real density{1.0};
};

// One triangle of a mesh: three indices into the vertex list, wound
// counter-clockwise as seen from OUTSIDE, so the face normal (v1-v0) x (v2-v0)
// points out of the solid. Winding is what carries the sign in the volume
// integral, so a mesh wound the wrong way describes negative volume and is
// rejected rather than silently giving a body inside out.
struct MeshTriangle
{
    std::size_t v0{};
    std::size_t v1{};
    std::size_t v2{};
};

// A solid uniform body bounded by a closed, outward-wound triangle mesh (M26's
// third mass primitive, and the general one: the sphere and the box are shapes
// this could also describe). Its mass properties are the volume integrals over
// the enclosed solid, computed as a signed sum of tetrahedra from the origin to
// each triangle, so the origin need not be inside the mesh and need not be the
// centre of mass.
struct SolidMesh
{
    std::vector<math::Vec3> vertices;
    std::vector<MeshTriangle> triangles;
    math::Real density{1.0};
};

// The mass distribution of a compound body, ready to drop into a RigidBody3D.
// The 3D sibling of MassProperties (ADR 0007's step 2), and the piece that lets
// a body's physical properties be COMPUTED from its geometry rather than typed
// in by hand: mass and inertia follow from shape and density.
//
// Unlike the 2D version, the inertia is not one number but three PRINCIPAL
// moments plus the orientation of the principal-axis frame in the frame the
// parts were given in. A compound body's inertia tensor is generally not
// diagonal in the lab frame; diagonalizing it is what produces these, and it is
// exactly what RigidBody3D stores (three moments and an orientation), so the
// result drops straight in.
struct MassProperties3D
{
    math::Real mass{0.0};
    math::Vec3 center_of_mass{};
    // Principal moments about the centre of mass, ascending, matching
    // RigidBody3D::inertia once the body carries `orientation`.
    math::Vec3 inertia{};
    // The principal-axis frame: a body with this orientation and these moments
    // has the compound body's full inertia tensor in the lab frame.
    math::Quat orientation{};
};

// Compute the mass properties of a body made of solid uniform spheres.
//
// Mass is density times volume per part; the centre of mass is the mass-weighted
// mean of the centres; the inertia tensor about that centre is the sum of each
// sphere's own tensor ((2/5) m r^2, isotropic) shifted by the parallel-axis
// theorem, m (|d|^2 I - d d^T); and the three principal moments and their frame
// are that tensor diagonalized.
//
// Returns a zero-mass result (nothing usable) for an empty list, or if any part
// has a non-positive or non-finite radius or density, rather than a meaningless
// number. This mirrors the 2D convention.
MassProperties3D mass_properties_3d(const std::vector<SolidSphere>& parts);

// The mass properties of a single solid box. Three distinct moments,
// (1/12) m (w_j^2 + w_k^2) about each of its own axes (w the full widths),
// carried into the lab frame by its orientation. Invalid input yields a
// zero-mass result. A single sphere goes through the list overload above (a
// one-element vector); only the box, a distinct type, needs its own.
MassProperties3D mass_properties_3d(const SolidBox& box);

// The mass properties of a solid uniform body bounded by a closed triangle
// mesh, the general primitive. Mass is density times the enclosed volume; the
// centre of mass and the inertia tensor are the volume integrals of position
// and of position-squared, each accumulated as a signed sum over the tetrahedra
// spanned by the origin and every triangle (the divergence theorem), then the
// tensor is shifted to the centre of mass and diagonalized. An axis-aligned box
// mesh reproduces `SolidBox` exactly.
//
// Returns a zero-mass result for a malformed mesh (fewer than four vertices or
// triangles, a triangle index out of range, a non-finite vertex, a
// non-positive or non-finite density) or one whose enclosed volume is not
// positive and finite, which is what an inside-out winding or a degenerate or
// non-closed mesh produces. This mirrors the sphere and box convention.
MassProperties3D mass_properties_3d(const SolidMesh& mesh);

// Merge two mass-property sets into one body, the way a compound is assembled
// from primitives of any kind. Masses add, the centre of mass is their
// mass-weighted mean, and each tensor is reconstructed (R diag(moments) R^T),
// shifted to the shared centre by the parallel-axis theorem, summed, and
// diagonalized again. Combining is commutative and associative up to the
// eigenvalue ordering, so the order parts are given in does not matter.
//
// A zero-mass operand is the identity: it returns the other unchanged, so a
// fold over parts can start from an empty body.
MassProperties3D combine(const MassProperties3D& a, const MassProperties3D& b);

// Build a RigidBody3D at rest (no velocity, no spin) from mass properties: the
// position is the centre of mass, and the inertia and orientation are the
// principal ones. The collision radius is left zero; a caller that wants the
// body to collide sets it.
RigidBody3D rigid_body_from(const MassProperties3D& properties);
} // namespace malloy::rigid
