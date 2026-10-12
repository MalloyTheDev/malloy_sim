#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include <malloy/math/vec3.hpp>
#include <malloy/view3d/view3d.hpp>

namespace malloy::view3d
{
// Triangle meshes for the frame's solid primitives (M41), so a viewer draws
// triangles and lines it is handed and computes no geometry of its own.
//
// Every triangle runs COUNTER-CLOCKWISE seen from the side its normal points
// to, and every closed mesh has its normals pointing out, which is what back-
// face culling with counter-clockwise front faces (OpenGL's and raylib's
// default) needs to show the near side. `normal` is the unit outward normal of
// the face, for flat shading.
struct Triangle
{
    Vec3 a{};
    Vec3 b{};
    Vec3 c{};
    Vec3 normal{};
};

// The twelve triangles of a box, two per face of box_faces; each face's
// normal is the body axis turned by the box's orientation.
std::array<Triangle, 12> box_triangles(const Box& box);

// How finely a sphere is cut: `rings` latitude bands from the +z pole and
// `slices` longitudes, giving slices * (2 rings - 2) triangles.
struct SphereDetail
{
    int rings{16};
    int slices{24};
};

// Level of detail by how many spheres a frame holds, so a dense gas stays
// cheap: up to 64 spheres 16 x 24 (720 triangles each), up to 512 10 x 16
// (288), up to 4096 6 x 10 (100), and beyond that 4 x 6 (36).
SphereDetail sphere_detail(std::size_t sphere_count);

// A UV sphere around `center`. Its vertices lie on the sphere (the poles
// exactly), neighbouring triangles share vertices bit for bit, so the mesh is
// closed, and it lies inside the true sphere. A detail below 2 x 3 is raised
// to it.
std::vector<Triangle> sphere_triangles(const Vec3& center, Real radius,
                                       const SphereDetail& detail);

// The same sphere placed from a unit one: append `unit` (sphere_triangles of
// the origin and radius 1) scaled by `radius` and moved to `center`. The
// vertices are bit for bit those of sphere_triangles(center, radius, detail)
// and the normals are the unit mesh's, so a viewer drawing many spheres makes
// the unit mesh once per level of detail and does no trigonometry per sphere.
void place_sphere(const std::vector<Triangle>& unit, const Vec3& center, Real radius,
                  std::vector<Triangle>& out);

// The two triangles of a ground patch, counter-clockwise from its free side.
std::array<Triangle, 2> patch_triangles(const Patch& patch);

// A "nice" grid spacing (1, 2 or 5 times a power of ten) near span / lines.
// Zero when span or lines is not positive and finite.
Real grid_spacing(Real span, int lines);

// The grid lines across a patch, about twenty each way. They sit at multiples
// of the spacing measured in the plane's own (u, v) coordinates from the point
// of the plane nearest the world origin (normal * offset), NOT from the patch
// centre, so the grid stays put in the world when the patch is re-centred, and
// motion along the ground stays visible however the camera follows. Empty when
// the patch is not finite or the numbers are too large for a grid to resolve.
std::vector<Segment> patch_grid(const Patch& patch);

// An arrow as five lines: the shaft and a four-line head at the tip, a fifth
// of the shaft long. Empty for a zero-length arrow.
std::vector<Segment> arrow_lines(const Arrow& arrow);

// Diffuse (Lambert) shade in [ambient, 1] for a unit face normal and a unit
// direction TOWARD the light: ambient + (1 - ambient) max(0, n . l). Computed
// per face on the CPU, so a viewer needs no shader.
Real shade(const Vec3& normal, const Vec3& to_light, Real ambient);
} // namespace malloy::view3d
