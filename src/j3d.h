/*
(c) 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001
for Los Alamos National Laboratory (LANL), which is operated by Triad National
Security, LLC for the U.S. Department of Energy/National Nuclear Security
Administration. All rights in the program are reserved by Triad National Security,
LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
The Government is granted for itself and others acting on its behalf a nonexclusive,
paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display
publicly, and to permit others to do so.
*/

#ifndef _J3D_H_
#define _J3D_H_

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <string.h>

#include "j3d-config.h"

#ifdef J3D_ENABLE_Kokkos
#include "Kokkos_Core.hpp"
#endif

#define J3D_MAX_NUM_MOMENTS ((J3D_MAX_ORDER + 1) * (J3D_MAX_ORDER + 2) * (J3D_MAX_ORDER + 3) / 6)

static const double ONE_SIXTH = 0.166666666666666666666666666666666666666666666666666667;
static const double ONE_THIRD = 0.333333333333333333333333333333333333333333333333333333;

/**
 * \brief Copy the contents of one int pointer to another.
 *        Assumes that memory regions dont overlap.
 * \param [in, out] dest The array we will be copying to.
 * \param [in] src The array we are copying from.
 * \param [in] n Number of bytes to copy.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_memcpy(void* const dest, const void* const src, const int n);

/**
 * \brief Set the elements of a pointer to (char) zero.
 * \param [in, out] ptr The array we are setting.
 * \param [in] n number of bytes to set.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_memset0(void* const ptr, const int n);

/**
 * \struct j3d_vec3
 * \brief a 3-vector
 */
struct j3d_vec3{
  double xyz[3];

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
  j3d_vec3(){}

KOKKOS_FUNCTION
  ~j3d_vec3(){}

KOKKOS_FUNCTION
  j3d_vec3& operator=(const j3d_vec3& v) {
    j3d_memcpy(this->xyz, v.xyz, 3 * sizeof(double));
    return *this;
  }
#endif

};

/**
 * \struct j3d_plane
 * \brief A plane.
 */
struct j3d_plane {
  j3d_vec3 n; /*!< Unit length vector normal to plane. */
  double d;   /*!< Signed perpendicular distance to the origin. */

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
  j3d_plane(){}

KOKKOS_FUNCTION
  ~j3d_plane(){}

KOKKOS_FUNCTION
  j3d_plane& operator=(const j3d_plane& pl) {
    this->n = pl.n;
    this->d = pl.d;
    return *this;
  }
#endif
};

/**
 * \struct j3d_poly
 * \brief  A polyhedron that can be convex, concave, or multi-connected.
 */
typedef struct {
  double verts[J3D_MAX_VERTS][3];
  int dir_edges[2 * J3D_MAX_EDGES];
  int indices[J3D_MAX_VERTS + 1];
  int num_dir_edges;
  int num_verts;
} j3d_poly;

/**
 * \brief Compute the dot product of two vectors.
 * \param [in] va One of the input vectors.
 * \param [in] vb One of the input vectors.
 * \return Dot product of va and vb.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
double j3d_dot(const double* const va, const double* const vb);

/**
 * \brief Compute the weighted average of two vertices.
 * \param [out] vr Resulting weighted average vertex.
 * \param [in] va One of the input vectors.
 * \param [in] wa The weight of vector va.
 * \param [in] vb One of the input vectors.
 * \param [in] wb The weight of vb.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_wav(double* const vr,
             const double* const va,
             const double wa,
             const double* const vb,
             const double wb);

/**
 * \brief Normalize a vertex
 * \param[in, out] vert The vertex to be normalized.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_norm(double* const vert);

/**
 * \brief Compute the number of moments for a given order.
 * \param [int] order The given order.
 * \return Number of moments.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_num_moments(const int order);

/**
 * \brief Shift the moments back after shifting the poly by its center
 * \param[in, out] moments Array containing the moments.  Must be atleast
 *                 (order+1)*(order+2)*(order+3)/6 long.
 * \param[in] order The order of the moments.
 * \param[in] vc Center of the poly
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_shift_moments(double* const moments, int order, const double* const vc);

/**
 * \brief Compute the center of a list of vertices.
 * \param [out] vc Center of the vertices.
 * \param [in] verts List of vertices.
 * \param [in] numverts Number of vertices.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_center(double* const vc, const double (*const verts)[3], const int numverts);

/**
 * \brief From a set of vertices defining a tetrahedra
 *        get the faces describing it.
 * \param [out] faces The four faces describing the tetrahedra.
 * \param [in] verts The four vertices describing the tetrahedra.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_tet_faces_from_verts(j3d_plane* faces, const j3d_vec3* const verts);

/** \brief From two vertices describing an axis aligned box,
 *         recover the faces describing the box.
 *  \param [out] faces The six faces describing the box.
 *  \param [in] bounds The two vertices describing the
 *              axis aligned box.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_box_faces_from_verts(j3d_plane* faces, const j3d_vec3* const bounds);

/**
 * \brief Get the signed volume of the tetrahedron defined by the
 *        input vertices.
 * \param [in] verts Four vertices defining a tetrahedron.
 * \return The signed volume of the input tetrahedron.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
double j3d_orient(const j3d_vec3* const verts);

/**
 * \brief Initialize a poly from a BREP representation.
 * \param [out] poly The polyhedron being initialized.
 * \param [in] vertices The array of length numverts containing the vertices.
 * \param [in] numverts Number of vertices in the polyhedra.
 * \param [in] faceinds Array of arrays where each internal array
 *             corresponds to a face and that array contains
 *             the node indices on the face in the outward normal direction.
 * \param [in] numvertsperface Number of vertices per face.
 * \param [in] numfaces Number of faces in the polyhedra.
 * \return Status code indicating if initialization was successful
 *         1 -> successful, 0 -> unsuccessful.
 *
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_init_poly(j3d_poly* poly,
                  j3d_vec3* vertices,
                  int numverts,
                  int** faceinds,
                  int* numvertsperface,
                  int numfaces);

/**
 * \brief Initialize a poly from a BREP representation.
 * \param [out] poly The polyhedron being initialized.
 * \param [in] vertices The array of length numverts containing the vertices.
 * \param [in] numverts Number of vertices in the polyhedra.
 * \param [in] facevertinds Flattened array, where each block
 *             corresponds to a face and that array contains
 *             the node indices on the face in the outward normal direction.
 * \param [in] facevertoffsets Offsets into facevertinds array to determine where 
 *             indices of a specific face's vertices begin
 * \param [in] numfaces Number of faces in the polyhedra.
 * \return Status code indicating if initialization was successful
 *         1 -> successful, 0 -> unsuccessful.
 *
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_brep_to_poly(j3d_poly* poly,
                     const j3d_vec3* const vertices,
                     const int numverts,
                     const int* const facevertinds,
                     const int* const facevertoffsets,
                     const int numfaces);

/**
 * \brief Initialize a j3d poloyhedron as a tetrahedron from
 *        a list of vectors.
 * \param[out] poly The poly to be initialized.
 * \param[in] verts An array of four vectors, giving
 *            the vertices of the tetrahedron.
 * \return Status code indicating if initialization was successful
 *         1 -> successful, 0 -> unsuccessful.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_init_tet(j3d_poly* poly, const j3d_vec3* const verts);

/**
 * \brief Initialize a poly as an axis aligned cube.
 * \param [out] poly The polyhedron to be initialized.
 * \param [in] rbounds An array of two vectors, giving the
 *             lower and upper corners of the cube.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_init_box(j3d_poly* poly, const j3d_vec3* const rbounds);

/**
 * \brief Copy the contents from one poly to another
 * \param [out] destination The poly we are copying to.
 * \param [in] source The poly we are copying from
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_copy_poly(j3d_poly* destination, const j3d_poly* const source);

/**
 * \brief Clip and cap polyhedra against an arbitrary number of planes.
 * \param [in, out] poly The polyhedron to be clipped.
 * \param [in] planes An array of planes to clip the polyhedra with.
 * \param [in] nplanes Number of planes.
 * \return status Code indicating if clip and cap was successful
 *         1 -> successful, 0 -> unsuccessful.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_clip(j3d_poly* poly, const j3d_plane* const planes, const int nplanes);

/**
 * \brief Split and cap a list of polyhedra against a single plane.
 * \param [in] inpolys Array of polys to split.
 * \param [in] npolys Number of polys we will be splitting.
 * \param [in] plane Plane we will be splitting with.
 * \param [out] out_pos The array of resulting polys on the positive side
 *              of the splitting plane.  Must be atleast npolys long.
 * \param [out] out_neg The array of resulting polys on the negative side
 *              of the splitting plane.  Must be atleast npolys long.
 * \return status Code indicating if split and cap was successful
 *         1 -> successful, 0 -> unsuccessful.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_split(j3d_poly* inpolys,
              const int npolys,
              const j3d_plane plane,
              j3d_poly* out_pos,
              j3d_poly* out_neg);

/**
 * \brief explicitly computes the zeroth order moment of a polyhedra
 * \param[in] poly Polyhedra to compute the zeroth order moment over
 * \param[out] moments Array containing the zeroth order moment.
 *             Must have a length of atleast 1.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_0(const j3d_poly* const poly, double* moments);

/**
 * \brief explicitly computes the first order moment of a polyhedra
 * \param[in] poly Polyhedra to compute the first order moments over
 * \param[out] moments Array containing the first order moments.
 *             Must have a length of atleast 4.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_1(const j3d_poly* const poly, double* moments);

/**
 * \brief explicitly computes the second order moment of a polyhedra
 * \param[in] poly Polyhedra to compute the second order moments over
 * \param[out] moments Array containing the second order moments.
 *             Must have a length of atleast 10.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_2(const j3d_poly* const poly, double* moments);

/**
 * \brief explicitly computes the third order moment of a polyhedra
 * \param[in] poly Polyhedra to compute the third order moments over
 * \param[out] moments Array containing the third order moments.
 *             Must have a length of atleast 20.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_3(const j3d_poly* const poly, double* moments);

/**
 * \brief computes the nth order moments of a polyhedra via
 *        the recursive method of Koehl.
 * \param [in] poly The polyhedra to compute the nth order moments over.
 * \param [out] moments Array containing the nth order moments.
 *              Must have a length of atleast (order+1)*(order+2)*(order+3)/6.
 * \param [out] order Order or the moments being computed.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_n(const j3d_poly* const poly, double* moments, int order);

/**
 * \brief computes the nth order moments of a polyhedra either explicitly
 *        or recursively depending on the order.
 * \param [in] poly The polyhedra to compute the nth order moments over.
 * \param [out] moments Array containing the nth order moments.
 *              Must have a length of atleast (order+1)*(order+2)*(order+3)/6.
 *              Orer of moments is row major, i.e. 1,x,y,z,x^2,xy,xz,y^2,yz,z^2,x^3.
 * \param [out] order Order or the moments being computed.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments(const j3d_poly* const poly, double* moments, int order);

/**
 * \brief print the information in a poly.
 * \param[in] poly Polyhedra to be print.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_print_poly(const j3d_poly* const poly);

#endif
