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

#ifndef _J2D_H_
#define _J2D_H_

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "j3d-config.h"

#ifdef J3D_ENABLE_Kokkos
#include "Kokkos_Core.hpp"
#endif

#define J2D_MAX_NUM_MOMENTS ((J2D_MAX_ORDER + 1) * (J2D_MAX_ORDER + 2) / 2)

/**
 * \brief Copy the contents of one int pointer to another.
 *        Assumes that memory regions dont overlap.
 * \param [in, out] dest The array we will be copying to.
 * \param [in] src The array we are copying from.
 * \param [in] num_elems Number of elements to copy.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_memcpy(void* const dest, const void* const src, const int n);

/**
 * \brief Set the elements of a pointer to (char) zero.
 * \param [in, out] ptr The array we are setting.
 * \param [in] n number of bytes to set.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_memset0(void* const ptr, const int n);

/**
 * \struct j2d_vec2
 * \brief a 2-vector
 */
struct j2d_vec2 {
  double xy[2];

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
  j2d_vec2(){}

KOKKOS_FUNCTION
  ~j2d_vec2(){}

KOKKOS_FUNCTION
  j2d_vec2& operator=(const j2d_vec2& v) {
    j2d_memcpy(this->xy, v.xy, 2 * sizeof(double));
    return *this;
  }
#endif

};

/**
 * \struct j2d_plane
 * \brief A plane
 */
struct j2d_plane {
  j2d_vec2 n; /*!< Unit length vector normal to the plane. */
  double d;   /*!< Signed perpendicular distance to the origin. */

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
  j2d_plane(){}

KOKKOS_FUNCTION
  ~j2d_plane(){}

KOKKOS_FUNCTION
  j2d_plane& operator=(const j2d_plane& pl) {
    this->d = pl.d;
    this->n = pl.n;
    return *this;
  }
#endif

};

/**
 * \struct j2d_poly
 * \brief A polygon.  Can be convex, concave, or multi-connected.
 */
typedef struct {
  double verts[J2D_MAX_VERTS][2];
  int dir_edges[2 * J2D_MAX_VERTS];
  int num_verts;
} j2d_poly;

/**
 * \brief Compute the dot product of two vectors.
 * \param [in] va One of the input vectors.
 * \param [in] vb One of the input vectors.
 * \return Dot product of va and vb.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
double j2d_dot(const double* const va, const double* const vb);

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
void j2d_wav(double* const vr,
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
void j2d_norm(double* const vert);

/**
 * \brief Compute the number of moments for a given order.
 * \param [int] order The given order.
 * \return Number of moments.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j2d_num_moments(const int order);

/**
 * \brief Shift the moments back after shifting the poly by its center
 * \param[in, out] moments Array containing the moments.  Must be atleast
 *                 (order+1)*(order+2)/2 long.
 * \param[in] order The order of the moments.
 * \param[in] vc Center of the poly
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_shift_moments(double* const moments, int order, const double* const vc);

/**
 * \brief Compute the center of a list of vertices.
 * \param [out] vc Center of the vertices.
 * \param [in] verts List of vertices.
 * \param [in] numverts Number of vertices.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_center(double* const vc, const double (*const verts)[2], const int numverts);

/**
 * \brief Compute the signed volume of a triangle described by three points
 * \param [in] va, vb, vcVertices describing a triangle.
 * \return The signed volume of the input triangle.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
double j2d_orient(const double* const va, const double* const vb, const double* const vc);

/**
 * \brief Get all the faces of a boundary description of a polygon.
 * \param [out] faces Array of planes of length numverts describing the faces of the  polygon.
 * \param [in] vertices Array of vertices of length numverts containing the vertices
 *             of the polygon in counter-clockwise order.
 * \param [in] numverts Number of vertices in the polygon
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_poly_faces_from_verts(j2d_plane* faces, const j2d_vec2* const vertices, const int numverts);

/**
 * \brief Initialize a simply-connected general polygeon from a list of vertices.
 * \param [out] poly The polygon we are initializing.
 * \param [in] vertices Array of length numverts containing the vertices
 *             of the polygon in counter-clockwise order.
 * \param [in] numverts Number of vertices in the polygon.
 * \return Status code indicating if initialization was successful
 *         1 -> successful, 0 -> unsuccessful.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j2d_init_poly(j2d_poly* poly, const j2d_vec2* const vertices, const int numverts);

/**
 * \brief Initialize a polygon as an axis aligned box.
 * \param [out] poly The polygon to initialize.
 * \param [in] rbounds An array of two vectors defining
 *             the lower and upper corners of the box.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j2d_init_box(j2d_poly* poly, const j2d_vec2* const rbounds);

/**
 * \brief Copy the contents from one poly to another
 * \param [out] destination The poly we are copying to.
 * \param [in] source The poly we are copying from
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_copy_poly(j2d_poly* destination, const j2d_poly* const source);

/**
 * \brief Clip and cap polygon against an arbitrary number of planes.
 * \param [in, out] poly The polygon to be clipped.
 * \param [in] planes An array of planes to clip the polygon with.
 * \param [in] nplanes Number of planes.
 * \return status Code indicating if clip and cap was successful
 *         1 -> successful, 0 -> unsuccessful.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j2d_clip(j2d_poly* poly, const j2d_plane* const planes, const int nplanes);

/**
 * \brief Split and cap a list of polygon against a single plane.
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
int j2d_split(j2d_poly* inpolys,
              const int npolys,
              const j2d_plane plane,
              j2d_poly* out_pos,
              j2d_poly* out_neg);

/**
 * \brief computes the nth order moments of a polygon via
 *        the recursive method of Koehl.
 * \param [in] poly The polygon to compute the nth order moments over.
 * \param [out] moments Array containing the nth order moments.
 *              Must have a length of atleast (order+1)*(order+2)/2.
 *              Order of moments is row major, i.e. 1,x,y,x^2,xy,y^2,x^3,x^2y.
 * \param [out] order Order or the moments being computed.
 */
#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_get_moments(const j2d_poly* const poly, double* moments, int order);

#endif
