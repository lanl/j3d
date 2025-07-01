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

#ifndef _TEST_HELPERS_H_
#define _TEST_HELPERS_H_

#include <stdlib.h>
#include <stdbool.h>
#include "j3d.h"
#include "j2d.h"

#define ASSERT_INT_EQ(a, b) { \
  if ((a) != (b)) { \
    fprintf(stderr, "FAILURE: %d != %d in file %s on line %d.\n", \
                    a, b, __FILE__, __LINE__); \
    return EXIT_FAILURE; \
  } \
} \

#define ASSERT_DOUBLE_EQ(a, b, tol) { \
  if (fabs((a) - (b)) > (tol)) { \
    fprintf(stderr, "FAILURE: %f != %f with a tolerance of %f in file %s on line %d.\n", \
                    a, b, tol, __FILE__, __LINE__); \
    return EXIT_FAILURE; \
  } \
} \

#define ASSERT_DOUBLE_LT(a, b) { \
  if (!((a) < (b))) { \
    fprintf(stderr, "FAILURE: %f !< %f in file %s on line %d.\n", \
                    a, b, __FILE__, __LINE__); \
    return EXIT_FAILURE; \
  } \
} \

////////////////////////////
/// COMPARISON UTILITIES ///
////////////////////////////
bool vec_eq_3d(j3d_vec3 a, j3d_vec3 b, double tol);
bool vec_eq_2d(j2d_vec2 a, j2d_vec2 b, double tol);

/////////////////////////////////
// RANDOM GENERATION UTILITIES //
/////////////////////////////////
int rand_int(int N);
double rand_uniform();
j3d_vec3 rand_uvec_3d();
double rand_tet_3d(j3d_vec3* verts, double minvol);
j2d_vec2 rand_uvec_2d();
double rand_tri_2d(j2d_vec2* verts, double minarea);

/////////////////////////////////
///// J3D TESTING UTILITIES /////
/////////////////////////////////
j3d_plane thru_cent_3d(j3d_poly* poly);
j3d_plane thru_face_3d(j3d_poly* poly);
j3d_plane thru_edge_cent_3d(j3d_poly* poly);
j3d_plane thru_edge_rand_3d(j3d_poly* poly);
j3d_plane thru_vert_cent_3d(j3d_poly* poly);
j3d_plane thru_vert_rand_3d(j3d_poly* poly);
j3d_plane point_plane_3d(j3d_vec3 p0, j3d_vec3 p1, j3d_vec3 p2);

/////////////////////////////////
///// J2D TESTING UTILITIES /////
/////////////////////////////////
j2d_plane thru_cent_2d(j2d_poly* poly);
j2d_plane thru_edge_2d(j2d_poly* poly);
j2d_plane thru_vert_cent_2d(j2d_poly* poly);
j2d_plane thru_vert_rand_2d(j2d_poly* poly);
j2d_plane point_plane_2d(j2d_vec2 p0, j2d_vec2 p1);

#endif
