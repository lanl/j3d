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

#include "test_helpers.h"


////////////////////////////
/// COMPARISON UTILITIES ///
////////////////////////////
bool vec_eq_3d(j3d_vec3 a, j3d_vec3 b, double tol) {
  return fabs(a.xyz[0] - b.xyz[0]) < tol &&
         fabs(a.xyz[1] - b.xyz[1]) < tol &&
         fabs(a.xyz[2] - b.xyz[2]) < tol;
}

bool vec_eq_2d(j2d_vec2 a, j2d_vec2 b, double tol) {
  return fabs(a.xy[0] - b.xy[0]) < tol &&
         fabs(a.xy[1] - b.xy[1]) < tol;
}

/////////////////////////////////
// RANDOM GENERATION UTILITIES //
/////////////////////////////////
int rand_int(int N) {
  // generates random number in the interval [0,N-1]
  return rand() % N;
}

double rand_uniform() {
  return ((double)rand()) / ((double)RAND_MAX);
}

j3d_vec3 rand_uvec_3d() {
  // generates a random unit length vector
  j3d_vec3 tmp;
  tmp.xyz[0] = rand_uniform();
  tmp.xyz[1] = rand_uniform();
  tmp.xyz[2] = rand_uniform();
  j3d_norm(tmp.xyz);
  return tmp;
}

double rand_tet_3d(j3d_vec3* verts, double minvol) {
  // generates a random tetrahedron with vertices on the unit sphere,
  // guaranteeing a volume of at least MIN_VOL (to avoid degenerate cases)
  int v;
  j3d_vec3 swp;
  double tetvol = 0.0;
  while (tetvol < minvol) {
    for (v = 0; v < 4; ++v) {
      verts[v] = rand_uvec_3d();
    }
    tetvol = j3d_orient(verts);
    if (tetvol < 0.0) {
      swp = verts[2];
      verts[2] = verts[3];
      verts[3] = swp;
      tetvol = -tetvol;
    }
  }
  return tetvol;
}

j2d_vec2 rand_uvec_2d() {
  // generates a random unit length vector
  j2d_vec2 tmp;
  tmp.xy[0] = rand_uniform();
  tmp.xy[1] = rand_uniform();
  j2d_norm(tmp.xy);
  return tmp;
}

double rand_tri_2d(j2d_vec2* verts, double minarea) {
  // generates a random triangle with vertices on the unit circle,
  // guaranteeing an area of at least minarea (to avoid degenerate cases)
  int v;
  j2d_vec2 swp;
  double tetvol = 0.0;
  while (tetvol < minarea) {
    for (v = 0; v < 3; ++v) {
      verts[v] = rand_uvec_2d();
    }
    tetvol = j2d_orient(verts[0].xy, verts[1].xy, verts[2].xy);
    if (tetvol < 0.0) {
      swp = verts[1];
      verts[1] = verts[2];
      verts[2] = swp;
      tetvol = -tetvol;
    }
  }
  return tetvol;
}

/////////////////////////////////
///// J3D TESTING UTILITIES /////
/////////////////////////////////
j3d_plane thru_cent_3d(j3d_poly* poly) {
  // make a random plane through the center
  // of the poly
  double cent[3];
  j3d_center(cent, poly->verts, poly->num_verts);

  j3d_plane plane;
  plane.n.xyz[0] = rand_uniform();
  plane.n.xyz[1] = rand_uniform();
  plane.n.xyz[2] = rand_uniform();
  j3d_norm(plane.n.xyz);

  plane.d = -j3d_dot(cent, plane.n.xyz);

  return plane;
}

j3d_plane thru_face_3d(j3d_poly* poly) {
  // make a plane coplanar with a face of the poly
  int v0 = rand_int(poly->num_verts);
  int n_nbrs_v0 = poly->indices[v0 + 1] - poly->indices[v0];
  int pnbr1 = rand_int(n_nbrs_v0);
  int pnbr2 = (pnbr1 + 1) % n_nbrs_v0;
  int v1 = poly->dir_edges[poly->indices[v0] + pnbr1];
  int v2 = poly->dir_edges[poly->indices[v0] + pnbr2];

  j3d_vec3 p0, p1, p2;
  j3d_memcpy(p0.xyz, poly->verts[v0], 3 * sizeof(double));
  j3d_memcpy(p1.xyz, poly->verts[v1], 3 * sizeof(double));
  j3d_memcpy(p2.xyz, poly->verts[v2], 3 * sizeof(double));
  return point_plane_3d(p0, p1, p2);
}

j3d_plane thru_edge_cent_3d(j3d_poly* poly) {
  // make a plane coplanar with an edge and the centroid of the poly
  int v0 = rand_int(poly->num_verts);
  int pnbr1 = rand_int(poly->indices[v0 + 1] - poly->indices[v0]);
  int v1 = poly->dir_edges[poly->indices[v0] + pnbr1];

  j3d_vec3 p0, p1, centroid;
  j3d_memcpy(p0.xyz, poly->verts[v0], 3 * sizeof(double));
  j3d_memcpy(p1.xyz, poly->verts[v1], 3 * sizeof(double));
  j3d_center(centroid.xyz, poly->verts, poly->num_verts);
  return point_plane_3d(p0, p1, centroid);
}

j3d_plane thru_edge_rand_3d(j3d_poly* poly) {
  // make a plane coplanar with an edge, but otherwise randomly oriented
  int v0 = rand_int(poly->num_verts);
  int pnbr1 = rand_int(poly->indices[v0 + 1] - poly->indices[v0]);
  int v1 = poly->dir_edges[poly->indices[v0] + pnbr1];

  j3d_vec3 p0, p1;
  j3d_memcpy(p0.xyz, poly->verts[v0], 3 * sizeof(double));
  j3d_memcpy(p1.xyz, poly->verts[v1], 3 * sizeof(double));
  return point_plane_3d(p0, p1, rand_uvec_3d());
}

j3d_plane thru_vert_cent_3d(j3d_poly* poly) {
  // make a plane coplanar with a vertex and the centroid, but otherwise randomly oriented
  int v0 = rand_int(poly->num_verts);
  j3d_vec3 p0, centroid;
  j3d_memcpy(p0.xyz, poly->verts[v0], 3 * sizeof(double));
  j3d_center(centroid.xyz, poly->verts, poly->num_verts);
  return point_plane_3d(p0, centroid, rand_uvec_3d());
}

j3d_plane thru_vert_rand_3d(j3d_poly* poly) {
  // make a plane coplanar with a vertex, but otherwise randomly oriented
  int v0 = rand_int(poly->num_verts);
  j3d_vec3 p0;
  j3d_memcpy(p0.xyz, poly->verts[v0], 3 * sizeof(double));
  j3d_plane plane;
  plane.n = rand_uvec_3d();
  plane.d = -j3d_dot(plane.n.xyz, p0.xyz);
  return plane;
}

j3d_plane point_plane_3d(j3d_vec3 p0, j3d_vec3 p1, j3d_vec3 p2) {
  // generate a plane passing through the three points given
  const double tol = 1e-4;
  j3d_plane plane;
  if (vec_eq_3d(p0, p1, tol)) {
    p0 = rand_uvec_3d();
  }
  if (vec_eq_3d(p1, p2, tol)) {
    p1 = rand_uvec_3d();
  }
  if (vec_eq_3d(p0, p2, tol)) {
    p0 = rand_uvec_3d();
  }
  plane.n.xyz[0] = (p1.xyz[1] - p0.xyz[1]) * (p2.xyz[2] - p0.xyz[2]) -
                   (p1.xyz[2] - p0.xyz[2]) * (p2.xyz[1] - p0.xyz[1]);
  plane.n.xyz[1] = (p1.xyz[2] - p0.xyz[2]) * (p2.xyz[0] - p0.xyz[0]) -
                   (p1.xyz[0] - p0.xyz[0]) * (p2.xyz[2] - p0.xyz[2]);
  plane.n.xyz[2] = (p1.xyz[0] - p0.xyz[0]) * (p2.xyz[1] - p0.xyz[1]) -
                   (p1.xyz[1] - p0.xyz[1]) * (p2.xyz[0] - p0.xyz[0]);
  j3d_norm(plane.n.xyz);
  plane.d = -j3d_dot(plane.n.xyz, p0.xyz);
  return plane;
}


/////////////////////////////////
///// J3D TESTING UTILITIES /////
/////////////////////////////////
j2d_plane thru_cent_2d(j2d_poly* poly) {
  double cent[2];
  j2d_center(cent, poly->verts, poly->num_verts);

  j2d_plane plane;
  plane.n.xy[0] = rand_uniform();
  plane.n.xy[1] = rand_uniform();
  j2d_norm(plane.n.xy);

  plane.d = -j2d_dot(cent, plane.n.xy);

  return plane;
}

j2d_plane thru_edge_2d(j2d_poly* poly) {
  // make a plane coplanar with an edge and the centroid of the poly
  int v0 = rand_int(poly->num_verts);
  int pnbr = rand_int(2);
  int v1 = poly->dir_edges[2 * v0 + pnbr];

  j2d_vec2 p0, p1;
  j2d_memcpy(p0.xy, poly->verts[v0], 2 * sizeof(double));
  j2d_memcpy(p1.xy, poly->verts[v1], 2 * sizeof(double));
  return point_plane_2d(p0, p1);
}

j2d_plane thru_vert_cent_2d(j2d_poly* poly) {
  // make a plane coplanar with a vertex and the centroid, but otherwise randomly oriented
  int v0 = rand_int(poly->num_verts);
  j2d_vec2 p0, centroid;
  j2d_memcpy(p0.xy, poly->verts[v0], 2 * sizeof(double));
  j2d_center(centroid.xy, poly->verts, poly->num_verts);
  return point_plane_2d(p0, centroid);
}

j2d_plane thru_vert_rand_2d(j2d_poly* poly) {
  // make a plane coplanar with a vertex and the centroid, but otherwise randomly oriented
  int v0 = rand_int(poly->num_verts);
  j2d_vec2 p0;
  j2d_memcpy(p0.xy, poly->verts[v0], 2 * sizeof(double));
  j2d_plane plane;
  plane.n = rand_uvec_2d();
  plane.d = -j2d_dot(plane.n.xy, p0.xy);
  return plane;
}

j2d_plane point_plane_2d(j2d_vec2 p0, j2d_vec2 p1) {
  // generate a plane passing through the three points given
  const double tol = 1e-4;
  j2d_plane plane;
  if (vec_eq_2d(p0, p1, tol)) {
    p1 = rand_uvec_2d();
  }
  plane.n.xy[0] = p0.xy[1] - p1.xy[1];
  plane.n.xy[1] = p1.xy[0] - p0.xy[0];
  j2d_norm(plane.n.xy);
  plane.d = -j2d_dot(plane.n.xy, p0.xy);
  return plane;
}
