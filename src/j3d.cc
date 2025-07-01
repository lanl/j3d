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

#include "j3d.h"

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_memcpy(void* const dest, const void* const src, const int n) {
#ifdef J3D_ENABLE_Kokkos
  char* const d = (char*)dest;
  const char* const s = (const char*)src;
  for (int i = 0; i < n; ++i) {
    d[i] = s[i];
  }
#else
  memcpy(dest, src, n);
#endif
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_memset0(void* const ptr, const int n) {
#ifdef J3D_ENABLE_Kokkos
  char* const p = (char*)ptr;
  for (int i = 0; i < n; ++i) {
    p[i] = (char)0;
  }
#else
  memset(ptr, 0, n);
#endif
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
double j3d_dot(const double* const va, const double* const vb) {
  return va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_wav(double* const vr,
             const double* const va,
             const double wa,
             const double* const vb,
             const double wb) {
  const double rec = 1.0 / (wa + wb);
  vr[0] = (wa * va[0] + wb * vb[0]) * rec;
  vr[1] = (wa * va[1] + wb * vb[1]) * rec;
  vr[2] = (wa * va[2] + wb * vb[2]) * rec;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_norm(double* const vert) {
  const double rec = 1.0 / sqrt(j3d_dot(vert, vert));
  vert[0] *= rec;
  vert[1] *= rec;
  vert[2] *= rec;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_num_moments(const int order) {
  return (order + 1) * (order + 2) * (order + 3) / 6;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_shift_moments(double* const moments, int order, const double* const vc) {
  // var declarations
  int m, i, j, k, corder;
  int mm, mi, mj, mk, mcorder;

  // store moments of a shifted polygon
  double moments2[J3D_MAX_NUM_MOMENTS];
  j3d_memset0(moments2, j3d_num_moments(order) * sizeof(double));

  // calculate and save Pascal's triangle
  double B[J3D_MAX_ORDER + 1][J3D_MAX_ORDER + 1];
  B[0][0] = 1.0;
  for (corder = 1, m = 1; corder <= order; ++corder) {
    for (i = corder; i >= 0; --i, ++m) {
      j = corder - i;
      B[i][corder] = 1.0;
      if (i > 0 && j > 0) {
        B[i][corder] = B[i][corder - 1] + B[i - 1][corder - 1];
      }
    }
  }

  // shift moments back to the original position using
  // \int_\Omega x^i y^j z^k d\vec r =
  // \int_\omega (x+\xi)^i (y+\eta)^j (z+\zeta)^k d\vec r =
  // \sum_{a,b,c=0}^{i,j,k} \binom{i}{a} \binom{j}{b} \binom{k}{c}
  // \xi^{i-a} \eta^{j-b} \zeta^{k-c} \int_\omega x^a y^b z^c d\vec r
  for (corder = 1, m = 1; corder <= order; ++corder) {
    for (i = corder; i >= 0; --i) {
      for (j = corder - i; j >= 0; --j, ++m) {
        k = corder - i - j;
        for (mcorder = 0, mm = 0; mcorder <= corder; ++mcorder) {
          for (mi = mcorder; mi >= 0; --mi) {
            for (mj = mcorder - mi; mj >= 0; --mj, ++mm) {
              mk = mcorder - mi - mj;
              if (mi <= i && mj <= j && mk <= k) {
                moments2[m] += B[mi][i] * B[mj][j] * B[mk][k] * pow(vc[0], (i - mi)) *
                               pow(vc[1], (j - mj)) * pow(vc[2], (k - mk)) * moments[mm];
              }
            }
          }
        }
      }
    }
  }

  // assign shifted moments
  for (m = 1; m < j3d_num_moments(order); ++m) {
    moments[m] = moments2[m];
  }
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_center(double* const vc, const double (*const verts)[3], const int numverts) {
  vc[0] = 0.0;
  vc[1] = 0.0;
  vc[2] = 0.0;
  for (int i = 0; i < numverts; ++i) {
    vc[0] += verts[i][0];
    vc[1] += verts[i][1];
    vc[2] += verts[i][2];
  }
  const double rec = 1.0 / numverts;
  vc[0] *= rec;
  vc[1] *= rec;
  vc[2] *= rec;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_tet_faces_from_verts(j3d_plane* faces, const j3d_vec3* const verts) {
  double tmpcent[3];
  faces[0].n.xyz[0] = ((verts[3].xyz[1] - verts[1].xyz[1]) * (verts[2].xyz[2] - verts[1].xyz[2]) -
                       (verts[2].xyz[1] - verts[1].xyz[1]) * (verts[3].xyz[2] - verts[1].xyz[2]));
  faces[0].n.xyz[1] = ((verts[2].xyz[0] - verts[1].xyz[0]) * (verts[3].xyz[2] - verts[1].xyz[2]) -
                       (verts[3].xyz[0] - verts[1].xyz[0]) * (verts[2].xyz[2] - verts[1].xyz[2]));
  faces[0].n.xyz[2] = ((verts[3].xyz[0] - verts[1].xyz[0]) * (verts[2].xyz[1] - verts[1].xyz[1]) -
                       (verts[2].xyz[0] - verts[1].xyz[0]) * (verts[3].xyz[1] - verts[1].xyz[1]));
  j3d_norm(faces[0].n.xyz);
  tmpcent[0] = ONE_THIRD * (verts[1].xyz[0] + verts[2].xyz[0] + verts[3].xyz[0]);
  tmpcent[1] = ONE_THIRD * (verts[1].xyz[1] + verts[2].xyz[1] + verts[3].xyz[1]);
  tmpcent[2] = ONE_THIRD * (verts[1].xyz[2] + verts[2].xyz[2] + verts[3].xyz[2]);
  faces[0].d = -j3d_dot(faces[0].n.xyz, tmpcent);

  faces[1].n.xyz[0] = ((verts[2].xyz[1] - verts[0].xyz[1]) * (verts[3].xyz[2] - verts[2].xyz[2]) -
                       (verts[2].xyz[1] - verts[3].xyz[1]) * (verts[0].xyz[2] - verts[2].xyz[2]));
  faces[1].n.xyz[1] = ((verts[3].xyz[0] - verts[2].xyz[0]) * (verts[2].xyz[2] - verts[0].xyz[2]) -
                       (verts[0].xyz[0] - verts[2].xyz[0]) * (verts[2].xyz[2] - verts[3].xyz[2]));
  faces[1].n.xyz[2] = ((verts[2].xyz[0] - verts[0].xyz[0]) * (verts[3].xyz[1] - verts[2].xyz[1]) -
                       (verts[2].xyz[0] - verts[3].xyz[0]) * (verts[0].xyz[1] - verts[2].xyz[1]));
  j3d_norm(faces[1].n.xyz);
  tmpcent[0] = ONE_THIRD * (verts[2].xyz[0] + verts[3].xyz[0] + verts[0].xyz[0]);
  tmpcent[1] = ONE_THIRD * (verts[2].xyz[1] + verts[3].xyz[1] + verts[0].xyz[1]);
  tmpcent[2] = ONE_THIRD * (verts[2].xyz[2] + verts[3].xyz[2] + verts[0].xyz[2]);
  faces[1].d = -j3d_dot(faces[1].n.xyz, tmpcent);

  faces[2].n.xyz[0] = ((verts[1].xyz[1] - verts[3].xyz[1]) * (verts[0].xyz[2] - verts[3].xyz[2]) -
                       (verts[0].xyz[1] - verts[3].xyz[1]) * (verts[1].xyz[2] - verts[3].xyz[2]));
  faces[2].n.xyz[1] = ((verts[0].xyz[0] - verts[3].xyz[0]) * (verts[1].xyz[2] - verts[3].xyz[2]) -
                       (verts[1].xyz[0] - verts[3].xyz[0]) * (verts[0].xyz[2] - verts[3].xyz[2]));
  faces[2].n.xyz[2] = ((verts[1].xyz[0] - verts[3].xyz[0]) * (verts[0].xyz[1] - verts[3].xyz[1]) -
                       (verts[0].xyz[0] - verts[3].xyz[0]) * (verts[1].xyz[1] - verts[3].xyz[1]));
  j3d_norm(faces[2].n.xyz);
  tmpcent[0] = ONE_THIRD * (verts[3].xyz[0] + verts[0].xyz[0] + verts[1].xyz[0]);
  tmpcent[1] = ONE_THIRD * (verts[3].xyz[1] + verts[0].xyz[1] + verts[1].xyz[1]);
  tmpcent[2] = ONE_THIRD * (verts[3].xyz[2] + verts[0].xyz[2] + verts[1].xyz[2]);
  faces[2].d = -j3d_dot(faces[2].n.xyz, tmpcent);

  faces[3].n.xyz[0] = ((verts[0].xyz[1] - verts[2].xyz[1]) * (verts[1].xyz[2] - verts[0].xyz[2]) -
                       (verts[0].xyz[1] - verts[1].xyz[1]) * (verts[2].xyz[2] - verts[0].xyz[2]));
  faces[3].n.xyz[1] = ((verts[1].xyz[0] - verts[0].xyz[0]) * (verts[0].xyz[2] - verts[2].xyz[2]) -
                       (verts[2].xyz[0] - verts[0].xyz[0]) * (verts[0].xyz[2] - verts[1].xyz[2]));
  faces[3].n.xyz[2] = ((verts[0].xyz[0] - verts[2].xyz[0]) * (verts[1].xyz[1] - verts[0].xyz[1]) -
                       (verts[0].xyz[0] - verts[1].xyz[0]) * (verts[2].xyz[1] - verts[0].xyz[1]));
  j3d_norm(faces[3].n.xyz);
  tmpcent[0] = ONE_THIRD * (verts[0].xyz[0] + verts[1].xyz[0] + verts[2].xyz[0]);
  tmpcent[1] = ONE_THIRD * (verts[0].xyz[1] + verts[1].xyz[1] + verts[2].xyz[1]);
  tmpcent[2] = ONE_THIRD * (verts[0].xyz[2] + verts[1].xyz[2] + verts[2].xyz[2]);
  faces[3].d = -j3d_dot(faces[3].n.xyz, tmpcent);
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_box_faces_from_verts(j3d_plane* faces, const j3d_vec3* const bounds) {
  faces[0].n.xyz[0] = 0.0;
  faces[0].n.xyz[1] = 0.0;
  faces[0].n.xyz[2] = 1.0;
  faces[0].d = -bounds[0].xyz[2];
  faces[2].n.xyz[0] = 0.0;
  faces[2].n.xyz[1] = 1.0;
  faces[2].n.xyz[2] = 0.0;
  faces[2].d = -bounds[0].xyz[1];
  faces[4].n.xyz[0] = 1.0;
  faces[4].n.xyz[1] = 0.0;
  faces[4].n.xyz[2] = 0.0;
  faces[4].d = -bounds[0].xyz[0];
  faces[1].n.xyz[0] = 0.0;
  faces[1].n.xyz[1] = 0.0;
  faces[1].n.xyz[2] = -1.0;
  faces[1].d = bounds[1].xyz[2];
  faces[3].n.xyz[0] = 0.0;
  faces[3].n.xyz[1] = -1.0;
  faces[3].n.xyz[2] = 0.0;
  faces[3].d = bounds[1].xyz[1];
  faces[5].n.xyz[0] = -1.0;
  faces[5].n.xyz[1] = 0.0;
  faces[5].n.xyz[2] = 0.0;
  faces[5].d = bounds[1].xyz[0];
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
double j3d_orient(const j3d_vec3* const verts) {
  double adx, bdx, cdx;
  double ady, bdy, cdy;
  double adz, bdz, cdz;
  adx = verts[0].xyz[0] - verts[3].xyz[0];
  bdx = verts[1].xyz[0] - verts[3].xyz[0];
  cdx = verts[2].xyz[0] - verts[3].xyz[0];
  ady = verts[0].xyz[1] - verts[3].xyz[1];
  bdy = verts[1].xyz[1] - verts[3].xyz[1];
  cdy = verts[2].xyz[1] - verts[3].xyz[1];
  adz = verts[0].xyz[2] - verts[3].xyz[2];
  bdz = verts[1].xyz[2] - verts[3].xyz[2];
  cdz = verts[2].xyz[2] - verts[3].xyz[2];
  return -ONE_SIXTH * (adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) +
                       cdx * (ady * bdz - adz * bdy));
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_init_poly(j3d_poly* poly,
                  j3d_vec3* vertices,
                  int numverts,
                  int** faceinds,
                  int* numvertsperface,
                  int numfaces) {
  int v, vprev, vcur, vnext, f, np, numdiredges;

  if (numverts > J3D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_init_poly: Max vertex buffer size exceeded.  Increase J3D_MAX_VERTS.\n");
#endif
    return 0;
  }

  // direct access to vertex buffer
  double(*const vertbuffer)[3] = poly->verts;
  int* const diredgebuffer = poly->dir_edges;
  int* const indexbuffer = poly->indices;
  poly->num_verts = numverts;

  // for finding a starting point to build connectivity
  // for storing degrees
  int n_buffer[J3D_MAX_VERTS];

  // for building a map to construct neighbors
  int nn_buffer[J3D_MAX_VERTS * J3D_MAX_VERTS];

  int* const degree = n_buffer;
  j3d_memset0(degree, numverts * sizeof(int));

  for (f = 0; f < numfaces; ++f) {
    for (v = 0; v < numvertsperface[f]; ++v) {
      ++degree[faceinds[f][v]];
    }
  }

  for (v = 0; v < numverts; ++v) {
    j3d_memcpy(vertbuffer[v], vertices[v].xyz, 3 * sizeof(double));
  }

  const int* const d_read = n_buffer;
  numdiredges = 0;
  for (v = 0; v < numverts; ++v) {
    indexbuffer[v] = numdiredges;
    numdiredges += d_read[v];
  }

  if (numdiredges > J3D_MAX_EDGES * 2) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_init_poly: Max edge buffer size exceeded.  Increase J3D_MAX_EDGES.\n");
#endif
    return 0;
  }

  indexbuffer[numverts] = numdiredges;
  poly->num_dir_edges = numdiredges;

  // return if we were given an invalid poly
  for (v = 0; v < numverts; ++v) {
    if (d_read[v] < 3) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
      fprintf(stderr, "j3d_init_poly: Min degree is too small.\n");
#endif
      return 0;
    }
  }

  // build graph connectivity by correctly orienting directed edges for each
  // vertex
  // build a mapping that for every node and its next node
  // what is its previous node
  // build a mapping from a node to one of its neighbors
  // it doesnt matter
  // we just need a neighbor to start at
  int* const curr_next_to_prev = nn_buffer;
  int* const start = n_buffer;
  for (f = 0; f < numfaces; ++f) {
    const int* const fverts = faceinds[f];
    const int nvperf = numvertsperface[f];
    vprev = fverts[nvperf - 1];
    for (v = 0; v < nvperf - 1; ++v) {
      vcur = fverts[v];
      vnext = fverts[(v + 1)];
      curr_next_to_prev[vcur * numverts + vnext] = vprev;
      start[vcur] = vnext;
      vprev = vcur;
    }
    vnext = fverts[0];
    vcur = fverts[v];
    curr_next_to_prev[vcur * numverts + vnext] = vprev;
    start[vcur] = vnext;
  }

  const int* const M_next = nn_buffer;
  const int* const s_read = n_buffer;
  np = 0;
  for (v = 0; v < numverts; ++v) {
    vnext = s_read[v];
    do {
      diredgebuffer[np] = vnext;
      vnext = M_next[v * numverts + vnext];
      ++np;
    } while (vnext != s_read[v]);
  }

  return 1;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_brep_to_poly(j3d_poly* poly,
                     const j3d_vec3* const vertices,
                     const int numverts,
                     const int* const facevertinds,
                     const int* const facevertoffsets,
                     const int numfaces) {
  int v, vcur, vnext, f, np, numdiredges;

  if (numverts > J3D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_brep_to_poly: Max vertex buffer size exceeded.  Increase J3D_MAX_VERTS.\n");
#endif
    return 0;
  }

  // direct access to vertex buffer
  double(*const vertbuffer)[3] = poly->verts;
  int* const diredgebuffer = poly->dir_edges;
  int* const indexbuffer = poly->indices;
  poly->num_verts = numverts;
  
  // for counting degrees
  // for finding a starting point for each vertex to build connectivity
  int n_buffer[J3D_MAX_VERTS];

  // for map to build connectivity
  int nn_buffer[J3D_MAX_VERTS * J3D_MAX_VERTS];

  int* const degree = n_buffer;
  j3d_memset0(degree, numverts * sizeof(int));

  for (np = 0; np < facevertoffsets[numfaces]; ++np) {
    ++degree[facevertinds[np]];
  }

  for (v = 0; v < numverts; ++v) {
    j3d_memcpy(vertbuffer[v], vertices[v].xyz, 3 * sizeof(double));
  }

  const int* const d_read = n_buffer;
  numdiredges = 0;
  for (v = 0; v < numverts; ++v) {
    indexbuffer[v] = numdiredges;
    numdiredges += d_read[v];
  }

  if (numdiredges > J3D_MAX_EDGES * 2) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_brep_to_poly: Max edge buffer size exceeded.  Increase J3D_MAX_EDGES.\n");
#endif
    return 0;
  }

  indexbuffer[numverts] = numdiredges;
  poly->num_dir_edges = numdiredges;

  // return if we were given an invalid poly
  for (v = 0; v < numverts; ++v) {
    if (d_read[v] < 3) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
      fprintf(stderr, "j3d_brep_to_poly: Min degree is too small.\n");
#endif
      return 0;
    }
  }

  // build graph connectivity by correctly orienting directed edges for each
  // vertex

  // build a mapping that for every node and its next node
  // what is its previous node
  // build a mapping from a node to one of its neighbors
  // it doesnt matter
  // we just need a neighbor to start at
  int* const curr_next_to_prev = nn_buffer;
  int* const start = n_buffer;
  f = 0;
  np = 1;
  while (np < facevertoffsets[numfaces]) {
    vcur = facevertinds[np-1];
    vnext = facevertinds[np];
    curr_next_to_prev[vcur*numverts+vnext] = facevertinds[facevertoffsets[f+1]-1];
    start[vcur] = vnext;
    ++np;
    while (np != facevertoffsets[f + 1]) {
      vcur = facevertinds[np-1];
      vnext = facevertinds[np];
      curr_next_to_prev[vcur*numverts+vnext] = facevertinds[np-2];
      start[vcur] = vnext;
      ++np;
    }
    vcur = facevertinds[np - 1];
    vnext = facevertinds[facevertoffsets[f]];
    curr_next_to_prev[vcur*numverts+vnext] = facevertinds[np-2];
    start[vcur] = vnext;
    ++np;
    ++f;
  }

  const int* const M_next = nn_buffer;
  const int* const s_read = n_buffer;
  np = 0;
  for (v = 0; v < numverts; ++v) {
    vnext = s_read[v];
    do {
      diredgebuffer[np] = vnext;
      vnext = M_next[v * numverts + vnext];
      ++np;
    } while (vnext != s_read[v]);
  }

  return 1;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_init_tet(j3d_poly* poly, const j3d_vec3* const verts) {
  if (J3D_MAX_VERTS < 4) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_init_tet: J3D_MAX_VERTS is less than 4.  Increase J3D_MAX_VERTS.");
#endif
    return 0;
  }
  if (J3D_MAX_EDGES < 6) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_init_tet: J3D_MAX_EDGES is less than 4.  Increase J3D_MAX_EDGES.");
#endif
    return 0;
  }

  // direct access to vertex buffer
  double(*const vertbuffer)[3] = poly->verts;
  int* const diredgebuffer = poly->dir_edges;
  int* const indexbuffer = poly->indices;

  poly->num_dir_edges = 12;
  poly->num_verts = 4;

  diredgebuffer[0] = 1;
  diredgebuffer[1] = 3;
  diredgebuffer[2] = 2;
  diredgebuffer[3] = 2;
  diredgebuffer[4] = 3;
  diredgebuffer[5] = 0;
  diredgebuffer[6] = 0;
  diredgebuffer[7] = 3;
  diredgebuffer[8] = 1;
  diredgebuffer[9] = 1;
  diredgebuffer[10] = 2;
  diredgebuffer[11] = 0;

  indexbuffer[0] = 0;
  indexbuffer[1] = 3;
  indexbuffer[2] = 6;
  indexbuffer[3] = 9;
  indexbuffer[4] = 12;

  for (int v = 0; v < 4; ++v) {
    j3d_memcpy(vertbuffer[v], verts[v].xyz, 3 * sizeof(double));
  }
  return 1;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_init_box(j3d_poly* poly, const j3d_vec3* const rbounds) {
  if (J3D_MAX_VERTS < 8) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_init_tet: J3D_MAX_VERTS is less than 4.  Increase J3D_MAX_VERTS.");
#endif
    return 0;
  }
  if (J3D_MAX_EDGES < 12) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_init_tet: J3D_MAX_EDGES is less than 4.  Increase J3D_MAX_EDGES.");
#endif
    return 0;
  }

  // direct access to vertex buffer
  double(*const vertbuffer)[3] = poly->verts;
  int* const diredgebuffer = poly->dir_edges;
  int* const indexbuffer = poly->indices;

  poly->num_dir_edges = 24;
  poly->num_verts = 8;

  diredgebuffer[0] = 1;
  diredgebuffer[1] = 4;
  diredgebuffer[2] = 3;
  diredgebuffer[3] = 2;
  diredgebuffer[4] = 5;
  diredgebuffer[5] = 0;
  diredgebuffer[6] = 3;
  diredgebuffer[7] = 6;
  diredgebuffer[8] = 1;
  diredgebuffer[9] = 0;
  diredgebuffer[10] = 7;
  diredgebuffer[11] = 2;
  diredgebuffer[12] = 7;
  diredgebuffer[13] = 0;
  diredgebuffer[14] = 5;
  diredgebuffer[15] = 4;
  diredgebuffer[16] = 1;
  diredgebuffer[17] = 6;
  diredgebuffer[18] = 5;
  diredgebuffer[19] = 2;
  diredgebuffer[20] = 7;
  diredgebuffer[21] = 6;
  diredgebuffer[22] = 3;
  diredgebuffer[23] = 4;

  indexbuffer[0] = 0;
  indexbuffer[1] = 3;
  indexbuffer[2] = 6;
  indexbuffer[3] = 9;
  indexbuffer[4] = 12;
  indexbuffer[5] = 15;
  indexbuffer[6] = 18;
  indexbuffer[7] = 21;
  indexbuffer[8] = 24;

  vertbuffer[0][0] = rbounds[0].xyz[0];
  vertbuffer[0][1] = rbounds[0].xyz[1];
  vertbuffer[0][2] = rbounds[0].xyz[2];
  vertbuffer[1][0] = rbounds[1].xyz[0];
  vertbuffer[1][1] = rbounds[0].xyz[1];
  vertbuffer[1][2] = rbounds[0].xyz[2];
  vertbuffer[2][0] = rbounds[1].xyz[0];
  vertbuffer[2][1] = rbounds[1].xyz[1];
  vertbuffer[2][2] = rbounds[0].xyz[2];
  vertbuffer[3][0] = rbounds[0].xyz[0];
  vertbuffer[3][1] = rbounds[1].xyz[1];
  vertbuffer[3][2] = rbounds[0].xyz[2];
  vertbuffer[4][0] = rbounds[0].xyz[0];
  vertbuffer[4][1] = rbounds[0].xyz[1];
  vertbuffer[4][2] = rbounds[1].xyz[2];
  vertbuffer[5][0] = rbounds[1].xyz[0];
  vertbuffer[5][1] = rbounds[0].xyz[1];
  vertbuffer[5][2] = rbounds[1].xyz[2];
  vertbuffer[6][0] = rbounds[1].xyz[0];
  vertbuffer[6][1] = rbounds[1].xyz[1];
  vertbuffer[6][2] = rbounds[1].xyz[2];
  vertbuffer[7][0] = rbounds[0].xyz[0];
  vertbuffer[7][1] = rbounds[1].xyz[1];
  vertbuffer[7][2] = rbounds[1].xyz[2];

  return 0;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_copy_poly(j3d_poly* destination, const j3d_poly* const source) {
  destination->num_dir_edges = source->num_dir_edges;
  destination->num_verts = source->num_verts;

  j3d_memcpy(destination->verts, source->verts[0], 3 * source->num_verts * sizeof(double));
  j3d_memcpy(destination->dir_edges, source->dir_edges, source->num_dir_edges * sizeof(int));
  j3d_memcpy(destination->indices, source->indices, (source->num_verts + 1) * sizeof(int));
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_clip(j3d_poly* poly, const j3d_plane* const planes, const int nplanes) {
  // direct access to poly buffers
  double(*const vertbuffer)[3] = poly->verts;
  int* const indexbuffer = poly->indices;
  int* const diredgebuffer = poly->dir_edges;
  int* const ndiredges = &poly->num_dir_edges;
  int* const nverts = &poly->num_verts;
  if (*nverts <= 0) return 0;

  // variable declarations
  int v, p, np, onv, onde, vprev, vcur, vnext, vstart, numunclipped, numclipped, numnew, vnext_idx, degree_v;

  // for temporarily storing vertices
  // for storing signed distances
  double d_buffer[J3D_MAX_VERTS * 3];

  // for map to cache next vertex in brep walk
  // for marking clipped verts
  // for compression mapping
  // for temporarily storing directed edges
  int buffer[J3D_MAX_VERTS * J3D_MAX_VERTS];

  // loop over each clip plane
  for (p = 0; p < nplanes; ++p) {
    // calculate signed distances to the clip plane
    onv = *nverts;
    onde = *ndiredges;
    numclipped = 0;
    double* const sdists = d_buffer;
    int* const clipped = buffer;
    for (v = 0; v < onv; ++v) {
      sdists[v] = planes[p].d + j3d_dot(vertbuffer[v], planes[p].n.xyz);
      np = (sdists[v] < 0.0);
      clipped[v] = np;
      numclipped += np;
    }

    // skip this face if the poly lies entirely on one side of it
    if (numclipped == 0) {
      continue;
    }
    if (numclipped == onv) {
      *nverts = 0;
      *ndiredges = 0;
      return 1;
    }

    // count number of new vertices
    v = 0;
    numnew = 0;
    int* const temp_diredgebuffer = &buffer[J3D_MAX_VERTS];
    for (np = 0; np < *ndiredges; ++np) {
      if (np == indexbuffer[v + 1]) { ++v; }
      temp_diredgebuffer[2 * numnew] = v;
      temp_diredgebuffer[2 * numnew + 1] = np;
      numnew += clipped[diredgebuffer[np]] * (1 - clipped[v]);
    }

    // check bounds are not exceeded
    if (*nverts + numnew > J3D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
      fprintf(stderr,
        "j3d_clip: Max vertex buffer size exceeded.  Increase J3D_MAX_VERTS.\n");
#endif
      return 0;
    }
    if (*ndiredges + (3 * numnew) > J3D_MAX_EDGES * 2) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
      fprintf(stderr,
        "j3d_clip: Max edge buffer size exceeded.  Increase J3D_MAX_EDGES.\n");
#endif
      return 0;
    }

    // check all edges and insert new vertices on the bisected edges
    j3d_memset0(&diredgebuffer[indexbuffer[*nverts]], numnew * 3 * sizeof(int));
    for (np = 0; np < numnew; ++np) {
      v = temp_diredgebuffer[2 * np];
      vnext_idx = temp_diredgebuffer[2 * np + 1];
      vnext = diredgebuffer[vnext_idx];

      j3d_wav(vertbuffer[onv + np], vertbuffer[v], -sdists[vnext], vertbuffer[vnext], sdists[v]);
      diredgebuffer[onde + 3 * np] = v;
      diredgebuffer[vnext_idx] = onv + np;
    }

    // go through and compress the vertex list, removing clipped verts and
    // re-indexing accordingly `clipped` will be used to map uncompressed vertex indices
    // start by compressiong old vertices
    numunclipped = 0;
    *ndiredges = 0;
    double* const temp_vertbuffer = d_buffer;
    for (v = 0; v < onv; ++v) {
      if (!clipped[v]) {
        degree_v = indexbuffer[v + 1] - indexbuffer[v];
        j3d_memcpy(&temp_diredgebuffer[*ndiredges], &diredgebuffer[indexbuffer[v]], degree_v * sizeof(int));
        indexbuffer[numunclipped] = *ndiredges;
        *ndiredges += degree_v;
        j3d_memcpy(&temp_vertbuffer[numunclipped * 3], vertbuffer[v], 3 * sizeof(double));
        clipped[v] = numunclipped++;
      }
    }

    // now compress new vertices
    // they all have degree 3
    j3d_memcpy(&temp_diredgebuffer[*ndiredges], &diredgebuffer[indexbuffer[onv]], 3 * numnew * sizeof(int));
    j3d_memcpy(&temp_vertbuffer[numunclipped * 3], vertbuffer[onv], 3 * numnew * sizeof(double));
    indexbuffer[numunclipped] = *ndiredges;
    for (v = onv; v < onv + numnew; ++v) {
      indexbuffer[numunclipped + 1] = indexbuffer[numunclipped] + 3;
      clipped[v] = numunclipped++;
    }
    *ndiredges += (3 * numnew);
    *nverts = numunclipped;
    // update neighbor indices after compression
    const int* const map = buffer;
    for (np = 0; np < *ndiredges; ++np) {
      temp_diredgebuffer[np] = map[temp_diredgebuffer[np]];
    }
    // copy back to original poly buffers
    j3d_memcpy(diredgebuffer, temp_diredgebuffer, (*ndiredges) * sizeof(int));
    j3d_memcpy(vertbuffer, temp_vertbuffer, 3 * (*nverts) * sizeof(double));

    // We will reset the value onv for cap
    onv -= numclipped;

    // build a map akin to what we see in reduce
    // to efficiently walk the open faces during cap.
    v = 0;
    np = 1;
    int* const prev_curr_to_next = buffer;
    while (np < indexbuffer[onv]) {
      while (np != indexbuffer[v + 1]) {
        prev_curr_to_next[diredgebuffer[np] * (*nverts) + v] = diredgebuffer[np - 1];
        ++np;
      }
      prev_curr_to_next[diredgebuffer[indexbuffer[v]] * (*nverts) + v] = diredgebuffer[indexbuffer[v + 1] - 1];
      ++v;
      ++np;
    }

    // for each new vert, search around the faces for its new neighbors and
    // doubly-link everything
    const int* const M_next = buffer;
    for (vstart = onv; vstart < *nverts; ++vstart) {
      vprev = vstart;
      vcur = diredgebuffer[indexbuffer[vprev]];
      vnext = prev_curr_to_next[vprev * (*nverts) + vcur];
      while (vnext < onv) {
        vprev = vcur;
        vcur = vnext;
        vnext = M_next[vprev * (*nverts) + vcur];
      }
      // all new vertices have degree 3 so this is fine
      // in n valence case
      diredgebuffer[indexbuffer[vstart] + 1] = vnext;
      diredgebuffer[indexbuffer[vnext] + 2] = vstart;
    }
  }

  return 1;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j3d_split(j3d_poly* inpolys,
              const int npolys,
              const j3d_plane plane,
              j3d_poly* out_pos,
              j3d_poly* out_neg) {
  int v, np, npnxt, onv, vcur, vnext, vstart, pnext, nright, cside, p;
  double newpos[3];

  int side[J3D_MAX_VERTS];
  double sdists[J3D_MAX_VERTS];

  double(*vertbuffer)[3];
  int* diredgebuffer;
  int* indexbuffer;
  int* ndiredges;
  int* nverts;

  j3d_poly* outpolys[2];

  for (p = 0; p < npolys; ++p) {
    vertbuffer = inpolys[p].verts;
    diredgebuffer = inpolys[p].dir_edges;
    indexbuffer = inpolys[p].indices;
    ndiredges = &inpolys[p].num_dir_edges;
    nverts = &inpolys[p].num_verts;

    outpolys[0] = &out_pos[p];
    outpolys[1] = &out_neg[p];
    if (*nverts <= 0) {
      out_pos[p].num_verts = 0;
      out_pos[p].num_dir_edges = 0;
      out_neg[p].num_verts = 0;
      out_neg[p].num_dir_edges = 0;
      continue;
    }

    // calculate signed distances to the clip plane
    nright = 0;
    j3d_memset0(side, J3D_MAX_VERTS * sizeof(int));
    for (v = 0; v < *nverts; ++v) {
      sdists[v] = plane.d + j3d_dot(vertbuffer[v], plane.n.xyz);
      if (sdists[v] < 0.0) {
        side[v] = 1;
        nright++;
      }
    }

    // return if the poly lies entirely on one side of it
    if (nright == 0) {
      out_pos[p] = inpolys[p];
      out_neg[p].num_verts = 0;
      out_neg[p].num_dir_edges = 0;
      continue;
    }
    if (nright == *nverts) {
      out_neg[p] = inpolys[p];
      out_pos[p].num_verts = 0;
      out_pos[p].num_dir_edges = 0;
      continue;
    }

    // check all edges and insert new vertices on the bisected edges
    onv = *nverts;
    for (vcur = 0; vcur < onv; ++vcur) {
      if (side[vcur]) continue;
      const int vcur_degree = indexbuffer[vcur + 1] - indexbuffer[vcur];
      for (np = 0; np < vcur_degree; ++np) {
        vnext = diredgebuffer[indexbuffer[vcur] + np];
        if (!side[vnext]) continue;
        j3d_wav(newpos, vertbuffer[vcur], -sdists[vnext], vertbuffer[vnext], sdists[vcur]);
        if (*nverts == J3D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
          fprintf(stderr, "j3d_split: Max vertex buffer size exceeded.  Increase J3D_MAX_VERTS.\n");
#endif
          return 0;
        }
        j3d_memcpy(vertbuffer[*nverts], newpos, 3 * sizeof(double));
        indexbuffer[*nverts] = *ndiredges;
        *ndiredges += 3;
        if (*ndiredges > J3D_MAX_EDGES * 2) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
          fprintf(stderr, "j3d_split: Max edge buffer exceeded.  Increase J3D_MAX_EDGES.\n");
#endif
          return 0;
        }
        diredgebuffer[indexbuffer[*nverts]] = vcur;
        diredgebuffer[indexbuffer[vcur] + np] = *nverts;
        ++(*nverts);

        if (*nverts == J3D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
          fprintf(stderr, "j3d_split: Max vertex buffer size exceeded.  Increase J3D_MAX_VERTS.\n");
#endif
          return 0;
        }
        j3d_memcpy(vertbuffer[*nverts], newpos, 3 * sizeof(double));
        side[*nverts] = 1;
        indexbuffer[*nverts] = *ndiredges;
        *ndiredges += 3;
        if (*ndiredges > J3D_MAX_EDGES) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
          fprintf(stderr, "j3d_split: Max edge buffer exceeded.  Increase J3D_MAX_EDGES.\n");
#endif
          return 0;
        }
        diredgebuffer[indexbuffer[*nverts]] = vnext;
        int* const vnext_nbrs = &diredgebuffer[indexbuffer[vnext]];
        npnxt = 0;
        while (vnext_nbrs[npnxt] != vcur) {
          ++npnxt;
        }
        vnext_nbrs[npnxt] = *nverts;
        ++(*nverts);
      }
    }
    indexbuffer[*nverts] = *ndiredges;

    // for each new vert, search around the faces for its new neighbors and
    // doubly-link everything
    for (vstart = onv; vstart < *nverts; ++vstart) {
      vcur = vstart;
      vnext = diredgebuffer[indexbuffer[vcur]];
      do {
        const int* const nbrs = &diredgebuffer[indexbuffer[vnext]];
        np = 0;
        while (nbrs[np] != vcur) {
          ++np;
        }
        vcur = vnext;
        pnext = (np + 1) % (indexbuffer[vnext + 1] - indexbuffer[vnext]);
        vnext = diredgebuffer[indexbuffer[vcur] + pnext];
      } while (vcur < onv);
      diredgebuffer[indexbuffer[vstart] + 2] = vcur;
      diredgebuffer[indexbuffer[vcur] + 1] = vstart;
    }

    // copy and compress vertices into their new buffers reusing side[] for
    // reindexing
    onv = *nverts;
    outpolys[0]->num_verts = 0;
    outpolys[1]->num_verts = 0;
    outpolys[0]->num_dir_edges = 0;
    outpolys[1]->num_dir_edges = 0;
    for (v = 0; v < onv; ++v) {
      cside = side[v];

      // rewrite indices and directed edges
      const int degree_v = indexbuffer[v + 1] - indexbuffer[v];

      j3d_memcpy(&(outpolys[cside]->dir_edges[outpolys[cside]->num_dir_edges]),
             &diredgebuffer[indexbuffer[v]],
             degree_v * sizeof(int));

      outpolys[cside]->indices[outpolys[cside]->num_verts] = outpolys[cside]->num_dir_edges;
      outpolys[cside]->num_dir_edges += degree_v;

      // rewrite vertices
      j3d_memcpy(outpolys[cside]->verts[outpolys[cside]->num_verts], vertbuffer[v], 3 * sizeof(double));

      // compression mapping
      side[v] = (outpolys[cside]->num_verts)++;
    }
    outpolys[0]->indices[outpolys[0]->num_verts] = outpolys[0]->num_dir_edges;
    outpolys[1]->indices[outpolys[1]->num_verts] = outpolys[1]->num_dir_edges;

    for (np = 0; np < outpolys[0]->num_dir_edges; ++np) {
      outpolys[0]->dir_edges[np] = side[outpolys[0]->dir_edges[np]];
    }
    for (np = 0; np < outpolys[1]->num_dir_edges; ++np) {
      outpolys[1]->dir_edges[np] = side[outpolys[1]->dir_edges[np]];
    }
  }

  return 1;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_0(const j3d_poly* const poly, double* moments) {
  // direct access to poly data
  const double(*const vertbuffer)[3] = poly->verts;
  const int* const indexbuffer = poly->indices;
  const int* const diredgebuffer = poly->dir_edges;
  const int ndiredges = poly->num_dir_edges;
  const int nverts = poly->num_verts;

  // variable declarations
  int vprev, vcurr, vnext, v, n;

  // some 3D arrays to reuse while fanning
  double v0[3];
  double v1[3];
  double v2[3];
  // center we will use for moments shifting
  double vc[3];

  // moments in the zeroth order case should be size 1
  j3d_memset0(moments, 1 * sizeof(double));

  // if empty, exit
  if (nverts <= 0) {
    return;
  }

  // compute center of poly for shifting
  j3d_center(vc, vertbuffer, nverts);

  // for tracking visited edges
  int emarks[J3D_MAX_VERTS * J3D_MAX_VERTS];
  j3d_memset0(emarks, nverts * nverts * sizeof(int));

  // for intelligent graph walk
  // akin to what we saw in init poly quadratic
  int prev_curr_to_next[J3D_MAX_VERTS * J3D_MAX_VERTS];
  v = 0;
  n = 1;
  while (n < ndiredges) {
    while (n != indexbuffer[v + 1]) {
      prev_curr_to_next[diredgebuffer[n] * nverts + v] = diredgebuffer[n - 1];
      ++n;
    }
    prev_curr_to_next[diredgebuffer[indexbuffer[v]] * nverts + v] = diredgebuffer[indexbuffer[v + 1] - 1];
    ++v;
    ++n;
  }

  // for each node, check each of its neighbors
  // and use the "adjacency matrix" emarks to
  // check if we have walked that face yet.
  // If not, walk it
  v = 0;
  const int* const M_next = prev_curr_to_next;
  for (n = 0; n < ndiredges; ++n) {
    if (n == indexbuffer[v + 1]) {
      ++v;
    }
    // check if we encountered this face
    if (emarks[v * nverts + diredgebuffer[n]]) {
      continue;
    }

    // if this face hasn't been visited then fan it out
    vprev = v;
    vcurr = diredgebuffer[n];
    vnext = M_next[vprev * nverts + vcurr];
    j3d_memcpy(v0, vertbuffer[vprev], 3 * sizeof(double));
    j3d_memcpy(v1, vertbuffer[vcurr], 3 * sizeof(double));
    emarks[v * nverts + vcurr] = 1;

    // shift point we are "fanning"
    // from
    v0[0] = v0[0] - vc[0];
    v0[1] = v0[1] - vc[1];
    v0[2] = v0[2] - vc[2];

    v1[0] = v1[0] - vc[0];
    v1[1] = v1[1] - vc[1];
    v1[2] = v1[2] - vc[2];

    while (vnext != v) {
      j3d_memcpy(v2, vertbuffer[vnext], 3 * sizeof(double));

      v2[0] = v2[0] - vc[0];
      v2[1] = v2[1] - vc[1];
      v2[2] = v2[2] - vc[2];

      const double v00 = v0[0];
      const double v01 = v0[1];
      const double v02 = v0[2];
      const double v10 = v1[0];
      const double v11 = v1[1];
      const double v12 = v1[2];
      const double v20 = v2[0];
      const double v21 = v2[1];
      const double v22 = v2[2];

      moments[0] += (-v20 * v11 * v02 + v10 * v21 * v02 + v20 * v01 * v12 - v00 * v21 * v12 -
                     v10 * v01 * v22 + v00 * v11 * v22);

      vprev = vcurr;
      vcurr = vnext;
      vnext = M_next[vprev * nverts + vcurr];
      emarks[vprev * nverts + vcurr] = 1;
      j3d_memcpy(v1, v2, 3 * sizeof(double));
    }
    emarks[vcurr * nverts + vnext] = 1;
  }
  moments[0] /= 6.0;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_1(const j3d_poly* const poly, double* moments) {
  // direct access to poly data
  const double(*const vertbuffer)[3] = poly->verts;
  const int* const indexbuffer = poly->indices;
  const int* const diredgebuffer = poly->dir_edges;
  const int ndiredges = poly->num_dir_edges;
  const int nverts = poly->num_verts;

  // variable declarations
  int vprev, vcurr, vnext, v, n;

  // some 3D arrays to reuse while fanning
  double v0[3];
  double v1[3];
  double v2[3];
  // center we will use for moments shifting
  double vc[3];

  // moments in the first order case should be size 4
  j3d_memset0(moments, 4 * sizeof(double));

  // if empty, exit
  if (nverts <= 0) {
    return;
  }

  // compute center of poly for shifting
  j3d_center(vc, vertbuffer, nverts);

  // for tracking visited edges
  int emarks[J3D_MAX_VERTS * J3D_MAX_VERTS];
  j3d_memset0(emarks, nverts * nverts * sizeof(int));

  // for intelligent graph walk
  // akin to what we saw in init poly quadratic
  int prev_curr_to_next[J3D_MAX_VERTS * J3D_MAX_VERTS];
  v = 0;
  n = 1;
  while (n < ndiredges) {
    while (n != indexbuffer[v + 1]) {
      prev_curr_to_next[diredgebuffer[n] * nverts + v] = diredgebuffer[n - 1];
      ++n;
    }
    prev_curr_to_next[diredgebuffer[indexbuffer[v]] * nverts + v] = diredgebuffer[indexbuffer[v + 1] - 1];
    ++v;
    ++n;
  }

  // for each node, check each of its neighbors
  // and use the "adjacency matrix" emarks to
  // check if we have walked that face yet.
  // If not, walk it
  v = 0;
  const int* const M_next = prev_curr_to_next;
  for (n = 0; n < ndiredges; ++n) {
    if (n == indexbuffer[v + 1]) {
      ++v;
    }
    // check if we encountered this face
    if (emarks[v * nverts + diredgebuffer[n]]) {
      continue;
    }

    // if this face hasn't been visited then fan it out
    vprev = v;
    vcurr = diredgebuffer[n];
    vnext = M_next[vprev * nverts + vcurr];
    j3d_memcpy(v0, vertbuffer[vprev], 3 * sizeof(double));
    j3d_memcpy(v1, vertbuffer[vcurr], 3 * sizeof(double));
    emarks[v * nverts + vcurr] = 1;

    // shift point we are "fanning"
    // from
    v0[0] = v0[0] - vc[0];
    v0[1] = v0[1] - vc[1];
    v0[2] = v0[2] - vc[2];

    v1[0] = v1[0] - vc[0];
    v1[1] = v1[1] - vc[1];
    v1[2] = v1[2] - vc[2];

    while (vnext != v) {
      j3d_memcpy(v2, vertbuffer[vnext], 3 * sizeof(double));

      v2[0] = v2[0] - vc[0];
      v2[1] = v2[1] - vc[1];
      v2[2] = v2[2] - vc[2];

      const double v00 = v0[0];
      const double v01 = v0[1];
      const double v02 = v0[2];
      const double v10 = v1[0];
      const double v11 = v1[1];
      const double v12 = v1[2];
      const double v20 = v2[0];
      const double v21 = v2[1];
      const double v22 = v2[2];

      const double sixv = (-v20 * v11 * v02 + v10 * v21 * v02 + v20 * v01 * v12 - v00 * v21 * v12 -
                           v10 * v01 * v22 + v00 * v11 * v22);

      moments[0] += sixv;
      moments[1] += sixv * (v00 + v10 + v20);
      moments[2] += sixv * (v01 + v11 + v21);
      moments[3] += sixv * (v02 + v12 + v22);

      vprev = vcurr;
      vcurr = vnext;
      vnext = M_next[vprev * nverts + vcurr];
      j3d_memcpy(v1, v2, 3 * sizeof(double));
      emarks[vprev * nverts + vcurr] = 1;
    }
    emarks[vcurr * nverts + vnext] = 1;
  }
  moments[0] /= 6.0;
  moments[1] /= 24.0;
  moments[2] /= 24.0;
  moments[3] /= 24.0;

  // shift moments since only zeroth order moment is invariant
  // to translation
  const double vc0 = vc[0];
  const double vc1 = vc[1];
  const double vc2 = vc[2];

  double moments2[4];

  moments2[1] = vc0 * moments[0] + moments[1];
  moments2[2] = vc1 * moments[0] + moments[2];
  moments2[3] = vc2 * moments[0] + moments[3];

  j3d_memcpy(&moments[1], &moments2[1], 3 * sizeof(double));
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_2(const j3d_poly* const poly, double* moments) {
  // direct access to poly data
  const double(*const vertbuffer)[3] = poly->verts;
  const int* const indexbuffer = poly->indices;
  const int* const diredgebuffer = poly->dir_edges;
  const int ndiredges = poly->num_dir_edges;
  const int nverts = poly->num_verts;

  // variable declarations
  int vprev, vcurr, vnext, v, n;

  // some 3D arrays to reuse while fanning
  double v0[3];
  double v1[3];
  double v2[3];
  // center we will use for moments shifting
  double vc[3];

  // moments in the first order case should be size 10
  j3d_memset0(moments, 10 * sizeof(double));

  // if empty, exit
  if (nverts <= 0) {
    return;
  }

  // compute center of poly for shifting
  j3d_center(vc, vertbuffer, nverts);

  // for tracking visited edges
  int emarks[J3D_MAX_VERTS * J3D_MAX_VERTS];
  j3d_memset0(emarks, nverts * nverts * sizeof(int));

  // for intelligent graph walk
  // akin to what we saw in init poly quadratic
  int prev_curr_to_next[J3D_MAX_VERTS * J3D_MAX_VERTS];
  v = 0;
  n = 1;
  while (n < ndiredges) {
    while (n != indexbuffer[v + 1]) {
      prev_curr_to_next[diredgebuffer[n] * nverts + v] = diredgebuffer[n - 1];
      ++n;
    }
    prev_curr_to_next[diredgebuffer[indexbuffer[v]] * nverts + v] = diredgebuffer[indexbuffer[v + 1] - 1];
    ++v;
    ++n;
  }

  // for each node, check each of its neighbors
  // and use the "adjacency matrix" emarks to
  // check if we have walked that face yet.
  // If not, walk it
  v = 0;
  const int* const M_next = prev_curr_to_next;
  for (n = 0; n < ndiredges; ++n) {
    if (n == indexbuffer[v + 1]) {
      ++v;
    }
    // check if we encountered this face
    if (emarks[v * nverts + diredgebuffer[n]]) {
      continue;
    }

    // if this face hasn't been visited then fan it out
    vprev = v;
    vcurr = diredgebuffer[n];
    vnext = M_next[vprev * nverts + vcurr];
    j3d_memcpy(v0, vertbuffer[vprev], 3 * sizeof(double));
    j3d_memcpy(v1, vertbuffer[vcurr], 3 * sizeof(double));
    emarks[v * nverts + vcurr] = 1;

    // shift point we are "fanning"
    // from
    v0[0] = v0[0] - vc[0];
    v0[1] = v0[1] - vc[1];
    v0[2] = v0[2] - vc[2];

    v1[0] = v1[0] - vc[0];
    v1[1] = v1[1] - vc[1];
    v1[2] = v1[2] - vc[2];

    while (vnext != v) {
      j3d_memcpy(v2, vertbuffer[vnext], 3 * sizeof(double));

      v2[0] = v2[0] - vc[0];
      v2[1] = v2[1] - vc[1];
      v2[2] = v2[2] - vc[2];

      const double v00 = v0[0];
      const double v01 = v0[1];
      const double v02 = v0[2];
      const double v10 = v1[0];
      const double v11 = v1[1];
      const double v12 = v1[2];
      const double v20 = v2[0];
      const double v21 = v2[1];
      const double v22 = v2[2];

      const double sixv = (-v20 * v11 * v02 + v10 * v21 * v02 + v20 * v01 * v12 - v00 * v21 * v12 -
                           v10 * v01 * v22 + v00 * v11 * v22);

      moments[0] += sixv;
      moments[1] += sixv * (v00 + v10 + v20);
      moments[2] += sixv * (v01 + v11 + v21);
      moments[3] += sixv * (v02 + v12 + v22);
      moments[4] += sixv * (v00 * v00 + v00 * v10 + v10 * v10 + v00 * v20 + v10 * v20 + v20 * v20);
      moments[5] += sixv * (2.0 * v00 * v01 + v01 * v10 + v00 * v11 + 2.0 * v10 * v11 + v01 * v20 +
                            v11 * v20 + v00 * v21 + v10 * v21 + 2.0 * v20 * v21);
      moments[6] += sixv * (2.0 * v00 * v02 + v02 * v10 + v00 * v12 + 2.0 * v10 * v12 + v02 * v20 +
                            v12 * v20 + v00 * v22 + v10 * v22 + 2.0 * v20 * v22);
      moments[7] += sixv * (v01 * v01 + v01 * v11 + v11 * v11 + v01 * v21 + v11 * v21 + v21 * v21);
      moments[8] += sixv * (2.0 * v01 * v02 + v02 * v11 + v01 * v12 + 2.0 * v11 * v12 + v02 * v21 +
                            v12 * v21 + v01 * v22 + v11 * v22 + 2.0 * v21 * v22);
      moments[9] += sixv * (v02 * v02 + v02 * v12 + v12 * v12 + v02 * v22 + v12 * v22 + v22 * v22);

      vprev = vcurr;
      vcurr = vnext;
      vnext = M_next[vprev * nverts + vcurr];
      j3d_memcpy(v1, v2, 3 * sizeof(double));
      emarks[vprev * nverts + vcurr] = 1;
    }
    emarks[vcurr * nverts + vnext] = 1;
  }
  moments[0] /= 6.0;
  moments[1] /= 24.0;
  moments[2] /= 24.0;
  moments[3] /= 24.0;
  moments[4] /= 60.0;
  moments[5] /= 120.0;
  moments[6] /= 120.0;
  moments[7] /= 60.0;
  moments[8] /= 120.0;
  moments[9] /= 60.0;

  // shift moments since only zeroth order moment is invariant
  // to translation
  const double vc0 = vc[0];
  const double vc1 = vc[1];
  const double vc2 = vc[2];

  double moments2[10];

  moments2[1] = vc0 * moments[0] + moments[1];
  moments2[2] = vc1 * moments[0] + moments[2];
  moments2[3] = vc2 * moments[0] + moments[3];
  moments2[4] = vc0 * vc0 * moments[0] + 2.0 * vc0 * moments[1] + moments[4];
  moments2[5] = vc0 * vc1 * moments[0] + vc0 * moments[2] + vc1 * moments[1] + moments[5];
  moments2[6] = vc0 * vc2 * moments[0] + vc0 * moments[3] + vc2 * moments[1] + moments[6];
  moments2[7] = vc1 * vc1 * moments[0] + 2.0 * vc1 * moments[2] + moments[7];
  moments2[8] = vc1 * vc2 * moments[0] + vc1 * moments[3] + vc2 * moments[2] + moments[8];
  moments2[9] = vc2 * vc2 * moments[0] + 2.0 * vc2 * moments[3] + moments[9];

  j3d_memcpy(&moments[1], &moments2[1], 9 * sizeof(double));
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_3(const j3d_poly* const poly, double* moments) {
  // direct access to poly data
  const double(*const vertbuffer)[3] = poly->verts;
  const int* const indexbuffer = poly->indices;
  const int* const diredgebuffer = poly->dir_edges;
  const int ndiredges = poly->num_dir_edges;
  const int nverts = poly->num_verts;

  // variable declarations
  int vprev, vcurr, vnext, v, n;

  // some 3D arrays to reuse while fanning
  double v0[3];
  double v1[3];
  double v2[3];
  // center we will use for moments shifting
  double vc[3];

  // moments in the first order case should be size 20
  j3d_memset0(moments, 20 * sizeof(double));

  // if empty, exit
  if (nverts <= 0) {
    return;
  }

  // compute center of poly for shifting
  j3d_center(vc, vertbuffer, nverts);

  // for tracking visited edges
  int emarks[J3D_MAX_VERTS * J3D_MAX_VERTS];
  j3d_memset0(emarks, nverts * nverts * sizeof(int));

  // for intelligent graph walk
  // akin to what we saw in init poly quadratic
  int prev_curr_to_next[J3D_MAX_VERTS * J3D_MAX_VERTS];
  v = 0;
  n = 1;
  while (n < ndiredges) {
    while (n != indexbuffer[v + 1]) {
      prev_curr_to_next[diredgebuffer[n] * nverts + v] = diredgebuffer[n - 1];
      ++n;
    }
    prev_curr_to_next[diredgebuffer[indexbuffer[v]] * nverts + v] = diredgebuffer[indexbuffer[v + 1] - 1];
    ++v;
    ++n;
  }

  // for each node, check each of its neighbors
  // and use the "adjacency matrix" emarks to
  // check if we have walked that face yet.
  // If not, walk it
  v = 0;
  const int* const M_next = prev_curr_to_next;
  for (n = 0; n < ndiredges; ++n) {
    if (n == indexbuffer[v + 1]) {
      ++v;
    }
    // check if we encountered this face
    if (emarks[v * nverts + diredgebuffer[n]]) {
      continue;
    }

    // if this face hasn't been visited then fan it out
    vprev = v;
    vcurr = diredgebuffer[n];
    vnext = M_next[vprev * nverts + vcurr];
    j3d_memcpy(v0, vertbuffer[vprev], 3 * sizeof(double));
    j3d_memcpy(v1, vertbuffer[vcurr], 3 * sizeof(double));
    emarks[v * nverts + vcurr] = 1;

    // shift point we are "fanning"
    // from
    v0[0] = v0[0] - vc[0];
    v0[1] = v0[1] - vc[1];
    v0[2] = v0[2] - vc[2];

    v1[0] = v1[0] - vc[0];
    v1[1] = v1[1] - vc[1];
    v1[2] = v1[2] - vc[2];

    while (vnext != v) {
      j3d_memcpy(v2, vertbuffer[vnext], 3 * sizeof(double));

      // shift other two points
      // on triangle fan
      v2[0] = v2[0] - vc[0];
      v2[1] = v2[1] - vc[1];
      v2[2] = v2[2] - vc[2];

      const double v00 = v0[0];
      const double v01 = v0[1];
      const double v02 = v0[2];
      const double v10 = v1[0];
      const double v11 = v1[1];
      const double v12 = v1[2];
      const double v20 = v2[0];
      const double v21 = v2[1];
      const double v22 = v2[2];

      const double sixv = (-v20 * v11 * v02 + v10 * v21 * v02 + v20 * v01 * v12 - v00 * v21 * v12 -
                           v10 * v01 * v22 + v00 * v11 * v22);

      moments[0] += sixv;
      moments[1] += sixv * (v00 + v10 + v20);
      moments[2] += sixv * (v01 + v11 + v21);
      moments[3] += sixv * (v02 + v12 + v22);
      moments[4] += sixv * (v00 * v00 + v00 * v10 + v10 * v10 + v00 * v20 + v10 * v20 + v20 * v20);
      moments[5] += sixv * (2.0 * v00 * v01 + v01 * v10 + v00 * v11 + 2.0 * v10 * v11 + v01 * v20 +
                            v11 * v20 + v00 * v21 + v10 * v21 + 2.0 * v20 * v21);
      moments[6] += sixv * (2.0 * v00 * v02 + v02 * v10 + v00 * v12 + 2.0 * v10 * v12 + v02 * v20 +
                            v12 * v20 + v00 * v22 + v10 * v22 + 2.0 * v20 * v22);
      moments[7] += sixv * (v01 * v01 + v01 * v11 + v11 * v11 + v01 * v21 + v11 * v21 + v21 * v21);
      moments[8] += sixv * (2.0 * v01 * v02 + v02 * v11 + v01 * v12 + 2.0 * v11 * v12 + v02 * v21 +
                            v12 * v21 + v01 * v22 + v11 * v22 + 2.0 * v21 * v22);
      moments[9] += sixv * (v02 * v02 + v02 * v12 + v12 * v12 + v02 * v22 + v12 * v22 + v22 * v22);
      moments[10] += sixv * (v00 * v00 * v00 + v00 * v00 * v10 + v00 * v10 * v10 + v10 * v10 * v10 +
                             v00 * v00 * v20 + v00 * v10 * v20 + v10 * v10 * v20 + v00 * v20 * v20 +
                             v10 * v20 * v20 + v20 * v20 * v20);
      moments[11] +=
        sixv * (3.0 * v00 * v00 * v01 + 2.0 * v00 * v01 * v10 + v01 * v10 * v10 + v00 * v00 * v11 +
                2.0 * v00 * v10 * v11 + 3.0 * v10 * v10 * v11 + 2.0 * v00 * v01 * v20 +
                v01 * v10 * v20 + v00 * v11 * v20 + 2.0 * v10 * v11 * v20 + v01 * v20 * v20 +
                v11 * v20 * v20 + v00 * v00 * v21 + v00 * v10 * v21 + v10 * v10 * v21 +
                2.0 * v00 * v20 * v21 + 2.0 * v10 * v20 * v21 + 3.0 * v20 * v20 * v21);
      moments[12] +=
        sixv * (3.0 * v00 * v00 * v02 + 2.0 * v00 * v02 * v10 + v02 * v10 * v10 + v00 * v00 * v12 +
                2.0 * v00 * v10 * v12 + 3.0 * v10 * v10 * v12 + 2.0 * v00 * v02 * v20 +
                v02 * v10 * v20 + v00 * v12 * v20 + 2.0 * v10 * v12 * v20 + v02 * v20 * v20 +
                v12 * v20 * v20 + v00 * v00 * v22 + v00 * v10 * v22 + v10 * v10 * v22 +
                2.0 * v00 * v20 * v22 + 2.0 * v10 * v20 * v22 + 3.0 * v20 * v20 * v22);
      moments[13] +=
        sixv * (3.0 * v00 * v01 * v01 + v01 * v01 * v10 + 2.0 * v00 * v01 * v11 +
                2.0 * v01 * v10 * v11 + v00 * v11 * v11 + 3.0 * v10 * v11 * v11 + v01 * v01 * v20 +
                v01 * v11 * v20 + v11 * v11 * v20 + 2.0 * v00 * v01 * v21 + v01 * v10 * v21 +
                v00 * v11 * v21 + 2.0 * v10 * v11 * v21 + 2.0 * v01 * v20 * v21 +
                2.0 * v11 * v20 * v21 + v00 * v21 * v21 + v10 * v21 * v21 + 3.0 * v20 * v21 * v21);
      moments[14] +=
        sixv * (6.0 * v00 * v01 * v02 + 2.0 * v01 * v02 * v10 + 2.0 * v00 * v02 * v11 +
                2.0 * v02 * v10 * v11 + 2.0 * v00 * v01 * v12 + 2.0 * v01 * v10 * v12 +
                2.0 * v00 * v11 * v12 + 6.0 * v10 * v11 * v12 + 2.0 * v01 * v02 * v20 +
                v02 * v11 * v20 + v01 * v12 * v20 + 2.0 * v11 * v12 * v20 + 2.0 * v00 * v02 * v21 +
                v02 * v10 * v21 + v00 * v12 * v21 + 2.0 * v10 * v12 * v21 + 2.0 * v02 * v20 * v21 +
                2.0 * v12 * v20 * v21 + 2.0 * v00 * v01 * v22 + v01 * v10 * v22 + v00 * v11 * v22 +
                2.0 * v10 * v11 * v22 + 2.0 * v01 * v20 * v22 + 2.0 * v11 * v20 * v22 +
                2.0 * v00 * v21 * v22 + 2.0 * v10 * v21 * v22 + 6.0 * v20 * v21 * v22);
      moments[15] +=
        sixv * (3.0 * v00 * v02 * v02 + v02 * v02 * v10 + 2.0 * v00 * v02 * v12 +
                2.0 * v02 * v10 * v12 + v00 * v12 * v12 + 3.0 * v10 * v12 * v12 + v02 * v02 * v20 +
                v02 * v12 * v20 + v12 * v12 * v20 + 2.0 * v00 * v02 * v22 + v02 * v10 * v22 +
                v00 * v12 * v22 + 2.0 * v10 * v12 * v22 + 2.0 * v02 * v20 * v22 +
                2.0 * v12 * v20 * v22 + v00 * v22 * v22 + v10 * v22 * v22 + 3.0 * v20 * v22 * v22);
      moments[16] += sixv * (v01 * v01 * v01 + v01 * v01 * v11 + v01 * v11 * v11 + v11 * v11 * v11 +
                             v01 * v01 * v21 + v01 * v11 * v21 + v11 * v11 * v21 + v01 * v21 * v21 +
                             v11 * v21 * v21 + v21 * v21 * v21);
      moments[17] +=
        sixv * (3.0 * v01 * v01 * v02 + 2.0 * v01 * v02 * v11 + v02 * v11 * v11 + v01 * v01 * v12 +
                2.0 * v01 * v11 * v12 + 3.0 * v11 * v11 * v12 + 2.0 * v01 * v02 * v21 +
                v02 * v11 * v21 + v01 * v12 * v21 + 2.0 * v11 * v12 * v21 + v02 * v21 * v21 +
                v12 * v21 * v21 + v01 * v01 * v22 + v01 * v11 * v22 + v11 * v11 * v22 +
                2.0 * v01 * v21 * v22 + 2.0 * v11 * v21 * v22 + 3.0 * v21 * v21 * v22);
      moments[18] +=
        sixv * (3.0 * v01 * v02 * v02 + v02 * v02 * v11 + 2.0 * v01 * v02 * v12 +
                2.0 * v02 * v11 * v12 + v01 * v12 * v12 + 3.0 * v11 * v12 * v12 + v02 * v02 * v21 +
                v02 * v12 * v21 + v12 * v12 * v21 + 2.0 * v01 * v02 * v22 + v02 * v11 * v22 +
                v01 * v12 * v22 + 2.0 * v11 * v12 * v22 + 2.0 * v02 * v21 * v22 +
                2.0 * v12 * v21 * v22 + v01 * v22 * v22 + v11 * v22 * v22 + 3.0 * v21 * v22 * v22);
      moments[19] += sixv * (v02 * v02 * v02 + v02 * v02 * v12 + v02 * v12 * v12 + v12 * v12 * v12 +
                             v02 * v02 * v22 + v02 * v12 * v22 + v12 * v12 * v22 + v02 * v22 * v22 +
                             v12 * v22 * v22 + v22 * v22 * v22);

      vprev = vcurr;
      vcurr = vnext;
      vnext = M_next[vprev * nverts + vcurr];
      j3d_memcpy(v1, v2, 3 * sizeof(double));
      emarks[vprev * nverts + vcurr] = 1;
    }
    emarks[vcurr * nverts + vnext] = 1;
  }
  moments[0] /= 6.0;
  moments[1] /= 24.0;
  moments[2] /= 24.0;
  moments[3] /= 24.0;
  moments[4] /= 60.0;
  moments[5] /= 120.0;
  moments[6] /= 120.0;
  moments[7] /= 60.0;
  moments[8] /= 120.0;
  moments[9] /= 60.0;
  moments[10] /= 120.0;
  moments[11] /= 360.0;
  moments[12] /= 360.0;
  moments[13] /= 360.0;
  moments[14] /= 720.0;
  moments[15] /= 360.0;
  moments[16] /= 120.0;
  moments[17] /= 360.0;
  moments[18] /= 360.0;
  moments[19] /= 120.0;

  // shift moments since only zeroth order moment is invariant
  // to translation
  const double vc0 = vc[0];
  const double vc1 = vc[1];
  const double vc2 = vc[2];

  double moments2[20];

  moments2[1] = vc0 * moments[0] + moments[1];
  moments2[2] = vc1 * moments[0] + moments[2];
  moments2[3] = vc2 * moments[0] + moments[3];
  moments2[4] = vc0 * vc0 * moments[0] + 2.0 * vc0 * moments[1] + moments[4];
  moments2[5] = vc0 * vc1 * moments[0] + vc0 * moments[2] + vc1 * moments[1] + moments[5];
  moments2[6] = vc0 * vc2 * moments[0] + vc0 * moments[3] + vc2 * moments[1] + moments[6];
  moments2[7] = vc1 * vc1 * moments[0] + 2.0 * vc1 * moments[2] + moments[7];
  moments2[8] = vc1 * vc2 * moments[0] + vc1 * moments[3] + vc2 * moments[2] + moments[8];
  moments2[9] = vc2 * vc2 * moments[0] + 2.0 * vc2 * moments[3] + moments[9];
  moments2[10] = vc0 * vc0 * vc0 * moments[0] + 3.0 * vc0 * vc0 * moments[1] +
                 3.0 * vc0 * moments[4] + moments[10];
  moments2[11] = vc0 * vc0 * vc1 * moments[0] + vc0 * vc0 * moments[2] +
                 2.0 * vc0 * vc1 * moments[1] + 2.0 * vc0 * moments[5] + vc1 * moments[4] +
                 moments[11];
  moments2[12] = vc0 * vc0 * vc2 * moments[0] + vc0 * vc0 * moments[3] +
                 2.0 * vc0 * vc2 * moments[1] + 2.0 * vc0 * moments[6] + vc2 * moments[4] +
                 moments[12];
  moments2[13] = vc0 * vc1 * vc1 * moments[0] + 2.0 * vc0 * vc1 * moments[2] + vc0 * moments[7] +
                 vc1 * vc1 * moments[1] + 2.0 * vc1 * moments[5] + moments[13];
  moments2[14] = vc0 * vc1 * vc2 * moments[0] + vc0 * vc1 * moments[3] + vc0 * vc2 * moments[2] +
                 vc0 * moments[8] + vc1 * vc2 * moments[1] + vc1 * moments[6] + vc2 * moments[5] +
                 moments[14];
  moments2[15] = vc0 * vc2 * vc2 * moments[0] + 2.0 * vc0 * vc2 * moments[3] + vc0 * moments[9] +
                 vc2 * vc2 * moments[1] + 2.0 * vc2 * moments[6] + moments[15];
  moments2[16] = vc1 * vc1 * vc1 * moments[0] + 3.0 * vc1 * vc1 * moments[2] +
                 3.0 * vc1 * moments[7] + moments[16];
  moments2[17] = vc1 * vc1 * vc2 * moments[0] + vc1 * vc1 * moments[3] +
                 2.0 * vc1 * vc2 * moments[2] + 2.0 * vc1 * moments[8] + vc2 * moments[7] +
                 moments[17];
  moments2[18] = vc1 * vc2 * vc2 * moments[0] + 2.0 * vc1 * vc2 * moments[3] + vc1 * moments[9] +
                 vc2 * vc2 * moments[2] + 2.0 * vc2 * moments[8] + moments[18];
  moments2[19] = vc2 * vc2 * vc2 * moments[0] + 3.0 * vc2 * vc2 * moments[3] +
                 3.0 * vc2 * moments[9] + moments[19];


  j3d_memcpy(&moments[1], &moments2[1], 19 * sizeof(double));
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments_n(const j3d_poly* const poly, double* moments, int order) {
  // direct access to poly data
  const double(*const vertbuffer)[3] = poly->verts;
  const int* const indexbuffer = poly->indices;
  const int* const diredgebuffer = poly->dir_edges;
  const int nverts = poly->num_verts;
  const int ndiredges = poly->num_dir_edges;

  // some declarations
  double sixv;
  int vprev, vcurr, vnext, v, n, m, i, j, k, corder;
  double v0[3];
  double v1[3];
  double v2[3];
  double vc[3];

  // moments in the first order case should be size
  // j3d_num_moments(order)
  j3d_memset0(moments, j3d_num_moments(order) * sizeof(double));

  // if empty, exit
  if (nverts <= 0) {
    return;
  }

  // compute center of poly for shifting
  j3d_center(vc, vertbuffer, nverts);

  // for tracking visited edges
  int emarks[J3D_MAX_VERTS * J3D_MAX_VERTS];
  j3d_memset0(emarks, nverts * nverts * sizeof(int));

  // Storage for coefficients keep two layers of the pyramid of coefficients
  // Note: Uses twice as much space as needed, but indexing is faster this way
  int prevlayer = 0;
  int curlayer = 1;
  int P = order + 1;
  double S[2 * (J3D_MAX_ORDER + 1) * (J3D_MAX_ORDER + 1)];
  double D[2 * (J3D_MAX_ORDER + 1) * (J3D_MAX_ORDER + 1)];
  double C[2 * (J3D_MAX_ORDER + 1) * (J3D_MAX_ORDER + 1)];

  // for intelligent graph walk
  // akin to what we saw in init poly quadratic
  int prev_curr_to_next[J3D_MAX_VERTS * J3D_MAX_VERTS];
  v = 0;
  n = 1;
  while (n < ndiredges) {
    while (n != indexbuffer[v + 1]) {
      prev_curr_to_next[diredgebuffer[n] * nverts + v] = diredgebuffer[n - 1];
      ++n;
    }
    prev_curr_to_next[diredgebuffer[indexbuffer[v]] * nverts + v] = diredgebuffer[indexbuffer[v + 1] - 1];
    ++v;
    ++n;
  }

  // for each node, check each of its neighbors
  // and use the "adjacency matrix" emarks to
  // check if we have walked that face yet.
  // If not, walk it
  v = 0;
  const int* const M_next = prev_curr_to_next;
  for (n = 0; n < ndiredges; ++n) {
    if (n == indexbuffer[v + 1]) {
      ++v;
    }
    // check if we encountered this face
    if (emarks[v * nverts + diredgebuffer[n]]) {
      continue;
    }

    // if this face hasn't been visited then
    // initialize face looping and fan it out
    vprev = v;
    vcurr = diredgebuffer[n];
    vnext = M_next[vprev * nverts + vcurr];
    j3d_memcpy(v0, vertbuffer[vprev], 3 * sizeof(double));
    j3d_memcpy(v1, vertbuffer[vcurr], 3 * sizeof(double));
    emarks[v * nverts + vcurr] = 1;

    // shift point we are "fanning" from
    v0[0] = v0[0] - vc[0];
    v0[1] = v0[1] - vc[1];
    v0[2] = v0[2] - vc[2];

    v1[0] = v1[0] - vc[0];
    v1[1] = v1[1] - vc[1];
    v1[2] = v1[2] - vc[2];

    while (vnext != v) {
      j3d_memcpy(v2, vertbuffer[vnext], 3 * sizeof(double));

      // shift other two points
      // on triangle fan
      v2[0] = v2[0] - vc[0];
      v2[1] = v2[1] - vc[1];
      v2[2] = v2[2] - vc[2];

      sixv = (-v2[0] * v1[1] * v0[2] + v1[0] * v2[1] * v0[2] + v2[0] * v0[1] * v1[2] -
              v0[0] * v2[1] * v1[2] - v1[0] * v0[1] * v2[2] + v0[0] * v1[1] * v2[2]);

      // calculate the moments using the fast recursive method of Koehl (2012)
      // essentially building a set of trinomial pyramids, one layer at a time

      // base case
      S[prevlayer * P * P] = 1.0;
      D[prevlayer * P * P] = 1.0;
      C[prevlayer * P * P] = 1.0;
      moments[0] += ONE_SIXTH * sixv;

      // build up successive polynomial orders
      for (corder = 1, m = 1; corder <= order; ++corder) {
        double* const C_curlayer = &C[curlayer * P * P];
        double* const D_curlayer = &D[curlayer * P * P];
        double* const S_curlayer = &S[curlayer * P * P];
        double* const C_prevlayer = &C[prevlayer * P * P];
        double* const D_prevlayer = &D[prevlayer * P * P];
        double* const S_prevlayer = &S[prevlayer * P * P];
        for (i = corder; i >= 0; --i) {
          double* const C_curlayer_i = &C_curlayer[i * P];
          double* const D_curlayer_i = &D_curlayer[i * P];
          double* const S_curlayer_i = &S_curlayer[i * P];
          double* const C_prevlayer_i = &C_prevlayer[i * P];
          double* const D_prevlayer_i = &D_prevlayer[i * P];
          double* const S_prevlayer_i = &S_prevlayer[i * P];
          for (j = corder - i; j >= 0; --j, ++m) {
            k = corder - i - j;
            C_curlayer_i[j] = 0;
            D_curlayer_i[j] = 0;
            S_curlayer_i[j] = 0;
            if (i > 0) {
              C_curlayer_i[j] += v2[0] * C_prevlayer[((i - 1) * P) + j];
              D_curlayer_i[j] += v1[0] * D_prevlayer[((i - 1) * P) + j];
              S_curlayer_i[j] += v0[0] * S_prevlayer[((i - 1) * P) + j];
            }
            if (j > 0) {
              C_curlayer_i[j] += v2[1] * C_prevlayer_i[j - 1];
              D_curlayer_i[j] += v1[1] * D_prevlayer_i[j - 1];
              S_curlayer_i[j] += v0[1] * S_prevlayer_i[j - 1];
            }
            if (k > 0) {
              C_curlayer_i[j] += v2[2] * C_prevlayer_i[j];
              D_curlayer_i[j] += v1[2] * D_prevlayer_i[j];
              S_curlayer_i[j] += v0[2] * S_prevlayer_i[j];
            }
            D_curlayer_i[j] += C_curlayer_i[j];
            S_curlayer_i[j] += D_curlayer_i[j];
            moments[m] += sixv * S_curlayer_i[j];
          }
        }
        curlayer = 1 - curlayer;
        prevlayer = 1 - prevlayer;
      }

      vprev = vcurr;
      vcurr = vnext;
      vnext = M_next[vprev * nverts + vcurr];
      j3d_memcpy(v1, v2, 3 * sizeof(double));
      emarks[vprev * nverts + vcurr] = 1;
    }
    emarks[vcurr * nverts + vnext] = 1;
  }

  // reuse C to recursively compute the leading multinomial coefficients
  C[prevlayer * P * P] = 1.0;
  for (corder = 1, m = 1; corder <= order; ++corder) {
    double* const C_curlayer = &C[curlayer * P * P];
    double* const C_prevlayer = &C[prevlayer * P * P];
    for (i = corder; i >= 0; --i) {
      double* const C_curlayer_i = &C_curlayer[i * P];
      double* const C_prevlayer_i = &C_prevlayer[i * P];
      for (j = corder - i; j >= 0; --j, ++m) {
        k = corder - i - j;
        C_curlayer_i[j] = 0.0;
        if (i > 0) C_curlayer_i[j] += C_prevlayer[((i - 1) * P) + j];
        if (j > 0) C_curlayer_i[j] += C_prevlayer_i[j - 1];
        if (k > 0) C_curlayer_i[j] += C_prevlayer_i[j];
        moments[m] /= C_curlayer_i[j] * (corder + 1) * (corder + 2) * (corder + 3);
      }
    }
    curlayer = 1 - curlayer;
    prevlayer = 1 - prevlayer;
  }

  // shift moments back to initial position
  j3d_shift_moments(moments, order, vc);
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_get_moments(const j3d_poly* const poly, double* moments, int order) {
  if (order > J3D_MAX_ORDER) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j3d_get_moments: order is too high.  Increase J3D_MAX_ORDER.\n");
#endif
    return;
  }

  switch (order) {
  case 0:
    j3d_get_moments_0(poly, moments);
    break;
  case 1:
    j3d_get_moments_1(poly, moments);
    break;
  case 2:
    j3d_get_moments_2(poly, moments);
    break;
  case 3:
    j3d_get_moments_3(poly, moments);
  default:
    j3d_get_moments_n(poly, moments, order);
    break;
  }
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j3d_print_poly(const j3d_poly* const poly) {
  printf("number of vertices: %d\n", poly->num_verts);
  printf("number of directed edges: %d\n", poly->num_dir_edges);
  printf("vertex coordinates:\n");
  for (int v = 0; v < poly->num_verts; ++v) {
    printf("v%d: (%f, %f, %f)\n", v, poly->verts[v][0], poly->verts[v][1], poly->verts[v][2]);
  }
  printf("vertex neighbors:\n");
  for (int v = 0; v < poly->num_verts; ++v) {
    printf("v%d: ", v);
    const int* const nbrs = &(poly->dir_edges[poly->indices[v]]);
    const int degree_v = poly->indices[v + 1] - poly->indices[v];
    for (int n = 0; n < degree_v; ++n) {
      printf("%d ", nbrs[n]);
    }
    printf("\n");
  }

  printf("vertex neigbhor indices:\n");
  for (int v = 0; v < poly->num_verts + 1; ++v) {
    printf("v%d: %d\n", v, poly->indices[v]);
  }
}
