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

#include "j2d.h"

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_memcpy(void* const dest, const void* const src, const int n) {
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
void j2d_memset0(void* const ptr, const int n) {
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
double j2d_dot(const double* const va, const double* const vb) {
  return va[0] * vb[0] + va[1] * vb[1];
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_wav(double* const vr,
             const double* const va,
             const double wa,
             const double* const vb,
             const double wb) {
  const double rec = 1.0 / (wa + wb);
  vr[0] = (wa * va[0] + wb * vb[0]) * rec;
  vr[1] = (wa * va[1] + wb * vb[1]) * rec;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_norm(double* const vert) {
  const double rec = 1.0 / sqrt(j2d_dot(vert, vert));
  vert[0] *= rec;
  vert[1] *= rec;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j2d_num_moments(const int order) {
  return (order + 1) * (order + 2) / 2;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_shift_moments(double* const moments, int order, const double* const vc) {
  // var declarations
  int m, i, j, corder;
  int mm, mi, mj, mcorder;

  // store moments of a shifted polygon
  double moments2[J2D_MAX_NUM_MOMENTS];
  for (m = 0; m < j2d_num_moments(order); ++m) {
    moments2[m] = 0.0;
  }

  // calculate and save Pascal's triangle
  double B[J2D_MAX_ORDER + 1][J2D_MAX_ORDER + 1];
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
  // \int_\Omega x^i y^j d\vec r =
  // \int_\omega (x+\xi)^i (y+\eta)^j d\vec r =
  // \sum_{a,b,c=0}^{i,j,k} \binom{i}{a} \binom{j}{b}
  // \xi^{i-a} \eta^{j-b} \int_\omega x^a y^b d\vec r
  for (corder = 1, m = 1; corder <= order; ++corder) {
    for (i = corder; i >= 0; --i, ++m) {
      j = corder - i;
      for (mcorder = 0, mm = 0; mcorder <= corder; ++mcorder) {
        for (mi = mcorder; mi >= 0; --mi, ++mm) {
          mj = mcorder - mi;
          if (mi <= i && mj <= j) {
            moments2[m] +=
              B[mi][i] * B[mj][j] * pow(vc[0], (i - mi)) * pow(vc[1], (j - mj)) * moments[mm];
          }
        }
      }
    }
  }

  // assign shifted moments
  for (m = 1; m < j2d_num_moments(order); ++m) {
    moments[m] = moments2[m];
  }
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_center(double* const vc, const double (*const verts)[2], const int numverts) {
  j2d_memset0(vc, 2 * sizeof(double));
  for (int i = 0; i < numverts; ++i) {
    vc[0] += verts[i][0];
    vc[1] += verts[i][1];
  }
  const double rec = 1.0 / numverts;
  vc[0] *= rec;
  vc[1] *= rec;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
double j2d_orient(const double* const va, const double* const vb, const double* const vc) {
  return 0.5 * ((va[0] - vc[0]) * (vb[1] - vc[1]) - (vb[0] - vc[0]) * (va[1] - vc[1]));
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_poly_faces_from_verts(j2d_plane* faces, const j2d_vec2* const vertices, const int numverts) {
  // dummy vars
  int f;
  double v0[2];
  double v1[2];

  // calculate a centroid and a unit normal for each face
  for (f = 0; f < numverts; ++f) {
    j2d_memcpy(v0, vertices[f].xy, 2 * sizeof(double));
    j2d_memcpy(v1, vertices[(f + 1) % numverts].xy, 2 * sizeof(double));

    // normal of the edge
    faces[f].n.xy[0] = v0[1] - v1[1];
    faces[f].n.xy[1] = v1[0] - v0[0];

    // normalize the normals and set the signed distance to origin
    j2d_norm(faces[f].n.xy);
    faces[f].d = -j2d_dot(faces[f].n.xy, v0);
  }
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j2d_init_poly(j2d_poly* poly, const j2d_vec2* const vertices, const int numverts) {
  if (numverts > J2D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j2d_init_poly: Max vertex buffer size exceeded.  Increase J2D_MAX_VERTS.\n");
#endif
    return 0;
  }

  // direct access to vertex buffer
  double(*const vertbuffer)[2] = poly->verts;
  int* const diredgebuffer = poly->dir_edges;
  int* nverts = &poly->num_verts;

  // init the poly
  *nverts = numverts;
  int v;
  for (v = 0; v < *nverts; ++v) {
    j2d_memcpy(vertbuffer[v], vertices[v].xy, 2 * sizeof(double));
    diredgebuffer[2 * v] = (v + 1) % (*nverts);
    diredgebuffer[(2 * v) + 1] = (*nverts + v - 1) % (*nverts);
  }

  return 1;
}

int j2d_init_box(j2d_poly* poly, const j2d_vec2* const rbounds) {
  if (4 > J2D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j2d_init_box: Max vertex buffer size exceeded.  Increase J2D_MAX_VERTS.\n");
#endif
    return 0;
  }

  poly->num_verts = 4;

  // direct access to vertex buffer
  double(*const vertbuffer)[2] = poly->verts;
  int* const diredgebuffer = poly->dir_edges;

  diredgebuffer[0] = 1;
  diredgebuffer[1] = 3;
  diredgebuffer[2] = 2;
  diredgebuffer[3] = 0;
  diredgebuffer[4] = 3;
  diredgebuffer[5] = 1;
  diredgebuffer[6] = 0;
  diredgebuffer[7] = 2;

  vertbuffer[0][0] = rbounds[0].xy[0];
  vertbuffer[0][1] = rbounds[0].xy[1];
  vertbuffer[1][0] = rbounds[1].xy[0];
  vertbuffer[1][1] = rbounds[0].xy[1];
  vertbuffer[2][0] = rbounds[1].xy[0];
  vertbuffer[2][1] = rbounds[1].xy[1];
  vertbuffer[3][0] = rbounds[0].xy[0];
  vertbuffer[3][1] = rbounds[1].xy[1];

  return 0;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_copy_poly(j2d_poly* destination, const j2d_poly* const source) {
  destination->num_verts = source->num_verts;

  j2d_memcpy(destination->verts, source->verts, 2 * source->num_verts * sizeof(double));
  j2d_memcpy(destination->dir_edges, source->dir_edges, 2 * source->num_verts * sizeof(int));
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j2d_clip(j2d_poly* poly, const j2d_plane* planes, const int nplanes) {
  // variable declarations
  int v, p, np, onv, vstart, vcur, vnext, numunclipped;

  // direct access to vertex buffer
  double(*const vertbuffer)[2] = poly->verts;
  int* const diredgebuffer = poly->dir_edges;
  int* const nverts = &poly->num_verts;
  if (*nverts <= 0) {
    return 0;
  }

  // signed distances to the clipping plane
  double sdists[J2D_MAX_VERTS];
  double smin, smax;

  // for marking clipped vertices
  int clipped[J2D_MAX_VERTS];

  // loop over each clip plane
  for (p = 0; p < nplanes; ++p) {
    // calculate signed distances to the clip plane
    onv = *nverts;
    smin = 1.0e30;
    smax = -1.0e30;
    j2d_memset0(clipped, J2D_MAX_VERTS * sizeof(int));
    for (v = 0; v < onv; ++v) {
      sdists[v] = planes[p].d + j2d_dot(vertbuffer[v], planes[p].n.xy);
      if (sdists[v] < smin) smin = sdists[v];
      if (sdists[v] > smax) smax = sdists[v];
      if (sdists[v] < 0.0) clipped[v] = 1;
    }

    // skip this face if the poly lies entirely on one side of it
    if (smin >= 0.0) continue;
    if (smax <= 0.0) {
      *nverts = 0;
      return 1;
    }

    // check all edges and insert new vertices on the bisected edges
    for (vcur = 0; vcur < onv; ++vcur) {
      if (clipped[vcur]) continue;
      for (np = 0; np < 2; ++np) {
        vnext = diredgebuffer[(2 * vcur) + np];
        if (!clipped[vnext]) continue;
        if (*nverts == J2D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
          fprintf(stderr, "j2d_clip: Max vertex buffer exceeded.  Increase J2D_MAX_VERTS.\n");
#endif
          return 0;
        }
        diredgebuffer[(2 * (*nverts)) + 1 - np] = vcur;
        diredgebuffer[(2 * (*nverts)) + np] = -1;
        diredgebuffer[(2 * vcur) + np] = *nverts;
        j2d_wav(
          vertbuffer[*nverts], vertbuffer[vcur], -sdists[vnext], vertbuffer[vnext], sdists[vcur]);
        ++(*nverts);
      }
    }

    // for each new vert, search around the poly for its new neighbors
    // and doubly-link everything
    for (vstart = onv; vstart < *nverts; ++vstart) {
      if (diredgebuffer[(2 * vstart) + 1] >= 0) continue;
      vcur = diredgebuffer[2 * vstart];
      do {
        vcur = diredgebuffer[2 * vcur];
      } while (vcur < onv);
      diredgebuffer[(2 * vstart) + 1] = vcur;
      diredgebuffer[2 * vcur] = vstart;
    }

    // go through and compress the vertex list, removing clipped verts
    // and re-indexing accordingly (reusing `clipped` to re-index everything)
    numunclipped = 0;
    for (v = 0; v < *nverts; ++v) {
      if (!clipped[v]) {
        // copy vertices
        j2d_memcpy(vertbuffer[numunclipped], vertbuffer[v], 2 * sizeof(double));

        // copy neighbors
        j2d_memcpy(&diredgebuffer[2 * numunclipped], &diredgebuffer[2 * v], 2 * sizeof(int));

        clipped[v] = numunclipped++;
      }
    }
    *nverts = numunclipped;
    for (v = 0; v < 2 * (*nverts); ++v) {
      diredgebuffer[v] = clipped[diredgebuffer[v]];
    }
  }
  return 1;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
int j2d_split(j2d_poly* inpolys,
              const int npolys,
              const j2d_plane plane,
              j2d_poly* out_pos,
              j2d_poly* out_neg) {
  // direct access to vertex buffer
  int* nverts;
  int p;
  double(*vertbuffer)[2];
  int* diredgebuffer;
  int v, np, onv, vcur, vnext, vstart, nright, cside;
  double newpos[2];
  int side[J2D_MAX_VERTS];
  double sdists[J2D_MAX_VERTS];
  j2d_poly* outpolys[2];

  for (p = 0; p < npolys; ++p) {
    nverts = &inpolys[p].num_verts;
    vertbuffer = inpolys[p].verts;
    diredgebuffer = inpolys[p].dir_edges;
    outpolys[0] = &out_pos[p];
    outpolys[1] = &out_neg[p];
    if (*nverts <= 0) {
      out_pos[p].num_verts = 0;
      out_neg[p].num_verts = 0;
      continue;
    }


    // calculate signed distances to the clip plane
    nright = 0;
    j2d_memset0(side, J2D_MAX_VERTS * sizeof(int));
    for (v = 0; v < *nverts; ++v) {
      sdists[v] = plane.d + j2d_dot(vertbuffer[v], plane.n.xy);
      sdists[v] *= -1;
      if (sdists[v] < 0.0) {
        side[v] = 1;
        nright++;
      }
    }

    // return if the poly lies entirely on one side of it
    if (nright == 0) {
      *(outpolys[0]) = inpolys[p];
      outpolys[1]->num_verts = 0;
      continue;
    }
    if (nright == *nverts) {
      *(outpolys[1]) = inpolys[p];
      outpolys[0]->num_verts = 0;
      continue;
    }

    // check all edges and insert new vertices on the bisected edges
    onv = *nverts;
    for (vcur = 0; vcur < onv; ++vcur) {
      if (side[vcur]) continue;
      for (np = 0; np < 2; ++np) {
        vnext = diredgebuffer[(2 * vcur) + np];
        if (!side[vnext]) continue;
        j2d_wav(newpos, vertbuffer[vcur], -sdists[vnext], vertbuffer[vnext], sdists[vcur]);
        if (*nverts == J2D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
          fprintf(stderr, "j2d_split: Max vertex buffer size exceeded.  Increase J2D_MAX_VERTS.\n");
#endif
          return 0;
        }
        j2d_memcpy(vertbuffer[*nverts], newpos, 2 * sizeof(double));
        diredgebuffer[(2 * vcur) + np] = *nverts;
        diredgebuffer[(2 * (*nverts)) + np] = -1;
        diredgebuffer[(2 * (*nverts)) + 1 - np] = vcur;
        ++(*nverts);
        if (*nverts == J2D_MAX_VERTS) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
          fprintf(stderr, "j2d_split: Max vertex buffer size exceeded.  Increase J2D_MAX_VERTS.\n");
#endif
          return 0;
        }
        side[*nverts] = 1;
        j2d_memcpy(vertbuffer[*nverts], newpos, 2 * sizeof(double));
        diredgebuffer[(2 * (*nverts)) + 1 - np] = -1;
        diredgebuffer[(2 * (*nverts)) + np] = vnext;
        diredgebuffer[(2 * vnext) + 1 - np] = *nverts;
        ++(*nverts);
      }
    }

    // for each new vert, search around the poly for its new neighbors
    // and doubly-link everything
    for (vstart = onv; vstart < *nverts; ++vstart) {
      if (diredgebuffer[(2 * vstart) + 1] >= 0) continue;
      vcur = diredgebuffer[2 * vstart];
      do {
        vcur = diredgebuffer[2 * vcur];
      } while (vcur < onv);
      diredgebuffer[(2 * vstart) + 1] = vcur;
      diredgebuffer[2 * vcur] = vstart;
    }

    // copy and compress vertices into their new buffers
    // reusing side[] for reindexing
    onv = *nverts;
    outpolys[0]->num_verts = 0;
    outpolys[1]->num_verts = 0;
    for (v = 0; v < onv; ++v) {
      cside = side[v];
      j2d_memcpy(outpolys[cside]->verts[outpolys[cside]->num_verts], vertbuffer[v], 2 * sizeof(double));
      j2d_memcpy(&(outpolys[cside]->dir_edges[2 * outpolys[cside]->num_verts]),
             &diredgebuffer[2 * v],
             2 * sizeof(int));
      side[v] = (outpolys[cside]->num_verts)++;
    }

    for (v = 0; v < 2 * outpolys[0]->num_verts; ++v) {
      outpolys[0]->dir_edges[v] = side[outpolys[0]->dir_edges[v]];
    }
    for (v = 0; v < 2 * outpolys[1]->num_verts; ++v) {
      outpolys[1]->dir_edges[v] = side[outpolys[1]->dir_edges[v]];
    }
  }
  return 1;
}

#ifdef J3D_ENABLE_Kokkos
KOKKOS_FUNCTION
#endif
void j2d_get_moments(const j2d_poly* const poly, double* moments, int order) {
  if (order > J2D_MAX_ORDER) {
#if !defined(NDEBUG) && !defined(J3D_ENABLE_Kokkos)
    fprintf(stderr, "j2d_get_moments: order is too high.  Increase J2D_MAX_ORDER.\n");
#endif
    return;
  }

  // var declarations
  int vcur, vnext, m, i, j, corder;
  double twoa;
  double v0[2];
  double v1[2];
  double vc[2];

  // direct access to vertex buffer
  const double(*const vertbuffer)[2] = poly->verts;
  const int* const diredgebuffer = poly->dir_edges;
  const int* const nverts = &poly->num_verts;

  // zero the moments
  for (m = 0; m < j2d_num_moments(order); ++m) {
    moments[m] = 0.0;
  }

  if (*nverts <= 0) return;


  j2d_center(vc, vertbuffer, *nverts);

  // Storage for coefficients
  // keep two layers of the triangle of coefficients
  int prevlayer = 0;
  int curlayer = 1;
  double D[J2D_MAX_ORDER + 1][2];
  double C[J2D_MAX_ORDER + 1][2];

  // iterate over edges and compute a sum over simplices
  for (vcur = 0; vcur < *nverts; ++vcur) {
    vnext = diredgebuffer[2 * vcur];
    j2d_memcpy(v0, vertbuffer[vcur], 2 * sizeof(double));
    j2d_memcpy(v1, vertbuffer[vnext], 2 * sizeof(double));

    v0[0] = v0[0] - vc[0];
    v0[1] = v0[1] - vc[1];
    v1[0] = v1[0] - vc[0];
    v1[1] = v1[1] - vc[1];

    twoa = (v0[0] * v1[1] - v0[1] * v1[0]);

    // calculate the moments
    // using the fast recursive method of Koehl (2012)
    // essentially building a set of Pascal's triangles, one layer at a time

    // base case
    D[0][prevlayer] = 1.0;
    C[0][prevlayer] = 1.0;
    moments[0] += 0.5 * twoa;

    // build up successive polynomial orders
    for (corder = 1, m = 1; corder <= order; ++corder) {
      for (i = corder; i >= 0; --i, ++m) {
        j = corder - i;
        C[i][curlayer] = 0;
        D[i][curlayer] = 0;
        if (i > 0) {
          C[i][curlayer] += v1[0] * C[i - 1][prevlayer];
          D[i][curlayer] += v0[0] * D[i - 1][prevlayer];
        }
        if (j > 0) {
          C[i][curlayer] += v1[1] * C[i][prevlayer];
          D[i][curlayer] += v0[1] * D[i][prevlayer];
        }
        D[i][curlayer] += C[i][curlayer];
        moments[m] += twoa * D[i][curlayer];
      }
      curlayer = 1 - curlayer;
      prevlayer = 1 - prevlayer;
    }
  }

  // reuse C to recursively compute the leading multinomial coefficients
  C[0][prevlayer] = 1.0;
  for (corder = 1, m = 1; corder <= order; ++corder) {
    for (i = corder; i >= 0; --i, ++m) {
      j = corder - i;
      C[i][curlayer] = 0.0;
      if (i > 0) C[i][curlayer] += C[i - 1][prevlayer];
      if (j > 0) C[i][curlayer] += C[i][prevlayer];
      moments[m] /= C[i][curlayer] * (corder + 1) * (corder + 2);
    }
    curlayer = 1 - curlayer;
    prevlayer = 1 - prevlayer;
  }

  j2d_shift_moments(moments, order, vc);
}
