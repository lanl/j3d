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

#undef POLY_ORDER
#undef NUM_MOMENTS
#undef NTETS
#undef STACK_SIZE
#undef NTHETA
#undef NPHI
#undef NVERTS
#undef NFACES

#define POLY_ORDER 0
#define NUM_MOMENTS ((POLY_ORDER + 1) * (POLY_ORDER + 2) * (POLY_ORDER + 3) / 6)
#define NTETS 128
#define STACK_SIZE 64
#define NTHETA 3
#define NPHI 3
#define NVERTS (NTHETA * NPHI)
#define NFACES (2 * NTHETA * NPHI)

const double TOL = 1e-4;

j3d_plane (*choptions_3d[6])(j3d_poly* poly);

int split_tets_thru_centroid() {
  printf("=== BEGIN SPLIT TET THRU CENTROID TEST ===\n");

  const double MIN_VOL = 1.0e-8;

  // a very basic sanity check. Splits a tet through its centroid
  // and checks to see whether the two resulting volumes add up to equal the original

  // variables: the polyhedra and their moments
  int m, i;
  j3d_vec3 verts[4];
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];
  j3d_plane splane;

  // generate a random tet and clip plane
  j3d_poly opoly[NTETS], poly1[NTETS], poly2[NTETS];
  for (i = 0; i < NTETS; ++i) {
    rand_tet_3d(verts, MIN_VOL);
    j3d_init_tet(&opoly[i], verts);
    if (i == 13) {
      splane = thru_cent_3d(&opoly[i]);
    }
  }

  // split them all about the same plane
  j3d_split(opoly, NTETS, splane, poly1, poly2);

  for (i = 0; i < NTETS; ++i) {
    // reduce the original and its two parts
    j3d_get_moments(&opoly[i], om, POLY_ORDER);
    j3d_get_moments(&poly1[i], m1, POLY_ORDER);
    j3d_get_moments(&poly2[i], m2, POLY_ORDER);

    // make sure the sum of moments equals the original
    for (m = 0; m < NUM_MOMENTS; ++m) {
      ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
    }

    // make sure neither of the two resulting volumes is larger than the original
    // (within some tolerance)
    ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
    ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));
  }
  printf("=== END SPLIT TET THRU CENTROID TEST ===\n");
  return EXIT_SUCCESS;
}

int split_nonconvex() {
  printf("=== BEGIN SPLIT NON CONVEX TEST ===\n");
  // a very basic sanity check. Splits a nonconvex poly (a zig-zag prism of sorts)
  // and checks to see whether the two resulting volumes add up to equal the original

  const int NZIGS = 3;
  const double ZOFF = 0.1;

  // variables: the polyhedra and their moments
  int m, v;
  j3d_plane splane;
  j3d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];

  // explicitly create the nonconvex poly
  opoly.num_verts = 4 * NZIGS;
  for (v = 0; v < NZIGS; ++v) {
    opoly.verts[v][0] = 1.0 * v;
    opoly.verts[v][1] = ZOFF + 1.0 * (v % 2);
    opoly.verts[v][2] = 0.0;
    opoly.verts[v + NZIGS][0] = 1.0 * (NZIGS - v - 1);
    opoly.verts[v + NZIGS][1] = 1.0 * ((NZIGS - v - 1) % 2);
    opoly.verts[v + NZIGS][2] = 0.0;
    opoly.verts[v + 2 * NZIGS][0] = 1.0 * v;
    opoly.verts[v + 2 * NZIGS][1] = ZOFF + 1.0 * (v % 2);
    opoly.verts[v + 2 * NZIGS][2] = 1.0;
    opoly.verts[v + 3 * NZIGS][0] = 1.0 * (NZIGS - v - 1);
    opoly.verts[v + 3 * NZIGS][1] = 1.0 * ((NZIGS - v - 1) % 2);
    opoly.verts[v + 3 * NZIGS][2] = 1.0;
  }

  opoly.indices[0] = 0;
  for (v = 0; v < 4 * NZIGS; ++v) {
    opoly.indices[v + 1] = opoly.indices[v] + 3;
  }
  opoly.num_dir_edges = 3 * 4 * NZIGS;

  for (v = 0; v < 2 * NZIGS; ++v) {
    opoly.dir_edges[opoly.indices[v] + 0] = (v + 1) % (2 * NZIGS);
    opoly.dir_edges[opoly.indices[v] + 1] = (v + 2 * NZIGS - 1) % (2 * NZIGS);
    opoly.dir_edges[opoly.indices[v] + 2] = (v + 2 * NZIGS);
    opoly.dir_edges[opoly.indices[v + 2 * NZIGS] + 0] =
      (v + 2 * NZIGS - 1) % (2 * NZIGS) + 2 * NZIGS;
    opoly.dir_edges[opoly.indices[v + 2 * NZIGS] + 1] = (v + 1) % (2 * NZIGS) + 2 * NZIGS;
    opoly.dir_edges[opoly.indices[v + 2 * NZIGS] + 2] = v;
  }

  // split along the x-axis (two single connected components this direction)
  poly1 = opoly;
  poly2 = opoly;
  splane.n.xyz[0] = 1.0;
  splane.n.xyz[1] = 0.0;
  splane.n.xyz[2] = 0.0;
  splane.d = -0.5 * (NZIGS - 1);
  j3d_clip(&poly1, &splane, 1);
  splane.n.xyz[0] *= -1;
  splane.n.xyz[1] *= -1;
  splane.n.xyz[2] *= -1;
  splane.d *= -1;
  j3d_clip(&poly2, &splane, 1);
  j3d_get_moments(&opoly, om, POLY_ORDER);
  j3d_get_moments(&poly1, m1, POLY_ORDER);
  j3d_get_moments(&poly2, m2, POLY_ORDER);
  for (m = 0; m < NUM_MOMENTS; ++m) {
    ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
  }
  ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
  ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));

  // split along the z-axis (two single connected components this direction)
  poly1 = opoly;
  poly2 = opoly;
  splane.n.xyz[0] = 0.0;
  splane.n.xyz[1] = 0.0;
  splane.n.xyz[2] = 1.0;
  splane.d = -0.5;
  j3d_clip(&poly1, &splane, 1);
  splane.n.xyz[0] *= -1;
  splane.n.xyz[1] *= -1;
  splane.n.xyz[2] *= -1;
  splane.d *= -1;
  j3d_clip(&poly2, &splane, 1);
  j3d_get_moments(&opoly, om, POLY_ORDER);
  j3d_get_moments(&poly1, m1, POLY_ORDER);
  j3d_get_moments(&poly2, m2, POLY_ORDER);
  for (m = 0; m < NUM_MOMENTS; ++m) {
    ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
  }
  ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
  ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));

  // split along the y-axis (multiple connected components this direction)
  poly1 = opoly;
  poly2 = opoly;
  splane.n.xyz[0] = 0.0;
  splane.n.xyz[1] = 1.0;
  splane.n.xyz[2] = 0.0;
  splane.d = -0.5 * (1.0 + ZOFF);
  j3d_clip(&poly1, &splane, 1);
  splane.n.xyz[0] *= -1;
  splane.n.xyz[1] *= -1;
  splane.n.xyz[2] *= -1;
  splane.d *= -1;
  j3d_clip(&poly2, &splane, 1);
  j3d_get_moments(&opoly, om, POLY_ORDER);
  j3d_get_moments(&poly1, m1, POLY_ORDER);
  j3d_get_moments(&poly2, m2, POLY_ORDER);
  for (m = 0; m < NUM_MOMENTS; ++m) {
    ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
  }
  ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
  ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));

  printf("=== END SPLIT NON CONVEX TEST ===\n");
  return EXIT_SUCCESS;
}

int recursive_splitting_nondegenerate() {
  printf("=== BEGIN RECURSIVE SPLITTING NONDEGENERATE TEST ===\n");
  // recursively splits a polyhedron (starting from a tet) through the centroid,
  // checking to see that the resulting volumes add up properly.

  const int MAX_DEPTH = 16;
  const int NUM_TRIALS = 10;
  const double MIN_VOL = 1.0e-8;

  // explicit stack-based implementation
  int nstack, depth, t, m;
  j3d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j3d_vec3 verts[4];
  j3d_plane splane;
  j3d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];

  // do many trials
  //printf("Recursively splitting %d tetrahedra, maximum splits per tet is %d.\n", NUM_TRIALS, MAX_DEPTH);
  for (t = 0; t < NUM_TRIALS; ++t) {
    // generate a random tet
    rand_tet_3d(verts, MIN_VOL);
    j3d_init_tet(&opoly, verts);

    // push the starting tet to the stack
    nstack = 0;
    polystack[nstack] = opoly;
    depthstack[nstack] = 0;
    ++nstack;

    // recursively split the poly
    while (nstack > 0) {
      // pop the stack
      --nstack;
      opoly = polystack[nstack];
      depth = depthstack[nstack];

      // generate a randomly oriented plane
      // through the centroid of the poly
      splane = thru_cent_3d(&opoly);

      // split the poly by making two copies of the original poly
      // and them clipping them against the same plane, with one
      // oriented oppositely
      poly1 = opoly;
      poly2 = opoly;
      j3d_clip(&poly1, &splane, 1);
      splane.n.xyz[0] *= -1;
      splane.n.xyz[1] *= -1;
      splane.n.xyz[2] *= -1;
      splane.d *= -1;
      j3d_clip(&poly2, &splane, 1);

      // reduce the original and its two parts
      j3d_get_moments(&opoly, om, POLY_ORDER);
      j3d_get_moments(&poly1, m1, POLY_ORDER);
      j3d_get_moments(&poly2, m2, POLY_ORDER);

      // make sure the sum of moments equals the original
      for (m = 0; m < NUM_MOMENTS; ++m) {
        ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
      }

      // make sure neither of the two resulting volumes is larger than the original
      // (within some tolerance)
      ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
      ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));

      //printf("nstack = %d, depth = %d, opoly = %.10e, p1 = %.10e, p2 = %.10e, err = %.10e\n",
      //nstack, depth, om[0], m1[0], m2[0], fabs(1.0 - om[0]/(m1[0] + m2[0])));

      // push the children to the stack if they have
      // an acceptably large volume
      if (depth < MAX_DEPTH) {
        if (m1[0] > MIN_VOL) {
          polystack[nstack] = poly1;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
        if (m2[0] > MIN_VOL) {
          polystack[nstack] = poly2;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
      }
    }
  }
  printf("=== END RECURSIVE SPLITTING NONDEGENERATE TEST ===\n");
  return EXIT_SUCCESS;
}

int recursive_splitting_degenerate() {
  printf("=== BEGIN RECURSIVE SPLITTING DEGENERATE TEST ===\n");
  // recursively splits a polyhedron (starting from a tet) with a cut plane
  // that is always degenerate with the tet in some way,
  // checking to see that the resulting volumes add up properly.

  const double MIN_VOL = 1.0e-8;
  const int MAX_DEPTH = 16;
  const int NUM_TRIALS = 10;

  // explicit stack-based implementation
  int nstack, depth, t, chopt, m;
  j3d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j3d_vec3 verts[4];
  j3d_plane splane;
  j3d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];

  // do many trials
  // printf("Recursively splitting %d tetrahedra, maximum splits per tet is %d.\n", NUM_TRIALS, MAX_DEPTH);
  for (t = 0; t < NUM_TRIALS; ++t) {
    // generate a random tet
    rand_tet_3d(verts, MIN_VOL);
    j3d_init_tet(&opoly, verts);
    // push the starting tet to the stack
    nstack = 0;
    polystack[nstack] = opoly;
    depthstack[nstack] = 0;
    ++nstack;

    // recursively split the poly
    while (nstack > 0) {
      // pop the stack
      --nstack;
      opoly = polystack[nstack];
      depth = depthstack[nstack];

      // generate a random plane from one of a few
      // possible degenerate configurations, ensuring that it
      // has a valid unit normal
      chopt = rand_int(6);
      do {
        splane = choptions_3d[chopt](&opoly);
      } while (splane.n.xyz[0] == 0.0 && splane.n.xyz[1] == 0.0 && splane.n.xyz[2] == 0.0);
      // split the poly by making two copies of the original poly
      // and them clipping them against the same plane, with one
      // oriented oppositely
      poly1 = opoly;
      poly2 = opoly;
      j3d_clip(&poly1, &splane, 1);
      splane.n.xyz[0] *= -1;
      splane.n.xyz[1] *= -1;
      splane.n.xyz[2] *= -1;
      splane.d *= -1;
      j3d_clip(&poly2, &splane, 1);
      // reduce the original and its two parts
      j3d_get_moments(&opoly, om, POLY_ORDER);
      j3d_get_moments(&poly1, m1, POLY_ORDER);
      j3d_get_moments(&poly2, m2, POLY_ORDER);
      // make sure the sum of moments equals the original
      for (m = 0; m < NUM_MOMENTS; ++m) {
        ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
      }

      // make sure neither of the two resulting volumes is larger than the original
      // (within some tolerance)
      ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
      ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));

      //printf("nstack = %d, depth = %d, opoly = %.10e, p1 = %.10e, p2 = %.10e, err = %.10e\n",
      //nstack, depth, om[0], m1[0], m2[0], fabs(1.0 - om[0]/(m1[0] + m2[0])));

      // push the children to the stack if they have
      // an acceptably large volume
      if (depth < MAX_DEPTH) {
        if (m1[0] > MIN_VOL) {
          polystack[nstack] = poly1;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
        if (m2[0] > MIN_VOL) {
          polystack[nstack] = poly2;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
      }
    }
  }
  printf("=== END RECURSIVE SPLITTING DEGENERATE TEST ===\n");
  return EXIT_SUCCESS;
}

int recursive_splitting_degenerate_perturbed() {
  printf("=== BEGIN RECURSIVE SPLITTING DEGENERATE PERTURBED ===\n");
  // recursively splits a polyhedron (starting from a tet) with a cut plane
  // that is always degenerate with the tet in some way,
  // checking to see that the resulting volumes add up properly.
  // In this one, the cut plane is perturbed by a small amount.

  const int MIN_PERTURB_ORDER = -17;
  const int MAX_PERTURB_ORDER = -1;
  const double MIN_VOL = 1.0e-8;
  const int MAX_DEPTH = 16;
  const int NUM_TRIALS = 1;

  // explicit stack-based implementation
  int nstack, depth, t, chopt, m;
  j3d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j3d_vec3 verts[4];
  j3d_plane splane;
  j3d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];
  double perturb;

  // do many trials
  // printf("Recursively splitting %d tetrahedra, maximum splits per tet is %d.\n", NUM_TRIALS, MAX_DEPTH);
  for (t = 0; t < NUM_TRIALS; ++t) {
    // compute the order of magnitude by which to perturb the clip plane,
    // determined by the trial number
    perturb = pow(10, MIN_PERTURB_ORDER + t % (MAX_PERTURB_ORDER - MIN_PERTURB_ORDER));
    //printf("omag = %d, pow = %.10e\n", MIN_PERTURB_ORDER + t%(MAX_PERTURB_ORDER - MIN_PERTURB_ORDER), perturb);

    // generate a random tet
    rand_tet_3d(verts, MIN_VOL);
    j3d_init_tet(&opoly, verts);

    // push the starting tet to the stack
    nstack = 0;
    polystack[nstack] = opoly;
    depthstack[nstack] = 0;
    ++nstack;

    // recursively split the poly
    while (nstack > 0) {
      // pop the stack
      --nstack;
      opoly = polystack[nstack];
      depth = depthstack[nstack];

      // generate a random plane from one of a few
      // possible degenerate configurations, ensuring that it
      // has a valid unit normal
      chopt = rand_int(6);
      do {
        splane = choptions_3d[chopt](&opoly);
      } while (splane.n.xyz[0] == 0.0 && splane.n.xyz[1] == 0.0 && splane.n.xyz[2] == 0.0);

      // randomly perturb the plane
      splane.n.xyz[0] *= 1.0 + perturb * (rand_uniform() - 0.5);
      splane.n.xyz[1] *= 1.0 + perturb * (rand_uniform() - 0.5);
      splane.n.xyz[2] *= 1.0 + perturb * (rand_uniform() - 0.5);
      splane.d *= 1.0 + perturb * (rand_uniform() - 0.5);

      // split the poly by making two copies of the original poly
      // and them clipping them against the same plane, with one
      // oriented oppositely
      poly1 = opoly;
      poly2 = opoly;
      j3d_clip(&poly1, &splane, 1);
      splane.n.xyz[0] *= -1;
      splane.n.xyz[1] *= -1;
      splane.n.xyz[2] *= -1;
      splane.d *= -1;
      j3d_clip(&poly2, &splane, 1);

      // reduce the original and its two parts
      j3d_get_moments(&opoly, om, POLY_ORDER);
      j3d_get_moments(&poly1, m1, POLY_ORDER);
      j3d_get_moments(&poly2, m2, POLY_ORDER);

      // make sure the sum of moments equals the original
      for (m = 0; m < NUM_MOMENTS; ++m) {
        ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
      }

      // make sure neither of the two resulting volumes is larger than the original
      // (within some tolerance)
      ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
      ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));

      //printf("nstack = %d, depth = %d, opoly = %.10e, p1 = %.10e, p2 = %.10e, err = %.10e\n",
      //nstack, depth, om[0], m1[0], m2[0], fabs(1.0 - om[0]/(m1[0] + m2[0])));

      // push the children to the stack if they have
      // an acceptably large volume
      if (depth < MAX_DEPTH) {
        if (m1[0] > MIN_VOL) {
          polystack[nstack] = poly1;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
        if (m2[0] > MIN_VOL) {
          polystack[nstack] = poly2;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
      }
    }
  }
  printf("=== END RECURSIVE SPLITTING DEGENERATE PERTURBED ===\n");
  return EXIT_SUCCESS;
}

int torus_load_and_chop() {
  printf("=== BEGIN TORUS LOAD AND CHOP TEST ===\n");

  const double MIN_VOL = 1.0e-8;
  const int MAX_DEPTH = 16;

  // recursively splits a polyhedron (starting from a torus) with a cut plane
  // that is always degenerate with the poly in some way,
  // checking to see that the resulting volumes add up properly.
  // This checks non-convex geometry, and also j3d_init_poly() when more than
  // three edges per vertex are needed.

  int i, j, f, m;

  // explicit stack-based implementation
  int nstack, depth, chopt;
  j3d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j3d_plane splane;
  j3d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];

  // torus parameters
  const double R_MAJ = 1.0;
  const double R_MIN = 0.1;
  const double TWOPI = 6.28318530718;

  j3d_vec3 vertices[NVERTS];

  // generate the vertex coordinates
  for (i = 0; i < NPHI; ++i) {
    for (j = 0; j < NTHETA; ++j) {
      vertices[i * NTHETA + j].xyz[0] =
        (R_MAJ + R_MIN * cos(j * TWOPI / NTHETA)) * cos(i * TWOPI / NPHI);
      vertices[i * NTHETA + j].xyz[1] =
        (R_MAJ + R_MIN * cos(j * TWOPI / NTHETA)) * sin(i * TWOPI / NPHI);
      vertices[i * NTHETA + j].xyz[2] = R_MIN * sin(j * TWOPI / NTHETA);

      // shift away from the origin so that none of the moments are identically zero
      vertices[i * NTHETA + j].xyz[0] += 0.5;
      vertices[i * NTHETA + j].xyz[1] += 0.7;
      vertices[i * NTHETA + j].xyz[2] += 0.9;
    }
  }

  // generate triangular faces
  int vertsperface[NFACES];
  int rawinds[NFACES][3];
  for (i = 0; i < NPHI; ++i)
    for (j = 0; j < NTHETA; ++j) {
      vertsperface[2 * (i * NTHETA + j)] = 3;
      rawinds[2 * (i * NTHETA + j)][0] = (i)*NTHETA + (j);
      rawinds[2 * (i * NTHETA + j)][1] = ((i + 1) % NPHI) * NTHETA + (j);
      rawinds[2 * (i * NTHETA + j)][2] = ((i + 1) % NPHI) * NTHETA + ((j + 1) % NTHETA);
      vertsperface[2 * (i * NTHETA + j) + 1] = 3;
      rawinds[2 * (i * NTHETA + j) + 1][0] = (i)*NTHETA + (j);
      rawinds[2 * (i * NTHETA + j) + 1][1] = ((i + 1) % NPHI) * NTHETA + ((j + 1) % NTHETA);
      rawinds[2 * (i * NTHETA + j) + 1][2] = (i)*NTHETA + ((j + 1) % NTHETA);
    }

  // make a double-pointer for the faces
  int* faceinds[NFACES];
  for (f = 0; f < NFACES; ++f) faceinds[f] = &rawinds[f][0];

  // initialize a general polyhedron
  j3d_init_poly(&opoly, vertices, NVERTS, faceinds, vertsperface, NFACES);

  // push the torus to the stack
  nstack = 0;
  polystack[nstack] = opoly;
  depthstack[nstack] = 0;
  ++nstack;

  // recursively split the poly
  while (nstack > 0) {
    // pop the stack
    --nstack;
    opoly = polystack[nstack];
    depth = depthstack[nstack];

    // generate a random plane from one of a few
    // possible degenerate configurations, ensuring that it
    // has a valid unit normal
    chopt = rand_int(6);
    do {
      splane = choptions_3d[chopt](&opoly);
    } while (splane.n.xyz[0] == 0.0 && splane.n.xyz[1] == 0.0 && splane.n.xyz[2] == 0.0);

    // split the poly by making two copies of the original poly
    // and them clipping them against the same plane, with one
    // oriented oppositely
    poly1 = opoly;
    poly2 = opoly;
    j3d_clip(&poly1, &splane, 1);
    splane.n.xyz[0] *= -1;
    splane.n.xyz[1] *= -1;
    splane.n.xyz[2] *= -1;
    splane.d *= -1;
    j3d_clip(&poly2, &splane, 1);

    // reduce the original and its two parts
    j3d_get_moments(&opoly, om, POLY_ORDER);
    j3d_get_moments(&poly1, m1, POLY_ORDER);
    j3d_get_moments(&poly2, m2, POLY_ORDER);

    // make sure the sum of moments equals the original
    for (m = 0; m < NUM_MOMENTS; ++m) {
      ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
    }

    // make sure neither of the two resulting volumes is larger than the original
    // (within some tolerance)
    ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
    ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));

    //printf("nstack = %d, depth = %d, opoly = %.10e, p1 = %.10e, p2 = %.10e, err = %.10e\n",
    //nstack, depth, om[0], m1[0], m2[0], fabs(1.0 - om[0]/(m1[0] + m2[0])));

    // push the children to the stack if they have
    // an acceptably large volume
    if (depth < MAX_DEPTH) {
      if (m1[0] > MIN_VOL) {
        polystack[nstack] = poly1;
        depthstack[nstack] = depth + 1;
        ++nstack;
      }
      if (m2[0] > MIN_VOL) {
        polystack[nstack] = poly2;
        depthstack[nstack] = depth + 1;
        ++nstack;
      }
    }
  }
  printf("=== END TORUS LOAD AND CHOP TEST ===\n");
  return EXIT_SUCCESS;
}

int alt_init_torus_load_and_chop() {
  printf("=== BEGIN TORUS LOAD AND CHOP TEST WITH THE ALTERNATIVE INITIALIZATION FROM BREP ===\n");

  const double MIN_VOL = 1.0e-8;
  const int MAX_DEPTH = 16;

  // recursively splits a polyhedron (starting from a torus) with a cut plane
  // that is always degenerate with the poly in some way,
  // checking to see that the resulting volumes add up properly.
  // This checks non-convex geometry, and also j3d_brep_to_poly() when more than
  // three edges per vertex are needed.

  int i, j, k, m;

  // explicit stack-based implementation
  int nstack, depth, chopt;
  j3d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j3d_plane splane;
  j3d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];

  // torus parameters
  const double R_MAJ = 1.0;
  const double R_MIN = 0.1;
  const double TWOPI = 6.28318530718;

  j3d_vec3 vertices[NVERTS];

  // generate the vertex coordinates
  for (i = 0; i < NPHI; ++i)
    for (j = 0; j < NTHETA; ++j) {
      vertices[i * NTHETA + j].xyz[0] =
        (R_MAJ + R_MIN * cos(j * TWOPI / NTHETA)) * cos(i * TWOPI / NPHI);
      vertices[i * NTHETA + j].xyz[1] =
        (R_MAJ + R_MIN * cos(j * TWOPI / NTHETA)) * sin(i * TWOPI / NPHI);
      vertices[i * NTHETA + j].xyz[2] = R_MIN * sin(j * TWOPI / NTHETA);

      // shift away from the origin so that none of the moments are identically zero
      vertices[i * NTHETA + j].xyz[0] += 0.5;
      vertices[i * NTHETA + j].xyz[1] += 0.7;
      vertices[i * NTHETA + j].xyz[2] += 0.9;
    }

  // generate triangular faces
  int facevertoffsets[NFACES + 1];
  int facevertinds[NFACES * 3];
  int face_id;
  for (i = 0; i < NPHI; ++i) {
    for (j = 0; j < NTHETA; ++j) {
      for (k = 0; k < 2; ++k) {
        face_id = 2 * (i * NTHETA + j) + k;
        facevertoffsets[face_id] = face_id * 3;
        facevertinds[face_id * 3] = (i)*NTHETA + (j);
        facevertinds[face_id * 3 + 1] = ((i + 1) % NPHI) * NTHETA + ((j + k) % NTHETA);
        facevertinds[face_id * 3 + 2] = ((i + 1 - k) % NPHI) * NTHETA + ((j + 1) % NTHETA);
      }
    }
  }
  facevertoffsets[NFACES] = NFACES * 3;

  // initialize a general polyhedron
  j3d_brep_to_poly(&opoly, vertices, NVERTS, facevertinds, facevertoffsets, NFACES);

  // push the torus to the stack
  nstack = 0;
  polystack[nstack] = opoly;
  depthstack[nstack] = 0;
  ++nstack;

  // recursively split the poly
  while (nstack > 0) {
    // pop the stack
    --nstack;
    opoly = polystack[nstack];
    depth = depthstack[nstack];

    // generate a random plane from one of a few
    // possible degenerate configurations, ensuring that it
    // has a valid unit normal
    chopt = rand_int(6);
    do {
      splane = choptions_3d[chopt](&opoly);
    } while (splane.n.xyz[0] == 0.0 && splane.n.xyz[1] == 0.0 && splane.n.xyz[2] == 0.0);

    // split the poly by making two copies of the original poly
    // and them clipping them against the same plane, with one
    // oriented oppositely
    poly1 = opoly;
    poly2 = opoly;
    j3d_clip(&poly1, &splane, 1);
    splane.n.xyz[0] *= -1;
    splane.n.xyz[1] *= -1;
    splane.n.xyz[2] *= -1;
    splane.d *= -1;
    j3d_clip(&poly2, &splane, 1);

    // reduce the original and its two parts
    j3d_get_moments(&opoly, om, POLY_ORDER);
    j3d_get_moments(&poly1, m1, POLY_ORDER);
    j3d_get_moments(&poly2, m2, POLY_ORDER);

    // make sure the sum of moments equals the original
    for (m = 0; m < NUM_MOMENTS; ++m) {
      ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
    }

    // make sure neither of the two resulting volumes is larger than the original
    // (within some tolerance)
    ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
    ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));

    //printf("nstack = %d, depth = %d, opoly = %.10e, p1 = %.10e, p2 = %.10e, err = %.10e\n",
    //nstack, depth, om[0], m1[0], m2[0], fabs(1.0 - om[0]/(m1[0] + m2[0])));

    // push the children to the stack if they have
    // an acceptably large volume
    if (depth < MAX_DEPTH) {
      if (m1[0] > MIN_VOL) {
        polystack[nstack] = poly1;
        depthstack[nstack] = depth + 1;
        ++nstack;
      }
      if (m2[0] > MIN_VOL) {
        polystack[nstack] = poly2;
        depthstack[nstack] = depth + 1;
        ++nstack;
      }
    }
  }
  printf("=== END TORUS LOAD AND CHOP TEST WITH THE ALTERNATIVE INITIALIZATION FROM BREP ===\n");
  return EXIT_SUCCESS;
}

#undef POLY_ORDER
#define POLY_ORDER 5
#undef NUM_MOMENTS
#define NUM_MOMENTS ((POLY_ORDER + 1) * (POLY_ORDER + 2) * (POLY_ORDER + 3) / 6)

int moments() {
  printf("=== BEGIN MOMENTS TEST ===\n");
  // check the moments against an analytic test case
  // (an axis-aligned box) up to some arbitrary order

  const int NUM_TRIALS = 1000;

  int i, j, k, mind, curorder;
  j3d_poly poly;
  j3d_vec3 box[2];
  double moments[NUM_MOMENTS];
  double exact;

  // initialize the box
  box[0].xyz[0] = 1.0;
  box[0].xyz[1] = 1.0;
  box[0].xyz[2] = 1.0;
  box[1].xyz[0] = 1.5;
  box[1].xyz[1] = 2.0;
  box[1].xyz[2] = 2.5;
  j3d_init_box(&poly, box);


  int trial;
  // printf("Computing moments of order %d, %d trials.\n", POLY_ORDER, NUM_TRIALS);
  for (trial = 0; trial < NUM_TRIALS; ++trial) {
    // get the moments from j3d_get_moments
    j3d_get_moments(&poly, moments, POLY_ORDER);
    // Check all the moments against the analytic solution
    // NOTE: This is the order that j3d_get_moments puts out!
    for (curorder = 0, mind = 0; curorder <= POLY_ORDER; ++curorder) {
      //printf("Order = %d\n", curorder);
      for (i = curorder; i >= 0; --i)
        for (j = curorder - i; j >= 0; --j, ++mind) {
          k = curorder - i - j;
          exact = 1.0 / ((i + 1) * (j + 1) * (k + 1)) *
                  (pow(box[1].xyz[0], i + 1) - pow(box[0].xyz[0], i + 1)) *
                  (pow(box[1].xyz[1], j + 1) - pow(box[0].xyz[1], j + 1)) *
                  (pow(box[1].xyz[2], k + 1) - pow(box[0].xyz[2], k + 1));
          // printf(" Int[ x^%d y^%d z^%d dV ] = %.10e, analytic = %.10e, frac = %f, error = %.10e\n",
          //        i, j, k, moments[mind], exact, moments[mind]/exact, fabs(1.0 - moments[mind]/exact));
          ASSERT_DOUBLE_EQ(moments[mind], exact, TOL);
        }
    }
  }
  printf("=== END MOMENTS TEST ===\n");
  return EXIT_SUCCESS;
}

int clip_and_split_tet() {
  printf("=== BEGIN CLIP AND SPLIT TET TEST ===\n");
  // number of trials
  const int NUM_TRIALS = 100;

  // initialize tet in the positive quadrant
  j3d_poly tet;
  // initialize number of edges and vertices
  tet.num_verts = 4;
  tet.num_dir_edges = 12;
  // initialize vertices
  tet.verts[0][0] = 0.0;
  tet.verts[0][1] = 0.0;
  tet.verts[0][2] = 0.0;
  tet.verts[1][0] = 1.0;
  tet.verts[1][1] = 0.0;
  tet.verts[1][2] = 0.0;
  tet.verts[2][0] = 0.0;
  tet.verts[2][1] = 1.0;
  tet.verts[2][2] = 0.0;
  tet.verts[3][0] = 0.0;
  tet.verts[3][1] = 0.0;
  tet.verts[3][2] = 1.0;
  // initialize indices
  tet.indices[0] = 0;
  tet.indices[1] = 3;
  tet.indices[2] = 6;
  tet.indices[3] = 9;
  tet.indices[4] = 12;
  // initialize neighbors
  tet.dir_edges[0] = 3;
  tet.dir_edges[1] = 2;
  tet.dir_edges[2] = 1;
  tet.dir_edges[3] = 2;
  tet.dir_edges[4] = 3;
  tet.dir_edges[5] = 0;
  tet.dir_edges[6] = 3;
  tet.dir_edges[7] = 1;
  tet.dir_edges[8] = 0;
  tet.dir_edges[9] = 1;
  tet.dir_edges[10] = 2;
  tet.dir_edges[11] = 0;

  // some variables for use later
  j3d_poly tet_clip_1;
  j3d_poly tet_clip_2;
  j3d_poly tet_split;
  j3d_plane rand_plane;

  for (int i = 0; i < NUM_TRIALS; ++i) {
    // copy polys
    j3d_copy_poly(&tet_clip_1, &tet);
    j3d_copy_poly(&tet_clip_2, &tet);
    j3d_copy_poly(&tet_split, &tet);
    // make random plane through centroid
    rand_plane = thru_cent_3d(&tet);
    // split
    j3d_poly pos_tet;
    j3d_poly neg_tet;
    j3d_split(&tet_split, 1, rand_plane, &pos_tet, &neg_tet);
    // clip
    j3d_clip(&tet_clip_1, &rand_plane, 1);
    rand_plane.d *= -1.0;
    rand_plane.n.xyz[0] *= -1.0;
    rand_plane.n.xyz[1] *= -1.0;
    rand_plane.n.xyz[2] *= -1.0;
    j3d_clip(&tet_clip_2, &rand_plane, 1);
    // variables to store moments
    double moments_pos_tet_1[20];
    double moments_pos_tet_2[20];
    double moments_neg_tet_1[20];
    double moments_neg_tet_2[20];
    double moments_tet_clip_1_1[20];
    double moments_tet_clip_1_2[20];
    double moments_tet_clip_2_1[20];
    double moments_tet_clip_2_2[20];
    // we expect that the moments of pos_tet == tet_clip_1
    // we expect that the moments of neg_tet == tet_clip_2
    // ZEROTH ORDER MOMENTS
    j3d_get_moments_0(&pos_tet, moments_pos_tet_1);
    j3d_get_moments_n(&pos_tet, moments_pos_tet_2, 0);
    j3d_get_moments_0(&neg_tet, moments_neg_tet_1);
    j3d_get_moments_n(&neg_tet, moments_neg_tet_2, 0);
    j3d_get_moments_0(&tet_clip_1, moments_tet_clip_1_1);
    j3d_get_moments_n(&tet_clip_1, moments_tet_clip_1_2, 0);
    j3d_get_moments_0(&tet_clip_2, moments_tet_clip_2_1);
    j3d_get_moments_n(&tet_clip_2, moments_tet_clip_2_2, 0);
    // check that recursive calculation is equivalent to explicit
    ASSERT_DOUBLE_EQ(moments_pos_tet_1[0], moments_pos_tet_2[0], TOL);
    ASSERT_DOUBLE_EQ(moments_neg_tet_1[0], moments_neg_tet_2[0], TOL);
    ASSERT_DOUBLE_EQ(moments_tet_clip_1_1[0], moments_tet_clip_1_2[0], TOL);
    ASSERT_DOUBLE_EQ(moments_tet_clip_2_1[0], moments_tet_clip_2_2[0], TOL);
    // check that moments of intersected poly on either side of plane
    // are equivalent
    ASSERT_DOUBLE_EQ(moments_pos_tet_1[0], moments_tet_clip_1_1[0], TOL);
    ASSERT_DOUBLE_EQ(moments_neg_tet_1[0], moments_tet_clip_2_1[0], TOL);
    // FIRST ORDER MOMENTS
    j3d_get_moments_1(&pos_tet, moments_pos_tet_1);
    j3d_get_moments_n(&pos_tet, moments_pos_tet_2, 1);
    j3d_get_moments_1(&neg_tet, moments_neg_tet_1);
    j3d_get_moments_n(&neg_tet, moments_neg_tet_2, 1);
    j3d_get_moments_1(&tet_clip_1, moments_tet_clip_1_1);
    j3d_get_moments_n(&tet_clip_1, moments_tet_clip_1_2, 1);
    j3d_get_moments_1(&tet_clip_2, moments_tet_clip_2_1);
    j3d_get_moments_n(&tet_clip_2, moments_tet_clip_2_2, 1);
    // check that recursive calculation is equivalent to explicit
    for (int j = 0; j < 4; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tet_1[j], moments_pos_tet_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tet_1[j], moments_neg_tet_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tet_clip_1_1[j], moments_tet_clip_1_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tet_clip_2_1[j], moments_tet_clip_2_2[j], TOL);
    }
    // check that moments of intersected poly on either side of plane
    // are equivalent
    for (int j = 0; j < 4; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tet_1[j], moments_tet_clip_1_1[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tet_1[j], moments_tet_clip_2_1[j], TOL);
    }
    // SECOND ORDER MOMENTS
    j3d_get_moments_2(&pos_tet, moments_pos_tet_1);
    j3d_get_moments_n(&pos_tet, moments_pos_tet_2, 2);
    j3d_get_moments_2(&neg_tet, moments_neg_tet_1);
    j3d_get_moments_n(&neg_tet, moments_neg_tet_2, 2);
    j3d_get_moments_2(&tet_clip_1, moments_tet_clip_1_1);
    j3d_get_moments_n(&tet_clip_1, moments_tet_clip_1_2, 2);
    j3d_get_moments_2(&tet_clip_2, moments_tet_clip_2_1);
    j3d_get_moments_n(&tet_clip_2, moments_tet_clip_2_2, 2);
    // check that recursive calculation is equivalent to explicit
    for (int j = 0; j < 10; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tet_1[j], moments_pos_tet_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tet_1[j], moments_neg_tet_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tet_clip_1_1[j], moments_tet_clip_1_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tet_clip_2_1[j], moments_tet_clip_2_2[j], TOL);
    }
    // check that moments of intersected poly on either side of plane
    // are equivalent
    for (int j = 0; j < 10; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tet_1[j], moments_tet_clip_1_1[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tet_1[j], moments_tet_clip_2_1[j], TOL);
    }
    // THIRD ORDER MOMENTS
    j3d_get_moments_3(&pos_tet, moments_pos_tet_1);
    j3d_get_moments_n(&pos_tet, moments_pos_tet_2, 3);
    j3d_get_moments_3(&neg_tet, moments_neg_tet_1);
    j3d_get_moments_n(&neg_tet, moments_neg_tet_2, 3);
    j3d_get_moments_3(&tet_clip_1, moments_tet_clip_1_1);
    j3d_get_moments_n(&tet_clip_1, moments_tet_clip_1_2, 3);
    j3d_get_moments_3(&tet_clip_2, moments_tet_clip_2_1);
    j3d_get_moments_n(&tet_clip_2, moments_tet_clip_2_2, 3);
    // check that recursive calculation is equivalent to explicit
    for (int j = 0; j < 20; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tet_1[j], moments_pos_tet_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tet_1[j], moments_neg_tet_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tet_clip_1_1[j], moments_tet_clip_1_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tet_clip_2_1[j], moments_tet_clip_2_2[j], TOL);
    }
    // check that moments of intersected poly on either side of plane
    // are equivalent
    for (int j = 0; j < 20; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tet_1[j], moments_tet_clip_1_1[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tet_1[j], moments_tet_clip_2_1[j], TOL);
    }
  }
  printf("=== END CLIP AND SPLIT TET TEST ===\n");
  return EXIT_SUCCESS;
}

int clip_and_split_concave_tri_bipyramid() {
  printf("=== BEGIN CLIP AND SPLIT CONCAVE TRI BIPYRAMID TEST ===\n");
  // number of trials
  const int NUM_TRIALS = 100;

  // initialize tri in the positive quadrant
  j3d_poly tri;
  // initialize number of edges and vertices
  tri.num_verts = 5;
  tri.num_dir_edges = 18;
  // initialize vertices
  tri.verts[0][0] = 0.0;
  tri.verts[0][1] = 0.0;
  tri.verts[0][2] = 0.0;
  tri.verts[1][0] = 0.0;
  tri.verts[1][1] = 3.0;
  tri.verts[1][2] = 0.0;
  tri.verts[2][0] = 0.0;
  tri.verts[2][1] = 0.0;
  tri.verts[2][2] = 3.0;
  tri.verts[3][0] = 5.0;
  tri.verts[3][1] = 1.0;
  tri.verts[3][2] = 1.0;
  tri.verts[4][0] = 8.0;
  tri.verts[4][1] = 1.0;
  tri.verts[4][2] = 1.0;

  // initialize indices
  tri.indices[0] = 0;
  tri.indices[1] = 4;
  tri.indices[2] = 8;
  tri.indices[3] = 12;
  tri.indices[4] = 15;
  tri.indices[5] = 18;
  // initialize neighbors
  tri.dir_edges[0] = 2;
  tri.dir_edges[1] = 3;
  tri.dir_edges[2] = 1;
  tri.dir_edges[3] = 4;
  tri.dir_edges[4] = 2;
  tri.dir_edges[5] = 4;
  tri.dir_edges[6] = 0;
  tri.dir_edges[7] = 3;
  tri.dir_edges[8] = 1;
  tri.dir_edges[9] = 3;
  tri.dir_edges[10] = 0;
  tri.dir_edges[11] = 4;
  tri.dir_edges[12] = 2;
  tri.dir_edges[13] = 1;
  tri.dir_edges[14] = 0;
  tri.dir_edges[15] = 2;
  tri.dir_edges[16] = 0;
  tri.dir_edges[17] = 1;

  // some variables for use later
  j3d_poly tri_clip_1;
  j3d_poly tri_clip_2;
  j3d_poly tri_split;
  j3d_plane rand_plane;

  for (int i = 0; i < NUM_TRIALS; ++i) {
    // copy polys
    j3d_copy_poly(&tri_clip_1, &tri);
    j3d_copy_poly(&tri_clip_2, &tri);
    j3d_copy_poly(&tri_split, &tri);
    // make random plane through centroid
    rand_plane = thru_cent_3d(&tri);
    // split
    j3d_poly pos_tri;
    j3d_poly neg_tri;
    j3d_split(&tri_split, 1, rand_plane, &pos_tri, &neg_tri);
    // clip
    j3d_clip(&tri_clip_1, &rand_plane, 1);
    rand_plane.d *= -1.0;
    rand_plane.n.xyz[0] *= -1.0;
    rand_plane.n.xyz[1] *= -1.0;
    rand_plane.n.xyz[2] *= -1.0;
    j3d_clip(&tri_clip_2, &rand_plane, 1);
    // variables to store moments
    double moments_pos_tri_1[20];
    double moments_pos_tri_2[20];
    double moments_neg_tri_1[20];
    double moments_neg_tri_2[20];
    double moments_tri_clip_1_1[20];
    double moments_tri_clip_1_2[20];
    double moments_tri_clip_2_1[20];
    double moments_tri_clip_2_2[20];
    // we expect that the moments of pos_tri == tri_clip_1
    // we expect that the moments of neg_tri == tri_clip_2
    // ZEROTH ORDER MOMENTS
    j3d_get_moments_0(&pos_tri, moments_pos_tri_1);
    j3d_get_moments_n(&pos_tri, moments_pos_tri_2, 0);
    j3d_get_moments_0(&neg_tri, moments_neg_tri_1);
    j3d_get_moments_n(&neg_tri, moments_neg_tri_2, 0);
    j3d_get_moments_0(&tri_clip_1, moments_tri_clip_1_1);
    j3d_get_moments_n(&tri_clip_1, moments_tri_clip_1_2, 0);
    j3d_get_moments_0(&tri_clip_2, moments_tri_clip_2_1);
    j3d_get_moments_n(&tri_clip_2, moments_tri_clip_2_2, 0);
    // check that recursive calculation is equivalent to explicit
    ASSERT_DOUBLE_EQ(moments_pos_tri_1[0], moments_pos_tri_2[0], TOL);
    ASSERT_DOUBLE_EQ(moments_neg_tri_1[0], moments_neg_tri_2[0], TOL);
    ASSERT_DOUBLE_EQ(moments_tri_clip_1_1[0], moments_tri_clip_1_2[0], TOL);
    ASSERT_DOUBLE_EQ(moments_tri_clip_2_1[0], moments_tri_clip_2_2[0], TOL);
    // check that moments of intersected poly on either side of plane
    // are equivalent
    ASSERT_DOUBLE_EQ(moments_pos_tri_1[0], moments_tri_clip_1_1[0], TOL);
    ASSERT_DOUBLE_EQ(moments_neg_tri_1[0], moments_tri_clip_2_1[0], TOL);
    // FIRST ORDER MOMENTS
    j3d_get_moments_1(&pos_tri, moments_pos_tri_1);
    j3d_get_moments_n(&pos_tri, moments_pos_tri_2, 1);
    j3d_get_moments_1(&neg_tri, moments_neg_tri_1);
    j3d_get_moments_n(&neg_tri, moments_neg_tri_2, 1);
    j3d_get_moments_1(&tri_clip_1, moments_tri_clip_1_1);
    j3d_get_moments_n(&tri_clip_1, moments_tri_clip_1_2, 1);
    j3d_get_moments_1(&tri_clip_2, moments_tri_clip_2_1);
    j3d_get_moments_n(&tri_clip_2, moments_tri_clip_2_2, 1);
    // check that recursive calculation is equivalent to explicit
    for (int j = 0; j < 4; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tri_1[j], moments_pos_tri_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tri_1[j], moments_neg_tri_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tri_clip_1_1[j], moments_tri_clip_1_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tri_clip_2_1[j], moments_tri_clip_2_2[j], TOL);
    }
    // check that moments of intersected poly on either side of plane
    // are equivalent
    for (int j = 0; j < 4; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tri_1[j], moments_tri_clip_1_1[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tri_1[j], moments_tri_clip_2_1[j], TOL);
    }
    // SECOND ORDER MOMENTS
    j3d_get_moments_2(&pos_tri, moments_pos_tri_1);
    j3d_get_moments_n(&pos_tri, moments_pos_tri_2, 2);
    j3d_get_moments_2(&neg_tri, moments_neg_tri_1);
    j3d_get_moments_n(&neg_tri, moments_neg_tri_2, 2);
    j3d_get_moments_2(&tri_clip_1, moments_tri_clip_1_1);
    j3d_get_moments_n(&tri_clip_1, moments_tri_clip_1_2, 2);
    j3d_get_moments_2(&tri_clip_2, moments_tri_clip_2_1);
    j3d_get_moments_n(&tri_clip_2, moments_tri_clip_2_2, 2);
    // check that recursive calculation is equivalent to explicit
    for (int j = 0; j < 10; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tri_1[j], moments_pos_tri_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tri_1[j], moments_neg_tri_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tri_clip_1_1[j], moments_tri_clip_1_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tri_clip_2_1[j], moments_tri_clip_2_2[j], TOL);
    }
    // check that moments of intersected poly on either side of plane
    // are equivalent
    for (int j = 0; j < 10; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tri_1[j], moments_tri_clip_1_1[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tri_1[j], moments_tri_clip_2_1[j], TOL);
    }
    // THIRD ORDER MOMENTS
    j3d_get_moments_3(&pos_tri, moments_pos_tri_1);
    j3d_get_moments_n(&pos_tri, moments_pos_tri_2, 3);
    j3d_get_moments_3(&neg_tri, moments_neg_tri_1);
    j3d_get_moments_n(&neg_tri, moments_neg_tri_2, 3);
    j3d_get_moments_3(&tri_clip_1, moments_tri_clip_1_1);
    j3d_get_moments_n(&tri_clip_1, moments_tri_clip_1_2, 3);
    j3d_get_moments_3(&tri_clip_2, moments_tri_clip_2_1);
    j3d_get_moments_n(&tri_clip_2, moments_tri_clip_2_2, 3);
    // check that recursive calculation is equivalent to explicit
    for (int j = 0; j < 20; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tri_1[j], moments_pos_tri_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tri_1[j], moments_neg_tri_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tri_clip_1_1[j], moments_tri_clip_1_2[j], TOL);
      ASSERT_DOUBLE_EQ(moments_tri_clip_2_1[j], moments_tri_clip_2_2[j], TOL);
    }
    // check that moments of intersected poly on either side of plane
    // are equivalent
    for (int j = 0; j < 20; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_tri_1[j], moments_tri_clip_1_1[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_tri_1[j], moments_tri_clip_2_1[j], TOL);
    }
  }

  printf("=== END CLIP AND SPLIT CONCAVE TRI BIPYRAMID TEST ===\n");
  return EXIT_SUCCESS;
}

int copy_tet() {
  printf("=== BEGIN COPY TET TEST ===\n");

  j3d_poly tet;
  j3d_vec3 vecs[4];

  vecs[0].xyz[0] = 0.0;
  vecs[0].xyz[1] = 0.0;
  vecs[0].xyz[2] = 0.0;
  vecs[1].xyz[0] = 1.0;
  vecs[1].xyz[1] = 0.0;
  vecs[1].xyz[2] = 0.0;
  vecs[2].xyz[0] = 0.0;
  vecs[2].xyz[1] = 1.0;
  vecs[2].xyz[2] = 0.0;
  vecs[3].xyz[0] = 0.0;
  vecs[3].xyz[1] = 0.0;
  vecs[3].xyz[2] = 1.0;

  j3d_init_tet(&tet, vecs);

  j3d_poly tet_copy;
  j3d_copy_poly(&tet_copy, &tet);

  ASSERT_INT_EQ(tet_copy.num_verts, 4);
  ASSERT_INT_EQ(tet_copy.num_dir_edges, 12);

  ASSERT_INT_EQ(tet_copy.indices[0], 0);
  ASSERT_INT_EQ(tet_copy.indices[1], 3);
  ASSERT_INT_EQ(tet_copy.indices[2], 6);
  ASSERT_INT_EQ(tet_copy.indices[3], 9);
  ASSERT_INT_EQ(tet_copy.indices[4], 12);

  ASSERT_INT_EQ(tet_copy.dir_edges[0], 1);
  ASSERT_INT_EQ(tet_copy.dir_edges[1], 3);
  ASSERT_INT_EQ(tet_copy.dir_edges[2], 2);
  ASSERT_INT_EQ(tet_copy.dir_edges[3], 2);
  ASSERT_INT_EQ(tet_copy.dir_edges[4], 3);
  ASSERT_INT_EQ(tet_copy.dir_edges[5], 0);
  ASSERT_INT_EQ(tet_copy.dir_edges[6], 0);
  ASSERT_INT_EQ(tet_copy.dir_edges[7], 3);
  ASSERT_INT_EQ(tet_copy.dir_edges[8], 1);
  ASSERT_INT_EQ(tet_copy.dir_edges[9], 1);
  ASSERT_INT_EQ(tet_copy.dir_edges[10], 2);
  ASSERT_INT_EQ(tet_copy.dir_edges[11], 0);

  ASSERT_DOUBLE_EQ(tet_copy.verts[0][0], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[0][1], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[0][2], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[1][0], 1.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[1][1], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[1][2], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[2][0], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[2][1], 1.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[2][2], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[3][0], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[3][1], 0.0, TOL);
  ASSERT_DOUBLE_EQ(tet_copy.verts[3][2], 1.0, TOL);

  printf("=== END COPY TET TEST ===\n");
  return EXIT_SUCCESS;
}

#undef NUM_TESTS
#define NUM_TESTS 11

int main() {
  choptions_3d[0] = thru_cent_3d;
  choptions_3d[1] = thru_face_3d;
  choptions_3d[2] = thru_edge_cent_3d;
  choptions_3d[3] = thru_edge_rand_3d;
  choptions_3d[4] = thru_vert_cent_3d;
  choptions_3d[5] = thru_vert_rand_3d;

  int tests[NUM_TESTS];

  tests[0] = split_tets_thru_centroid();
  tests[1] = split_nonconvex();
  tests[2] = recursive_splitting_nondegenerate();
  tests[3] = recursive_splitting_degenerate();
  tests[4] = recursive_splitting_degenerate_perturbed();
  tests[5] = torus_load_and_chop();
  tests[6] = alt_init_torus_load_and_chop();
  tests[7] = moments();
  tests[8] = clip_and_split_tet();
  tests[9] = clip_and_split_concave_tri_bipyramid();
  tests[10] = copy_tet();

  for (int i = 0; i < NUM_TESTS; ++i) {
    if (tests[i] == EXIT_FAILURE) {
      return 1;
    }
  }

  return 0;
}

#undef POLY_ORDER
#undef NUM_MOMENTS
#undef NTETS
#undef STACK_SIZE
#undef NTHETA
#undef NPHI
#undef NVERTS
#undef NFACES
#undef NUM_TESTS
