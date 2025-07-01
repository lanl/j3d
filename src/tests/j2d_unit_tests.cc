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
#undef NTRIS
#undef STACK_SIZE
#undef NVERTS

#define POLY_ORDER 2
#define NUM_MOMENTS ((POLY_ORDER + 1) * (POLY_ORDER + 2) / 2)
#define NTRIS 128
#define STACK_SIZE 64
#define NVERTS 20

const double TOL = 1e-4;

j2d_plane (*choptions_2d[4])(j2d_poly* poly);

int split_tris_thru_centroid() {
  printf("=== BEGIN SPLIT TRI THRU CENTROID TESTS ===\n");

  const double MIN_AREA = 1e-8;

  // a very basic sanity check. Splits a tet through its centroid
  // and checks to see whether the two resulting volumes add up to equal the original

  // variables: the polyhedra and their moments
  int m, i;
  j2d_vec2 verts[3];
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];
  j2d_plane splane;

  // generate a random tet and clip plane
  j2d_poly opoly[NTRIS], poly1[NTRIS], poly2[NTRIS];
  for (i = 0; i < NTRIS; ++i) {
    rand_tri_2d(verts, MIN_AREA);
    j2d_init_poly(&opoly[i], verts, 3);
    if (i == 13) {
      splane = thru_cent_2d(&opoly[i]);
    }
  }

  // split them all about the same plane
  j2d_split(opoly, NTRIS, splane, poly1, poly2);

  for (i = 0; i < NTRIS; ++i) {
    // reduce the original and its two parts
    j2d_get_moments(&opoly[i], om, POLY_ORDER);
    j2d_get_moments(&poly1[i], m1, POLY_ORDER);
    j2d_get_moments(&poly2[i], m2, POLY_ORDER);

    // make sure the sum of moments equals the original
    for (m = 0; m < NUM_MOMENTS; ++m) {
      ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
    }

    //printf(" original = %f, parts = %f %f, sum = %f\n", om[0], m1[0], m2[0], m1[0]+m2[0]);

    // make sure neither of the two resulting volumes is larger than the original
    // (within some tolerance)
    ASSERT_DOUBLE_LT(m1[0], om[0] * (1.0 + TOL));
    ASSERT_DOUBLE_LT(m2[0], om[0] * (1.0 + TOL));
  }
  printf("=== END SPLIT TRI THRU CENTROID TESTS ===\n");
  return EXIT_SUCCESS;
}

int recursive_splitting_nondegenerate() {
  printf("=== BEGIN RECURSIVE SPLITTING NONDEGENERATE TEST ===\n");

  const int MAX_DEPTH = 16;
  const int NUM_TRIALS = 100;
  const double MIN_AREA = 1e-8;

  // recursively splits a polygon (starting from a tet) through the centroid,
  // checking to see that the resulting volumes add up properly.

  // explicit stack-based implementation
  int nstack, depth, t, m;
  j2d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j2d_vec2 verts[3];
  j2d_plane splane;
  j2d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];

  // do many trials
  // printf("Recursively splitting %d triangles, maximum splits per tri is %d.\n", NUM_TRIALS, MAX_DEPTH);
  for (t = 0; t < NUM_TRIALS; ++t) {
    // generate a random tet
    rand_tri_2d(verts, MIN_AREA);
    j2d_init_poly(&opoly, verts, 3);

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
      splane = thru_cent_2d(&opoly);

      // split the poly by making two copies of the original poly
      // and them clipping them against the same plane, with one
      // oriented oppositely
      poly1 = opoly;
      poly2 = opoly;
      j2d_clip(&poly1, &splane, 1);
      splane.n.xy[0] *= -1;
      splane.n.xy[1] *= -1;
      splane.d *= -1;
      j2d_clip(&poly2, &splane, 1);

      // reduce the original and its two parts
      j2d_get_moments(&opoly, om, POLY_ORDER);
      j2d_get_moments(&poly1, m1, POLY_ORDER);
      j2d_get_moments(&poly2, m2, POLY_ORDER);

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
        if (m1[0] > MIN_AREA) {
          polystack[nstack] = poly1;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
        if (m2[0] > MIN_AREA) {
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

  const int MAX_DEPTH = 16;
  const int NUM_TRIALS = 100;
  const double MIN_AREA = 1e-8;

  // recursively splits a polygon (starting from a tet) with a cut plane
  // that is always degenerate with the tet in some way,
  // checking to see that the resulting volumes add up properly.

  // explicit stack-based implementation
  int nstack, depth, t, chopt, m;
  j2d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j2d_vec2 verts[3];
  j2d_plane splane;
  j2d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];

  // do many trials
  // printf("Recursively splitting %d triangles, maximum splits per tri is %d.\n", NUM_TRIALS, MAX_DEPTH);
  for (t = 0; t < NUM_TRIALS; ++t) {
    // generate a random tet
    rand_tri_2d(verts, MIN_AREA);
    j2d_init_poly(&opoly, verts, 3);

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
      chopt = rand_int(4);
      do {
        splane = choptions_2d[chopt](&opoly);
      } while (splane.n.xy[0] == 0.0 && splane.n.xy[1] == 0.0);

      // split the poly by making two copies of the original poly
      // and them clipping them against the same plane, with one
      // oriented oppositely
      poly1 = opoly;
      poly2 = opoly;
      j2d_clip(&poly1, &splane, 1);
      splane.n.xy[0] *= -1;
      splane.n.xy[1] *= -1;
      splane.d *= -1;
      j2d_clip(&poly2, &splane, 1);

      // reduce the original and its two parts
      j2d_get_moments(&opoly, om, POLY_ORDER);
      j2d_get_moments(&poly1, m1, POLY_ORDER);
      j2d_get_moments(&poly2, m2, POLY_ORDER);

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
        if (m1[0] > MIN_AREA) {
          polystack[nstack] = poly1;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
        if (m2[0] > MIN_AREA) {
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
  printf("=== BEGIN RECURSIVE SPLITTING DEGENERATE PERTURBED TEST ===\n");

  const int MIN_PERTURB_ORDER = -17;
  const int MAX_PERTURB_ORDER = -1;
  const int MAX_DEPTH = 16;
  const int NUM_TRIALS = 100;
  const double MIN_AREA = 1e-8;

  // recursively splits a polygon (starting from a tet) with a cut plane
  // that is always degenerate with the tet in some way,
  // checking to see that the resulting volumes add up properly.
  // In this one, the cut plane is perturbed by a small amount.

  // explicit stack-based implementation
  int nstack, depth, t, chopt, m;
  j2d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j2d_vec2 verts[3];
  j2d_plane splane;
  j2d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];
  double perturb;

  // do many trials
  // printf("Recursively splitting %d triangles, maximum splits per tri is %d.\n", NUM_TRIALS, MAX_DEPTH);
  for (t = 0; t < NUM_TRIALS; ++t) {
    // compute the order of magnitude by which to perturb the clip plane,
    // determined by the trial number
    perturb = pow(10, MIN_PERTURB_ORDER + t % (MAX_PERTURB_ORDER - MIN_PERTURB_ORDER));
    //printf("omag = %d, pow = %.10e\n", MIN_PERTURB_ORDER + t%(MAX_PERTURB_ORDER - MIN_PERTURB_ORDER), perturb);

    // generate a random tet
    rand_tri_2d(verts, MIN_AREA);
    j2d_init_poly(&opoly, verts, 3);

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
      chopt = rand_int(4);
      do {
        splane = choptions_2d[chopt](&opoly);
      } while (splane.n.xy[0] == 0.0 && splane.n.xy[1] == 0.0);

      // randomly perturb the plane
      splane.n.xy[0] *= 1.0 + perturb * (rand_uniform() - 0.5);
      splane.n.xy[1] *= 1.0 + perturb * (rand_uniform() - 0.5);
      splane.d *= 1.0 + perturb * (rand_uniform() - 0.5);

      // split the poly by making two copies of the original poly
      // and them clipping them against the same plane, with one
      // oriented oppositely
      poly1 = opoly;
      poly2 = opoly;
      j2d_clip(&poly1, &splane, 1);
      splane.n.xy[0] *= -1;
      splane.n.xy[1] *= -1;
      splane.d *= -1;
      j2d_clip(&poly2, &splane, 1);

      // reduce the original and its two parts
      j2d_get_moments(&opoly, om, POLY_ORDER);
      j2d_get_moments(&poly1, m1, POLY_ORDER);
      j2d_get_moments(&poly2, m2, POLY_ORDER);

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
        if (m1[0] > MIN_AREA) {
          polystack[nstack] = poly1;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
        if (m2[0] > MIN_AREA) {
          polystack[nstack] = poly2;
          depthstack[nstack] = depth + 1;
          ++nstack;
        }
      }
    }
  }
  printf("=== END RECURSIVE SPLITTING DEGENERATE PERTURBED TEST ===\n");
  return EXIT_SUCCESS;
}

int random_verts() {
  printf("=== BEGIN RANDOM VERTS TEST ===\n");

  const int MAX_DEPTH = 16;
  const double MIN_AREA = 1e-8;

  // recursively splits a polygon with randomly generated verts with a cut plane
  // that is always degenerate with the poly in some way,
  // checking to see that the resulting volumes add up properly.

  int v, m;

  // explicit stack-based implementation
  int nstack, depth, chopt;
  j2d_poly polystack[STACK_SIZE];
  int depthstack[STACK_SIZE];

  // variables: the polyhedra and their moments
  j2d_plane splane;
  j2d_vec2 verts[NVERTS];
  j2d_poly opoly, poly1, poly2;
  double om[NUM_MOMENTS];
  double m1[NUM_MOMENTS];
  double m2[NUM_MOMENTS];

  // randomly generate vertices
  // Note: This will be nasty, tangled, and nonconvex
  for (v = 0; v < NVERTS; ++v) {
    verts[v] = rand_uvec_2d();
  }

  // initialize a general polygon
  j2d_init_poly(&opoly, verts, NVERTS);

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
    chopt = rand_int(4);
    do {
      splane = choptions_2d[chopt](&opoly);
    } while (splane.n.xy[0] == 0.0 && splane.n.xy[1] == 0.0);
    //splane = thru_cent_2d(&opoly);

    // split the poly by making two copies of the original poly
    // and them clipping them against the same plane, with one
    // oriented oppositely
    poly1 = opoly;
    poly2 = opoly;
    j2d_clip(&poly1, &splane, 1);
    splane.n.xy[0] *= -1;
    splane.n.xy[1] *= -1;
    splane.d *= -1;
    j2d_clip(&poly2, &splane, 1);

    // reduce the original and its two parts
    j2d_get_moments(&opoly, om, POLY_ORDER);
    j2d_get_moments(&poly1, m1, POLY_ORDER);
    j2d_get_moments(&poly2, m2, POLY_ORDER);

    // make sure the sum of moments equals the original
    for (m = 0; m < NUM_MOMENTS; ++m) {
      ASSERT_DOUBLE_EQ(om[m], m1[m] + m2[m], TOL);
    }

    // make sure neither of the two resulting volumes is larger than the original
    // (within some tolerance)
    // Note: Not here! We allow negative volumes and moments for inside-out polygons...
    // ASSERT_DOUBLE_LT(m1[0], om[0]*(1.0 + TOL));
    // ASSERT_DOUBLE_LT(m2[0], om[0]*(1.0 + TOL));

    // printf("nstack = %d, depth = %d, opoly = %.10e, p1 = %.10e, p2 = %.10e, err = %.10e\n",
    // nstack, depth, om[0], m1[0], m2[0], fabs(1.0 - om[0]/(m1[0] + m2[0])));

    // push the children to the stack if they have
    // an acceptably large volume
    if (depth < MAX_DEPTH) {
      if (fabs(m1[0]) > MIN_AREA) {
        polystack[nstack] = poly1;
        depthstack[nstack] = depth + 1;
        ++nstack;
      }
      if (fabs(m2[0]) > MIN_AREA) {
        polystack[nstack] = poly2;
        depthstack[nstack] = depth + 1;
        ++nstack;
      }
    }
  }
  printf("=== END RANDOM VERTS TEST ===\n");
  return EXIT_SUCCESS;
}

#undef POLY_ORDER
#define POLY_ORDER 5
#undef NUM_MOMENTS
#define NUM_MOMENTS ((POLY_ORDER + 1) * (POLY_ORDER + 2) / 2)

int moments() {
  printf("=== BEGIN MOMENTS TEST ===\n");

  // check the moments against an analytic test case
  // (an axis-aligned box) up to some arbitrary order

  int i, j, mind, curorder;
  j2d_poly poly;
  j2d_vec2 box[2];
  double moments[NUM_MOMENTS];
  double exact;

  // initialize the box
  box[0].xy[0] = -1.0;
  box[0].xy[1] = 1.0;
  box[1].xy[0] = 1.5;
  box[1].xy[1] = 2.0;
  j2d_init_box(&poly, box);

  // get the moments from j2d_get_moments
  j2d_get_moments(&poly, moments, POLY_ORDER);

  // Check all the moments against the analytic solution
  // NOTE: This is the order that j2d_get_moments puts out!
  for (curorder = 0, mind = 0; curorder <= POLY_ORDER; ++curorder) {
    //printf("Order = %d\n", curorder);
    for (i = curorder; i >= 0; --i, ++mind) {
      j = curorder - i;
      exact = 1.0 / ((i + 1) * (j + 1)) * (pow(box[1].xy[0], i + 1) - pow(box[0].xy[0], i + 1)) *
              (pow(box[1].xy[1], j + 1) - pow(box[0].xy[1], j + 1));
      // printf(" Int[ x^%d y^%d dA ] = %.10e, analytic = %.10e, frac = %f, error = %.10e\n",
      // i, j, moments[mind], exact, moments[mind]/exact, fabs(1.0 - moments[mind]/exact));
      ASSERT_DOUBLE_EQ(moments[mind], exact, TOL);
    }
  }
  printf("=== END MOMENTS TEST ===\n");
  return EXIT_SUCCESS;
}


int split_and_clip_square() {
  printf("=== BEGIN SPLIT AND CLIP SQUARE TEST ===\n");

  const int NUM_TRIALS = 100;

  // initialize square in the positive quadrant
  j2d_poly sqr;
  // initialize number of vertices
  sqr.num_verts = 4;
  // initialize vertices
  sqr.verts[0][0] = 0.0;
  sqr.verts[0][1] = 0.0;
  sqr.verts[1][0] = 1.0;
  sqr.verts[1][1] = 0.0;
  sqr.verts[2][0] = 1.0;
  sqr.verts[2][1] = 0.0;
  sqr.verts[3][0] = 1.0;
  sqr.verts[3][1] = 1.0;
  // initialize neighbors
  sqr.dir_edges[0] = 1;
  sqr.dir_edges[1] = 2;
  sqr.dir_edges[2] = 3;
  sqr.dir_edges[3] = 0;
  sqr.dir_edges[4] = 0;
  sqr.dir_edges[5] = 3;
  sqr.dir_edges[6] = 2;
  sqr.dir_edges[7] = 1;

  // some variables for use later
  j2d_poly sqr_clip_1;
  j2d_poly sqr_clip_2;
  j2d_poly sqr_split;
  j2d_plane rand_plane;
  for (int i = 0; i < NUM_TRIALS; ++i) {
    // copy polys
    j2d_copy_poly(&sqr_clip_1, &sqr);
    j2d_copy_poly(&sqr_clip_2, &sqr);
    j2d_copy_poly(&sqr_split, &sqr);
    // make random plane through centroid
    rand_plane = thru_cent_2d(&sqr);
    // split
    j2d_poly pos_sqr;
    j2d_poly neg_sqr;
    j2d_split(&sqr_split, 1, rand_plane, &pos_sqr, &neg_sqr);
    // clip
    j2d_clip(&sqr_clip_1, &rand_plane, 1);
    rand_plane.d *= -1.0;
    rand_plane.n.xy[0] *= -1.0;
    rand_plane.n.xy[1] *= -1.0;
    j2d_clip(&sqr_clip_2, &rand_plane, 1);
    // arrays for storing moments
    double moments_pos_sqr[NUM_MOMENTS];
    double moments_neg_sqr[NUM_MOMENTS];
    double moments_sqr_clip_1[NUM_MOMENTS];
    double moments_sqr_clip_2[NUM_MOMENTS];
    // we expect that the moments of pos_sqr == sqr_clip_1
    // we expect that the moments of neg_sqr == sqr_clip_2
    // ELEVENTH ORDER MOMENTS
    j2d_get_moments(&pos_sqr, moments_pos_sqr, POLY_ORDER);
    j2d_get_moments(&neg_sqr, moments_neg_sqr, POLY_ORDER);
    j2d_get_moments(&sqr_clip_1, moments_sqr_clip_1, POLY_ORDER);
    j2d_get_moments(&sqr_clip_2, moments_sqr_clip_2, POLY_ORDER);
    // check that moments of intersected poly on either side of plane
    // are equivalent
    for (int j = 0; j < NUM_MOMENTS; ++j) {
      ASSERT_DOUBLE_EQ(moments_pos_sqr[j], moments_sqr_clip_1[j], TOL);
      ASSERT_DOUBLE_EQ(moments_neg_sqr[j], moments_sqr_clip_2[j], TOL);
    }
  }
  printf("=== END SPLIT AND CLIP SQUARE TEST ===\n");
  return EXIT_SUCCESS;
}

int copy_box() {
  printf("=== BEGIN COPY BOX ===\n");

  j2d_poly box;
  j2d_vec2 vecs[2];

  vecs[0].xy[0] = 0.0;
  vecs[0].xy[1] = 0.0;
  vecs[1].xy[0] = 1.0;
  vecs[1].xy[1] = 1.0;

  j2d_init_box(&box, vecs);

  j2d_poly box_copy;
  j2d_copy_poly(&box_copy, &box);

  ASSERT_INT_EQ(box_copy.num_verts, 4);

  ASSERT_INT_EQ(box_copy.dir_edges[0], 1);
  ASSERT_INT_EQ(box_copy.dir_edges[1], 3);
  ASSERT_INT_EQ(box_copy.dir_edges[2], 2);
  ASSERT_INT_EQ(box_copy.dir_edges[3], 0);
  ASSERT_INT_EQ(box_copy.dir_edges[4], 3);
  ASSERT_INT_EQ(box_copy.dir_edges[5], 1);
  ASSERT_INT_EQ(box_copy.dir_edges[6], 0);
  ASSERT_INT_EQ(box_copy.dir_edges[7], 2);

  ASSERT_DOUBLE_EQ(box_copy.verts[0][0], 0.0, TOL);
  ASSERT_DOUBLE_EQ(box_copy.verts[0][1], 0.0, TOL);
  ASSERT_DOUBLE_EQ(box_copy.verts[1][0], 1.0, TOL);
  ASSERT_DOUBLE_EQ(box_copy.verts[1][1], 0.0, TOL);
  ASSERT_DOUBLE_EQ(box_copy.verts[2][0], 1.0, TOL);
  ASSERT_DOUBLE_EQ(box_copy.verts[2][1], 1.0, TOL);
  ASSERT_DOUBLE_EQ(box_copy.verts[3][0], 0.0, TOL);
  ASSERT_DOUBLE_EQ(box_copy.verts[3][1], 1.0, TOL);

  printf("=== END COPY BOX ===\n");
  return EXIT_SUCCESS;
}

#undef NUM_TESTS
#define NUM_TESTS 8

int main() {
  choptions_2d[0] = thru_cent_2d;
  choptions_2d[1] = thru_edge_2d;
  choptions_2d[2] = thru_vert_cent_2d;
  choptions_2d[3] = thru_vert_rand_2d;

  int tests[NUM_TESTS];

  tests[0] = split_tris_thru_centroid();
  tests[1] = recursive_splitting_nondegenerate();
  tests[2] = recursive_splitting_degenerate();
  tests[3] = recursive_splitting_degenerate_perturbed();
  tests[4] = random_verts();
  tests[5] = moments();
  tests[6] = split_and_clip_square();
  tests[7] = copy_box();

  for (int i = 0; i < NUM_TESTS; ++i) {
    if (tests[i] == EXIT_FAILURE) {
      return 1;
    }
  }

  return 0;
}

#undef POLY_ORDER
#undef NUM_MOMENTS
#undef NTRIS
#undef STACK_SIZE
#undef NVERTS
#undef NUM_TESTS
