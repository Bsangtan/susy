// -----------------------------------------------------------------
// Modified rectangular Wilson loops of fundamental links
// Evaluate in different spatial dirs to check rotational invariance
// Checked that results are gauge invariant
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Walk around path of fundamental links specified by dir, sign and kind
// dir lists the directions in the path
// sign lists whether to go forward (1) or backwards (-1)
// kind tells us whether to use linkf (1) or mom = (linkf^{-1})^dag (-1)
// length is the length of the path, and of each array
// Use tempmat1 to accumulate linkf product along path
// Use tempmat2 for temporary storage
void rsymm_path(int *dir, int *sign, int *kind, int length) {
  register int i;
  register site *s;
  int j;
  su3_matrix_f *mat;
  msg_tag *mtag = NULL;

  // Initialize tempmat1 with first link in path
  if (sign[0] > 0) {    // Gather from site - dir[0], no adjoint
    if (kind[0] > 0) {
      mtag = start_gather_site(F_OFFSET(linkf[dir[0]]), sizeof(su3_matrix_f),
                               goffset[dir[0]] + 1, EVENANDODD, gen_pt[0]);
    }
    else if (kind[0] < 0) {
      mtag = start_gather_site(F_OFFSET(mom[dir[0]]), sizeof(su3_matrix_f),
                               goffset[dir[0]] + 1, EVENANDODD, gen_pt[0]);
    }
    else {
      node0_printf("rsymm_path: unrecognized kind[0] = %d\n", kind[0]);
      terminate(1);
    }

    wait_gather(mtag);
    FORALLSITES(i, s)
      su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(tempmat1[i]));
    cleanup_gather(mtag);
  }

  else if (sign[0] < 0) {    // Take adjoint, no gather
    FORALLSITES(i, s) {
      if (kind[0] > 0)
        su3_adjoint_f(&(s->linkf[dir[0]]), &(tempmat1[i]));
      else if (kind[0] < 0)
        su3_adjoint_f(&(s->mom[dir[0]]), &(tempmat1[i]));
      else {
        node0_printf("rsymm_path: unrecognized kind[0] = %d\n", kind[0]);
        terminate(1);
      }
    }
  }
  else {
    node0_printf("rsymm_path: unrecognized sign[0] = %d\n", sign[0]);
    terminate(1);
  }

  // Accumulate subsequent links in product in tempmat1
  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {    // mult_su3_nn_f then gather from site - dir[j]
      FORALLSITES(i, s) {
        if (kind[j] > 0)
          mult_su3_nn_f(&(tempmat1[i]), &(s->linkf[dir[j]]), &(tempmat2[i]));
        else if (kind[j] < 0)
          mult_su3_nn_f(&(tempmat1[i]), &(s->mom[dir[j]]), &(tempmat2[i]));
        else {
          node0_printf("rsymm_path: unrecognized kind[%d] = %d\n", j, kind[j]);
          terminate(1);
        }
      }
      mtag = start_gather_field(tempmat2, sizeof(su3_matrix_f),
                                goffset[dir[j]] + 1, EVENANDODD, gen_pt[0]);

      wait_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(tempmat1[i]));
      cleanup_gather(mtag);
    }

    else if (sign[j] < 0) {    // Gather from site + dir[j] then mult_su3_na_f
      mtag = start_gather_field(tempmat1, sizeof(su3_matrix_f),
                                goffset[dir[j]], EVENANDODD, gen_pt[1]);

      wait_gather(mtag);
      FORALLSITES(i, s) {
          mat = (su3_matrix_f *)(gen_pt[1][i]);
        if (kind[j] > 0)
          mult_su3_na_f(mat, &(s->linkf[dir[j]]), &(tempmat2[i]));
        else if (kind[j] < 0)
          mult_su3_na_f(mat, &(s->mom[dir[j]]), &(tempmat2[i]));
        else {
          node0_printf("rsymm_path: unrecognized kind[%d] = %d\n", j, kind[j]);
          terminate(1);
        }
      }
      FORALLSITES(i, s)   // Don't want to overwrite tempmat1 too soon
        su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
      cleanup_gather(mtag);
    }
    else {
      node0_printf("rsymm_path: unrecognized sign[%d] = %d\n", j, sign[j]);
      terminate(1);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print both usual and transformed Wilson loops
// Use tempmat1 for temporary storage
void rsymm(int project) {
  register int i;
  register site *s;
  int dir_normal, dir_inv, dist, dist_inv, mu, length, max = MAX_X + 1;
  int dir[4 * max], sign[4 * max], kind[4 * max];
  double rsymm_loop, wloop, invlink[NUMLINK], invlink_sum = 0.0, td;
  double invlinkSq = 0.0;
  complex tc;
  su3_matrix_f tmat;

  node0_printf("rsymm: MAX = %d\n", max);

  // Optionally consider only hermitian part of links
  // Save links in f_U for future use (mom will be used for inverses)
  if (project == 1) {
    node0_printf("rsymm considering polar-projected links\n");
    FORALLSITES(i, s) {
      for (mu = XUP; mu < NUMLINK; mu++) {
        su3mat_copy_f(&(s->linkf[mu]), &(s->f_U[mu]));
        polar(&(s->f_U[mu]), &tmat, &(s->linkf[mu]));
      }
    }
  }

  // Compute and optionally check inverse matrices
  // Temporarily store the adjoint of the inverse in momentum matrices,
  // since it transforms like the original link
  for (mu = XUP; mu < NUMLINK; mu++) {
    FORALLSITES(i, s) {
      invert(&(s->linkf[mu]), &tmat);
      su3_adjoint_f(&tmat, &(s->mom[mu]));

#ifdef DEBUG_CHECK
#define INV_TOL 1e-12
#define INV_TOL_SQ 1e-24
      // Check inversion -- tmat should be unit matrix
      int j, k;
      mult_su3_nn_f(&(s->mom[mu]), &(s->linkf[mu]), &tmat);
      for (j = 0; j < NCOL; j++) {
        if (fabs(1 - tmat.e[j][j].real) > INV_TOL
            || fabs(tmat.e[j][j].imag) > INV_TOL) {
          printf("Link inversion fails on node%d:\n", this_node);
          dumpmat_f(&tmat);
        }
        for (k = j + 1; k < NCOL; k++) {
          if (cabs_sq(&(tmat.e[j][k])) > INV_TOL_SQ
              || cabs_sq(&(tmat.e[k][j])) > INV_TOL_SQ) {
            printf("Link inversion fails on node%d:\n", this_node);
            dumpmat_f(&tmat);
          }
        }
      }
      // Check left multiplication in addition to right multiplication
      mult_su3_nn_f(&(s->linkf[mu]), &(s->mom[mu]), &tmat);
      for (j = 0; j < NCOL; j++) {
        if (fabs(1 - tmat.e[j][j].real) > INV_TOL
            || fabs(tmat.e[j][j].imag) > INV_TOL) {
          printf("Link inversion fails on node%d:\n", this_node);
          dumpmat_f(&tmat);
        }
        for (k = j + 1; k < NCOL; k++) {
          if (cabs_sq(&(tmat.e[j][k])) > INV_TOL_SQ
              || cabs_sq(&(tmat.e[k][j])) > INV_TOL_SQ) {
            printf("Link inversion fails on node%d:\n", this_node);
            dumpmat_f(&tmat);
          }
        }
      }
#endif
    }
  }

  // First check average value of the inverted link
  // Tr[U^{-1} (U^{-1})^dag] / N and corresponding width
  // Just like link_trace() but use s->mom instead of s->linkf
  for (dir_inv = XUP; dir_inv < NUMLINK; dir_inv++) {
    invlink[dir_inv] = 0.0;
    FORALLSITES(i, s) {
      td = realtrace_su3_f(&(s->mom[dir_inv]), &(s->mom[dir_inv]));
      invlink[dir_inv] += td;
      invlinkSq += td * td;
    }
    invlink[dir_inv] /= ((double)volume * NCOL);
    g_doublesum(&(invlink[dir_inv]));
    invlink_sum += invlink[dir_inv];
  }
  invlink_sum /= ((double)NUMLINK);
  invlinkSq /= ((double)volume * NCOL * NCOL * NUMLINK);
  g_doublesum(&(invlinkSq));

  node0_printf("INVLINK");
  for (dir_inv = XUP; dir_inv < NUMLINK; dir_inv++)
    node0_printf(" %.6g", invlink[dir_inv]);
  td = sqrt(invlinkSq - invlink_sum * invlink_sum);
  node0_printf(" %.6g %.6g\n", invlink_sum, td);

  // Construct and print all loops up to max x max
  // in all NUMLINK * (NUMLINK - 1) directions
  // Invert all links in the second direction in each loop
  for (dir_normal = XUP; dir_normal < NUMLINK; dir_normal++) {
    for (dir_inv = XUP; dir_inv < NUMLINK; dir_inv++) {
      if (dir_inv == dir_normal)
        continue;

      for (dist = 1; dist <= max; dist++) {
        for (dist_inv = 1; dist_inv <= max; dist_inv++) {
          // Set up rectangular Wilson loop path as list of dir, sign * kind
          length = 2 * (dist + dist_inv);
          for (i = 0; i < dist; i++) {
            dir[i] = dir_normal;
            sign[i] = 1;
            kind[i] = 1;
          }
          for (i = dist; i < dist + dist_inv; i++) {
            dir[i] = dir_inv;
            sign[i] = 1;
            kind[i] = -1;
          }
          for (i = dist + dist_inv; i < 2 * dist + dist_inv; i++) {
            dir[i] = dir_normal;
            sign[i] = -1;
            kind[i] = 1;
          }
          for (i = 2 * dist + dist_inv; i < length; i++) {
            dir[i] = dir_inv;
            sign[i] = -1;
            kind[i] = -1;
          }
#ifdef DEBUG_CHECK
          node0_printf("path %d [%d] %d [%d] length %d: ",
                       dir_normal, dist, dir_inv, dist_inv, length);
          for (i = 0; i < length; i++)
            node0_printf(" (%d)*%d(%d) ", sign[i], dir[i], kind[i]);
          node0_printf("\n");
#endif

          // path and rsymm_path accumulate the product in tempmat1
          path(dir, sign, length);
          wloop = 0.0;
          FORALLSITES(i, s) {
            tc = trace_su3_f(&(tempmat1[i]));
            wloop += tc.real;
          }
          g_doublesum(&wloop);

          rsymm_path(dir, sign, kind, length);
          rsymm_loop = 0.0;
          FORALLSITES(i, s) {
            tc = trace_su3_f(&(tempmat1[i]));
            rsymm_loop += tc.real;
          }
          g_doublesum(&rsymm_loop);
          // Format: normal [dir] inverted [dir] usual transformed
          node0_printf("RSYMM %d [%d] %d [%d] %.8g %.8g\n",
                       dist, dir_normal, dist_inv, dir_inv,
                       wloop / volume, rsymm_loop / volume);
        } // dist_inv
      } // dist
    } // dir_inv
  } // dir_normal

  // Restore original links from f_U
  if (project == 1) {
    FORALLSITES(i, s) {
      for (mu = XUP; mu < NUMLINK; mu++)
        su3mat_copy_f(&(s->f_U[mu]), &(s->linkf[mu]));
    }
  }
}
// -----------------------------------------------------------------
