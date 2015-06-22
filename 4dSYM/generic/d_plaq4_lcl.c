// -----------------------------------------------------------------
// Measure average space--space and space--time plaquettes,
// Use tempmat1 for temporary storage

// #define LOCAL_PLAQ turns on plaquette measurement on hypersurfaces
// plaq_prll[xx] gives the average of the in a fixed MY_X hypersurface,
// as a function of the MY_X coordinate
// Similarly the other plaquettes are given by plaq_perp[xx]

// #define MIN_PLAQ turns on measurement of minimum plaquette per config
// for tuning nHYP as in Hasenfratz & Knechtli, hep-lat/0103029

// #define PLAQ_DIST prints out all plaquettes for plotting distribution
// CAUTION: Do not run PLAQ_DIST with MPI!

#define MIN_PLAQ
#ifdef LOCAL_PLAQ
#define MY_X x
#define MY_DIR XUP          // Printed out first time the routine is called
#define MY_N nx
static int print_dir = 0;   // Controls printing out MY_DIR
#endif

//#define PLAQ_DIST
#include "generic_includes.h"

void d_plaquette_lcl(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  register su3_matrix_f *m1, *m4;
  double ss_sum = 0.0, st_sum = 0.0, cur_plaq;
#ifdef MIN_PLAQ
  double min_plaq = 200.0 * NCOL;
#endif
  msg_tag *mtag0, *mtag1;
  su3_matrix_f mtmp;

#ifdef LOCAL_PLAQ
  int xx;
  double *plaq_perp = malloc(MY_N * sizeof(*plaq_perp));
  double *plaq_prll = malloc(MY_N * sizeof(*plaq_prll));

  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] = 0.0;
    plaq_prll[xx] = 0.0;
  }
#endif

#ifdef PLAQ_DIST
  if (this_node != 0) {
    printf("plaquette: don't run PLAQ_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  // We can exploit a symmetry under dir<-->dir2
  for (dir = YUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(linkf[dir2]), sizeof(su3_matrix_f),
                                goffset[dir], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[dir]), sizeof(su3_matrix_f),
                                goffset[dir2], EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->linkf[dir]);
        m4 = &(s->linkf[dir2]);
        mult_su3_an_f(m4, m1, &(tempmat1[i]));
      }
      wait_gather(mtag0);
      wait_gather(mtag1);

      if (dir == TUP || dir2 == TUP) {
        FORALLSITES(i, s) {
          m1 = (su3_matrix_f *)(gen_pt[0][i]);
          m4 = (su3_matrix_f *)(gen_pt[1][i]);
          mult_su3_nn_f(&(tempmat1[i]), m1, &mtmp);
          cur_plaq = (double)realtrace_su3_f(m4, &mtmp);
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
#ifdef PLAQ_DIST
          printf("PLAQ_DIST %d %d %d %d %d %d %.4g\n",
                 s->x, s->y, s->z, s->t, dir, dir2, cur_plaq);
#endif
          st_sum += cur_plaq;
#ifdef LOCAL_PLAQ
          if (dir == MY_DIR || dir2 == MY_DIR)
            plaq_perp[s->MY_X] += cur_plaq;
          else
            plaq_prll[s->MY_X] += cur_plaq;
#endif
        }
      }
      else {
        FORALLSITES(i, s) {
          m1 = (su3_matrix_f *)(gen_pt[0][i]);
          m4 = (su3_matrix_f *)(gen_pt[1][i]);
          mult_su3_nn_f(&(tempmat1[i]), m1, &mtmp);
          cur_plaq = (double)realtrace_su3_f(m4, &mtmp);
#ifdef MIN_PLAQ
          if (cur_plaq < min_plaq)
            min_plaq = cur_plaq;
#endif
#ifdef PLAQ_DIST
          printf("PLAQ_DIST %d %d %d %d %d %d %.4g\n",
                 s->x, s->y, s->z, s->t, dir, dir2, cur_plaq);
#endif
          ss_sum += cur_plaq;
#ifdef LOCAL_PLAQ
          if (dir == MY_DIR || dir2 == MY_DIR)
            plaq_perp[s->MY_X] += cur_plaq;
          else
            plaq_prll[s->MY_X] += cur_plaq;
#endif
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

  // Average over four plaquettes that involve the temporal link
  // and six that do not
  *ss_plaq = ss_sum / ((double)(6.0 * volume));
  *st_plaq = st_sum / ((double)(4.0 * volume));

#ifdef LOCAL_PLAQ
  for (xx = 0; xx < MY_N; xx++) {
    g_doublesum(&(plaq_perp[xx]));
    g_doublesum(&(plaq_prll[xx]));
  }
  // Normalization
  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] *= ((double)(MY_N)) / (3.0 * volume):
    plaq_prll[xx] *= ((double)(MY_N)) / (3.0 * volume);
  }

  // Print out
  if (this_node == 0) {
    if (print_dir == 0) {
      printf("LOCAL_PLAQ [0=XUP,..., 3=TUP] dir=%d\n", MY_DIR);
      print_dir = 1;
    }
    printf("THIN_PLAQ_PERP");
    for (xx = 0; xx < MY_N; xx++)
      printf(" %e", (double)plaq_perp[xx]);
    printf("\n");
    printf("THIN_PLAQ_PRLL");
    for (xx = 0; xx < MY_N; xx++)
      printf(" %e", (double)plaq_prll[xx]);
    printf("\n");
  }
  free(plaq_perp);
  free(plaq_prll);
#endif

#ifdef MIN_PLAQ
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("MIN_PLAQ %.8g", min_plaq);
#endif
}
// -----------------------------------------------------------------
