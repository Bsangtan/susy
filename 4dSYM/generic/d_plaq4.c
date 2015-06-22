// -----------------------------------------------------------------
// Measure average space--space and space--time plaquettes
// Use tempmat1 for temporary storage
#include "generic_includes.h"

void d_plaquette(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  register su3_matrix_f *m1, *m4;
  double ss_sum = 0.0, st_sum = 0.0;
  msg_tag *mtag0, *mtag1;
  su3_matrix_f tmat;

  // We can exploit a symmetry under dir<-->dir2
  for (dir = YUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      // gen_pt[0] is U_b(x+a), gen_pt[1] is U_a(x+b)
      mtag0 = start_gather_site(F_OFFSET(linkf[dir2]), sizeof(su3_matrix_f),
                                goffset[dir], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[dir]), sizeof(su3_matrix_f),
                                goffset[dir2], EVENANDODD, gen_pt[1]);

      // tempmat1 = Udag_b(x) U_a(x)
      FORALLSITES(i, s) {
        m1 = &(s->linkf[dir]);
        m4 = &(s->linkf[dir2]);
        mult_su3_an_f(m4, m1, &(tempmat1)[i]);
      }
      wait_gather(mtag0);
      wait_gather(mtag1);

      // Compute tr[Udag_a(x+b) Udag_b(x) U_a(x) U_b(x+a)]
      if (dir == TUP || dir2 == TUP) {
        FORALLSITES(i, s) {
          m1 = (su3_matrix_f *)(gen_pt[0][i]);
          m4 = (su3_matrix_f *)(gen_pt[1][i]);
          mult_su3_nn_f(&(tempmat1[i]), m1, &tmat);
          st_sum += (double)realtrace_su3_f(m4, &tmat);
        }
      }
      else {
        FORALLSITES(i, s) {
          m1 = (su3_matrix_f *)(gen_pt[0][i]);
          m4 = (su3_matrix_f *)(gen_pt[1][i]);
          mult_su3_nn_f(&(tempmat1[i]), m1, &tmat);
          ss_sum += (double)realtrace_su3_f(m4, &tmat);
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
}
// -----------------------------------------------------------------
