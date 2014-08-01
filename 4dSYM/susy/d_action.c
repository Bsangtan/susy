// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge momentum contribution to the action
double d_hmom_action() {
  register int i, mu;
  register site *s;
  double sum = 0;

  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      sum += (double)realtrace_su3_f(&(s->mom[mu]), &(s->mom[mu]));
  }
  g_doublesum(&sum);
//  node0_printf("gauge momentum %e\n",sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Include tunable coefficient C2 in the d^2 term
double d_gauge_action() {
  register int i;
  register site *s;
  int mu, nu;
  double g_action = 0.0, norm = 0.5 * C2;
  complex cg_action;
  su3_matrix_f tmat;

  compute_DmuUmu();
  FORALLSITES(i, s) {
    mult_su3_nn_f(&(s->DmuUmu), &(s->DmuUmu), &tmat);
    cg_action = trace_su3_f(&tmat);
    g_action += norm * cg_action.real;
  }

  compute_Fmunu();
  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      FORALLSITES(i, s)
        g_action += 2 * realtrace_su3_f(&(s->Fmunu[mu][nu]),
                                        &(s->Fmunu[mu][nu]));
    }
  }
  g_action *= kappa;
  g_doublesum(&g_action);
  return g_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Mass term for U(1) mode -- note factor of kappa
double d_bmass_action() {
  register int i;
  register site *s;
  int mu;
  double sum = 0.0, dum;

  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s) {
      dum = 1.0 / (double)NCOL;
      dum *= realtrace_su3_f(&(s->linkf[mu]), &(s->linkf[mu]));
      dum -= 1.0;
      sum += kappa * bmass * bmass * dum * dum;
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion contribution to the action
// Ignore the ampdeg term -- since the pseudofermion src is fixed
// throughout the trajectory, it has no effect on Delta S (checked)
// sol, however, depends on the gauge fields through the CG
double d_fermion_action(Twist_Fermion *src, Twist_Fermion **sol) {
  register int i, j;
  register site *s;
  double sum = 0;
  complex ctmp;

#ifdef DEBUG_CHECK  // Check ampdeg4 term
  FORALLSITES(i, s)
    sum += ampdeg4 * (double)magsq_TF(&(src[i]));
  g_doublesum(&sum);
  node0_printf("ampdeg|chi|^2 = %.4g\n", sum);
#endif

  FORALLSITES(i, s) {
    for (j = 0; j < Norder; j++) {
      ctmp = TF_dot(&(src[i]), &(sol[j][i]));   // src^dag.sol[j]
      sum += (double)(amp4[j] * ctmp.real);
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant contribution to the action
double d_det_action() {
  register int i, dir1, dir2;
  register site *s;
  double det_action = 0.0;
  complex det;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;
  su3_matrix_f tmat;

  for (dir1 = YUP; dir1 < NUMLINK; dir1++) {
    for (dir2 = XUP; dir2 < dir1; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(linkf[dir2]), sizeof(su3_matrix_f),
                                goffset[dir1], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[dir1]), sizeof(su3_matrix_f),
                                goffset[dir2], EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s)
        mult_su3_an_f(&(s->linkf[dir2]), &(s->linkf[dir1]), &(s->tempmat1));

      wait_gather(mtag0);
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(s->tempmat1), (su3_matrix_f *)(gen_pt[0][i]),
                      &(s->staple));
      }
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mult_su3_na_f((su3_matrix_f *)(gen_pt[1][i]), &(s->staple), &tmat);
        det = find_det(&tmat);
        det_action += (double)((det.real - 1.0) * (det.real - 1.0));
        det_action += (double)(det.imag * det.imag);
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  det_action *= kappa_u1;
  g_doublesum(&det_action);
  return det_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA operators' contributions to the action
// Note factor of kappa on both
// B_a is now included in the site
void d_ops_action(double *K_act, double *S_act) {
  register int i;
  register site *s;
  int a, b;
  Real tS[NUMLINK][NUMLINK];
  complex ctmp;
  su3_matrix_f tmat;

  compute_Bmu();
  *K_act = 0.0;
  for (a = 0; a < NUMLINK; a++) {
    for (b = 0; b < NUMLINK; b++)
      tS[a][b] = 0.0;
  }

  // Sum operators over volume, averaging over all 25 SUGRA components
  // In measurements we check that ctmp is purely real; omit here for speed
  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++) {
      for (b = 0; b < NUMLINK; b++) {
        mult_su3_nn_f(&(s->B[a]), &(s->B[b]), &tmat);
        ctmp = trace_su3_f(&tmat);
        tS[a][b] += ctmp.real;
        if (a == b)
          *K_act += ctmp.real;          // Sum, not average
      }
    }
  }
  g_doublesum(K_act);
  for (a = 0; a < NUMLINK; a++) {
    for (b = 0; b < NUMLINK; b++)
      g_doublesum(&(tS[a][b]));
  }

  // SUGRA trace subtraction
  for (a = 0; a < NUMLINK; a++)
    tS[a][a] -= *K_act / (Real)NUMLINK;

  // For now just average over all 25 SUGRA components,
  // hit with couplings and return
  *S_act = 0.0;
  for (a = 0; a < NUMLINK; a++) {
    for (b = 0; b < NUMLINK; b++)
      *S_act += tS[a][b];
  }
  *K_act *= CK;
  *S_act *= CS / (Real)(NUMLINK * NUMLINK);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print out zeros if fermion and determinant actions not included
double d_action(Twist_Fermion *src, Twist_Fermion **sol) {
  double g_act, bmass_act, h_act, f_act = 0.0, det_act = 0.0;
  double S_act = 0.0, K_act = 0.0, total;
  g_act = d_gauge_action();
  bmass_act = d_bmass_action();
  h_act = d_hmom_action();
  det_act = d_det_action();

  // Only compute operator contributions if non-zero
  if (CK > 1.0e-8 || CS > 1.0e-8)
    d_ops_action(&K_act, &S_act);

#ifndef PUREGAUGE
  f_act = d_fermion_action(src, sol);
#endif
  node0_printf("action: ");
  node0_printf("gauge %.8g bmass %.8g ", g_act, bmass_act);
  node0_printf("det %.8g ", det_act);
  node0_printf("fermion %.8g hmom %.8g ", f_act, h_act);
  node0_printf("konishi %.8g ", K_act);
  node0_printf("sugra %.8g ", S_act);
  total = g_act + bmass_act + det_act + f_act + h_act + K_act + S_act;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
