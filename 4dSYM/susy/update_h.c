// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
//#define FORCE_DEBUG
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update mom[NUMLINK] with the gauge force
// Include tunable coefficient C2 in the d^2 term of the action
// Use tr_dest and tempmat1 for temporary storage
double gauge_force(Real eps) {
  register int i, mu, nu, index;
  register site *s;
  register Real eb3, eb3UdU;
  double returnit = 0.0;
  complex tc, tc2;
  msg_tag *tag[NUMLINK], *tag0, *tag1;
  su3_matrix_f tmat1, tmat2, tmat3, UdU, *mat;

  // Two contributions from d^2 term
  compute_DmuUmu();   // Computes plaqdet if G is non-zero
  // If G is non-zero then below we will need G Tr[DmuUmu] plaqdet[mu][nu]
  // Save this in a modified plaqdet[mu][nu]
  if (G > IMAG_TOL && DIMF == NCOL * NCOL) {
    FORALLSITES(i, s) {
      tc = trace_su3_f(&DmuUmu[i]);
      for (mu = XUP; mu < NUMLINK; mu++) {
        for (nu = mu + 1; nu < NUMLINK; nu++) {
          CMUL(tc, plaqdet[mu][nu][i], tc2);
          plaqdet[mu][nu][i] = tc2;

          CMUL(tc, plaqdet[nu][mu][i], tc2);
          plaqdet[nu][mu][i] = tc2;
        }
      }
    }
  }

  // First we have the finite difference operator derivative times DmuUmu
  // Ubar_a(x) DmuUmu(x) - DmuUmu(x + a) Ubar_a(x)
  tag[0] = start_gather_field(DmuUmu, sizeof(su3_matrix_f),
                              goffset[0], EVENANDODD, gen_pt[0]);
  for (mu = XUP; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(DmuUmu, sizeof(su3_matrix_f),
                                       goffset[mu + 1], EVENANDODD,
                                       gen_pt[mu + 1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mat = (su3_matrix_f *)gen_pt[mu][i];
      mult_su3_an_f(&(s->linkf[mu]), &(DmuUmu[i]), &tmat1);
      mult_su3_na_f(mat, &(s->linkf[mu]), &tmat2);
      sub_su3_matrix_f(&tmat1, &tmat2, &tmat3);
      scalar_mult_su3_matrix_f(&tmat3, C2, &(s->f_U[mu]));    // Initialize
    }
    cleanup_gather(tag[mu]);
  }

  // Second we have DmuUmu times the plaquette determinant derivative
  // Sum_nu {Tr[DmuUmu(x)] plaqdet[nu][mu](x)
  //       + Tr[DmuUmu(x-nu)] plaqdet[mu][nu](x-nu)} U_mu^{-1}(x)
  // Only compute if G is non-zero
  // Use tr_dest and tempmat1 for temporary storage
  if (G > IMAG_TOL && DIMF == NCOL * NCOL) {
    for (mu = XUP; mu < NUMLINK; mu++) {
      // Invert U_a(x) and store in tempmat1
      // Zero tr_dest to hold sum
      FORALLSITES(i, s) {
        invert(&(s->linkf[mu]), &(tempmat1[i]));
        tr_dest[i] = cmplx(0.0, 0.0);
      }

      for (nu = XUP; nu < NUMLINK; nu++) {
        if (mu == nu)
          continue;

        // Gather Tr[DmuUmu] plaqdet[mu][nu] from x - nu
        tag0 = start_gather_field(plaqdet[mu][nu], sizeof(complex),
                                  goffset[nu] + 1, EVENANDODD, gen_pt[0]);

        // Add Tr[DmuUmu(x)] plaqdet[nu][mu](x) to sum while gather runs
        FORALLSITES(i, s)
          CSUM(tr_dest[i], plaqdet[nu][mu][i]);

        // Add Tr[DmuUmu(x-nu)] plaqdet[mu][nu](x-nu)} to sum
        wait_gather(tag0);
        FORALLSITES(i, s)
          CSUM(tr_dest[i], *((complex *)(gen_pt[0][i])));
        cleanup_gather(tag0);
      }

      // Now add to force
      FORALLSITES(i, s) {
        CMULREAL(tr_dest[i], G, tc);
        c_scalar_mult_add_su3mat_f(&(s->f_U[mu]), &(tempmat1[i]), &tc,
                                   &(s->f_U[mu]));
      }
    }
  }

  // Contribution from Fbar_{ab} F_{ab} term
  compute_Fmunu();
  for (mu = XUP; mu < NUMLINK; mu++) {
    for (nu = XUP; nu < NUMLINK; nu++) {
      if (mu == nu)
        continue;

      index = plaq_index[mu][nu];
      FORALLSITES(i, s) {
        if (mu > nu)
          scalar_mult_su3_matrix_f(&(Fmunu[index][i]), -1.0, &tmat2);
        else
          su3mat_copy_f(&(Fmunu[index][i]), &tmat2);
        mult_su3_an_f(&tmat2, &(s->linkf[nu]), &(tempmat2[i]));
      }

      tag0 = start_gather_site(F_OFFSET(linkf[nu]), sizeof(su3_matrix_f),
                               goffset[mu], EVENANDODD, gen_pt[0]);
      tag1 = start_gather_field(tempmat2, sizeof(su3_matrix_f),
                                goffset[nu] + 1, EVENANDODD, gen_pt[1]);

      // Set up tempmat1 while gathers run
      FORALLSITES(i, s) {
        if (mu > nu)
          scalar_mult_su3_matrix_f(&(Fmunu[index][i]), -1.0, &(tempmat1[i]));
        else
          su3mat_copy_f(&(Fmunu[index][i]), &(tempmat1[i]));
      }

      wait_gather(tag0);
      wait_gather(tag1);
      FORALLSITES(i, s) {
        mult_su3_na_f((su3_matrix_f *)gen_pt[0][i], &(tempmat1[i]), &tmat1);
        sub_su3_matrix_f(&tmat1, (su3_matrix_f *)gen_pt[1][i], &tmat2);
        scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmat2, 2.0, &(s->f_U[mu]));
      }
      cleanup_gather(tag0);
      cleanup_gather(tag1);
    }
  }

  // Factor of kappa = N / (2lambda) on both (D_a U_a)^2 and F^2 terms
  for (mu = XUP; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      scalar_mult_su3_matrix_f(&(s->f_U[mu]), kappa, &(s->f_U[mu]));
  }

  // Only compute U(1) mass term if non-zero -- note factor of kappa
  eb3 = 2.0 * kappa * bmass * bmass / (Real)(NCOL * NCOL);
  if (eb3 > IMAG_TOL) {
    for (mu = XUP; mu < NUMLINK; mu++) {
      FORALLSITES(i, s) {
        mult_su3_an_f(&(s->linkf[mu]), &(s->linkf[mu]), &UdU);
        tc = trace_su3_f(&UdU);
        eb3UdU = eb3 * (tc.real -(Real)NCOL);

        su3_adjoint_f(&(s->linkf[mu]), &tmat1);
        scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmat1, eb3UdU,
                                     &(s->f_U[mu]));
      }
    }
  }

  // Finally take adjoint and update the momentum
  // Subtract to reproduce -Adj(f_U)
  for (mu = XUP; mu < NUMLINK; mu++) {
    FORALLSITES(i, s) {
      su3_adjoint_f(&(s->f_U[mu]), &tmat2);
      scalar_mult_sub_su3_matrix_f(&(s->mom[mu]), &tmat2, eps, &(s->mom[mu]));
    }
  }

  // Compute average gauge force
  FORALLSITES(i, s) {
    for (mu = XUP; mu < NUMLINK; mu++)
      returnit += realtrace_su3_f(&(s->f_U[mu]), &(s->f_U[mu]));
  }
  g_doublesum(&returnit);

  // Add in force from separate determinant term
  returnit += det_force(eps);
  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion contribution to gauge link force,
//   f_U = Adj(Ms).D_U M(U, Ub).s - Adj[Adj(Ms).D_Ub M(U, Ub).s]
// "s" is sol while "Ms" is psol
// Copy these into persistent su3_vectors for easier gathering
// Use tempmat1, Tr_Uinv, tr_dest and Ddet[012] for temporary storage
void assemble_fermion_force(Twist_Fermion *sol, Twist_Fermion *psol) {
  register int i;
  register site *s;
  char **local_pt[2][4];
  int mu, nu, a, b, c, d, e, l, m, gather, next, flip = 0;
  int index, i_ab, i_de, j;
  Real permm, BC;
  complex tc, tc2, trace;
  msg_tag *mtag[NUMLINK], *tag0[2], *tag1[2], *tag2[2], *tag3[2];
  msg_tag *mtag0 = NULL, *mtag1 = NULL, *Trtag = NULL;
  su3_vector tvec, *vec, *vec0, *vec1, *vec2, *vec3;
  su3_matrix_f tmat, tmat2, *mat;

  // For gathering it is convenient to copy the input Twist_Fermions
  // into persistent site, link and plaquette su3_vectors
  // We can use "src" and "dest" vectors for this storage,
  // though we want to call them "sol" and "psol" for clarity
  su3_vector *site_sol = site_src, *site_psol = site_dest;
  su3_vector *link_sol[NUMLINK], *link_psol[NUMLINK];
  su3_vector *plaq_sol[NPLAQ], *plaq_psol[NPLAQ];
  for (mu = XUP; mu < NUMLINK; mu++) {
    link_sol[mu] = link_src[mu];
    link_psol[mu] = link_dest[mu];
  }
  for (mu = 0; mu < NPLAQ; mu++) {
    plaq_sol[mu] = plaq_src[mu];
    plaq_psol[mu] = plaq_dest[mu];
  }
  FORALLSITES(i, s) {
    su3vec_copy(&(sol[i].Fsite), &(site_sol[i]));
    su3vec_copy(&(psol[i].Fsite), &(site_psol[i]));
    for (mu = XUP; mu < NUMLINK; mu++) {
      su3vec_copy(&(sol[i].Flink[mu]), &(link_sol[mu][i]));
      su3vec_copy(&(psol[i].Flink[mu]), &(link_psol[mu][i]));
    }
    for (mu = 0; mu < NPLAQ; mu++) {
      su3vec_copy(&(sol[i].Fplaq[mu]), &(plaq_sol[mu][i]));
      su3vec_copy(&(psol[i].Fplaq[mu]), &(plaq_psol[mu][i]));
    }
  }

  // Clear the force collectors
  FORALLSITES(i, s) {
    for (mu = XUP; mu < NUMLINK; mu++)
      clear_su3mat_f(&(s->f_U[mu]));
  }

  // First calculate DU on chi_{munu} D_mu(U) psi_nu
#ifdef VP
  for (mu = XUP; mu < NUMLINK; mu++) {
    for (nu = XUP; nu < NUMLINK; nu++) {
      if (mu == nu)
        continue;

      index = plaq_index[mu][nu];
      mtag0 = start_gather_field(link_sol[nu], sizeof(su3_vector),
                                 goffset[mu], EVENANDODD, gen_pt[0]);
      wait_gather(mtag0);
      FORALLSITES(i, s) {
        vec = (su3_vector *)(gen_pt[0][i]);
        if (mu > nu) {    // plaq_psol is anti-symmetric under mu <--> nu
          scalar_mult_su3_vector(&(plaq_psol[index][i]), -1.0, &tvec);
        }                 // Suppress compiler error
        else
          su3vec_copy(&(plaq_psol[index][i]), &tvec);
        for (a = 0; a < DIMF; a++) {
          for (b = 0; b < DIMF; b++) {
            CMULJ_(tvec.c[a], vec->c[b], tc);
            CNEGATE(tc, tc);
            CMULREAL(tc, s->bc1[mu], tc);
            c_scalar_mult_add_su3mat_f(&(s->f_U[mu]), &(Lambda_prod[b][a]),
                                       &tc, &(s->f_U[mu]));
          }
        }
      }
      cleanup_gather(mtag0);

      FORALLSITES(i, s) {
        clear_su3mat_f(&(tempmat1[i]));
        if (mu > nu) {    // plaq_psol is anti-symmetric under mu <--> nu
          scalar_mult_su3_vector(&(plaq_psol[index][i]), -1.0, &tvec);
        }                 // Suppress compiler error
        else
          su3vec_copy(&(plaq_psol[index][i]), &tvec);
        for (a = 0; a < DIMF; a++) {
          for (b = 0; b < DIMF; b++) {
            CMULJ_(tvec.c[a], (link_sol[nu][i]).c[b], tc);
            c_scalar_mult_add_su3mat_f(&(tempmat1[i]), &(Lambda_prod[a][b]),
                                       &tc, &(tempmat1[i]));
          }
        }
      }
      mtag1 = start_gather_field(tempmat1, sizeof(su3_matrix_f),
                                 goffset[nu] + 1, EVENANDODD, gen_pt[1]);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mat = (su3_matrix_f *)(gen_pt[1][i]);
        add_su3_matrix_f(&(s->f_U[mu]), mat, &(s->f_U[mu]));
      }
      cleanup_gather(mtag1);
    }
  }

  // 2nd term
  for (mu = XUP; mu < NUMLINK; mu++) {
    for (nu = XUP; nu < NUMLINK; nu++) {
      if (mu == nu)
        continue;

      index = plaq_index[mu][nu];
      mtag0 = start_gather_field(link_psol[nu], sizeof(su3_vector),
                                 goffset[mu], EVENANDODD, gen_pt[0]);

      FORALLSITES(i, s) {
        clear_su3mat_f(&(tempmat1[i]));
        if (mu > nu) {    // plaq_sol is anti-symmetric under mu <--> nu
          scalar_mult_su3_vector(&(plaq_sol[index][i]), -1.0, &tvec);
        }                 // Suppress compiler error
        else
          su3vec_copy(&(plaq_sol[index][i]), &tvec);
        for (a = 0; a < DIMF; a++) {
          for (b = 0; b < DIMF; b++) {
            CMULJ_((link_psol[nu][i]).c[a], tvec.c[b], tc);
            c_scalar_mult_add_su3mat_f(&(tempmat1[i]), &(Lambda_prod[b][a]),
                                       &tc, &(tempmat1[i]));
          }
        }
      }
      mtag1 = start_gather_field(tempmat1, sizeof(su3_matrix_f),
                                 goffset[nu] + 1, EVENANDODD, gen_pt[1]);

      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mat = (su3_matrix_f *)(gen_pt[1][i]);
        sub_su3_matrix_f(&(s->f_U[mu]), mat, &(s->f_U[mu]));
      }
      cleanup_gather(mtag1);

      wait_gather(mtag0);
      FORALLSITES(i, s) {
        vec = (su3_vector *)(gen_pt[0][i]);
        if (mu > nu) {    // plaq_sol is anti-symmetric under mu <--> nu
          scalar_mult_su3_vector(&(plaq_sol[index][i]), -1.0, &tvec);
        }                 // Suppress compiler error
        else
          su3vec_copy(&(plaq_sol[index][i]), &tvec);
        for (a = 0; a < DIMF; a++) {
          for (b = 0; b < DIMF; b++) {
            CMULJ_(vec->c[a], tvec.c[b], tc);
            CMULREAL(tc, s->bc1[mu], tc);
            c_scalar_mult_add_su3mat_f(&(s->f_U[mu]), &(Lambda_prod[a][b]),
                                       &tc, &(s->f_U[mu]));
          }
        }
      }
      cleanup_gather(mtag0);
    }
  }
#endif
#ifdef SV
  // 3rd term, DUbar on eta Dbar_mu psi_mu (LtoS)
  mtag[0] = start_gather_field(site_psol, sizeof(su3_vector),
                               goffset[0], EVENANDODD, gen_pt[0]);
  for (mu = XUP; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(site_psol, sizeof(su3_vector),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      clear_su3mat_f(&tmat);
      vec = (su3_vector *)(gen_pt[mu][i]);    // site_psol(x + mu)
      for (a = 0; a < DIMF; a++) {
        for (b = 0; b < DIMF; b++) {
          CMULJ_((site_psol[i]).c[a], (link_sol[mu][i]).c[b], tc);
          CNEGATE(tc, tc);
          c_scalar_mult_add_su3mat_f(&tmat, &(Lambda_prod[a][b]), &tc, &tmat);
          CMULJ_((vec->c[a]), (link_sol[mu][i]).c[b], tc);
          CMULREAL(tc, s->bc1[mu], tc);
          c_scalar_mult_add_su3mat_f(&tmat, &(Lambda_prod[b][a]), &tc, &tmat);
        }
      }
      su3_adjoint_f(&tmat, &tmat2);
      scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmat2, 0.5, &(s->f_U[mu]));
    }
    cleanup_gather(mtag[mu]);
  }

  // 4th term, DUbar on psi_mu Dbar_mu eta
  mtag[0] = start_gather_field(site_sol, sizeof(su3_vector),
                               goffset[0], EVENANDODD, gen_pt[0]);
  for (mu = XUP; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(site_sol, sizeof(su3_vector),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      clear_su3mat_f(&tmat);
      vec = (su3_vector *)(gen_pt[mu][i]);
      for (a = 0; a < DIMF; a++) {
        for (b = 0; b < DIMF; b++) {
          CMULJ_((link_psol[mu][i]).c[a], vec->c[b], tc);
          CMULREAL(tc, s->bc1[mu], tc);
          CNEGATE(tc, tc);
          c_scalar_mult_add_su3mat_f(&tmat, &(Lambda_prod[a][b]), &tc, &tmat);
          CMULJ_((link_psol[mu][i]).c[a], (site_sol[i]).c[b], tc);
          c_scalar_mult_add_su3mat_f(&tmat, &(Lambda_prod[b][a]), &tc, &tmat);
        }
      }
      su3_adjoint_f(&tmat, &tmat2);
      scalar_mult_add_su3_matrix_f(&(s->f_U[mu]), &tmat2, 0.5, &(s->f_U[mu]));
    }
    cleanup_gather(mtag[mu]);
  }

  // Plaquette determinant contributions if G is non-zero
  // Use tempmat1, Tr_Uinv, tr_dest and Ddet[012] for temporary storage
  if (G > IMAG_TOL && DIMF == NCOL * NCOL) {
    complex *tempc = malloc(sites_on_node * sizeof(*tempc));
    complex *inv_term = malloc(sites_on_node * sizeof(*inv_term));
    complex *adj_term = malloc(sites_on_node * sizeof(*adj_term));
    complex detG = cmplx(0.0, -0.5 * G * sqrt((Real)NCOL));
    compute_plaqdet();

    // Save U_a(x)^{-1} and Udag_a(x)^{-1}
    // Use Ddet[0][a] and Ddet[1][a] for temporary storage, respectively
    su3_matrix_f *Uinv[NUMLINK], *Udag_inv[NUMLINK];
    for (a = XUP; a < NUMLINK; a++) {
      Uinv[a] = Ddet[0][a];
      Udag_inv[a] = Ddet[1][a];
      FORALLSITES(i, s) {
        invert(&(s->linkf[a]), &(Uinv[a][i]));
        su3_adjoint_f(&(s->linkf[a]), &tmat);
        invert(&tmat, &(Udag_inv[a][i]));
      }
    }

    // First connect link_sol with site_psol[DIMF - 1]^dag (LtoS)
    for (a = XUP; a < NUMLINK; a++) {
      // Save sum_j U_a(x)^{-1} Lambda^j psi_a^j(x) U_a(x)^{-1} in Ddet[2][a]
      // Save sum_j Tr[U_a(x)^{-1} Lambda^j] psi_a^j(x) in Tr_Uinv[a]
      FORALLSITES(i, s) {
        // Initialize accumulators for sum over j
        Tr_Uinv[a][i] = cmplx(0.0, 0.0);
        clear_su3mat_f(&(Ddet[2][a][i]));
        for (j = 0; j < DIMF; j++) {
          mult_su3_nn_f(&(Uinv[a][i]), &(Lambda[j]), &tmat);
          tc = trace_su3_f(&tmat);                // Accumulate trace
          CMUL(tc, (link_sol[a][i]).c[j], tc2);
          CSUM(Tr_Uinv[a][i], tc2);

          mult_su3_nn_f(&tmat, &(Uinv[a][i]), &tmat2);
          tc = (link_sol[a][i]).c[j];             // Accumulate product
          c_scalar_mult_add_su3mat_f(&(Ddet[2][a][i]), &tmat2, &tc,
                                     &(Ddet[2][a][i]));
        }

        // Save eta^D(x) plaqdet[a][b](x) in modified plaqdet[a][b]
        CONJG((site_psol[i]).c[DIMF - 1], tc);
        for (b = a + 1; b < NUMLINK; b++) {
          CMUL(tc, plaqdet[a][b][i], tc2);
          plaqdet[a][b][i] = tc2;

          CMUL(tc, plaqdet[b][a][i], tc2);
          plaqdet[b][a][i] = tc2;
        }
      }
    }
    for (a = XUP; a < NUMLINK; a++) {
      // Initialize accumulators for sums over b
      FORALLSITES(i, s) {
        tr_dest[i] = cmplx(0.0, 0.0);
        inv_term[i] = cmplx(0.0, 0.0);
        adj_term[i] = cmplx(0.0, 0.0);
      }
      for (b = XUP; b < NUMLINK; b++) {
        if (a == b)
          continue;

        // Gathers galore!    TODO: Can this be simplified?
        // 0) Tr_Uinv[b](x + a - b) in two steps  (Trtag)
        // 1) plaqdet[a][b](x - b)                (mtag0)
        // 2) plaqdet[b][a](x - b)                (mtag1)
        // 3) Tr_Uinv[b](x + a)                   (mtag[0])
        // 4) Tr_Uinv[b](x - b)                   (mtag[1])
        // 5) Tr_Uinv[a](x + b)                   (mtag[2])
        // 6) Tr_Uinv[a](x - b)                   (mtag[3])
        Trtag = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[a], EVENANDODD, gen_pt[0]);
        mtag0 = start_gather_field(plaqdet[a][b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[1]);
        mtag1 = start_gather_field(plaqdet[b][a], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[2]);
        mtag[0] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                     goffset[a], EVENANDODD, gen_pt[3]);
        mtag[1] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                     goffset[b] + 1, EVENANDODD, gen_pt[4]);
        mtag[2] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                     goffset[b], EVENANDODD, gen_pt[5]);
        mtag[3] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                     goffset[b] + 1, EVENANDODD, gen_pt[6]);

        // Stage two of Tr_Uinv[b](x + a - b) gather
        // Use tempc for temporary storage
        wait_gather(Trtag);
        FORALLSITES(i, s) {
          tc = *((complex *)(gen_pt[0][i]));
          CMULREAL(tc, s->bc1[a], tempc[i]);
        }
        cleanup_gather(Trtag);
        Trtag = start_gather_field(tempc, sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[0]);

        // Shorthand: D[a][b](x) is eta^D(x) plaqdet[a][b](x)
        //            T[a](x) is Tr[U_a(x)^{-1} psi_a(x)]
        wait_gather(mtag0);         // 1) D[a][b](x - b)
        wait_gather(mtag1);         // 2) D[b][a](x - b)
        wait_gather(mtag[0]);       // 3) T[b](x + a)
        wait_gather(mtag[1]);       // 4) T[b](x - b)
        wait_gather(mtag[2]);       // 5) T[a](x + b)
        wait_gather(mtag[3]);       // 6) T[a](x - b)
        wait_gather(Trtag);         // 0) T[b](x + a - b)
        FORALLSITES(i, s) {
          // Accumulate tr_dest
          // D[b][a](x) {T[a](x) + T[b](x + a)}
          tc = *((complex *)(gen_pt[3][i]));    // T[b](x + a)
          CMULREAL(tc, s->bc1[a], tc);
          CADD(Tr_Uinv[a][i], tc, tc2);
          CMUL(plaqdet[b][a][i], tc2, tc);
          CSUM(tr_dest[i], tc);

          // D[a][b](x - b) {T[a](x) + T[b](x - b)}
          CMULREAL(Tr_Uinv[a][i], s->bc1[OPP_LDIR(b)], trace);
          tc = *((complex *)(gen_pt[4][i]));    // T[b](x - b)
          CADD(trace, tc, tc2);
          tc = *((complex *)(gen_pt[1][i]));    // D[a][b](x - b)
          CMUL(tc, tc2, trace);
          CSUM(tr_dest[i], trace);

          // Accumulate adj_term
          // D[a][b](x) {T[b](x) + T[a](x + b)}
          tc = *((complex *)(gen_pt[5][i]));    // T[a](x + b)
          CMULREAL(tc, s->bc1[b], tc);
          CADD(Tr_Uinv[b][i], tc, tc2);
          CMUL(plaqdet[a][b][i], tc2, tc);
          CSUM(adj_term[i], tc);

          // D[b][a](x - b) {T[a](x - b) + T[b](x + a - b) bc1[a](x - b)}
          tc = *((complex *)(gen_pt[0][i]));    // T[b](x + a - b)
          tc2 = *((complex *)(gen_pt[6][i]));   // T[a](x - b)
          CADD(tc, tc2, trace);
          tc = *((complex *)(gen_pt[2][i]));    // D[b][a](x - b)
          CMUL(tc, trace, tc2);
          CSUM(adj_term[i], tc2);

          // Accumulate inv_term = sum_b D[b][a](x) + D[a][b](x - b)
          tc = *((complex *)(gen_pt[1][i]));    // D[a][b](x - b)
          CMULREAL(tc, s->bc1[OPP_LDIR(b)], tc);
          CADD(plaqdet[b][a][i], tc, tc2);
          CSUM(inv_term[i], tc2);
        }
        cleanup_gather(Trtag);
        cleanup_gather(mtag0);
        cleanup_gather(mtag1);
        cleanup_gather(mtag[0]);
        cleanup_gather(mtag[1]);
        cleanup_gather(mtag[2]);
        cleanup_gather(mtag[3]);
      }

      // Now add to force
      // Be sure detG coupling is included in the adjoint
      FORALLSITES(i, s) {
        // Initialize with tr_dest hitting U_a(x)^{-1}
        CMUL(tr_dest[i], detG, tc);
        c_scalar_mult_add_su3mat_f(&(s->f_U[a]), &(Uinv[a][i]), &tc,
                                   &(s->f_U[a]));

        // Subtract inv_term hitting U_a(x)^{-1} psi_a(x) U_a(x)^{-1}
        CMUL(inv_term[i], detG, tc);
        CNEGATE(tc, tc);
        c_scalar_mult_add_su3mat_f(&(s->f_U[a]), &(Ddet[2][a][i]), &tc,
                                   &(s->f_U[a]));

        // Finally add adj_term hitting Udag_a(x)^{-1} followed by adjoint
        CMUL(adj_term[i], detG, tc);
        c_scalar_mult_su3mat_f(&(Udag_inv[a][i]), &tc, &tmat);
        su3_adjoint_f(&tmat, &tmat2);
        add_su3_matrix_f(&(s->f_U[a]), &tmat2, &(s->f_U[a]));
      }
    }

    // Second connect site_sol[DIMF - 1] with link_psol^dag (StoL)
    compute_plaqdet();      // Reset plaqdet
    CNEGATE(detG, detG);    // As in action, anti-commute site and link
    for (a = XUP; a < NUMLINK; a++) {
      // Save sum_j U_a(x)^{-1} Lambda^j psi_a^j(x) U_a(x)^{-1} in Ddet[2][a]
      // Save sum_j Tr[U_a(x)^{-1} Lambda^j] psi_a^j(x) in Tr_Uinv[a]
      FORALLSITES(i, s) {
        // Initialize accumulators for sum over j
        Tr_Uinv[a][i] = cmplx(0.0, 0.0);
        clear_su3mat_f(&(Ddet[2][a][i]));
        for (j = 0; j < DIMF; j++) {
          mult_su3_nn_f(&(Uinv[a][i]), &(Lambda[j]), &tmat);
          tc = trace_su3_f(&tmat);                // Accumulate trace
          CMUL_J(tc, (link_psol[a][i]).c[j], tc2);
          CSUM(Tr_Uinv[a][i], tc2);

          mult_su3_nn_f(&tmat, &(Uinv[a][i]), &tmat2);
          CONJG((link_psol[a][i]).c[j], tc);      // Accumulate product
          c_scalar_mult_add_su3mat_f(&(Ddet[2][a][i]), &tmat2, &tc,
                                     &(Ddet[2][a][i]));
        }

        // Save eta^D(x) plaqdet[a][b](x) in modified plaqdet[a][b]
        tc = (site_sol[i]).c[DIMF - 1];
        for (b = a + 1; b < NUMLINK; b++) {
          CMUL(tc, plaqdet[a][b][i], tc2);
          plaqdet[a][b][i] = tc2;

          CMUL(tc, plaqdet[b][a][i], tc2);
          plaqdet[b][a][i] = tc2;
        }
      }
    }
    for (a = XUP; a < NUMLINK; a++) {
      // Initialize accumulators for sums over b
      FORALLSITES(i, s) {
        tr_dest[i] = cmplx(0.0, 0.0);
        inv_term[i] = cmplx(0.0, 0.0);
        adj_term[i] = cmplx(0.0, 0.0);
      }
      for (b = XUP; b < NUMLINK; b++) {
        if (a == b)
          continue;

        // Gathers again!  Structure unchanged, all data new
        // 0) Tr_Uinv[b](x + a - b) in two steps  (Trtag)
        // 1) plaqdet[a][b](x - b)                (mtag0)
        // 2) plaqdet[b][a](x - b)                (mtag1)
        // 3) Tr_Uinv[b](x + a)                   (mtag[0])
        // 4) Tr_Uinv[b](x - b)                   (mtag[1])
        // 5) Tr_Uinv[a](x + b)                   (mtag[2])
        // 6) Tr_Uinv[a](x - b)                   (mtag[3])
        Trtag = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[a], EVENANDODD, gen_pt[0]);
        mtag0 = start_gather_field(plaqdet[a][b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[1]);
        mtag1 = start_gather_field(plaqdet[b][a], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[2]);
        mtag[0] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                     goffset[a], EVENANDODD, gen_pt[3]);
        mtag[1] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                     goffset[b] + 1, EVENANDODD, gen_pt[4]);
        mtag[2] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                     goffset[b], EVENANDODD, gen_pt[5]);
        mtag[3] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                     goffset[b] + 1, EVENANDODD, gen_pt[6]);

        // Stage two of Tr_Uinv[b](x + a - b) gather
        // Use tempc for temporary storage
        wait_gather(Trtag);
        FORALLSITES(i, s) {
          tc = *((complex *)(gen_pt[0][i]));
          CMULREAL(tc, s->bc1[a], tempc[i]);
        }
        cleanup_gather(Trtag);
        Trtag = start_gather_field(tempc, sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[0]);

        // Shorthand: D[a][b](x) is eta^D(x) plaqdet[a][b](x)
        //            T[a](x) is Tr[U_a(x)^{-1} psi_a(x)]
        wait_gather(mtag0);         // 1) D[a][b](x - b)
        wait_gather(mtag1);         // 2) D[b][a](x - b)
        wait_gather(mtag[0]);       // 3) T[b](x + a)
        wait_gather(mtag[1]);       // 4) T[b](x - b)
        wait_gather(mtag[2]);       // 5) T[a](x + b)
        wait_gather(mtag[3]);       // 6) T[a](x - b)
        wait_gather(Trtag);         // 0) T[b](x + a - b)
        FORALLSITES(i, s) {
          // Accumulate tr_dest
          // D[b][a](x) {T[a](x) + T[b](x + a)}
          tc = *((complex *)(gen_pt[3][i]));    // T[b](x + a)
          CMULREAL(tc, s->bc1[a], tc);
          CADD(Tr_Uinv[a][i], tc, tc2);
          CMUL(plaqdet[b][a][i], tc2, tc);
          CSUM(tr_dest[i], tc);

          // D[a][b](x - b) {T[a](x) + T[b](x - b)}
          CMULREAL(Tr_Uinv[a][i], s->bc1[OPP_LDIR(b)], trace);
          tc = *((complex *)(gen_pt[4][i]));    // T[b](x - b)
          CADD(trace, tc, tc2);
          tc = *((complex *)(gen_pt[1][i]));    // D[a][b](x - b)
          CMUL(tc, tc2, trace);
          CSUM(tr_dest[i], trace);

          // Accumulate adj_term
          // D[a][b](x) {T[b](x) + T[a](x + b)}
          tc = *((complex *)(gen_pt[5][i]));    // T[a](x + b)
          CMULREAL(tc, s->bc1[b], tc);
          CADD(Tr_Uinv[b][i], tc, tc2);
          CMUL(plaqdet[a][b][i], tc2, tc);
          CSUM(adj_term[i], tc);

          // D[b][a](x - b) {T[a](x - b) + T[b](x + a - b) bc1[a](x - b)}
          tc = *((complex *)(gen_pt[0][i]));    // T[b](x + a - b)
          tc2 = *((complex *)(gen_pt[6][i]));   // T[a](x - b)
          CADD(tc, tc2, trace);
          tc = *((complex *)(gen_pt[2][i]));    // D[b][a](x - b)
          CMUL(tc, trace, tc2);
          CSUM(adj_term[i], tc2);

          // Accumulate inv_term = sum_b D[b][a](x) + D[a][b](x - b)
          tc = *((complex *)(gen_pt[1][i]));    // D[a][b](x - b)
          CMULREAL(tc, s->bc1[OPP_LDIR(b)], tc);
          CADD(plaqdet[b][a][i], tc, tc2);
          CSUM(inv_term[i], tc2);
        }
        cleanup_gather(Trtag);
        cleanup_gather(mtag0);
        cleanup_gather(mtag1);
        cleanup_gather(mtag[0]);
        cleanup_gather(mtag[1]);
        cleanup_gather(mtag[2]);
        cleanup_gather(mtag[3]);
      }

      // Now add to force
      // Be sure detG coupling is included in the adjoint
      FORALLSITES(i, s) {
        // Initialize with tr_dest hitting U_a(x)^{-1}
        CMUL(tr_dest[i], detG, tc);
        c_scalar_mult_add_su3mat_f(&(s->f_U[a]), &(Uinv[a][i]), &tc,
                                   &(s->f_U[a]));

        // Subtract inv_term hitting U_a(x)^{-1} psi_a(x) U_a(x)^{-1}
        CMUL(inv_term[i], detG, tc);
        CNEGATE(tc, tc);
        c_scalar_mult_add_su3mat_f(&(s->f_U[a]), &(Ddet[2][a][i]), &tc,
                                   &(s->f_U[a]));

        // Finally add adj_term hitting Udag_a(x)^{-1} followed by adjoint
        CMUL(adj_term[i], detG, tc);
        c_scalar_mult_su3mat_f(&(Udag_inv[a][i]), &tc, &tmat);
        su3_adjoint_f(&tmat, &tmat2);
        add_su3_matrix_f(&(s->f_U[a]), &tmat2, &(s->f_U[a]));
      }
    }
    free(tempc);
    free(inv_term);
    free(adj_term);
  }
#endif
#ifdef QCLOSED
  if (NUMLINK != 5) {
    node0_printf("ERROR: NUMLINK IS %d != 5\n", NUMLINK);
    terminate(1);
  }

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }

  // First Q-closed piece: chi_ab D_c chi_de epsilon_{abcde}
  // From setup_lamba.c, we see b > a and e > d
  // Start first set of gathers
  a = FQ_lookup[0][0];
  b = FQ_lookup[0][1];
  c = FQ_lookup[0][2];
  d = FQ_lookup[0][3];
  e = FQ_lookup[0][4];
  i_ab = plaq_index[a][b];
  i_de = plaq_index[d][e];

  tag0[0] = start_gather_field(plaq_psol[i_de], sizeof(su3_vector),
                               F1Q_d2[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(plaq_sol[i_ab], sizeof(su3_vector),
                               goffset[c], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(plaq_psol[i_de], sizeof(su3_vector),
                               goffset[c], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_field(plaq_sol[i_ab], sizeof(su3_vector),
                               F1Q_d1[0], EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {         // Start next set of gathers
      next = j + 1;
      a = FQ_lookup[next][0];
      b = FQ_lookup[next][1];
      c = FQ_lookup[next][2];
      d = FQ_lookup[next][3];
      e = FQ_lookup[next][4];
      i_ab = plaq_index[a][b];
      i_de = plaq_index[d][e];

      tag0[gather] = start_gather_field(plaq_psol[i_de], sizeof(su3_vector),
                                        F1Q_d2[next], EVENANDODD,
                                        local_pt[gather][0]);
      tag1[gather] = start_gather_field(plaq_sol[i_ab], sizeof(su3_vector),
                                        goffset[c], EVENANDODD,
                                        local_pt[gather][1]);
      tag2[gather] = start_gather_field(plaq_psol[i_de], sizeof(su3_vector),
                                        goffset[c], EVENANDODD,
                                        local_pt[gather][2]);
      tag3[gather] = start_gather_field(plaq_sol[i_ab], sizeof(su3_vector),
                                        F1Q_d1[next], EVENANDODD,
                                        local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = FQ_lookup[j][0];
    b = FQ_lookup[j][1];
    c = FQ_lookup[j][2];
    d = FQ_lookup[j][3];
    e = FQ_lookup[j][4];
    permm = perm[d][e][c][a][b];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      clear_su3mat_f(&(tempmat1[i]));
      vec0 = (su3_vector *)(local_pt[flip][0][i]);
      vec1 = (su3_vector *)(local_pt[flip][1][i]);
      vec2 = (su3_vector *)(local_pt[flip][2][i]);
      vec3 = (su3_vector *)(local_pt[flip][3][i]);
      for (l = 0; l < DIMF; l++) {
        for (m = 0; m < DIMF; m++) {
          CMULJ_(vec0->c[l], vec1->c[m], tc);
          CMULREAL(tc, permm, tc);
          BC = (s->bc3[a][b][c]) * (s->bc1[c]);
          CMULREAL(tc, BC, tc);
          c_scalar_mult_add_su3mat_f(&(tempmat1[i]), &(Lambda_prod[l][m]),
                                     &tc, &(tempmat1[i]));

          CMULJ_(vec2->c[l], vec3->c[m], tc);
          CMULREAL(tc, -permm, tc);
          BC = (s->bc2[OPP_LDIR(a)][OPP_LDIR(b)]) * (s->bc1[c]);
          CMULREAL(tc, BC, tc);
          c_scalar_mult_add_su3mat_f(&(tempmat1[i]), &(Lambda_prod[m][l]),
                                     &tc, &(tempmat1[i]));
        }
      }
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = (flip + 1) % 2;

    FORALLSITES(i, s) {
      su3_adjoint_f(&(tempmat1[i]), &tmat);
      scalar_mult_add_su3_matrix_f(&(s->f_U[c]), &tmat, -0.5, &(s->f_U[c]));
    }
  }

  // Second Q-closed piece
  // From setup_lamba.c, we see b > a and e > d
  // Start first set of gathers
  flip = 0;                       // Reset
  a = FQ_lookup[0][0];
  b = FQ_lookup[0][1];
  c = FQ_lookup[0][2];
  d = FQ_lookup[0][3];
  e = FQ_lookup[0][4];
  i_ab = plaq_index[a][b];
  i_de = plaq_index[d][e];

  tag0[0] = start_gather_field(plaq_psol[i_ab], sizeof(su3_vector),
                               F2Q_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(plaq_sol[i_de], sizeof(su3_vector),
                               goffset[c], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(plaq_psol[i_ab], sizeof(su3_vector),
                               goffset[c], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_field(plaq_sol[i_de], sizeof(su3_vector),
                               F2Q_d2[0], EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {         // Start next set of gathers
      next = j + 1;
      a = FQ_lookup[next][0];
      b = FQ_lookup[next][1];
      c = FQ_lookup[next][2];
      d = FQ_lookup[next][3];
      e = FQ_lookup[next][4];
      i_ab = plaq_index[a][b];
      i_de = plaq_index[d][e];

      tag0[gather] = start_gather_field(plaq_psol[i_ab], sizeof(su3_vector),
                                        F2Q_d1[next], EVENANDODD,
                                        local_pt[gather][0]);
      tag1[gather] = start_gather_field(plaq_sol[i_de], sizeof(su3_vector),
                                        goffset[c], EVENANDODD,
                                        local_pt[gather][1]);
      tag2[gather] = start_gather_field(plaq_psol[i_ab], sizeof(su3_vector),
                                        goffset[c], EVENANDODD,
                                        local_pt[gather][2]);
      tag3[gather] = start_gather_field(plaq_sol[i_de], sizeof(su3_vector),
                                        F2Q_d2[next], EVENANDODD,
                                        local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = FQ_lookup[j][0];
    b = FQ_lookup[j][1];
    c = FQ_lookup[j][2];
    d = FQ_lookup[j][3];
    e = FQ_lookup[j][4];
    permm = perm[a][b][c][d][e];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      clear_su3mat_f(&(tempmat1[i]));
      vec0 = (su3_vector *)(local_pt[flip][0][i]);
      vec1 = (su3_vector *)(local_pt[flip][1][i]);
      vec2 = (su3_vector *)(local_pt[flip][2][i]);
      vec3 = (su3_vector *)(local_pt[flip][3][i]);
      for (l = 0; l < DIMF; l++) {
        for (m = 0; m < DIMF; m++) {
          CMULJ_(vec0->c[l], vec1->c[m], tc);
          CMULREAL(tc, permm, tc);
          BC = (s->bc2[OPP_LDIR(a)][OPP_LDIR(b)]) * (s->bc1[c]);
          CMULREAL(tc, BC, tc);
          c_scalar_mult_add_su3mat_f(&(tempmat1[i]), &(Lambda_prod[l][m]),
                                     &tc, &(tempmat1[i]));

          CMULJ_(vec2->c[l], vec3->c[m], tc);
          CMULREAL(tc, -permm, tc);
          BC = (s->bc3[a][b][c]) * (s->bc1[c]);
          CMULREAL(tc, BC, tc);
          c_scalar_mult_add_su3mat_f(&(tempmat1[i]), &(Lambda_prod[m][l]),
                                     &tc, &(tempmat1[i]));
        }
      }
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = (flip + 1) % 2;

    FORALLSITES(i, s) {
      su3_adjoint_f(&(tempmat1[i]), &tmat);
      scalar_mult_add_su3_matrix_f(&(s->f_U[c]), &tmat, -0.5, &(s->f_U[c]));
    }
  }
#endif

  // Final adjoint and minus sign
  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      scalar_mult_su3_matrix_f(&(s->f_U[mu]), -1.0, &tmat);
      su3_adjoint_f(&tmat, &(s->f_U[mu]));
    }
  }
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Assume that the multiCG has been run, with the answers in sol[j]
// Accumulate f_U for each pole into fullforce, add to momenta
// Use fullforce-->Fmunu and tempTF for temporary storage
double fermion_force(Real eps, Twist_Fermion *src, Twist_Fermion **sol) {
  register int i;
  register site *s;
  int mu, n;
  double returnit = 0.0;
  su3_matrix_f **fullforce = malloc(NUMLINK * sizeof(*fullforce));

#ifdef FORCE_DEBUG
  int kick, ii, jj, iters = 0;
  Real final_rsq;
  double individ_force, old_action, new_action = 0.0;
  su3_matrix_f tmat1, tprint, tprint2;
  clear_su3mat_f(&tprint);
  clear_su3mat_f(&tmat1);
#endif

  for (mu = XUP; mu < NUMLINK; mu++) {
    fullforce[mu] = Fmunu[mu];    // Use Fmunu for temporary storage
    FORALLSITES(i, s)
      clear_su3mat_f(&(fullforce[mu][i]));
  }

  for (n = 0; n < Norder; n++) {
    fermion_op(sol[n], tempTF, PLUS);
    // Makes sense to multiply here by amp4[n]...
    FORALLSITES(i, s)
      scalar_mult_TF(&(tempTF[i]), amp4[n], &(tempTF[i]));

    assemble_fermion_force(sol[n], tempTF);
#ifdef FORCE_DEBUG
    individ_force = 0.0;
#endif
    for (mu = XUP; mu < NUMLINK; mu++) {
      FORALLSITES(i, s) {
        add_su3_matrix_f(&(fullforce[mu][i]), &(s->f_U[mu]),
                         &(fullforce[mu][i]));
#ifdef FORCE_DEBUG
//      if (s->x == 0 && s->y == 0 && s->z == 0 && s->t == 0 && mu == 3) {
//        printf("Fermion force mu=%d on site (%d, %d, %d, %d)\n",
//               mu, s->x, s->y, s->z ,s->t);
//        dumpmat_f(&(s->f_U[mu]));
//      }
      // Compute average gauge force
      individ_force += realtrace_su3_f(&(s->f_U[mu]), &(s->f_U[mu]));
#endif
      }
    }
#ifdef FORCE_DEBUG
    g_doublesum(&individ_force);
    node0_printf("Individ_force %d %.4g\n",
                 n, eps * sqrt(individ_force) / volume);

    // Check that force syncs with fermion action
    old_action = d_fermion_action(src, sol);
    fermion_rep();
    iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
    new_action = d_fermion_action(src, sol);
    node0_printf("EXITING  %.4g\n", new_action - old_action);
    if (fabs(new_action - old_action) > 1e-3)
      terminate(1);                             // Don't go further for now

#if 0
    // Do a scan of the fermion action
    for (mu = XUP; mu < NUMLINK; mu++) {
      FORALLSITES(i, s) {
        node0_printf("mu=%d on site (%d, %d, %d, %d)\n",
                     mu, s->x, s->y, s->z, s->t);
        tmat1 = s->linkf[mu];
        dumpmat_f(&(s->f_U[mu]));

        for (ii = 0; ii < NCOL; ii++) {
          for (jj = 0; jj < NCOL; jj++) {
            for (kick = -1; kick <= 1; kick += 2) {
              s->linkf[mu] = tmat1;
              s->linkf[mu].e[ii][jj].real += 0.001 * (Real)kick;

              fermion_rep();
              iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
              if (kick == -1)
                new_action -= d_fermion_action(src, sol);
              if (kick == 1) {
                new_action += d_fermion_action(src, sol);
                tprint.e[ii][jj].real = -250.0 * new_action;
              }
            }

            for (kick = -1; kick <= 1; kick += 2) {
              s->linkf[mu] = tmat1;
              s->linkf[mu].e[ii][jj].imag += 0.001 * (Real)kick;

              fermion_rep();
              iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
              if (kick == -1)
                new_action -= d_fermion_action(src, sol);
              if (kick == 1) {
                new_action += d_fermion_action(src, sol);
                node0_printf("XXXG%d%dI %.4g %.4g\n",
                             ii, jj, 0.001 * (Real)kick, 500 * new_action);
                tprint.e[ii][jj].imag = -250 * new_action;
              }
            }
          }
        }
        sub_su3_matrix_f(&tprint, &(s->f_U[mu]), &tprint2);
        node0_printf("mu=%d on site (%d, %d, %d, %d): %.4g\n",
                     mu, s->x, s->y, s->z, s->t,
                     realtrace_su3_f(&tprint2, &tprint2));
        dumpmat_f(&tprint);
        s->linkf[mu] = tmat1;

        fermion_rep();
        iters += congrad_multi_field(src, sol, niter, rsqmin, &final_rsq);
      }
    }   // End scan of the fermion action
#endif
#endif
  }

  // Update the momentum from the fermion force -- sum or eps
  // Opposite sign as to gauge force,
  // because dS_G / dU = 2F_g while ds_F / dU = -2F_f
  FORALLSITES(i, s) {
    for (mu = XUP; mu < NUMLINK; mu++) {
      scalar_mult_add_su3_matrix_f(&(s->mom[mu]), &(fullforce[mu][i]), eps,
                                   &(s->mom[mu]));
      returnit += realtrace_su3_f(&(fullforce[mu][i]), &(fullforce[mu][i]));
    }
  }
  g_doublesum(&returnit);

  free(fullforce);
  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------
