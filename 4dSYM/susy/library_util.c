// -----------------------------------------------------------------
// Mostly routines on individual Twist_Fermions,
// which could be moved into the libraries
// The last two are exceptions that loop over all sites
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void dump_TF(Twist_Fermion *in) {
  int mu;
  node0_printf("Fsite:   ");
  dumpmat(&(in->Fsite));
  FORALLDIR(mu) {
    node0_printf("Flink %d: ", mu);
    dumpmat(&(in->Flink[mu]));
  }
  for (mu = 0; mu < NPLAQ; mu++) {
    node0_printf("Fplaq %d: ", mu);
    dumpmat(&(in->Fplaq[mu]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a Twist_Fermion (hardly worth a function)
void copy_TF(Twist_Fermion *src, Twist_Fermion *dest) {
  *dest = *src;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clear a Twist_Fermion
void clear_TF(Twist_Fermion *in) {
  register int i;
  clear_mat(&(in->Fsite));
  FORALLDIR(i)
    clear_mat(&(in->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    clear_mat(&(in->Fplaq[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the squared magnitude of a Twist_Fermion, ReTr[adag.a]
Real magsq_TF(Twist_Fermion *in) {
  register int i;
  register Real sum;
  sum = realtrace(&(in->Fsite), &(in->Fsite));
  FORALLDIR(i)
    sum += realtrace(&(in->Flink[i]), &(in->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    sum += realtrace(&(in->Fplaq[i]), &(in->Fplaq[i]));
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the dot product of two Twist_Fermions, Tr[adag.b]
complex TF_dot(Twist_Fermion *a, Twist_Fermion *b) {
  register int i;
  complex sum, tc;
  sum = complextrace_an(&(a->Fsite), &(b->Fsite));
  FORALLDIR(i) {
    tc = complextrace_an(&(a->Flink[i]), &(b->Flink[i]));
    CSUM(sum, tc);
  }
  for (i = 0; i < NPLAQ; i++) {
    tc = complextrace_an(&(a->Fplaq[i]), &(b->Fplaq[i]));
    CSUM(sum, tc);
  }
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// c <-- c + s * b
void scalar_mult_sum_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_sum_matrix(&(b->Fsite), s, &(c->Fsite));
  FORALLDIR(i)
    scalar_mult_sum_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_sum_matrix(&(b->Fplaq[i]), s, &(c->Fplaq[i]));
}

// c <-- a + s * b
void scalar_mult_add_TF(Twist_Fermion *a, Twist_Fermion *b,
                        Real s, Twist_Fermion *c) {

  register int i;
  scalar_mult_add_matrix(&(a->Fsite), &(b->Fsite), s, &(c->Fsite));
  FORALLDIR(i)
    scalar_mult_add_matrix(&(a->Flink[i]), &(b->Flink[i]), s, &(c->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_add_matrix(&(a->Fplaq[i]), &(b->Fplaq[i]), s, &(c->Fplaq[i]));
}

// c <-- s1 * a + s2 * b
void scalar_mult_mult_add_TF(Twist_Fermion *a, Real s1, Twist_Fermion *b,
		Real s2, Twist_Fermion *c) {

	register int i;
	scalar_mult_mult_add_matrix(&(a->Fsite), s1, &(b->Fsite), s2, &(c->Fsite));
	FORALLDIR(i)
		scalar_mult_mult_add_matrix(&(a->Flink[i]), s1, &(b->Flink[i]), s2,
				&(c->Flink[i]));
	for (i = 0; i < NPLAQ; i++)
		scalar_mult_mult_add_matrix(&(a->Fplaq[i]), s1, &(b->Fplaq[i]), s2,
				&(c->Fplaq[i]));
}

void scalar_mult_TF(Twist_Fermion *src, Real s, Twist_Fermion *dest) {
  register int i;
  scalar_mult_matrix(&(src->Fsite), s, &(dest->Fsite));
  FORALLDIR(i)
    scalar_mult_matrix(&(src->Flink[i]), s, &(dest->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_matrix(&(src->Fplaq[i]), s, &(dest->Fplaq[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a gauge field as an array of NUMLINK matrices
void gauge_field_copy(field_offset src, field_offset dest) {
  register int i, dir, src2, dest2;
  register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    FORALLDIR(dir) {
      mat_copy((matrix *)F_PT(s, src2), (matrix *)F_PT(s, dest2));
      src2 += sizeof(matrix);
      dest2 += sizeof(matrix);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Shift a matrix without parallel transport
// The dir should come from goffset
void shiftmat(matrix *dat, matrix *temp, int dir) {
  register int i;
  register site *s;
  msg_tag *mtag;

  mtag = start_gather_field(dat, sizeof(matrix),
                            dir, EVENANDODD, gen_pt[0]);
  wait_gather(mtag);
  FORALLSITES(i, s)
    mat_copy((matrix *)gen_pt[0][i], &(temp[i]));
  cleanup_gather(mtag);
  FORALLSITES(i, s)
    mat_copy(&(temp[i]), &(dat[i]));
}
// -----------------------------------------------------------------
