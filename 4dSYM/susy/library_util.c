// -----------------------------------------------------------------
// Mostly routines on individual Twist_Fermions,
// which could be moved into the libraries
// The last two are exceptions that loop over all sites
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void conjTF(Twist_Fermion *src, Twist_Fermion *dest) {
  int mu, j;
  for (j = 0; j < DIMF; j++) {
    CONJG(src->Fsite.c[j], dest->Fsite.c[j]);
    for (mu = 0; mu < NUMLINK; mu++)
      CONJG(src->Flink[mu].c[j], dest->Flink[mu].c[j]);
    for (mu = 0; mu < NPLAQ; mu++)
      CONJG(src->Fplaq[mu].c[j], dest->Fplaq[mu].c[j]);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void dump_TF(Twist_Fermion *source) {
  int mu;
  node0_printf("Fsite:   ");
  dumpvec(&(source->Fsite));
  for (mu = 0; mu < NUMLINK; mu++) {
    node0_printf("Flink %d: ", mu);
    dumpvec(&(source->Flink[mu]));
  }
  for (mu = 0; mu < NPLAQ; mu++) {
    node0_printf("Fplaq %d: ", mu);
    dumpvec(&(source->Fplaq[mu]));
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
void clear_TF(Twist_Fermion *vec) {
  register int i;
  clearvec(&(vec->Fsite));
  for (i = 0; i < NUMLINK; i++)
    clearvec(&(vec->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    clearvec(&(vec->Fplaq[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the squared magnitude of a Twist_Fermion
Real magsq_TF(Twist_Fermion *vec) {
  register int i;
  register Real sum;
  sum = magsq_vec(&(vec->Fsite));
  for (i = 0; i < NUMLINK; i++)
    sum += magsq_vec(&(vec->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    sum += magsq_vec(&(vec->Fplaq[i]));
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the dot product of two Twist_Fermions, adag.b
complex TF_dot(Twist_Fermion *a, Twist_Fermion *b) {
  register int i;
  complex sum, tc;
  sum = dot(&(a->Fsite), &(b->Fsite));
  for (i = 0; i < NUMLINK; i++) {
    tc = dot(&(a->Flink[i]), &(b->Flink[i]));
    CSUM(sum, tc);
  }
  for (i = 0; i < NPLAQ; i++) {
    tc = dot(&(a->Fplaq[i]), &(b->Fplaq[i]));
    CSUM(sum, tc);
  }
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest <-- src1 + s * src2
void scalar_mult_add_TF(Twist_Fermion *src1, Twist_Fermion *src2,
                        Real s, Twist_Fermion *dest) {

  register int i;
  scalar_mult_add_vector(&(src1->Fsite), &(src2->Fsite), s, &(dest->Fsite));
  for (i = 0; i < NUMLINK; i++) {
    scalar_mult_add_vector(&(src1->Flink[i]), &(src2->Flink[i]), s,
                           &(dest->Flink[i]));
  }
  for (i = 0; i < NPLAQ; i++) {
    scalar_mult_add_vector(&(src1->Fplaq[i]), &(src2->Fplaq[i]), s,
                           &(dest->Fplaq[i]));

  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void scalar_mult_TF(Twist_Fermion *src, Real s, Twist_Fermion *dest) {
  register int i;
  scalar_mult_vector(&(src->Fsite), s, &(dest->Fsite));
  for (i = 0; i < NUMLINK; i++)
    scalar_mult_vector(&(src->Flink[i]), s, &(dest->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_vector(&(src->Fplaq[i]), s, &(dest->Fplaq[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a gauge field as an array of NUMLINK matrices
void gauge_field_copy_f(field_offset src, field_offset dest) {
  register int i, dir, src2, dest2;
  register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    FORALLDIR(dir) {
      mat_copy_f((matrix_f *)F_PT(s, src2), (matrix_f *)F_PT(s, dest2));
      src2 += sizeof(matrix_f);
      dest2 += sizeof(matrix_f);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Shift a matrix without parallel transport
// The dir should come from goffset
void shiftmat(matrix_f *dat, matrix_f *temp, int dir) {
  register int i;
  register site *s;
  msg_tag *mtag;

  mtag = start_gather_field(dat, sizeof(matrix_f),
                            dir, EVENANDODD, gen_pt[0]);
  wait_gather(mtag);
  FORALLSITES(i, s)
    mat_copy_f((matrix_f *)gen_pt[0][i], &(temp[i]));
  cleanup_gather(mtag);
  FORALLSITES(i, s)
    mat_copy_f(&(temp[i]), &(dat[i]));
}
// -----------------------------------------------------------------
