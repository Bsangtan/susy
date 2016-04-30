// -----------------------------------------------------------------
// Subtract two irrep matrices
// c <-- a - b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void sub_matrix(matrix *a, matrix *b, matrix *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++)
      CSUB(a->e[i][j], b->e[i][j], c->e[i][j]);
  }
}
// -----------------------------------------------------------------
