// -----------------------------------------------------------------
// Add result of complex scalar multiplication on irrep matrix
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void c_scalar_mult_add_mat(matrix *a, matrix *b, complex *s, matrix *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      c->e[i][j].real = a->e[i][j].real + b->e[i][j].real * s->real \
                                        - b->e[i][j].imag * s->imag;
      c->e[i][j].imag = a->e[i][j].imag + b->e[i][j].imag * s->real \
                                        + b->e[i][j].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------
