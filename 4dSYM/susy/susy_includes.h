// -----------------------------------------------------------------
// Include files for supersymmetric evolution
#include "../include/config.h"  // Keep this first
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>             // For print_var.c, setup.c, gauge_info.c
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);
int update();
void update_h(Real eps);
void update_u(Real eps);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Susy routines
// Lots of things to initialize and set up
void setup_lambda();
void setup_PtoP();
void setup_FQ();
void setup_offset();
void setup_qclosed_offset();
void setup_rhmc();
void fermion_rep();

// Helper routines for action and force computations
void compute_plaqdet();
void compute_DmuUmu();    // Computes plaqdet if G is non-zero
void compute_Fmunu();

// Gaussian random source
int grsource(Twist_Fermion *source);

// Action routines
double d_action(Twist_Fermion **source, Twist_Fermion ***sol);
double d_gauge_action(int do_det);
double d_fermion_action();

// Force routines
double gauge_force(Real eps);
double fermion_force(Real eps, Twist_Fermion *source, Twist_Fermion **psim);
double det_force(Real eps);

// Fermion matrix--vector operators (D & D^2) and multi-mass CG
void fermion_op(Twist_Fermion *src, Twist_Fermion *dest, int sign);
void hdelta0_field(Twist_Fermion *src, Twist_Fermion *dest);
int congrad_multi_field(Twist_Fermion *src, Twist_Fermion **psim,
                        int MaxCG, Real RsdCG, Real *size_r);

// Compute average link Tr[Udag U] / N_c
// Number of blocking steps only affects output formatting
void d_link(int bl);
void d_link_frep(int bl);

Real order(int i, int j, int k, int l, int m);
void epsilon();

// Basic Twist_Fermion and gauge field manipulations
// May eventually move to libraries
void dump_TF(Twist_Fermion *vec);
void conjTF(Twist_Fermion *src, Twist_Fermion *dest);
void copy_TF(Twist_Fermion *src, Twist_Fermion *dest);
void clear_TF(Twist_Fermion *dest);
Real magsq_TF(Twist_Fermion *vec);
complex TF_dot(Twist_Fermion *a, Twist_Fermion *b);
void scalar_mult_add_TF(Twist_Fermion *src1, Twist_Fermion *src2, Real s,
                        Twist_Fermion *dest);
void scalar_mult_TF(Twist_Fermion *src, Real s, Twist_Fermion *dest);
void gauge_field_copy_f(field_offset src, field_offset dest);

// Random gauge transformation for testing gauge invariance
void random_gauge_trans(Twist_Fermion *TF);

// Determinant-related routines
void measure_det();
void widths();        // Widths of plaquette and det distributions
complex find_det(su3_matrix_f *Q);
void det_project(su3_matrix_f *in, su3_matrix_f *out);

// Adjugate matrix needed by det_force
void adjugate(su3_matrix_f *in, su3_matrix_f *out);

// Matrix inverse is just adjugate divided by determinant
void invert(su3_matrix_f *in, su3_matrix_f *out);

// Modified Wilson loops use invert and path
void path(int *dir, int *sign, int length);
void rsymm();
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// More measurements
#ifdef CORR
// Konishi and SUGRA correlators
void setup_P();
void compute_Ba();
void d_correlator();    // Projected to zero spatial momentum

// Map (x, y, z, t) to scalar displacements r
Real A4map(int x_in, int y_in, int z_in, int t_in);
void d_correlator_r();  // Functions of r
#endif
#ifdef BILIN
// vevs to explore susy breaking
int d_bilinear();       // Fermion bilinear
int d_susyTrans();      // Supersymmetry transformation
#endif
#ifdef PL_CORR
// Polyakov loop correlator -- NOT CURRENTLY IN USE
void ploop_c();         // Computes Polyakov loop at each spatial site
void print_var3(char *label);
#endif
#ifdef WLOOP
// Wilson loops -- including determinant division and polar projection
// These look at correlators of products of temporal links,
// which requires gauge fixing
Real A4map_slice(int x, int y, int z);
void hvy_pot(int do_det);
void hvy_pot_polar();

// These construct explicit paths along lattice principal axes, for checking
void hvy_pot_loop(int do_det);
void hvy_pot_polar_loop();

// Use LAPACK in the polar projection
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/zheev.html
// First argument turns on eigenvector computations
// Second argument chooses between storing upper or lower triangle
// Third and fifth arguments are the dimensions of the matrix
// Fourth argument is that matrix, overwritten by the eigenvectors
// Sixth argument holds the computed eigenvalues
// Seventh argument is real workspace of size given by the eighth argument
// Ninth argument is real workspace of size 3 * NCOL - 2
// Final argument reports success or information about failure
void zheev_(char *doV, char *uplo, int *N1, double *store, int *N2,
            double *eigs, double *work, int *Nwork, double *Rwork, int *stat);
void polar(su3_matrix_f *a, su3_matrix_f *b);
#endif

// Monopole computation uses find_det
void monopole();

#ifdef MCRG
void block_mcrg(int bl);
void blocked_plaq(int Nstout, int bl);    // Also monitors det and widths
void blocked_ops(int Nstout, int bl);
void blocked_ploop(int Nstout, int bl);
void blocked_rsymm(int Nstout, int bl);
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Eigenvalue routines
#ifdef EIG
#include "primme.h"
int make_evs(int Nvec, Twist_Fermion **eigVec, double *eigVal, int flag);
void check_Dmat(int Nvec, Twist_Fermion **eigVec);

// Use LAPACK to diagonalize <psi_j | D | psi_i>
// on the subspace of Ddag.D eigenvalues psi
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/zgeev.html
// First two arguments turn off eigenvector computations
// Third and fifth arguments are the dimensions of the matrix
// Fourth argument is that matrix, which will be overwritten
// Sixth argument holds the computed eigenvalues
// Seventh and ninth arguments are eigenvectors
// Eighth and tenth arguments are the dimensions of the eigenvectors
// Eleventh argument is real workspace, of size given by the twelfth argument
// Thirteenth argument is real workspace, of size given by the third argument
// Final argument reports success or information about failure
void zgeev_(char *doL, char *doR, int *N1, double *store, int *N2, double *eigs,
            double *dumL, int *NL, double *dumR, int *NR,
            double *work, int *Nwork, double *Rwork, int *stat);
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Pfaffian phase
#ifdef PHASE
void d_phase();
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Stochastic mode number
#ifdef MODE
void coefficients();
void step(Twist_Fermion *src, Twist_Fermion *res);
#endif
// -----------------------------------------------------------------
