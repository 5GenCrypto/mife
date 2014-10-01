#ifndef DGSL__H
#define DGSL__H

#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include "gso.h"

/**
   Sampling algorithms
*/

typedef enum {
  DGSL_DETECT       = 0x0, //< detect which algorithm to use
  DGSL_IDENTITY     = 0x2, //< identity lattice
  DGSL_INLATTICE    = 0x3, //< c is in the lattice
  DGSL_COSET        = 0x4, //< c is not in the lattice
} dgsl_alg_t;

struct _dgsl_mp_t;

typedef struct _dgsl_rot_mp_t{
  fmpz_mat_t B; //< basis matrix
  mpfr_t *c; //< center
  fmpz *c_z; //< center
  mpfr_t sigma; //< Gaussian parameter
  dgs_disc_gauss_mp_t **D; //< storage for internal samplers
  int (*call)(fmpz *rop,  const struct _dgsl_rot_mp_t *self, gmp_randstate_t state); //< call this function

} dgsl_rot_mp_t;


typedef struct _dgsl_mp_t{
  fmpz_mat_t B; //< basis matrix
  mpfr_mat_t G; //< Gram-Schmidt matrix
  mpfr_t *c; //< center
  fmpz *c_z; //< center
  mpfr_t sigma; //< Gaussian parameter
  dgs_disc_gauss_mp_t **D; //< storage for internal samplers
  int (*call)(fmpz *rop,  const struct _dgsl_mp_t *self, gmp_randstate_t state); //< call this function

} dgsl_mp_t;

/**
   @param B basis matrix (copied), if matrix is 1 x n it is assumed that it represents a rotational basis mod X^n + 1
   @param sigma Gaussian width parameter (copied)
   @param c center
   @param algorithm
*/

dgsl_mp_t *dgsl_mp_init(const fmpz_mat_t B, mpfr_t sigma, mpfr_t *c, const dgsl_alg_t algorithm);

dgsl_rot_mp_t *dgsl_rot_mp_init(const fmpz_mat_t B, mpfr_t sigma, mpfr_t *c, const dgsl_alg_t algorithm);

/**
   Simple sampling when B is the identity

   @param rop return value (pre-allocated)
   @param self sampler
   @param state entropy source

*/
int dgsl_mp_call_inlattice(fmpz *rop,  const dgsl_mp_t *self, gmp_randstate_t state);
int dgsl_rot_mp_call_inlattice(fmpz *rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state);

int dgsl_mp_call_identity(fmpz *rop,  const dgsl_mp_t *self, gmp_randstate_t state);
int dgsl_rot_mp_call_identity(fmpz *rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state);

void dgsl_mp_clear(dgsl_mp_t *self);

void dgsl_rot_mp_clear(dgsl_rot_mp_t *self);

static inline void fmpz_mod_poly_sample_D(fmpz_mod_poly_t f, dgsl_mp_t *D, flint_rand_t randstate) {
  assert(D);
  assert(randstate->gmp_init);

  const long n = fmpz_mat_ncols(D->B);
  fmpz_mod_poly_realloc(f, n);
  D->call(f->coeffs, D, randstate->gmp_state);
  _fmpz_mod_poly_set_length(f, n);
  _fmpz_mod_poly_normalise(f);
}

static inline void fmpz_poly_sample_D(fmpz_poly_t f, dgsl_rot_mp_t *D, flint_rand_t randstate) {
  assert(D);
  assert(randstate->gmp_init);

  const long n = fmpz_mat_ncols(D->B);
  fmpz_poly_realloc(f, n);
  D->call(f->coeffs, D, randstate->gmp_state);
  _fmpz_poly_set_length(f, n);
  _fmpz_poly_normalise(f);
}

static inline void fmpz_mod_poly_sample_sigma(fmpz_mod_poly_t f, long len, mpfr_t sigma, flint_rand_t randstate) {
  fmpz_mat_t I;
  fmpz_mat_init(I, 1, len);
  fmpz_mat_one(I); dgsl_mp_t *D = dgsl_mp_init(I, sigma, NULL, DGSL_IDENTITY);

  fmpz_mod_poly_sample_D(f, D, randstate);

  dgsl_mp_clear(D);
  fmpz_mat_clear(I);
}

static inline void fmpz_poly_sample_sigma(fmpz_poly_t f, long len, mpfr_t sigma, flint_rand_t randstate) {
  fmpz_mat_t I;
  fmpz_mat_init(I, 1, len);
  fmpz_mat_one(I); dgsl_rot_mp_t *D = dgsl_rot_mp_init(I, sigma, NULL, DGSL_IDENTITY);

  fmpz_poly_sample_D(f, D, randstate);

  dgsl_rot_mp_clear(D);
  fmpz_mat_clear(I);
}

void _dgsl_rot_mp_sqrt_sigma_2(fmpq_poly_t rop, const fmpq_poly_t g, const mpfr_t sigma, const int r, const long n, const mpfr_prec_t prec);

#endif
