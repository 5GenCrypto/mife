#ifndef GPV__H
#define GPV__H

#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include "gso.h"

/**
   Sampling algorithms
*/

typedef enum {
  GPV_DETECT       = 0x0, //< detect which algorithm to use
  GPV_IDENTITY     = 0x2, //< identity lattice
  GPV_INLATTICE    = 0x3, //< c is in the lattice
  GPV_COSET        = 0x4, //< c is not in the lattice
} gpv_alg_t;

struct _gpv_mp_t;

typedef struct _gpv_mp_t{
  fmpz_mat_t B; //< basis matrix
  mpfr_mat_t G; //< Gram-Schmidt matrix
  mpfr_t *c; //< center
  fmpz *c_z; //< center
  mpfr_t sigma; //< Gaussian parameter
  dgs_disc_gauss_mp_t **D; //< storage for internal samplers
  int (*call)(fmpz *rop,  const struct _gpv_mp_t *self, gmp_randstate_t state); //< call this function

} gpv_mp_t;

/**
   @param B basis matrix (copied), if matrix is 1 x n it is assumed that it represents a rotational basis mod X^n + 1
   @param sigma Gaussian width parameter (copied)
   @param c center
   @param algorithm
*/

gpv_mp_t *gpv_mp_init(const fmpz_mat_t B, mpfr_t sigma, mpfr_t *c, const gpv_alg_t algorithm);

/**
   Simple sampling when B is the identity

   @param rop return value (pre-allocated)
   @param self sampler
   @param state entropy source

*/

int gpv_mp_call_simple(fmpz *rop,  const gpv_mp_t *self, gmp_randstate_t state);

void gpv_mp_clear(gpv_mp_t *self);

static inline void fmpz_poly_sample_D(fmpz_poly_t f, gpv_mp_t *D, flint_rand_t randstate) {
  assert(D);
  assert(randstate->gmp_init);

  const long n = fmpz_mat_ncols(D->B);
  fmpz_poly_realloc(f, n);
  do {
    D->call(f->coeffs, D, randstate->gmp_state);
  } while (fmpz_is_zero(f->coeffs +(n-1)));
  _fmpz_poly_set_length(f, n);
}

static inline void fmpz_poly_sample_sigma(fmpz_poly_t f, long len, mpfr_t sigma, flint_rand_t randstate) {
  fmpz_mat_t I;
  fmpz_mat_init(I, 1, len);
  fmpz_mat_one(I); gpv_mp_t *D = gpv_mp_init(I, sigma, NULL, GPV_IDENTITY);

  fmpz_poly_sample_D(f, D, randstate);

  gpv_mp_clear(D);
  fmpz_mat_clear(I);
}

#endif
