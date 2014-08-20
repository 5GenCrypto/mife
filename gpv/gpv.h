#ifndef GPV__H
#define GPV__H

#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>

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
  fmpz *c; //< center
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

gpv_mp_t *gpv_mp_init(const fmpz_mat_t B, const mpfr_t sigma,const fmpz *c, const gpv_alg_t algorithm);

/**
   Simple sampling when B is the identity

   @param rop return value (pre-allocated)
   @param self sampler
   @param state entropy source

*/

int gpv_mp_call_simple(fmpz *rop,  const gpv_mp_t *self, gmp_randstate_t state);


void gpv_mp_clear(gpv_mp_t *self);

#endif
