#ifndef DGSL__H
#define DGSL__H

#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <oz/oz.h>
#include "gso.h"

/**
   Sampling algorithms
*/

typedef enum {
  DGSL_DETECT        = 0x0, //< detect which algorithm to use
  DGSL_IDENTITY      = 0x2, //< identity lattice
  DGSL_INLATTICE     = 0x3, //< c is in the lattice
  DGSL_COSET         = 0x4, //< c is not necessarily in the lattice
  DGSL_GPV_INLATTICE = 0x8, //< c is in the lattice, use GPV
} dgsl_alg_t;

struct _dgsl_mp_t;

typedef struct _dgsl_rot_mp_t{
  long n;        //< dimension
  fmpz_poly_t B; //< basis matrix
  fmpz_poly_t c_z; //< center
  fmpq_poly_t c; //< center
  mpfr_t sigma; //< Gaussian parameter
  dgs_disc_gauss_mp_t **D; //< storage for internal samplers

  int (*call) (fmpz_poly_t rop,  const struct _dgsl_rot_mp_t *self, gmp_randstate_t state); //< call this function

  long   r;
  mpfr_t r_f;

  fmpq_poly_t sigma_sqrt;
  mpfr_prec_t prec;

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
   @param B basis matrix (copied), if matrix is 1 x n it is assumed that it
          represents a rotational basis mod X^n + 1
   @param sigma Gaussian width parameter (copied)
   @param c center
   @param algorithm
*/

dgsl_mp_t *dgsl_mp_init(const fmpz_mat_t B, mpfr_t sigma, mpfr_t *c, const dgsl_alg_t algorithm);

dgsl_rot_mp_t *dgsl_rot_mp_init(const long n, const fmpz_poly_t B, mpfr_t sigma, fmpq_poly_t c, const dgsl_alg_t algorithm, const oz_flag_t flags);

int dgsl_mp_call_inlattice(fmpz *rop,  const dgsl_mp_t *self, gmp_randstate_t state);
int dgsl_rot_mp_call_gpv_inlattice(fmpz_poly_t rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state);

int _dgsl_rot_mp_call_inlattice_multiplier(fmpz_poly_t rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state);
int dgsl_rot_mp_call_inlattice(fmpz_poly_t rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state);
int dgsl_rot_mp_call_plus1(fmpz_poly_t rop, const dgsl_rot_mp_t *self, gmp_randstate_t state, uint64_t flags);

/**
   Return a fresh sample when B is the identity

   @param rop return value (pre-allocated)
   @param self sampler with self->call == dgsl_mp_call_identity
   @param state entropy source
*/

int dgsl_mp_call_identity(fmpz *rop,  const dgsl_mp_t *self, gmp_randstate_t state);

/**
   Return a fresh sample when B is the identity

   @param rop return value (pre-allocated)
   @param self sampler with self->call == dgsl_mp_call_identity
   @param state entropy source
*/

int dgsl_rot_mp_call_identity(fmpz_poly_t rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state);

void dgsl_mp_clear(dgsl_mp_t *self);

void dgsl_rot_mp_clear(dgsl_rot_mp_t *self);

static inline void fmpz_mod_poly_sample_D(fmpz_mod_poly_t f, dgsl_rot_mp_t *D, flint_rand_t randstate) {
  assert(D);
  assert(randstate->gmp_init);

  fmpz_poly_t tmp;
  fmpz_poly_init(tmp);
  fmpz_poly_set_fmpz_mod_poly(tmp, f);
  D->call(tmp, D, randstate->gmp_state);
  fmpz_mod_poly_set_fmpz_poly(f, tmp);
  fmpz_poly_clear(tmp);
}

static inline void fmpz_poly_sample_D(fmpz_poly_t f, dgsl_rot_mp_t *D, flint_rand_t randstate) {
  assert(D); assert(randstate->gmp_init);
  D->call(f, D, randstate->gmp_state);
}

static inline void fmpz_poly_sample_D_plus1(fmpz_poly_t f, dgsl_rot_mp_t *D, flint_rand_t randstate) {
  assert(D); assert(randstate->gmp_init);
  assert(D->call == dgsl_rot_mp_call_inlattice);
  dgsl_rot_mp_call_plus1(f, D, randstate->gmp_state, OZ_VERBOSE);
}

static inline void fmpz_mod_poly_sample_sigma(fmpz_mod_poly_t f, long len, mpfr_t sigma, flint_rand_t randstate) {
  fmpz_poly_t I;
  fmpz_poly_init(I);
  fmpz_poly_one(I);
  dgsl_rot_mp_t *D = dgsl_rot_mp_init(len, I, sigma, NULL, DGSL_IDENTITY, OZ_VERBOSE);

  fmpz_mod_poly_sample_D(f, D, randstate);

  dgsl_rot_mp_clear(D);
  fmpz_poly_clear(I);
}

static inline void fmpz_poly_sample_sigma(fmpz_poly_t f, long len, mpfr_t sigma, flint_rand_t randstate) {
  fmpz_poly_t I;
  fmpz_poly_init(I);
  fmpz_poly_one(I);
  dgsl_rot_mp_t *D = dgsl_rot_mp_init(len, I, sigma, NULL, DGSL_IDENTITY, OZ_VERBOSE);

  fmpz_poly_sample_D(f, D, randstate);

  dgsl_rot_mp_clear(D);
  fmpz_poly_clear(I);
}

void _dgsl_rot_mp_sqrt_sigma_2(fmpq_poly_t rop, const fmpz_poly_t g, const mpfr_t sigma,
                              const int r, const long n, const mpfr_prec_t prec, const uint64_t flags);

void fmpz_poly_disc_gauss_rounding(fmpz_poly_t rop, const fmpq_poly_t x, const mpfr_t r_f, gmp_randstate_t randstate);

void fmpq_poly_sample_D1(fmpq_poly_t f, const int n, const mpfr_prec_t prec, gmp_randstate_t state);

#endif
