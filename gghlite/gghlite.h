#ifndef _GGHLITE_H_
#define _GGHLITE_H_

#include <stdint.h>
#include <math.h>

#include <gpv/gpv.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint-addons/flint-addons.h>

#include <gghlite/gghlite-defs.h>

/**
   Generate parameters for GGHLite instance requiring no randomness.
*/

void gghlite_pk_init_params(gghlite_pk_t self, size_t lambda, size_t kappa, uint64_t rerand_mask);

/**
   Generate parameters requiring randomness.
*/

void gghlite_init_instance(gghlite_t self, flint_rand_t randstate);

/**
   Generate a new GGHLite instance

   @param self        GGHLite instance.
   @param lambda      security parameter λ > 0 for 2^λ bit security.
   @param kappa       multi-linearity parameter 0< κ < KAPPA
   @param rerand_mask set i-th bit to generate rerandomisation elements for level i+1
   @param randstate   source of entropy
*/

void gghlite_init(gghlite_t self, const size_t lambda, const size_t kappa,
                  const uint64_t rerand_mask, flint_rand_t randstate);


void gghlite_pk_clear(gghlite_pk_t self);

void gghlite_clear(gghlite_t self, int clear_pk);

/**
   Print sizes of involved parameters to STDOUT
*/

void gghlite_print_params(const gghlite_pk_t self);

void fmpz_mod_poly_init_gghlite(fmpz_mod_poly_t op, gghlite_pk_t self);

void gghlite_rerand(fmpz_mod_poly_t rop, gghlite_pk_t self, long k, flint_rand_t randstate);

void gghlite_lift(fmpz_mod_poly_t rop, gghlite_pk_t self, long k, long kprime, int rerand, flint_rand_t randstate);

void gghlite_sample(fmpz_mod_poly_t rop, gghlite_pk_t self, long k, flint_rand_t randstate);

#endif //_GGHLITE_H_
