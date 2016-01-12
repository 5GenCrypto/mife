/**
   @file gghlite.h
   @brief GGHLite public API

   Parties who wish to use graded encoding schemes as a building block should only need this API to
   realise all functionality required.

   Both symmetric and asymmetric graded encodings schemes are supported.

   NAMING CONVENTION:

   - Functions starting with `gghlite_params_` should be thought of as methods on `gghlite_params_t`
   - Functions starting with `gghlite_sk_` should be thought of as methods on `gghlite_sk_t`
   - Functions starting with `gghlite_enc` should be thought of as methods on `gghlite_enc_t`
 */

#ifndef _GGHLITE_H_
#define _GGHLITE_H_

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#include <dgsl/dgsl.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mod_poly.h>

#include <gghlite/gghlite-defs.h>
#include <gghlite/misc.h>

/**
   @defgroup params Instance Generation & Parameters
*/


void gghlite_params_initzero(gghlite_params_t self, size_t lambda, size_t kappa, size_t gamma);

/**
   @brief Generate parameters for GGHLite instance requiring no randomness.

   @param self        GGHLite `params`, all fields are overwritten
   @param lambda      security parameter $λ > 0$
   @param kappa       multi-linearity parameter $κ > 0$
   @param rerand_mask generate re-randomisation elements for level $i$ if ``1<<(i-1) & rerand_mask``
   @param flags       flags controlling verbosity etc.

   @ingroup params
*/
void gghlite_params_init_gamma(gghlite_params_t self, size_t lambda, size_t kappa, size_t gamma, uint64_t rerand_mask, gghlite_flag_t flags);

static inline void gghlite_params_init(gghlite_params_t self, size_t lambda, size_t kappa, uint64_t rerand_mask, gghlite_flag_t flags) {
	gghlite_params_init_gamma(self, lambda, kappa, kappa, rerand_mask, flags);
}



/**
   @brief Generate parameters for GGHLite jigsaw puzzle instance requiring no randomness.

   This light wrapper sets `GGHLITE_FLAGS_ASYMMETRIC | GGHLITE_FLAGS_GOOD_G_INV` and sets rerand
   mask to 0x0.

   @param self        GGHLite `params`, all fields are overwritten
   @param lambda      security parameter $λ > 0$
   @param kappa       multi-linearity parameter $κ > 0$
   @param flags       flags controlling verbosity etc.

   @ingroup params
*/


static inline void gghlite_jigsaw_params_init(gghlite_params_t self, size_t lambda, size_t kappa, size_t gamma, gghlite_flag_t flags) {
  gghlite_params_init_gamma(self, lambda, kappa, gamma, 0x0, flags | GGHLITE_FLAGS_ASYMMETRIC | GGHLITE_FLAGS_GOOD_G_INV);
}

/**
   @brief Generate fields requiring randomness.

   @param self       GGHLite secret key, all fields but `params` are overwritten
   @param randstate  entropy source, assumes `flint_randinit(randstate)` and
                     `_flint_rand_init_gmp(randstate)` was called

   @ingroup params
*/

void gghlite_sk_init(gghlite_sk_t self, aes_randstate_t randstate);


void gghlite_sk_set_D_g(gghlite_sk_t self);
void gghlite_params_set_D_sigmas(gghlite_params_t params);

/**
   @brief Initialise a new GGHLite instance.

   @param self        GGHLite secret key, all fields are overwritten
   @param lambda      security parameter $λ > 0$
   @param kappa       multi-linearity parameter $κ > 0$
   @param rerand_mask generate re-randomisation elements for level $i$ if ``1<<(i-1) & rerand_mask``
   @param flags       flags controlling the behaviour of the algorithms such as verbosity etc.
   @param randstate   entropy source, assumes `flint_randinit(randstate) and
                      `_flint_rand_init_gmp(randstate)` was called

   @ingroup params
*/

void gghlite_init(gghlite_sk_t self, const size_t lambda, const size_t kappa, const size_t gamma,
                  const uint64_t rerand_mask, const gghlite_flag_t flags, aes_randstate_t randstate);


/**
   @brief Initialise a new GGHLite jigsaw puzzle instance.

   @param self        GGHLite secret key, all fields are overwritten
   @param lambda      security parameter $λ > 0$
   @param kappa       multi-linearity parameter $κ > 0$
   @param flags       flags controlling the behaviour of the algorithms such as verbosity etc.
   @param randstate   entropy source, assumes `flint_randinit(randstate) and
                      `_flint_rand_init_gmp(randstate)` was called

   @ingroup params
*/

static inline void gghlite_jigsaw_init_gamma(gghlite_sk_t self, size_t lambda, size_t kappa, size_t gamma,
                                       gghlite_flag_t flags, aes_randstate_t randstate) {
  gghlite_init(self, lambda, kappa, gamma, 0x0, flags | GGHLITE_FLAGS_ASYMMETRIC | GGHLITE_FLAGS_GOOD_G_INV,
               randstate);
}

static inline void gghlite_jigsaw_init(gghlite_sk_t self, size_t lambda, size_t kappa,
                                       gghlite_flag_t flags, aes_randstate_t randstate) {
	gghlite_jigsaw_init_gamma(self, lambda, kappa, kappa, flags, randstate);
}



/**
   @brief Get a shallow copy of `params` of `op`.

   @param rop  GGHlite `params`, all fields are overwritten
   @param  op  initialised GGHLite Instance

   @ingroup params
*/

void gghlite_params_ref(gghlite_params_t rop, gghlite_sk_t op);

/**
   @brief Clear GGHLite `params`.

   @param self all fields are cleared

   @ingroup params
*/

void gghlite_params_clear(gghlite_params_t self);

/**
   @brief Clear GGHLite instance.

   The parameter `clear_params` is useful to clear all secret data that was required to produce the
   public set of parameters stored in `self->params`. After `gghlite_sk_clear(self, 0)` was called,
   it is safe to forget about `self` and to use `self->params` independently (which can later be
   cleared by calling `gghlite_params_clear(self->params)`. A reference to `self->params` can be
   obtained by calling `gghlite_params_ref(·, self->params)`.

   @param self          all fields are cleared, except perhaps pk
   @param clear_params  if set, `params` are also cleared

   @ingroup params
*/

void gghlite_sk_clear(gghlite_sk_t self, int clear_params);

/* parameter estimation functions */
double gghlite_params_get_enc(const gghlite_params_t self);
void gghlite_params_test_kappa_enc_size(size_t lambda, size_t max_kappa, FILE *fp);

/**
   @brief Print sizes of parameters to STDOUT.

   @param self assumes `gghlite_params_init_params(self,·,·,·)` was called

   @ingroup params
*/

void gghlite_params_print(const gghlite_params_t self);

/**
   @brief Return true if params are for a symmetric graded encoding scheme.

   @param self      initialised GGHLite `params`

   @ingroup params
*/

static inline int gghlite_params_is_symmetric(const gghlite_params_t self) {
  return !(self->flags & GGHLITE_FLAGS_ASYMMETRIC);
}

/**
   @brief Return root-Hermite factor $δ_0$ required to break the scheme

   @param self      initialised GGHLite `params`


   @ingroup params
*/

double gghlite_params_get_delta_0(const gghlite_params_t self);


/**
   @brief Check if security constraints are satisfied.

   @param self      initialised GGHLite `params`


   @ingroup params
*/

int gghlite_params_check_sec(const gghlite_params_t self);

/**
   @brief Return true if rerandomisation elements are available for $i+1$ level.

   @param self      initialised GGHLite `params`

   @ingroup params
*/

static inline int gghlite_params_have_rerand(const gghlite_params_t self, const size_t i) {
  return (self->rerand_mask & (1ULL)<<i);
}

/**
   @brief Initialise a new clear text element
*/

#define gghlite_clr_init  fmpz_poly_init

/**
   @brief Clear a new clear text element
*/

#define gghlite_clr_clear fmpz_poly_clear

/**
   @brief Return true if two clear text elements are equal
*/

#define gghlite_clr_equal fmpz_poly_equal

/**
   @defgroup encodings Manipulating Encodings
*/


/**
   @brief Initialise encoding to zero at level 0.

   @param op   uninitialised polynomial
   @param self initialised GGHLite `params`

   @ingroup encodings
*/

void gghlite_enc_init(gghlite_enc_t op, const gghlite_params_t self);

#define gghlite_enc_clear fmpz_mod_poly_clear

/**
   @brief Rerandomise encoding at level $k$ in group $G_i$.

   Computes $f = f + ρ_0·b_{k,0} + ρ_1·b_{k,1}$ where $ρ_i ← D_{R,σ^*}$.

   @param f         initialised encoding at level $k$
   @param self      initialised GGHLite `params`
   @param k         level `1 ≤ k ≤ κ`
   @param i         group index (use zero in symmetric setting)
   @param randstate entropy source, assumes `flint_randinit(randstate)` and
                    `_flint_rand_init_gmp(randstate)` was called

   @note Note that we have no means to check that `rop` is indeed a level-$k$ encoding.

   @ingroup encodings
*/

void gghlite_enc_rerand(gghlite_enc_t rop, const gghlite_params_t self, const gghlite_enc_t op,
                        size_t k, size_t i, aes_randstate_t randstate);

/**
   @brief Raise encoding at level $k$ to level $l$ and re-randomise if requested.

   @param rop       initialised encoding
   @param self      initialised GGHLite `params`
   @param op        initialised encoding at level $l$
   @param l         targer level `1 ≤ l ≤ κ`
   @param k         current level of ``op`` satisfying `0 ≤ k ≤ l`
   @param i         group index (zero in symmetric setting)
   @param rerand    flag controlling if re-randomisation is run after raising
   @param randstate entropy source, assumes `flint_randinit(randstate)` and
                    `_flint_rand_init_gmp(randstate)` was called

   @note Note that we have no means to check that `op` is indeed a level-$k$ encoding.

   @ingroup encodings
*/

void gghlite_enc_raise(gghlite_enc_t rop, const gghlite_params_t self, const gghlite_enc_t op,
                       size_t l, size_t k, size_t i,
                       int rerand, aes_randstate_t randstate);

/**
   @brief Raise an encoding at level $0$ to level $l$ in group $G_0$.

   @param rop       initialised encoding
   @param self      initialised GGHLite `params`
   @param op        initialised encoding at level $l$
   @param l         targer level $1 ≤ l ≤ κ$
   @param randstate entropy source, assumes `flint_randinit(randstate)` and
                    `_flint_rand_init_gmp(randstate)` was called

   @note Note that we have no means to check that `op` is indeed a level-$k$ encoding.

   @ingroup encodings
*/

static inline void gghlite_enc_raise0(gghlite_enc_t rop, gghlite_params_t self, gghlite_enc_t op,
                                      size_t l, aes_randstate_t randstate) {

  // TODO: fix have_rerand API
  int rerand = (gghlite_params_have_rerand(self, l-1)) ? 1 : 0;
  gghlite_enc_raise(rop, self, op, l, 0, 0, rerand, randstate);
}


/**
   @brief Set `op` to an encoding of $c$ at level $k$ in group $G_i$.

   @param op        uninitialised polynomial
   @param c         a small integer
   @param self      initialised GGHLite `params`
   @param k         targer level $1 ≤ k ≤ κ$
   @param i         group index (zero in symmetric setting)
   @param rerand    flag controlling if re-randomisation is run after raising
   @param randstate entropy source, assumes `flint_randinit(randstate)` and
                    `_flint_rand_init_gmp(randstate)` was called

   @ingroup encodings
*/

static inline void gghlite_enc_set_ui(gghlite_enc_t op, unsigned long c, const gghlite_params_t self,
                                      const size_t k, const size_t i, const int rerand,
                                      aes_randstate_t randstate) {
  fmpz_mod_poly_oz_ntt_set_ui(op, c, self->n);
  if(k>0)
    gghlite_enc_raise(op, self, op, k, 0, i, rerand, randstate);
}

/**
   @brief Set `op` to an encoding of `c` at level 0.

   @param op        uninitialised polynomial
   @param c         a small integer
   @param self      initialised GGHLite `params`

   @ingroup encodings
*/

static inline void gghlite_enc_set_ui0(gghlite_enc_t op, unsigned long c, const gghlite_params_t self) {
  fmpz_mod_poly_oz_ntt_set_ui(op, c, self->n);
}

/**
   @brief Copy encoding

   @ingroup encodings
*/

#define gghlite_enc_set fmpz_mod_poly_set

/**
   @brief Sample a new random encoding at levek $k$ in group $i$.

   @param rop       initialised encoding, return value
   @param self      initialised GGHLite `params`
   @param k         targer level `0 ≤ k ≤ κ`
   @param i         group $G_i$ (zero in the symmetric case)
   @param randstate entropy source, assumes `flint_randinit(randstate) and
                    `_flint_rand_init_gmp(randstate)` was called

   @ingroup encodings
*/

void gghlite_enc_sample(gghlite_enc_t rop, gghlite_params_t self, size_t k, size_t i, aes_randstate_t randstate);

/**
   @brief Encode $f$ at level-$k$ in group $G_i$.

   @param rop       initialised encoding, return value
   @param self      initialised GGHLite instance
   @param f         an element in $\\ZZ[x]/(x^n+1)$
   @param k         targer level $0 ≤ k ≤ κ$
   @param i         group $G_i$ (zero in the symmetric case)
   @param rerand    flag controlling if re-randomisation is run after raising
   @param randstate entropy source, assumes `flint_randinit(randstate)` and
                    `_flint_rand_init_gmp(randstate)` was called

   @note If `self` is an asymmetric map only, then $k ≤ 1$ is required.

   @ingroup encodings
*/

void gghlite_enc_set_gghlite_clr(gghlite_enc_t rop, const gghlite_sk_t self, const gghlite_clr_t f,
                                 const size_t k, int *group, const int rerand,
                                 aes_randstate_t randstate);


/**
   @brief Encode $f$ at level-$0$.

   @param rop       initialised encoding, return value
   @param self      initialised GGHLite instance
   @param f         an element in $\\ZZ[x]/(x^n+1)$
   @param randstate entropy source, assumes `flint_randinit(randstate)` and
                    `_flint_rand_init_gmp(randstate)` was called

   @note If `self` is an asymmetric map only, then $k ≤ 1$ is required.

   @ingroup encodings
*/

static inline void gghlite_enc_set_gghlite_clr0(gghlite_enc_t rop, const gghlite_sk_t self, const gghlite_clr_t f,
                                                aes_randstate_t randstate) {
	int *group = malloc(self->params->gamma * sizeof(int));
	memset(group, 0, self->params->gamma * sizeof(int));
	group[0] = 1;
  gghlite_enc_set_gghlite_clr(rop, self, f, 0, group, 0, randstate);
  free(group);
}

/**
   @brief Compute $h = f·g$.

   @param h         initialised encoding, return value
   @param self      initialised GGHLite `params`
   @param f         valid encoding
   @param g         valid encoding
*/

static inline void gghlite_enc_mul(gghlite_enc_t h, const gghlite_params_t self, const gghlite_enc_t f, const gghlite_enc_t g) {
  fmpz_mod_poly_oz_ntt_mul(h, f, g, self->n);
}

/**
   @brief Compute $h = f+g$.

   @param h         initialised encoding, return value
   @param self      initialised GGHLite `params`
   @param f         valid encoding
   @param g         valid encoding

   @ingroup encodings
*/

static inline void gghlite_enc_add(gghlite_enc_t h, const gghlite_params_t self, const gghlite_enc_t f, const gghlite_enc_t g) {
  fmpz_mod_poly_add(h, f, g);
}

/**
   @brief Compute $h = f-g$.

   @param h         initialised encoding, return value
   @param self      initialised GGHLite `params`
   @param f         valid encoding
   @param g         valid encoding

   @ingroup encodings
*/

static inline void gghlite_enc_sub(gghlite_enc_t h, const gghlite_params_t self, const gghlite_enc_t f, const gghlite_enc_t g) {
  fmpz_mod_poly_sub(h, f, g);
}

/**
   @brief Extract canonical string from $f$

   @param rop       initialised encoding, return value
   @param self      initialised GGHLite `params`
   @param f        valid encoding at level-$k$

   @ingroup encodings
*/

void gghlite_enc_extract(fmpz_poly_t rop, const gghlite_params_t self, const gghlite_enc_t f);

/**
   @brief Return 1 if $f$ is an encoding of zero at level $κ$

   @param self      initialised GGHLite `params`
   @param f         valid encoding at level-$k$

   @ingroup encodings
*/

int gghlite_enc_is_zero(const gghlite_params_t self, const gghlite_enc_t op);


#endif //_GGHLITE_H_
