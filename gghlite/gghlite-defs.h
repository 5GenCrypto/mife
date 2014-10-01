#ifndef _DEFS_H_
#define _DEFS_H_

#include <gghlite/config.h>
#include <gghlite/misc.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <dgsl/dgsl.h>

typedef fmpz_mod_poly_t gghlite_enc_t;
typedef fmpz_poly_t     gghlite_clr_t;

/**
   Maximum supported multi-linearity level.
*/

#define KAPPA (sizeof(uint64_t)<<3)

/**
   GGHLite Public key.

   FIELDS:

  - the security parameter `λ`
  - the multi-linearity parameter `κ ≤ KAPPA`
  - the rerandomisation mask such that if the `i`-th bit is set
    then rerandomisation elements for level `i+1` are generated
  - the dimension of the lattice `n`
  - the number of bits `ℓ` in each coefficient used for extraction
  - the modulus `q`
  - Gaussian width parameter `σ`   for sampling `g`
  - Gaussian width parameter `σ'`  for sampling `b_{k,i}`
  - Gaussian width parameter `σ^*` for sampling `ρ_i`
  - the bound `ℓ_b` on `σ_n(rot(B^(k)))`
  - the bound `ℓ_g` on `|g^-1|`
  - the modulus `x^n + 1` to specify the cyclotomic ring `ZZ[x]/(x^n+1)`
  - the zero-testing parameter `p_{zt}`
  - level-`k` encodings of zero `x_{k,i}` for each level `k` specified by rerandomisation mask
  - one level-1 encoding of 1
  - the discrete Gaussian distribution `D_{R,σ'}`
  - the discrete Gaussian distribution `D_{R,σ^*}`
   
*/

struct _gghlite_pk_struct {
  size_t lambda;
  size_t kappa;
  uint64_t rerand_mask;

  long n;
  long ell;
  fmpz_t q;
  mpfr_t sigma;
  mpfr_t sigma_p;
  mpfr_t sigma_s;
  mpfr_t ell_b;
  mpfr_t ell_g;
  mpfr_t ell_q;

  gghlite_enc_t modulus;
  gghlite_enc_t pzt;
  gghlite_enc_t x[KAPPA][2];
  gghlite_enc_t y;
  dgsl_rot_mp_t *D_sigma_p;
  dgsl_rot_mp_t *D_sigma_s;
};

typedef struct _gghlite_pk_struct gghlite_pk_t[1];

/**
   GGHLite Secret Key.

   FIELDS:

  - the corresponding GGHLite public key
  - an ideal generator `<g>`
  - `g^-1` in `K[x]/(x^n+1)`
  - the masking element `z`
  - the masking element `h`
  - an element `a mod <g> = 1`
  - an element `b mod <g> = 0`
*/

struct _gghlite_struct {
  gghlite_pk_t pk;

  gghlite_clr_t g;
  fmpq_poly_t g_inv;
  gghlite_enc_t z;
  gghlite_clr_t h;
  gghlite_clr_t a;
  gghlite_clr_t b[KAPPA][2];


  uint64_t t_is_prime;
  uint64_t t_is_subideal;
  uint64_t t_sample;
};

typedef struct _gghlite_struct gghlite_t[1];

#endif /* _DEFS_H_ */
