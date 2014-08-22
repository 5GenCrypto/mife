#ifndef _DEFS_H_
#define _DEFS_H_

#include <gghlite/config.h>
#include <gghlite/misc.h>

/**
   Maximum supported multi-linearity level
*/

#define KAPPA (sizeof(uint64_t)<<3)

/**
   GGHLite Public key

*/

struct _gghlite_pk_struct {
  size_t lambda;        //< security parameter `λ`
  size_t kappa;         //< multi-linearity parameter `κ ≤ KAPPA`
  uint64_t rerand_mask; //< set i-th bit to generate rerandomisation elements for level i+1

  long n;         //< dimension of the lattice
  fmpz_t q;       //< modulus
  mpfr_t sigma;   //< Gaussian width `σ` for `g`
  mpfr_t sigma_p; //< Gaussian width `σ'` for `b_{k,i}`
  mpfr_t sigma_s; //< Gaussian width `σ^*` for `r_i`
  mpfr_t ell_b;   //< bound `ℓ_b` on `σ_n(rot(B^(k)))`
  mpfr_t ell_g;   //< bound `ℓ_g` on `|g^-1| `

  fmpz_mod_poly_t modulus;        //< modulus `x^n + 1`
  fmpz_mod_poly_t pzt;            //< zero-testing parameter `p_{zt}`
  fmpz_mod_poly_t x[KAPPA][2]; //< re-randomisers `x_{k,i}` for level `k`
  fmpz_mod_poly_t y;              //< encoding of 1
  gpv_mp_t *D_sigma_p;            //< discrete Gaussian distribution `D_{R,σ'}`
  gpv_mp_t *D_sigma_s;            //< discrete Gaussian distribution `D_{R,σ^*}`
};

typedef struct _gghlite_pk_struct gghlite_pk_t[1];

/**
   GGHLite Secret Key
*/

struct _gghlite_struct {
  gghlite_pk_t pk;          //< GGHLite public key

  fmpz_poly_t g;            //< ideal generator `<g>`
  fmpq_poly_t g_inv;        //< `g^-1` in `K[x]/(x^n+1)`
  fmpz_mod_poly_t z;        //< mask `z`
  fmpz_poly_t h;            //< mask `h`
  fmpz_poly_t a;            //< `a mod <g> = 1`
  fmpz_poly_t b[KAPPA][2];  //< `b mod <g> = 0`
};

typedef struct _gghlite_struct gghlite_t[1];

#endif /* _DEFS_H_ */
