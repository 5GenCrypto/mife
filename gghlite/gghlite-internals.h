#ifndef _GGHLITE_INTERNALS_H_
#define _GGHLITE_INTERNALS_H_

#include <math.h>
#include <gghlite/gghlite-defs.h>
#include "oz/oz.h"

/**
   Check if `|g^-1| ≤ ℓ_g`
*/

static inline int _gghlite_g_inv_check(const gghlite_pk_t self, fmpq_poly_t g_inv) {
  mpfr_t g_inv_norm;
  mpfr_init2(g_inv_norm, mpfr_get_prec(self->ell_g));
  fmpq_poly_eucl_norm_mpfr(g_inv_norm, g_inv, MPFR_RNDN);
  int r = (mpfr_cmp(g_inv_norm, self->ell_g) <= 0);
  mpfr_clear(g_inv_norm);
  return r;
}

/**
   Return precision used for floating point computations
*/

static inline mpfr_prec_t _gghlite_prec(const gghlite_pk_t self) {
  const mpfr_prec_t prec = 2*self->lambda;
  if (prec < 53)
    return 53;
  else
    return prec;
}

#ifndef GGHLITE_HEURISTICS

/**
   Compute `σ` in double precision.
*/

static inline double _gghlite_sigma(double n) {
  const double e  = 2.71828182845905;
  const double pi = 3.14159265358979;
  const double sigma = 4*pi*n * sqrt(e*log(8*n)/pi);
  return sigma;
}

/**
   Compute `σ`.

   CONSTRAINTS:

   #. `σ = 4·π·n·\sqrt{e·log(8n)/π}/p_g`, cf. [LSS14]_ p.16.
*/

void _gghlite_pk_set_sigma(gghlite_pk_t self);


/**
  Compute `ℓ_g` in double precision.
*/

static inline double _gghlite_ell_g(double n) {
  const double e  = 2.71828182845905;
  const double pi = 3.14159265358979;
  const double sigma = _gghlite_sigma(n);
  const double ell_g = 4*sqrt(pi*e*n)/(sigma);
  return ell_g;
}

#else
/**
   Compute `σ`.

   CONSTRAINTS:

   #. `σ = \sqrt{n}
*/

/**
   Compute `σ` in double precision.
*/

static inline double _gghlite_sigma(double n) {
  const double pi = 3.14159265358979;
  const double sigma = sqrt(2*pi*n);
  return sigma;
}

void _gghlite_pk_set_sigma(gghlite_pk_t self);


/**
  Compute `ℓ_g` in double precision.
*/

static inline double _gghlite_ell_g(double n) {
  const double ell_g = 1/sqrt(n);
  return ell_g;
}

#endif


/**
   Number of small primes to test in sloppy variant
*/

static inline int _gghlite_nsmall_primes(const gghlite_pk_t self) {
  /* we try about 1% small primes first, where 1% relates to the total number of primes needed for
     multi-modular result */
  const long n = self->n;
  int nsp = ceil((log2(_gghlite_sigma(n)) + log2(n)/2.0) * n/100.0/(FLINT_BITS -1));
  if (nsp < 20)
    nsp = 20;
  return nsp;
}

/**
  Compute `ℓ_g`.

  CONSTRAINTS:

   #. `ℓ_g = 4·sqrt(π·e·n)/(p_g·σ)`, cf. [LSS14]_ p.16`

   @note: We assume `p_b = 1`
*/

void _gghlite_pk_set_ell_g(gghlite_pk_t self);

/**
   Compute `σ'` in double precision.
*/

static inline double _gghlite_sigma_p(double n) {
  const double sigma = _gghlite_sigma(n);
  const double e  = 2.71828182845905;
  const double pi = 3.14159265358979;
  const double sigma_p0 = 2.0 * pow(n, 1.5) * sigma * sqrt(e*log(8*n)/pi);
  const double sigma_p1 = 7.0 * pow(n, 2.5) * pow(log(n), 1.5) * sigma;

  if (sigma_p1 > sigma_p0)
    return sigma_p1;
  else
    return sigma_p0;
}

/**
   Compute `σ'`.

   CONSTRAINTS:

   #. `σ' ≥ 2n^{1.5}·σ\sqrt{e·log(8n)/π}/p_b`, cf. [LSS14]_, Eq. (5), p.16
   #. `σ' ≥ 7n^{2.5}·ln(n)^{1.5}·σ`, cf. [LSS14]_, p.17
*/

void _gghlite_pk_set_sigma_p(gghlite_pk_t self);

void _gghlite_pk_set_D_sigma_p(gghlite_pk_t self);

/**
   Compute `ℓ_b` in double precision.
*/

static inline double _gghlite_ell_b(double n) {
  const double sigma_p = _gghlite_sigma_p(n);
  const double e  = 2.71828182845905;
  const double pi = 3.14159265358979;
  const double ell_b = 1.0/(2.0*sqrt(pi*e*n)) * sigma_p;
  return ell_b;
}

/**
   Compute `ℓ_b`.

   CONSTRAINTS:

   #. `ℓ_b = p_b/(2\sqrt{π·e·n})·σ'`, cf. [LSS14]_, p.17

   @note: We assume `p_b = 1`
*/

void _gghlite_pk_set_ell_b(gghlite_pk_t self);

static inline double _gghlite_sigma_s(double n, double lambda, double kappa, const uint64_t rerand_mask) {
  if(rerand_mask == 0)
    return 1;
  if (rerand_mask > 1)
    ggh_die("Re-randomisation at higher levels is not implemented yet.");
  const double sigma_p = _gghlite_sigma_p(n);
  const double pi = 3.14159265358979;
  const double ell_g = _gghlite_ell_g(n);
  const double ell_b = _gghlite_ell_b(n);
  const double eps = log(lambda)/kappa;

  const double sigma_s0 = pow(n, 1.5) * ell_g * sigma_p * sqrt(2*log(4*n/eps)/pi);
  const double sigma_s1 = pow(n, 1.5) * pow(sigma_p, 2.0) * sqrt(8*pi/eps)/ell_b;

  if (sigma_s0 > sigma_s1)
    return sigma_s0;
  else
    return sigma_s1;
}

/**
   Compute `σ^*`.

   CONSTRAINTS:

   #. `σ^* ≥ n^{1.5}·ℓ_g·σ'·\sqrt{2·log(4nε_ρ^{-1})/π}`, cf. [LSS14]_, p.17, Eq. (8)
   #. `σ^* ≥ n^{1.5}·(σ')²\sqrt{8πε_d^{-1}}/ℓ_b`, cf. [LSS14]_, p.17, Eq. (9) with
   `εₑ^{-1} = O(log λ/κ)`.
*/

void _gghlite_pk_set_sigma_s(gghlite_pk_t self);

/**
   Init `D_{σ^*}`.
*/

void _gghlite_pk_set_D_sigma_s(gghlite_pk_t self);

void _gghlite_pk_set_q(gghlite_pk_t self);

static inline void _gghlite_get_q_mpfr(mpfr_t q, const gghlite_pk_t self, mpfr_rnd_t rnd) {
  assert(!fmpz_is_zero(self->q));
  mpz_t qz;
  mpz_init(qz);
  fmpz_get_mpz(qz, self->q);
  mpfr_set_z(q, qz, rnd);
  mpz_clear(qz);
}

void _gghlite_pk_set_ell(gghlite_pk_t self);

void _gghlite_sample_g(gghlite_t self, flint_rand_t randstate);

void _gghlite_sample_z(gghlite_t self, flint_rand_t randstate);

void _gghlite_sample_h(gghlite_t self, flint_rand_t randstate);

void _gghlite_sample_b(gghlite_t self, flint_rand_t randstate);

void _gghlite_sample_a(gghlite_t self, flint_rand_t randstate);

void _gghlite_set_pzt(gghlite_t self);

void _gghlite_set_x(gghlite_t self);

void _gghlite_set_y(gghlite_t self);

void gghlite_print_norms(const gghlite_t self);

void gghlite_print_times(const gghlite_t self);

#define S_TO_SIGMA 0.398942280401433 /* 1/sqrt(2*pi) */

dgsl_rot_mp_t *_gghlite_dgsl_from_poly(fmpz_poly_t g, mpfr_t sigma, fmpq_poly_t c, dgsl_alg_t algorithm);

dgsl_rot_mp_t *_gghlite_dgsl_from_n(const long n, mpfr_t sigma);


#define MAX_K 1024

extern double delta_from_k[MAX_K];

/**
   Return suitable k for given δ_0.

   :param delta_0: root-Hermite factor δ_0
*/

static inline long _gghlite_k_from_delta(const double delta_0) {
  long k;
  for(k=40; k<MAX_K; k++) {
    if (delta_from_k[k] <= delta_0)
      break;
  }
  if (k == MAX_K)
    ggh_die("Cannot establish required block size");
  return k;
}

/**
   Return expected number of BKZ rounds.

   :param n: lattice dimension
   :param k: block size

   See Theorem 1 in *Analyzing Blockwise Lattice Algorithms using Dynamical Systems* by Guillaume
   Hanrot, Xavier Pujol, and Damien Stehlé.
*/

static inline double _gghlite_repeat_from_n_k(const long n, const long k) {
  return 3*log2(n)  - 2*log2(k) + log2(log2(n));
}

/**
   Return expected cost of BKZ with SVP oracle implemented by enumeration.

   :param self: GGHLite public key

   See *On the concrete hardness of Learning with Errors* by Martin R. Albrecht, Rachel Player and
   Sam Scott, Cryptology ePrint Archive, Report 2015/046
*/

double gghlite_pk_cost_bkz_enum(const gghlite_pk_t self);

/**
   Return expected cost of BKZ with SVP oracle implemented by sieving.

   :param self: GGHLite public key

   See *On the concrete hardness of Learning with Errors* by Martin R. Albrecht, Rachel Player and
   Sam Scott, Cryptology ePrint Archive, Report 2015/046
*/

double gghlite_pk_cost_bkz_sieve(const gghlite_pk_t self);


#endif /* _GGHLITE_INTERNALS_H_ */
