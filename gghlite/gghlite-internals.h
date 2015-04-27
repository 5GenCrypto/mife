/**
   @file gghlite-internals.h
   @brief GGHLite internal API

   @note Most users should not have to rely on this interface.
*/

#ifndef _GGHLITE_INTERNALS_H_
#define _GGHLITE_INTERNALS_H_

#include <math.h>
#include <gghlite/gghlite-defs.h>
#include <gghlite/misc.h>
#include <oz/oz.h>

/**
   @brief Check if $|g^{-1}| ≤ ℓ_g@
*/

static inline int _gghlite_g_inv_check(const gghlite_params_t self, fmpq_poly_t g_inv) {
  mpfr_t g_inv_norm;
  mpfr_init2(g_inv_norm, mpfr_get_prec(self->ell_g));
  fmpq_poly_eucl_norm_mpfr(g_inv_norm, g_inv, MPFR_RNDN);
  int r = (mpfr_cmp(g_inv_norm, self->ell_g) <= 0);
  mpfr_clear(g_inv_norm);
  return r;
}

/**
   @brief Return precision used for floating point computations.
*/

static inline mpfr_prec_t _gghlite_prec(const gghlite_params_t self) {
  const mpfr_prec_t prec = 2*self->lambda;
  if (prec < 53)
    return 53;
  else
    return prec;
}

#ifndef GGHLITE_HEURISTICS

/**
   Compute $σ$ in double precision.
*/

static inline double _gghlite_sigma(double n) {
  const double e  = 2.71828182845905;
  const double pi = 3.14159265358979;
  const double sigma = 4*pi*n * sqrt(e*log(8*n)/pi);
  return sigma;
}

/**
   Compute $σ$.

   CONSTRAINTS:

   - @f$σ = 4·π·n·\sqrt{e·\log(8n)/π}/p_g@f$, cf. [LSS14]_ p.16.
*/

void _gghlite_params_set_sigma(gghlite_params_t self);


/**
  Compute $ℓ_g$ in double precision.
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
   Compute @f$σ@f$.

   CONSTRAINTS:

   - $σ = \sqrt{n}$
*/

/**
   Compute @f$σ@f$ in double precision.
*/

static inline double _gghlite_sigma(double n) {
  const double pi = 3.14159265358979;
  const double sigma = sqrt(2*pi*n);
  return sigma;
}

void _gghlite_params_set_sigma(gghlite_params_t self);


/**
  Compute @f$ℓ_g@f$ in double precision.
*/

static inline double _gghlite_ell_g(double n) {
  const double ell_g = 1/sqrt(n);
  return ell_g;
}

#endif


/**
   Number of small primes to test in sloppy variant
*/

static inline int _gghlite_nsmall_primes(const gghlite_params_t self) {
  /* we try about 1% small primes first, where 1% relates to the total number of primes needed for
     multi-modular result */
  const long n = self->n;
  int nsp = ceil((log2(_gghlite_sigma(n)) + log2(n)/2.0) * n/100.0/(FLINT_BITS -1));
  if (nsp < 20)
    nsp = 20;
  return nsp;
}

/**
  Compute @f$ℓ_g@f$.

  CONSTRAINTS:

  - @f$ℓ_g = 4·\sqrt(π·e·n)/(p_g·σ)@f$, cf. [LSS14]_ p.16

  @note We assume @f$p_g = 1@f$
*/

void _gghlite_params_set_ell_g(gghlite_params_t self);

/**
   Compute @f$σ'@f$ in double precision.
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
   Compute @f$σ'@f$.

   CONSTRAINTS:

   - @f$σ' ≥ 2n^{1.5}·σ\sqrt{e·\log(8n)/π}/p_b@f$, cf. [LSS14]_, Eq. (5), p.16
   - @f$σ' ≥ 7n^{2.5}·ln(n)^{1.5}·σ@f$, cf. [LSS14]_, p.17
*/

void _gghlite_params_set_sigma_p(gghlite_params_t self);

void _gghlite_params_set_D_sigma_p(gghlite_params_t self);

/**
   Compute $ℓ_b$ in double precision.
*/

static inline double _gghlite_ell_b(double n) {
  const double sigma_p = _gghlite_sigma_p(n);
  const double e  = 2.71828182845905;
  const double pi = 3.14159265358979;
  const double ell_b = 1.0/(2.0*sqrt(pi*e*n)) * sigma_p;
  return ell_b;
}

/**
   Compute $ℓ_b$.

   CONSTRAINTS

   - @f$ℓ_b = p_b/(2\sqrt{π·e·n})·σ'@f$, cf. [LSS14]_, p.17

   @note We assume $p_b = 1$
*/

void _gghlite_params_set_ell_b(gghlite_params_t self);

static inline double _gghlite_sigma_s(double n, double lambda, double kappa, const uint64_t rerand_mask) {
  if(rerand_mask == 0)
    return 1;
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
   @brief Compute $σ^*$.

   CONSTRAINTS:

   - @f$σ^* ≥ n^{1.5}·ℓ_g·σ'·\sqrt{2·\log(4nε_ρ^{-1})/π}@f$, cf. [LSS14]_, p.17, Eq. (8)
   - @f$σ^* ≥ n^{1.5}·(σ')²\sqrt{8πε_d^{-1}}/ℓ_b@f$, cf. [LSS14]_, p.17, Eq. (9) with
   @f$εₑ^{-1} = O(\log λ/κ)@f$.
*/

void _gghlite_params_set_sigma_s(gghlite_params_t self);

/**
   @brief Init $D_{σ^*}$.
*/

void _gghlite_params_set_D_sigma_s(gghlite_params_t self);

void _gghlite_params_set_q(gghlite_params_t self);

static inline void _gghlite_params_get_q_mpfr(mpfr_t q, const gghlite_params_t self, mpfr_rnd_t rnd) {
  assert(!fmpz_is_zero(self->q));
  mpz_t qz;
  mpz_init(qz);
  fmpz_get_mpz(qz, self->q);
  mpfr_set_z(q, qz, rnd);
  mpz_clear(qz);
}

void _gghlite_params_set_ell(gghlite_params_t self);

void _gghlite_sk_sample_g(gghlite_sk_t self, flint_rand_t randstate);

/**
   @brief Sample $z_i$ and $z_i^{-1}$.
 */

void _gghlite_sk_sample_z(gghlite_sk_t self, flint_rand_t randstate);

void _gghlite_sk_sample_h(gghlite_sk_t self, flint_rand_t randstate);

void _gghlite_sk_sample_b(gghlite_sk_t self, flint_rand_t randstate);

void _gghlite_sk_sample_a(gghlite_sk_t self, flint_rand_t randstate);

void _gghlite_sk_set_pzt(gghlite_sk_t self);

/**
   @brief Set $c = (c⋅g)/z$ for some small $c$.
 */

void _gghlite_sk_set_x(gghlite_sk_t self);

/**
   @brief Set $y = (1 + c⋅g)/z$ for some small $c$ for each source group in rerand_mask.
 */

void _gghlite_sk_set_y(gghlite_sk_t self);

void gghlite_sk_print_norms(const gghlite_sk_t self);

void gghlite_sk_print_times(const gghlite_sk_t self);

/**
   @brief dgsl samples proportionally to $\exp(-(x-c)²/(2σ²))$ but GGHLite is specifiied with respect to
   $\exp(-π(x-c)²/σ²)$. So we divide by $\sqrt{2π}$ in several places which we store in this macro.
*/

#define S_TO_SIGMA 0.398942280401433

dgsl_rot_mp_t *_gghlite_dgsl_from_poly(fmpz_poly_t g, mpfr_t sigma, fmpq_poly_t c, dgsl_alg_t algorithm, const oz_flag_t flags);

dgsl_rot_mp_t *_gghlite_dgsl_from_n(const long n, mpfr_t sigma, const oz_flag_t flags);

#define MAX_K 1024 //<! maximum block size for BKZ estimation

extern double delta_from_k[MAX_K];

/**
   @brief Return suitable $k$ for given $δ_0$.

   @param delta_0 root-Hermite factor $δ_0$
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

   @param n lattice dimension
   @param k block size

   See Theorem 1 in *Analyzing Blockwise Lattice Algorithms using Dynamical Systems* by Guillaume
   Hanrot, Xavier Pujol, and Damien Stehlé.
*/

static inline double _gghlite_repeat_from_n_k(const long n, const long k) {
  return 3*log2(n)  - 2*log2(k) + log2(log2(n));
}

/**
   @brief Return expected cost of BKZ with SVP oracle implemented by enumeration.

   @param self GGHLite public key

   See *On the concrete hardness of Learning with Errors* by Martin R. Albrecht, Rachel Player and
   Sam Scott, Cryptology ePrint Archive, Report 2015/046
*/

double gghlite_params_cost_bkz_enum(const gghlite_params_t self);

/**
   @brief Return expected cost of BKZ with SVP oracle implemented by sieving.

   @param self GGHLite public key

   See *On the concrete hardness of Learning with Errors* by Martin R. Albrecht, Rachel Player and
   Sam Scott, Cryptology ePrint Archive, Report 2015/046
*/

double gghlite_params_cost_bkz_sieve(const gghlite_params_t self);

/**
   @brief Return true if self represents a symmetric graded encoding scheme.
*/

static inline int gghlite_sk_is_symmetric(const gghlite_sk_t self) {
  return !(self->params->flags & GGHLITE_FLAGS_ASYMMETRIC);
}

/**
   Return log2() of the Euclidean norm of the numerator of ``op``.

   @param self      initialised GGHLite public key
   @param op        initialised
*/

double gghlite_log2_eucl_norm(const gghlite_sk_t self, const gghlite_enc_t op,
                              const size_t level, const size_t group);

/**
   @brief Multiply $f$ by zero-testing parameter $p_{zt}$.

   @param rop       initialised encoding, return value
   @param self      initialised GGHLite `params`
   @param f         valid encoding at level-$k$

   @ingroup internal-encodings
*/

void _gghlite_enc_extract_raw(gghlite_clr_t rop, const gghlite_params_t self, const gghlite_enc_t f);

#endif /* _GGHLITE_INTERNALS_H_ */
