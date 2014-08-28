#ifndef _GGHLITE_INTERNALS_H_
#define _GGHLITE_INTERNALS_H_

#include <math.h>

#include <gghlite/gghlite-defs.h>

/**
   Return candidate `log(q)` for a given `log(n)` and `κ`.
*/

static inline int64_t _gghlite_log_q(const long log_n, const long kappa) {
#ifndef GGHLITE_HEURISTICS
  long log_q = (10.5*log_n + log2(kappa)/2.0 + 10.5*log2(log_n)) * (8*kappa);
#else
  long log_q = ( 9.0*log_n + log2(kappa)/2.0 +  8*log2(log_n)) * (8*kappa);
#endif
  return log_q;
}

/**
   Check if `|g^-1| ≤ ℓ_g`
*/

static inline int _gghlite_g_inv_check(const gghlite_pk_t self, fmpq_poly_t g_inv) {
  mpfr_t g_inv_norm;
  mpfr_init2(g_inv_norm, mpfr_get_prec(self->ell_g));
  _fmpq_vec_2norm_mpfr(g_inv_norm, fmpq_poly_numref(g_inv), fmpq_poly_denref(g_inv), self->n);
  int r = (mpfr_cmp(g_inv_norm, self->ell_g) <= 0);
  mpfr_clear(g_inv_norm);
  return r;
}



/**
   Return precision used for floating point computations
*/

static inline long _gghlite_prec(const gghlite_pk_t self) {
  if (self->lambda*self->lambda < 53)
    return 53;
  else
    return self->lambda*self->lambda;
}
/**
   Compute `n` and `q`.
*/

void _gghlite_pk_set_n_q(gghlite_pk_t self);


static inline void _gghlite_get_q_mpfr(mpfr_t q, const gghlite_pk_t self, mpfr_rnd_t rnd) {
  assert(!fmpz_is_zero(self->q));
  mpz_t qz;
  mpz_init(qz);
  fmpz_get_mpz(qz, self->q);
  mpfr_set_z(q, qz, rnd);
  mpz_clear(qz);
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
  const double sigma = sqrt(n);
  return sigma;
}

void _gghlite_pk_set_sigma(gghlite_pk_t self);


/**
  Compute `ℓ_g` in double precision.
*/

static inline double _gghlite_ell_g(double n) {
  const double ell_g = 1/sqrt(n*log(n));
  return ell_g;
}

#endif


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

static inline double _gghlite_sigma_s(double n, double lambda, double kappa) {
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

/**
   Set `x^n + 1`.
*/

static inline void _gghlite_pk_set_modulus(gghlite_pk_t self) {
  assert(self->n);
  assert(!fmpz_is_zero(self->q));
  fmpz_mod_poly_init(self->modulus, self->q);
  fmpz_mod_poly_set_coeff_ui(self->modulus, self->n, 1);
  fmpz_mod_poly_set_coeff_ui(self->modulus,       0, 1);
}

/**
   Check if security constraints are satisfied.

   @todo: refine this function to reflect actual lattice attacks
 */

static inline int _gghlite_check_sec(int64_t log_q, size_t n, size_t lambda) {
  /* cost by the lattice rule of thumb is >= lambda */
  int rt  = (3.0/8.0 * log_q < n/(double)lambda);
  /* heuristic cost of one LLLL: n^3·logB^2 */
  int lll = (lambda <= (3*log2(n) + 2*log2(log_q))) && (lambda <= (4*log2(n) + 2*log2(log_q)));
  return rt || lll;
}

/**
   Check security constraints are satisfied
*/

static inline int gghlite_check_sec(const gghlite_pk_t self) {
  return _gghlite_check_sec(fmpz_sizeinbase(self->q, 2), self->n, self->lambda);
}

static inline int _gghlite_check_func(int64_t log_q, size_t n, size_t lambda, size_t kappa) {
  double sigma_p = _gghlite_sigma_p(n);
  double sigma_s = _gghlite_sigma_s(n, lambda, kappa);
  double lhs = 8*kappa * (log2(3.0) + 1.5*log2(n) + log2(sigma_p) + log2(sigma_s));
  if (lhs < log_q)
    return 1;
  else
    return 0;
}

void _gghlite_sample_g(gghlite_t self, flint_rand_t randstate);

void _gghlite_sample_z(gghlite_t self, flint_rand_t randstate);

void _gghlite_sample_h(gghlite_t self, flint_rand_t randstate);

void _gghlite_sample_b(gghlite_t self, flint_rand_t randstate);

void _gghlite_sample_a(gghlite_t self, flint_rand_t randstate);

void _gghlite_set_pzt(gghlite_t self);

void _gghlite_set_x(gghlite_t self);

void _gghlite_set_y(gghlite_t self);

void gghlite_print_norms(const gghlite_t self);

#endif /* _GGHLITE_INTERNALS_H_ */
