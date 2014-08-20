#ifndef GGHLITE__H
#define GGHLITE__H

#include <stdint.h>
#include <math.h>

#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint-addons/flint-addons.h>

/**
   Maximum supported multi-linearity level
*/

#define MAXKAPPA (sizeof(uint64_t)<<3)

/**
   GGHLite Public key

   @todo should cyclotomic_polynomial be of type fmpz_mod_poly_t?
*/

struct _gghlite_pk_struct {
  size_t lambda;   //< security parameter `λ`
  size_t kappa;    //< multi-linearity parameter `κ ≤ MAXKAPPA`

  long n;         //< dimension of the lattice
  fmpz_t q;       //< modulus
  mpfr_t sigma;   //< Gaussian width `σ` for `g`
  mpfr_t sigma_p; //< Gaussian width `σ'` for `b_{k,i}`
  mpfr_t sigma_s; //< Gaussian width `σ^*` for `r_i`
  mpfr_t ell_b;   //< bound `ℓ_b` on `σ_n(rot(B^(k)))`
  mpfr_t ell_g;   //< bound `ℓ_g` on `|g^-1| `

  fmpz_poly_t cyclotomic_polynomial; //< modulus `x^n + 1`
  fmpz_mod_poly_t pzt;               //< zero-testing parameter `p_{zt}`
  fmpz_mod_poly_t x[MAXKAPPA][2];    //< re-randomisers `x_{k,i}` for level `k`
  fmpz_mod_poly_t y;                 //< encoding of 1
};

typedef struct _gghlite_pk_struct gghlite_pk_t[1];

/**
   GGHLite Secret Key
*/

struct _gghlite_struct {
  gghlite_pk_t pk;              //< GGHLite public key

  fmpz_poly_t g;                //< ideal generator `<g>`
  fmpq_poly_t g_inv;            //< `g^-1` in `K[x]/(x^n+1)`
  fmpz_mod_poly_t z;            //< mask `z`
  fmpz_poly_t h;                //< mask `h`
  fmpz_poly_t a;                //< `a mod <g> = 1`
  fmpz_poly_t b[MAXKAPPA][2];  //< `b mod <g> = 0`
};

typedef struct _gghlite_struct gghlite_t[1];

/**
   Return candidate `log(q)` for a given `log(n)` and `κ`.
*/

static inline int64_t _gghlite_log_q(const long log_n, const long kappa) {
  long log_q = (10.5*log_n + log2(kappa)/2.0 + log2(log_n)) * (8*kappa);
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
   Check if security constraints are satisfied.

   @todo: refine this function to reflect actual lattice attacks
 */

static inline int _gghlite_check_sec(int64_t log_q, size_t n, size_t lambda) {
  /* cost by the lattice rule of thumb is >= lambda */
  int rt  = (3.0/8.0 * log_q < n/(double)lambda);
  /* cost of one L2: d^(3+ε) * n * (d + logB) * logB */
  int lll = (lambda <= (5*log2(n) + log2(log_q))) && (lambda <= (4*log2(n) + 2*log2(log_q)));
  return rt || lll;
}


/**
   Check security constraints are satisfied
*/

static inline int gghlite_check_sec(const gghlite_pk_t self) {
  return _gghlite_check_sec(fmpz_sizeinbase(self->q, 2), self->n, self->lambda);
}

/**
   Return precision used for floating point computations
*/

static inline long _gghlite_prec(const gghlite_pk_t self) {
  if (2*self->lambda < 53)
    return 53;
  else
    return 2*self->lambda;
}

/**
   Compute `n` and `q`.
*/

void _gghlite_set_n_q(gghlite_pk_t self);

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

void _gghlite_set_sigma(gghlite_pk_t self);

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

/**
  Compute `ℓ_g`.

  CONSTRAINTS:

   #. `ℓ_g = 4·sqrt(π·e·n)/(p_g·σ)`, cf. [LSS14]_ p.16`

   @note: We assume `p_b = 1`
*/

void _gghlite_set_ell_g(gghlite_pk_t self);

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

void _gghlite_set_sigma_p(gghlite_pk_t self);

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

void _gghlite_set_ell_b(gghlite_pk_t self);

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

void _gghlite_set_sigma_s(gghlite_pk_t self);

/**
   Set `x^n + 1`.
*/

static inline void _gghlite_set_cyclotomic_polynomial(gghlite_pk_t self) {
  fmpz_poly_init(self->cyclotomic_polynomial);
  fmpz_poly_set_coeff_si(self->cyclotomic_polynomial, self->n, 1);
  fmpz_poly_set_coeff_si(self->cyclotomic_polynomial,       0, 1);
}

/**
   Generate parameters for GGHLite instance requiring no randomness.
*/

static inline void gghlite_init_step1(gghlite_t self, size_t lambda, size_t kappa) {
  assert(lambda > 0);
  assert((kappa > 0) && (kappa <= MAXKAPPA));

  self->pk->lambda = lambda;
  self->pk->kappa = kappa;

  _gghlite_set_n_q(self->pk);
  _gghlite_set_sigma(self->pk);
  _gghlite_set_ell_g(self->pk);
  _gghlite_set_sigma_p(self->pk);
  _gghlite_set_ell_b(self->pk);
  _gghlite_set_sigma_s(self->pk);
  _gghlite_set_cyclotomic_polynomial(self->pk);
}

void _gghlite_sample_g(gghlite_t self, flint_rand_t randstate);
void _gghlite_sample_z(gghlite_t self, flint_rand_t randstate);
void _gghlite_sample_h(gghlite_t self, flint_rand_t randstate);

/**
 */

void _gghlite_sample_b(gghlite_t self, const uint64_t rerand_mask, flint_rand_t randstate);

void _gghlite_set_pzt(gghlite_t self);

void _gghlite_set_x(gghlite_t self);

/**
   Generate parameters requiring randomness.
*/

static inline void gghlite_init_step2(gghlite_t self, const uint64_t rerand_mask, flint_rand_t randstate) {
  assert(self->pk->lambda);
  assert(self->pk->kappa);

  _gghlite_sample_g(self, randstate);
  _gghlite_sample_z(self, randstate);
  _gghlite_sample_h(self, randstate);
  _gghlite_sample_b(self, rerand_mask, randstate);

  _gghlite_set_pzt(self);
  _gghlite_set_x(self);
}

/**
   Generate a new GGHLite instance

   @param self        GGHLite instance.
   @param lambda      security parameter λ > 0 for 2^λ bit security.
   @param kappa       multi-linearity parameter 0< κ < MAXKAPPA
   @param rerand_mask set i-th bit to generate rerandomisation elements for level i+1
   @param randstate   source of entropy
*/

static inline void gghlite_init(gghlite_t self, const size_t lambda, const size_t kappa,
                                const uint64_t rerand_mask, flint_rand_t randstate) {
  gghlite_init_step1(self, lambda, kappa);
  gghlite_init_step2(self, rerand_mask, randstate);
}

static inline void gghlite_pk_clear(gghlite_pk_t self) {
  fmpz_mod_poly_clear(self->pzt);

  for(long k=0; k<self->kappa; k++) {
    fmpz_mod_poly_clear(self->x[k][0]);
    fmpz_mod_poly_clear(self->x[k][1]);
  }

  fmpz_poly_clear(self->cyclotomic_polynomial);
  mpfr_clear(self->sigma_s);
  mpfr_clear(self->ell_b);
  mpfr_clear(self->sigma_p);
  mpfr_clear(self->ell_g);
  mpfr_clear(self->sigma);
  fmpz_clear(self->q);
}

static inline void gghlite_clear(gghlite_t self, int clear_pk) {
  for(long k=0; k<self->pk->kappa; k++) {
    fmpz_poly_clear(self->b[k][0]);
    fmpz_poly_clear(self->b[k][1]);
  }
  fmpz_poly_clear(self->h);
  fmpz_mod_poly_clear(self->z);
  fmpz_poly_clear(self->g);
  if (clear_pk)
    gghlite_pk_clear(self->pk);
}

/**
   Print sizes of involved parameters to STDOUT
*/

void gghlite_print_params(const gghlite_pk_t self);

void gghlite_pk_clear(gghlite_pk_t self);
void gghlite_clear(gghlite_t self, int clear_pk);

#endif //GGHLITE__H
