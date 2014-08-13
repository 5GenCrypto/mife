#include "config.h"
#include "misc.h"
#include "gghlite.h"

#include <mpfr.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint-addons/flint-addons.h>
#include <dgs/dgs.h>
#include <gpv/gpv.h>

#include <math.h>
#include <stdint.h>




void _gghlite_sample_g(gghlite_t self, flint_rand_t randstate) {
  assert(self->pk);
  assert(self->pk->n);

  fmpz_poly_init(self->g);
  fmpq_poly_t g_inv;
  fmpq_poly_init(g_inv);


  mpfr_t g_inv_norm;
  mpfr_init2(g_inv_norm, fmpz_sizeinbase(self->pk->q,2));

  while(1) {
    fflush(0);
    fmpz_poly_sample_sigma(self->g, self->pk->n, self->pk->sigma, randstate);

#ifdef GGHLITE_CHECK_PRIMALITY    
    /* 1. check if prime */
    if (!fmpz_poly_ideal_is_probaprime(self->g)) {
      printf("!p");
      continue;
    }
#endif
    /* 2. check norm of inverse */
    fmpz_poly_invert_mod_fmpq(g_inv, self->g, self->pk->cyclotomic_polynomial);
    if (!_gghlite_g_inv_check(self->pk, g_inv)) {
      printf("!n");
      continue;
    }
    break;
  }
  printf("\n");
  mpfr_clear(g_inv_norm);
  fmpq_poly_clear(g_inv);
}

void _gghlite_sample_h(gghlite_t self, flint_rand_t randstate) {
  assert(self->pk);
  assert(self->pk->n);
  assert(fmpz_cmp_ui(self->pk->q,0)>0);

  mpz_t q;
  mpz_init(q);
  fmpz_get_mpz(q, self->pk->q);
  mpfr_t sqrt_q;
  mpfr_init2(sqrt_q, fmpz_sizeinbase(self->pk->q,2));
  mpfr_set_z(sqrt_q, q, MPFR_RNDN);
  mpz_clear(q);
  mpfr_sqrt(sqrt_q, sqrt_q, MPFR_RNDN);
  assert(mpfr_cmp_ui(sqrt_q,0)>0);

  fmpz_poly_init(self->h);
  fmpz_poly_sample_sigma(self->h, self->pk->n, sqrt_q, randstate);

  mpfr_clear(sqrt_q);
}

void _gghlite_set_cyclotomic_polynomial(gghlite_pk_t self) {
  fmpz_poly_init(self->cyclotomic_polynomial);
  fmpz_poly_set_coeff_si(self->cyclotomic_polynomial, self->n, 1);
  fmpz_poly_set_coeff_si(self->cyclotomic_polynomial,       0, 1);
}

void _gghlite_set_n_q(gghlite_pk_t self) {
  const int64_t lambda = self->lambda;
  const int64_t kappa  = self->kappa;
  int64_t n_base = kappa * lambda * ceil(log2(lambda));

  int64_t c     = 1;
  int64_t log_n = ceil(log2(c*n_base));
  int64_t n     = ((int64_t)(1))<<log_n;
  int64_t log_q = _gghlite_log_q(log_n, kappa);

  while(1) {
    if (_gghlite_check_sec(log_q, n, lambda))
      break;
    c*=2;
    log_n = ceil(log2(c*n_base));
    n = ((int64_t)(1))<<log_n;
    log_q = _gghlite_log_q(log_n, kappa);
  }

  self->n = n;
  fmpz_init_set_ui(self->q, 2);
  fmpz_pow_ui(self->q, self->q, log_q);

  fmpz_t zeta;
  fmpz_sub_ui(self->q, self->q, 1);
  fmpz_divexact_si(zeta, self->q, self->n);

  fmpz_t tmp;
  while(1) {
    fmpz_mul_ui(tmp, zeta, self->n);
    fmpz_add_ui(tmp, tmp, 1);
    if (fmpz_is_probabprime(tmp)) {
      fmpz_set(self->q, tmp);
      break;
    }
    fmpz_add_ui(zeta, zeta, 1);
  }
  fmpz_clear(tmp);
  fmpz_clear(zeta);
}


/**
   `σ = 2·n·\sqrt{e·log(8n)/π}`, cf. [LSS14]_ p.16.
*/

int _gghlite_set_sigma(gghlite_pk_t self) {
  assert(self->n > 0);

  mpfr_init2(self->sigma, _gghlite_prec(self));
  mpfr_set_ui(self->sigma, self->n, MPFR_RNDN);
  mpfr_mul_ui(self->sigma, self->sigma, 2, MPFR_RNDN);

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 8, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);

  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_div(tmp, tmp, pi, MPFR_RNDN);
  mpfr_clear(pi);

  mpfr_t e;
  mpfr_init2(e, _gghlite_prec(self));
  mpfr_set_ui(e, 1, MPFR_RNDN);
  mpfr_exp(e, e, MPFR_RNDN);
  mpfr_mul(tmp, e, tmp, MPFR_RNDN);
  mpfr_clear(e);

  mpfr_sqrt(tmp, tmp, MPFR_RNDN);

  mpfr_mul(self->sigma, self->sigma, tmp, MPFR_RNDN);
  mpfr_clear(tmp);

  return 0;
}

/**
   `σ' ≥ 2n^{3/2}σ\sqrt{e·log(8n)/π}`, cf. [LSS14]_, p.17
*/

int _gghlite_set_sigma_p(gghlite_pk_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma,0)>0);

  mpfr_init2(self->sigma_p, _gghlite_prec(self));
  mpfr_set(self->sigma_p, self->sigma, MPFR_RNDN);

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 8, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);

  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_div(tmp, tmp, pi, MPFR_RNDN);
  mpfr_clear(pi);

  mpfr_t e;
  mpfr_init2(e, _gghlite_prec(self));
  mpfr_set_ui(e, 1, MPFR_RNDN);
  mpfr_exp(e, e, MPFR_RNDN);
  mpfr_mul(tmp, e, tmp, MPFR_RNDN);
  mpfr_clear(e);

  mpfr_sqrt(tmp, tmp, MPFR_RNDN);

  mpfr_mul(self->sigma_p, self->sigma_p, tmp, MPFR_RNDN);

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);

  mpfr_t pow;
  mpfr_init2(pow, _gghlite_prec(self));
  mpfr_set_d(pow, 1.5, MPFR_RNDN);
  mpfr_pow(tmp, tmp, pow, MPFR_RNDN);
  mpfr_clear(pow);

  mpfr_mul(self->sigma_p, self->sigma_p, tmp, MPFR_RNDN);
  mpfr_clear(tmp);

  return 0;
}

/**
  `ℓ_b = p_b/(2\sqrt{π·e·n})·σ'`, cf. [LSS14]_, p.17

  We assume p_b = 1
*/

void _gghlite_set_ell_b(gghlite_pk_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma_p, 0)>0);

  mpfr_init2(self->ell_b, _gghlite_prec(self));

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_set_ui(tmp, self->n, MPFR_RNDN);

  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_mul(tmp, tmp, pi, MPFR_RNDN);
  mpfr_clear(pi);

  mpfr_t e;
  mpfr_init2(e, _gghlite_prec(self));
  mpfr_set_ui(e, 1, MPFR_RNDN);
  mpfr_exp(e, e, MPFR_RNDN);
  mpfr_mul(tmp, e, tmp, MPFR_RNDN);
  mpfr_clear(e);

  mpfr_sqrt(tmp, tmp, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
  mpfr_ui_div(tmp, 1, tmp, MPFR_RNDN);
  mpfr_mul(self->ell_b, tmp, self->sigma_p, MPFR_RNDN);
  mpfr_clear(tmp);

}

/*
  `ℓ_g 4·sqrt(π·e·n)/(p_g·σ)`, cf. [LSS14]_ p.16`

    We assume p_b = 1
*/

void _gghlite_set_ell_g(gghlite_pk_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma, 0)>0);

  mpfr_init2(self->ell_g, _gghlite_prec(self));

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_set_ui(tmp, self->n, MPFR_RNDN);

  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_mul(tmp, tmp, pi, MPFR_RNDN);
  mpfr_clear(pi);

  mpfr_t e;
  mpfr_init2(e, _gghlite_prec(self));
  mpfr_set_ui(e, 1, MPFR_RNDN);
  mpfr_exp(e, e, MPFR_RNDN);
  mpfr_mul(tmp, e, tmp, MPFR_RNDN);
  mpfr_clear(e);

  mpfr_sqrt(tmp, tmp, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 4, MPFR_RNDN);
  mpfr_div(self->ell_g, tmp, self->sigma, MPFR_RNDN);
  mpfr_clear(tmp);
}

/**
 `σ^* ≥ n^{1.5}·(σ')²\sqrt{2πε_d^{-1}}/ℓ_b`, cf. [LSS14]_, p.18, eq. (8) with `εₑ^{-1} = O(log λ/κ)`.
*/

int _gghlite_set_sigma_s(gghlite_pk_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma_p, 0)>0);
  assert(mpfr_cmp_ui(self->ell_b, 0)>0);

  mpfr_init2(self->sigma_s, _gghlite_prec(self));
  mpfr_pow_ui(self->sigma_s, self->sigma_p, 2, MPFR_RNDN); // σ^* := (σ')^2

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_set_ui(tmp, self->lambda, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);
  mpfr_div_ui(tmp, tmp, self->kappa, MPFR_RNDN); // ε_d

  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_div(tmp, pi, tmp, MPFR_RNDN);  // π/ε_d
  mpfr_clear(pi);

  mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN); // 2π/ε_d

  mpfr_sqrt(tmp, tmp, MPFR_RNDN); // sqrt(2π/ε_d)

  mpfr_mul(self->sigma_s, self->sigma_s, tmp, MPFR_RNDN); // σ^* := (σ')^2 · sqrt(2π/ε_d)

  mpfr_div(self->sigma_s, self->sigma_s, self->ell_b, MPFR_RNDN); // σ^* := (σ')^2 · sqrt(2π/ε_d)/ℓ_b

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);

  mpfr_t pow;
  mpfr_init2(pow, _gghlite_prec(self));
  mpfr_set_d(pow, 1.5, MPFR_RNDN);
  mpfr_pow(tmp, tmp, pow, MPFR_RNDN); // n^(3/2)
  mpfr_clear(pow);

  mpfr_mul(self->sigma_s, self->sigma_s, tmp, MPFR_RNDN); // σ^* := n^(3/2) · (σ')^2 · sqrt(2π/ε_d)/ℓ_b

  mpfr_clear(tmp);

  return 0;
}

void gghlite_init(gghlite_t self, const int64_t lambda, const int64_t kappa, flint_rand_t randstate) {
  gghlite_init_step1(self, lambda, kappa);
  gghlite_init_step2(self, randstate);
}

void gghlite_init_step1(gghlite_t self, int64_t lambda, int64_t kappa) {
  if (lambda < 1)
    ggh_die("λ ≥ 1 required.");
  self->pk->lambda = lambda;

  if (kappa < 1)
    ggh_die("κ ≥ 1 required.");
  self->pk->kappa = kappa;

  _gghlite_set_n_q(self->pk);
  _gghlite_set_sigma(self->pk);
  _gghlite_set_ell_g(self->pk);
  _gghlite_set_sigma_p(self->pk);
  _gghlite_set_ell_b(self->pk);
  _gghlite_set_sigma_s(self->pk);
  _gghlite_set_cyclotomic_polynomial(self->pk);
}

void gghlite_init_step2(gghlite_t self, flint_rand_t randstate) {
  assert(self->pk->lambda);
  assert(self->pk->kappa);
  _gghlite_sample_g(self, randstate);
  _gghlite_sample_h(self, randstate);
}

void gghlite_pk_clear(gghlite_pk_t self) {
  fmpz_poly_clear(self->cyclotomic_polynomial);
  mpfr_clear(self->sigma_s);
  mpfr_clear(self->ell_b);
  mpfr_clear(self->sigma_p);
  mpfr_clear(self->ell_g);
  mpfr_clear(self->sigma);
  fmpz_clear(self->q);

}

void gghlite_clear(gghlite_t self, int clear_pk) {
  if (clear_pk)
    gghlite_pk_clear(self->pk);
  fmpz_poly_clear(self->g);
  fmpz_poly_clear(self->h);
}

void gghlite_print(const gghlite_t self) {
  printf("GGHLite Instance:\n");
  printf("         λ: %7ld\n",self->pk->lambda);
  printf("         k: %7ld\n",self->pk->kappa);
  printf("         n: %7ld\n",self->pk->n);
  printf("   log₂(q): %7ld (check: %d)\n", fmpz_sizeinbase(self->pk->q, 2), gghlite_check_sec(self->pk));
  printf("   log₂(σ): %7.1f\n",  log2(mpfr_get_d(self->pk->sigma, MPFR_RNDN)));
  printf(" log₂(ℓ_g): %7.1f\n",  log2(mpfr_get_d(self->pk->ell_g, MPFR_RNDN)));
  printf("  log₂(σ'): %7.1f\n",  log2(mpfr_get_d(self->pk->sigma_p, MPFR_RNDN)));
  printf(" log₂(ℓ_b): %7.1f\n",  log2(mpfr_get_d(self->pk->ell_b, MPFR_RNDN)));
  printf(" log₂(σ^*): %7.1f\n", log2(mpfr_get_d(self->pk->sigma_s, MPFR_RNDN)));
};
