#include <string.h>

#include "gghlite-internals.h"
#include "gghlite.h"

void _gghlite_params_initzero(gghlite_params_t self, size_t lambda, size_t kappa, size_t gamma) {
  memset(self, 0, sizeof(struct _gghlite_params_struct));

  self->lambda = lambda;
  self->kappa = kappa;
	self->gamma = gamma;

  mpfr_init2(self->sigma, _gghlite_prec(self));
  mpfr_init2(self->ell_g, _gghlite_prec(self));
  mpfr_init2(self->sigma_p, _gghlite_prec(self));
  mpfr_init2(self->ell_b, _gghlite_prec(self));
  mpfr_init2(self->sigma_s, _gghlite_prec(self));
  mpfr_init2(self->xi, _gghlite_prec(self));
}

void _gghlite_params_set_ell(gghlite_params_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma, 0)>0);

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_mul(tmp, tmp, self->sigma, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 8, MPFR_RNDN);
  mpfr_log2(tmp, tmp, MPFR_RNDN);
  self->ell = ceil(mpfr_get_d(tmp, MPFR_RNDN));
  mpfr_clear(tmp);
}

void _gghlite_params_set_q(gghlite_params_t self) {
  const size_t kappa  = self->kappa;

  mpfr_t q_base;
  mpfr_init2(q_base, _gghlite_prec(self));

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_sqrt(tmp, tmp, MPFR_RNDN);
  mpfr_pow_ui(tmp, tmp, 3, MPFR_RNDN);
  mpfr_mul(tmp, tmp, self->sigma_p, MPFR_RNDN);
  mpfr_mul(tmp, tmp, self->sigma_p, MPFR_RNDN);

  mpfr_set(q_base, tmp, MPFR_RNDN); //(σ')^2 * n^(1.5)

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_sqrt(tmp, tmp, MPFR_RNDN);
  mpfr_pow_ui(tmp, tmp, 3, MPFR_RNDN);

  mpfr_mul(tmp, tmp, self->sigma_p, MPFR_RNDN);
  mpfr_mul(tmp, tmp, self->sigma_s, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);

  mpfr_add(q_base, q_base, tmp, MPFR_RNDN); // (σ')^2 · n^(1.5) + 2σ^* · σ' · n^(1.5)

  mpfr_pow_ui(q_base, q_base, kappa, MPFR_RNDN); // ((σ')^2 · n^(1.5) + 2σ^* · σ' · n^(1.5))^κ

  mpfr_mul_ui(q_base, q_base, self->n, MPFR_RNDN);
  mpfr_mul(q_base, q_base, self->ell_g, MPFR_RNDN); // n · ℓ_g · ((σ')^2 · n^(1.5) + 2σ^* · σ' · n^(1.5))^κ

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_sqrt(tmp, tmp, MPFR_RNDN);
  mpfr_pow_ui(tmp, tmp, kappa-1, MPFR_RNDN);
  mpfr_mul(q_base, tmp, q_base, MPFR_RNDN); // n · ℓ_g · sqrt(n)^(κ-1) · ((σ')^2 · n^(1.5) + 2σ^* · σ' · n^(1.5))^κ

  mpfr_t log_q_base;
  mpfr_init2(log_q_base, _gghlite_prec(self));
  mpfr_log2(log_q_base, q_base, MPFR_RNDN);

  mpfr_set_ui(tmp, self->ell, MPFR_RNDN);
  mpfr_add_ui(tmp, tmp, self->lambda, MPFR_RNDN);
  mpfr_div_ui(tmp, tmp, 2, MPFR_RNDN);

  mpfr_set(self->xi, tmp, MPFR_RNDN);

  mpfr_set_ui(tmp, self->ell, MPFR_RNDN);
  mpfr_add_ui(tmp, tmp, self->lambda, MPFR_RNDN);
  mpfr_add(tmp, tmp, log_q_base, MPFR_RNDN);
  mpfr_div(self->xi, self->xi, tmp, MPFR_RNDN);

  mpfr_set(tmp, self->xi, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
  mpfr_neg(tmp, tmp, MPFR_RNDN);
  mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
  mpfr_ui_div(tmp, 2, tmp, MPFR_RNDN);

  mpfr_pow(tmp, q_base, tmp, MPFR_RNDN);

  mpz_t q;
  mpz_init(q);
  mpfr_get_z(q, tmp, MPFR_RNDN);
  fmpz_set_mpz(self->q, q);
  mpz_clear(q);

  mpfr_clear(tmp);
  mpfr_clear(log_q_base);
  mpfr_clear(q_base);

  fmpz_fdiv_q_2exp(self->q, self->q, n_flog(self->n,2)+1);
  fmpz_mul_2exp(self->q, self->q, n_flog(self->n,2)+1);
  fmpz_add_ui(self->q, self->q, 1);
  while(1) {
    if(fmpz_is_probabprime(self->q))
      break;
    fmpz_add_ui(self->q, self->q, 2*self->n);
  }
}

#ifndef GGHLITE_HEURISTICS

void _gghlite_params_set_sigma(gghlite_params_t self) {
  assert(self->n > 0);

  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);

  mpfr_set_ui(self->sigma, self->n, MPFR_RNDN);
  mpfr_mul_ui(self->sigma, self->sigma, 4, MPFR_RNDN);
  mpfr_mul(self->sigma, self->sigma, pi, MPFR_RNDN);

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 8, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);
  mpfr_div(tmp, tmp, pi, MPFR_RNDN);


  mpfr_t e;
  mpfr_init2(e, _gghlite_prec(self));
  mpfr_set_ui(e, 1, MPFR_RNDN);
  mpfr_exp(e, e, MPFR_RNDN);
  mpfr_mul(tmp, e, tmp, MPFR_RNDN);
  mpfr_clear(e);

  mpfr_sqrt(tmp, tmp, MPFR_RNDN);

  mpfr_mul(self->sigma, self->sigma, tmp, MPFR_RNDN);
  mpfr_clear(tmp);
  mpfr_clear(pi);
}

void _gghlite_params_set_ell_g(gghlite_params_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma, 0)>0);

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
#else

void _gghlite_params_set_sigma(gghlite_params_t self) {
  assert(self->n > 0);

  mpfr_set_ui(self->sigma, self->n, MPFR_RNDN);

  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);

  mpfr_mul_ui(self->sigma, self->sigma, 2, MPFR_RNDN);
  mpfr_mul(self->sigma, self->sigma, pi, MPFR_RNDN);

  mpfr_sqrt(self->sigma, self->sigma, MPFR_RNDN);

  mpfr_clear(pi);
}

void _gghlite_params_set_ell_g(gghlite_params_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma, 0)>0);

  mpfr_set_ui(self->ell_g, self->n, MPFR_RNDN);
  mpfr_sqrt(self->ell_g, self->ell_g, MPFR_RNDN);
  mpfr_ui_div(self->ell_g, 1, self->ell_g, MPFR_RNDN);
}

#endif


void _gghlite_params_set_sigma_p(gghlite_params_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma,0)>0);

  mpfr_t pow;
  mpfr_init2(pow, _gghlite_prec(self));
  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_t e;
  mpfr_init2(e, _gghlite_prec(self));
  mpfr_set_ui(e, 1, MPFR_RNDN);
  mpfr_exp(e, e, MPFR_RNDN);

  mpfr_t sigma_p0;
  mpfr_init2(sigma_p0, _gghlite_prec(self));
  mpfr_t sigma_p1;
  mpfr_init2(sigma_p1, _gghlite_prec(self));

  /* `σ' ≥ 2n^{1.5}·σ\sqrt{e·log(8n)/π}` */
  mpfr_set(sigma_p0, self->sigma, MPFR_RNDN); // `σ`

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 8, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);  // `log(8n)`

  mpfr_div(tmp, tmp, pi, MPFR_RNDN); // `log(8n)/π`

  mpfr_mul(tmp, e, tmp, MPFR_RNDN); // `e·log(8n)/π`

  mpfr_sqrt(tmp, tmp, MPFR_RNDN); // `sqrt(e·log(8n)/π)`

  mpfr_mul(sigma_p0, sigma_p0, tmp, MPFR_RNDN); // `σ·sqrt(e·log(8n)/π)`

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_set_d(pow, 1.5, MPFR_RNDN);
  mpfr_pow(tmp, tmp, pow, MPFR_RNDN); // `n^(3/2)`
  mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN); // `2n^(3/2)`

  mpfr_mul(sigma_p0, sigma_p0, tmp, MPFR_RNDN); // `2n^(3/2)·σ·sqrt(e·log(8n)/π)`

  /*  `σ' ≥ 7n^{2.5}·ln(n)^{1.5}·σ` */
  mpfr_set(sigma_p1, self->sigma, MPFR_RNDN);

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_set_d(pow, 2.5, MPFR_RNDN);
  mpfr_pow(tmp, tmp, pow, MPFR_RNDN); // `n^(2.5)`
  mpfr_mul_ui(tmp, tmp, 7, MPFR_RNDN); // `7n^(2.5)`
  mpfr_mul(sigma_p1, sigma_p1, tmp, MPFR_RNDN);  // `7n^(5/2)·σ`

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);
  mpfr_set_d(pow, 1.5, MPFR_RNDN);
  mpfr_pow(tmp, tmp, pow, MPFR_RNDN); // `log(n)^(3/2)`

  mpfr_mul(sigma_p1, sigma_p1, tmp, MPFR_RNDN); // `7n^{5/2}·σ · log(n)^{3/2}`


  if (mpfr_cmp(sigma_p1, sigma_p0) >= 0)
    mpfr_set(self->sigma_p, sigma_p1, MPFR_RNDN);
  else
    mpfr_set(self->sigma_p, sigma_p0, MPFR_RNDN);

  mpfr_clear(tmp);
  mpfr_clear(e);
  mpfr_clear(pi);
  mpfr_clear(pow);
  mpfr_clear(sigma_p0);
  mpfr_clear(sigma_p1);
}

void _gghlite_params_set_D_sigma_p(gghlite_params_t self) {
  assert(self->n);
  assert(!mpfr_zero_p(self->sigma_p));
  const oz_flag_t flags = (self->flags & GGHLITE_FLAGS_QUIET) ? 0 : OZ_VERBOSE;
  self->D_sigma_p = _gghlite_dgsl_from_n(self->n, self->sigma_p, flags);
}


void _gghlite_params_set_ell_b(gghlite_params_t self) {
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma_p, 0)>0);

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


void _gghlite_params_set_sigma_s(gghlite_params_t self) {
  if (self->rerand_mask == 0) {
    /* if there is no re-randomisation there is not σ^* */
    mpfr_set_d(self->sigma_s, 1.0, MPFR_RNDN);
    return;
  } else if (self->rerand_mask != 1) {
    ggh_die("Re-randomisation at higher levels is not implemented yet.");
  }

  assert(self->kappa > 0);
	assert(self->gamma > 0);
  assert(self->lambda > 0);
  assert(self->n > 0);
  assert(mpfr_cmp_ui(self->sigma_p, 0)>0);
  assert(mpfr_cmp_ui(self->ell_b, 0)>0);
  assert(mpfr_cmp_ui(self->ell_g, 0)>0);

  mpfr_t sigma_s0, sigma_s1, tmp, pi, pow, eps;
  mpfr_init2(sigma_s0, _gghlite_prec(self));
  mpfr_init2(sigma_s1, _gghlite_prec(self));
  mpfr_init2(tmp, _gghlite_prec(self));
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_init2(pow, _gghlite_prec(self));
  mpfr_init2(eps, _gghlite_prec(self));
  mpfr_set_ui(eps, self->lambda, MPFR_RNDN);
  mpfr_log(eps, eps, MPFR_RNDN);
  mpfr_div_ui(eps, eps, self->kappa, MPFR_RNDN); // ε_d

  /* σ^* ≥ n^(3/2) · (σ')^2 · sqrt(8π/ε_d)/ℓ_b */

  mpfr_pow_ui(sigma_s0, self->sigma_p, 2, MPFR_RNDN); // σ^* := (σ')^2

  mpfr_div(tmp, pi, eps, MPFR_RNDN);  // π/ε_d
  mpfr_mul_ui(tmp, tmp, 8, MPFR_RNDN); // 8π/ε_d
  mpfr_sqrt(tmp, tmp, MPFR_RNDN); // sqrt(8π/ε_d)
  mpfr_mul(sigma_s0, sigma_s0, tmp, MPFR_RNDN); // σ^* := (σ')^2 · sqrt(8π/ε_d)
  mpfr_div(sigma_s0, sigma_s0, self->ell_b, MPFR_RNDN); // σ^* := (σ')^2 · sqrt(8π/ε_d)/ℓ_b

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_set_d(pow, 1.5, MPFR_RNDN);
  mpfr_pow(tmp, tmp, pow, MPFR_RNDN); // n^(3/2)

  mpfr_mul(sigma_s0, sigma_s0, tmp, MPFR_RNDN); // σ^* := n^(3/2) · (σ')^2 · sqrt(8π/ε_d)/ℓ_b

  /* `σ^* ≥ n^{1.5}·ℓ_g·σ'·\sqrt{2·log(4nε_ρ^{-1})/π}` */

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_div(tmp, tmp, eps, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 4, MPFR_RNDN);
  mpfr_log(tmp, tmp, MPFR_RNDN);
  mpfr_div(tmp, tmp, pi, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
  mpfr_sqrt(tmp, tmp, MPFR_RNDN);
  mpfr_mul(sigma_s1, tmp, self->sigma_p, MPFR_RNDN);
  mpfr_mul(sigma_s1, sigma_s1, self->ell_g, MPFR_RNDN);

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_set_d(pow, 1.5, MPFR_RNDN);
  mpfr_pow(tmp, tmp, pow, MPFR_RNDN); // n^(3/2)

  mpfr_mul(sigma_s1, sigma_s1, tmp,MPFR_RNDN);

  if (mpfr_cmp(sigma_s0, sigma_s1) >= 0)
    mpfr_set(self->sigma_s, sigma_s0, MPFR_RNDN);
  else
    mpfr_set(self->sigma_s, sigma_s1, MPFR_RNDN);

  if(self->flags & GGHLITE_FLAGS_GDDH_HARD) {
    /* GGH: 2^λ λ n^{4.5} κ vs. GGHLite n^{5.5} sqrt(κ) → 2^λ λ sqrt(κ)/n */
    mpfr_set_ui(tmp, 2, MPFR_RNDN);
    mpfr_pow_ui(tmp, tmp, self->lambda, MPFR_RNDN);
    mpfr_mul(self->sigma_s, self->sigma_s, tmp, MPFR_RNDN); // · 2^λ
    mpfr_mul_ui(self->sigma_s, self->sigma_s, self->lambda, MPFR_RNDN); // · λ 2^λ
    mpfr_sqrt_ui(tmp, self->kappa, MPFR_RNDN);
    mpfr_mul(self->sigma_s, self->sigma_s, tmp, MPFR_RNDN);  // · λ 2^λ sqrt(κ)
    mpfr_div_ui(self->sigma_s, self->sigma_s, self->n, MPFR_RNDN);
  }


  mpfr_clear(eps);
  mpfr_clear(pow);
  mpfr_clear(pi);
  mpfr_clear(tmp);
  mpfr_clear(sigma_s1);
  mpfr_clear(sigma_s0);
}

void _gghlite_params_set_D_sigma_s(gghlite_params_t self) {
  assert(self->n);
  assert(!mpfr_zero_p(self->sigma_s));
  const oz_flag_t flags = (self->flags & GGHLITE_FLAGS_QUIET) ? 0 : OZ_VERBOSE;
  self->D_sigma_s = _gghlite_dgsl_from_n(self->n, self->sigma_s, flags);
}

double gghlite_params_get_delta_0(const gghlite_params_t self) {

  /* Finding a short d·g in <g> */

  mpfr_t target_norm;
  mpfr_init2(target_norm, _gghlite_prec(self));

  mpfr_t tmp;
  mpfr_init2(tmp, _gghlite_prec(self));

  double delta_0_g = 0.0;

  if(self->rerand_mask) {
    mpfr_set(tmp, self->sigma_p, MPFR_RNDN);
    mpfr_pow_ui(tmp, tmp, self->kappa, MPFR_RNDN);
    mpfr_mul_ui(target_norm, tmp, 6, MPFR_RNDN); // 6⋅(σ')^κ
    mpfr_mul(target_norm, target_norm, self->sigma_s, MPFR_RNDN); // 6 ⋅ (σ')^κ ⋅ σ^*

    mpfr_set_ui(tmp, self->n, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    mpfr_pow_ui(tmp, tmp, 2*self->kappa + 3, MPFR_RNDN);
    mpfr_mul(target_norm, target_norm, tmp, MPFR_RNDN);
    mpfr_mul_ui(target_norm, target_norm, 2, MPFR_RNDN); // 6 ⋅ (σ')^κ ⋅ σ^* ⋅ sqrt(n)^(2κ+3)

    _gghlite_params_get_q_mpfr(tmp, self, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);

    mpfr_div(target_norm, tmp, target_norm, MPFR_RNDN);

    mpfr_set_ui(tmp, self->n, MPFR_RNDN);
    mpfr_ui_div(tmp, 1, tmp, MPFR_RNDN);

    mpfr_pow(tmp, target_norm, tmp, MPFR_RNDN);
    delta_0_g = mpfr_get_d(tmp, MPFR_RNDN);
  }

  /* finding (~b0,~b1) from (x0/x1) */

  _gghlite_params_get_q_mpfr(tmp, self, MPFR_RNDN);
  mpfr_sqrt(tmp, tmp, MPFR_RNDN);

  mpfr_set_ui(target_norm, 2, MPFR_RNDN);
  mpfr_mul_ui(target_norm, target_norm, self->n, MPFR_RNDN);
  mpfr_sqrt(target_norm, target_norm, MPFR_RNDN); //< sqrt(2n)
  mpfr_mul(target_norm, target_norm, self->sigma_p, MPFR_RNDN); //< sqrt(2n) · σ'

  if (self->rerand_mask) // if there are no rerandomisers we target random encodings
    mpfr_div(target_norm, target_norm, self->sigma, MPFR_RNDN); //< sqrt(2n) · σ'/σ

  mpfr_mul_d(target_norm, target_norm, 0.3, MPFR_RNDN); // TODO: Don't hardcode this
  mpfr_ui_div(target_norm, 1, target_norm, MPFR_RNDN); //< 1/(sqrt(2n) · σ'/σ · τ)
  mpfr_mul(target_norm, tmp, target_norm, MPFR_RNDN); //< sqrt(q)/(sqrt(2n) · σ'/σ · τ)

  mpfr_set_ui(tmp, self->n, MPFR_RNDN);
  mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
  mpfr_ui_div(tmp, 1, tmp, MPFR_RNDN);

  mpfr_pow(tmp, target_norm, tmp, MPFR_RNDN);
  double delta_0_ntru = mpfr_get_d(tmp, MPFR_RNDN);

  mpfr_clear(target_norm);
  mpfr_clear(tmp);

  if (self->rerand_mask) {
    if (delta_0_g > delta_0_ntru)
      return delta_0_g;
    else
      return delta_0_ntru;
  } else
    return delta_0_ntru;
}

double gghlite_params_cost_bkz_enum(const gghlite_params_t self) {
  const double delta_0 = gghlite_params_get_delta_0(self);
  if (delta_0 >= 1.0219) // LLL
    return 3*log2(self->n);
  const int k = _gghlite_k_from_delta(delta_0);
  const double r = _gghlite_repeat_from_n_k(self->n, k);
  return 0.270188776350190*k*log(k) - 1.0192050451318417*k + 16.10253135200765 + r;
}

double gghlite_params_cost_bkz_sieve(const gghlite_params_t self) {
  const double delta_0 = gghlite_params_get_delta_0(self);
  if (delta_0 >= 1.0219) // LLL
    return 3*log2(self->n);
  const int k = _gghlite_k_from_delta(delta_0);
  const double r = _gghlite_repeat_from_n_k(self->n, k);
  return 0.3366*k + 12.31 + r;
}


int gghlite_params_check_sec(const gghlite_params_t self) {
  double rt0 = gghlite_params_cost_bkz_enum(self);
  double rt1 = gghlite_params_cost_bkz_sieve(self);

  return ((rt0 >= self->lambda) && (rt1 >= self->lambda));
}


void gghlite_params_init_gamma(gghlite_params_t self, size_t lambda, size_t kappa, size_t gamma, uint64_t rerand_mask, gghlite_flag_t flags) {
  assert(lambda > 0);
  assert((kappa > 0) && (kappa <= KAPPA));
  assert((gamma > 0) && (gamma <= GAMMA));


  _gghlite_params_initzero(self, lambda, kappa, gamma);

  self->rerand_mask = rerand_mask;
  self->flags = flags;

  for(int log_n = 7; ; log_n++) {
    self->n = ((long)1)<<log_n;
    _gghlite_params_set_sigma(self);
    _gghlite_params_set_ell_g(self);
    _gghlite_params_set_ell(self);
    _gghlite_params_set_sigma_p(self);
    _gghlite_params_set_ell_b(self);
    _gghlite_params_set_sigma_s(self);
    _gghlite_params_set_q(self);

    if (gghlite_params_check_sec(self))
      break;
  }
}

void gghlite_params_clear(gghlite_params_t self) {
  fmpz_mod_poly_clear(self->pzt);

  const size_t bound = (gghlite_params_is_symmetric(self)) ? 1 : self->gamma;

  for(size_t i=0; i<bound; i++) {
    for(size_t k=0; k<self->kappa; k++) {
      if (gghlite_params_have_rerand(self, k)) {
        fmpz_mod_poly_clear(self->x[i][k][0]);
        fmpz_mod_poly_clear(self->x[i][k][1]);
      }
    }
    fmpz_mod_poly_clear(self->y[i]);
  }

  mpfr_clear(self->xi);
  dgsl_rot_mp_clear(self->D_sigma_s);
  mpfr_clear(self->sigma_s);
  mpfr_clear(self->ell_b);
  dgsl_rot_mp_clear(self->D_sigma_p);
  mpfr_clear(self->sigma_p);
  mpfr_clear(self->ell_g);
  mpfr_clear(self->sigma);
  fmpz_mod_poly_oz_ntt_precomp_clear(self->ntt);
  fmpz_clear(self->q);
}

void gghlite_params_ref(gghlite_params_t rop, gghlite_sk_t op) {
  memcpy(rop, op->params, sizeof(struct _gghlite_params_struct));
}
