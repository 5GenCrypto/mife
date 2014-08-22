#include <string.h>

#include "gghlite-internals.h"
#include "gghlite.h"

void _gghlite_pk_zero(gghlite_pk_t self) {
  memset(self, 0, sizeof(struct _gghlite_struct));
}

void _gghlite_pk_set_n_q(gghlite_pk_t self) {
  const int64_t lambda = self->lambda;
  const int64_t kappa  = self->kappa;
  int64_t n_base = kappa * lambda * ceil(log2(lambda));

  int64_t c     = 1;
  int64_t log_n = ceil(log2(c*n_base));
  int64_t n     = ((int64_t)(1))<<log_n;
  int64_t log_q = _gghlite_log_q(log_n, kappa);

  /* increase n until security check passes */
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

  /* increase q until it is probably prime and q % n == 1*/

  fmpz_t zeta;
  fmpz_init(zeta);
  fmpz_sub_ui(self->q, self->q, 1);
  fmpz_divexact_si(zeta, self->q, self->n);

  fmpz_t tmp;
  fmpz_init(tmp);
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


void _gghlite_pk_set_sigma(gghlite_pk_t self) {
  assert(self->n > 0);

  mpfr_t pi;
  mpfr_init2(pi, _gghlite_prec(self));
  mpfr_const_pi(pi, MPFR_RNDN);

  mpfr_init2(self->sigma, _gghlite_prec(self));
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

void _gghlite_pk_set_ell_g(gghlite_pk_t self) {
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

void _gghlite_pk_set_sigma_p(gghlite_pk_t self) {
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


  mpfr_init2(self->sigma_p, _gghlite_prec(self));
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

void _gghlite_pk_set_D_sigma_p(gghlite_pk_t self) {
  assert(!self->n);
  assert(!mpfr_zero_p(self->sigma_p));
  self->D_sigma_p = _gghlite_gpv_from_n(self->n, self->sigma_p);
}


void _gghlite_pk_set_ell_b(gghlite_pk_t self) {
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


/**
   CONSTRAINTS:

   #. `σ^* ≥ n^{1.5}·ℓ_g·σ'·\sqrt{2·log(4nε_ρ^{-1})/π}`, cf. [LSS14]_, p.17, Eq. (8)
   #. `σ^* ≥ n^{1.5}·(σ')²\sqrt{8πε_d^{-1}}/ℓ_b`, cf. [LSS14]_, p.17, Eq. (9) with
   `εₑ^{-1} = O(log λ/κ)`.
*/

void _gghlite_pk_set_sigma_s(gghlite_pk_t self) {
  assert(self->kappa > 0);
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

  mpfr_init2(self->sigma_s, _gghlite_prec(self));

  /* `σ^* ≥ n^{1.5}·ℓ_g·σ'·\sqrt{2·log(4nε_ρ^{-1})/π}` */
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

  mpfr_clear(eps);
  mpfr_clear(pow);
  mpfr_clear(pi);
  mpfr_clear(tmp);
  mpfr_clear(sigma_s1);
  mpfr_clear(sigma_s0);
}

void _gghlite_pk_set_D_sigma_s(gghlite_pk_t self) {
  assert(!self->n);
  assert(!mpfr_zero_p(self->sigma_s));
  self->D_sigma_s = _gghlite_gpv_from_n(self->n, self->sigma_s);
}

void gghlite_pk_init_params(gghlite_pk_t self, size_t lambda, size_t kappa, uint64_t rerand_mask) {
  assert(lambda > 0);
  assert((kappa > 0) && (kappa <= KAPPA));

  _gghlite_pk_zero(self);
  
  self->lambda = lambda;
  self->kappa = kappa;
  self->rerand_mask = rerand_mask;

  _gghlite_pk_set_n_q(self);
  _gghlite_pk_set_sigma(self);
  _gghlite_pk_set_ell_g(self);
  _gghlite_pk_set_sigma_p(self);
  _gghlite_pk_set_ell_b(self);
  _gghlite_pk_set_sigma_s(self);
  _gghlite_pk_set_modulus(self);
}

void gghlite_pk_clear(gghlite_pk_t self) {
  fmpz_mod_poly_clear(self->pzt);

  fmpz_mod_poly_clear(self->y);

  for(long k=0; k<self->kappa; k++) {
    if(self->rerand_mask && (1ULL)<<k) {
      fmpz_mod_poly_clear(self->x[k][0]);
      fmpz_mod_poly_clear(self->x[k][1]);
    }
  }

  fmpz_mod_poly_clear(self->modulus);
  gpv_mp_clear(self->D_sigma_s);
  mpfr_clear(self->sigma_s);
  mpfr_clear(self->ell_b);
  gpv_mp_clear(self->D_sigma_p);
  mpfr_clear(self->sigma_p);
  mpfr_clear(self->ell_g);
  mpfr_clear(self->sigma);
  fmpz_clear(self->q);
}

void gghlite_print_params(const gghlite_pk_t self) {
  const long lambda = self->lambda;
  const long kappa = self->kappa;
  const long n = self->n;
  printf("         λ: %7ld\n",lambda);
  printf("         k: %7ld\n",kappa);
  printf("         n: %7ld\n",n);
  printf("   log₂(q): %7ld (check: %d)\n", fmpz_sizeinbase(self->q, 2), gghlite_check_sec(self));
  printf("   log₂(σ): %7.1f dp: (%7.1f)\n", log2(mpfr_get_d(self->sigma,   MPFR_RNDN)), log2(_gghlite_sigma(n)));
  printf(" log₂(ℓ_g): %7.1f dp: (%7.1f)\n", log2(mpfr_get_d(self->ell_g,   MPFR_RNDN)), log2(_gghlite_ell_g(n)));
  printf("  log₂(σ'): %7.1f dp: (%7.1f)\n", log2(mpfr_get_d(self->sigma_p, MPFR_RNDN)), log2(_gghlite_sigma_p(n)));
  printf(" log₂(ℓ_b): %7.1f dp: (%7.1f)\n", log2(mpfr_get_d(self->ell_b,   MPFR_RNDN)), log2(_gghlite_ell_b(n)));
  printf(" log₂(σ^*): %7.1f dp: (%7.1f)\n", log2(mpfr_get_d(self->sigma_s, MPFR_RNDN)), log2(_gghlite_sigma_s(n, lambda, kappa)));
}
