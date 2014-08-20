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

void _gghlite_set_x(gghlite_t self) {
  fmpz_mod_poly_t z_inv;
  fmpz_mod_poly_init(z_inv, self->pk->q);

  fmpz_mod_poly_t acc;
  fmpz_mod_poly_init(acc, self->pk->q);

  fmpz_mod_poly_t modulus;
  fmpz_mod_poly_init(modulus, self->pk->q);
  fmpz_mod_poly_set_fmpz_poly(modulus, self->pk->cyclotomic_polynomial);

  fmpz_mod_poly_invert_mod(z_inv, self->z, modulus);

  fmpz_mod_poly_t tmp;
  fmpz_mod_poly_init(tmp, self->pk->q);

  for(long k=0; k<self->pk->kappa; k++) {
    if(fmpz_poly_is_zero(self->b[k][0]) || fmpz_poly_is_zero(self->b[k][1]))
      continue;

    fmpz_mod_poly_powmod_ui_binexp(acc, z_inv, k+1, modulus);

    fmpz_mod_poly_init(self->pk->x[k][0], self->pk->q);
    fmpz_mod_poly_set_fmpz_poly(tmp, self->b[k][0]);
    fmpz_mod_poly_mulmod(self->pk->x[k][0], tmp, acc, modulus);

    fmpz_mod_poly_init(self->pk->x[k][1], self->pk->q);
    fmpz_mod_poly_set_fmpz_poly(tmp, self->b[k][1]);
    fmpz_mod_poly_mulmod(self->pk->x[k][1], tmp, acc, modulus);
  }
  fmpz_mod_poly_clear(tmp);
  fmpz_mod_poly_clear(modulus);
  fmpz_mod_poly_clear(acc);
  fmpz_mod_poly_clear(z_inv);
}

void _gghlite_sample_b(gghlite_t self, uint64_t rerand_mask, flint_rand_t randstate) {
  assert(self->pk);
  assert(self->pk->n);
  assert(!fmpz_poly_is_zero(self->g));

  const long n = self->pk->n;
  gpv_mp_t *D = _gghlite_gpv_from_poly(self->g, self->pk->sigma_p, NULL, GPV_INLATTICE);

  mpfr_t sigma_n;
  mpfr_init2(sigma_n, _gghlite_prec(self->pk));

  mpfr_t norm;
  mpfr_init2(norm, _gghlite_prec(self->pk));

  mpfr_t sqrtn_sigma_p;
  mpfr_init2(sqrtn_sigma_p, _gghlite_prec(self->pk));
  mpfr_set_si(sqrtn_sigma_p, n, MPFR_RNDN);
  mpfr_sqrt(sqrtn_sigma_p, sqrtn_sigma_p, MPFR_RNDN);
  mpfr_mul(sqrtn_sigma_p, sqrtn_sigma_p, self->pk->sigma_p, MPFR_RNDN);

  fmpz_poly_t B;
  fmpz_poly_init2(B, 2*n);

  long fail[3] = {0,0,0};

#ifndef GGHLITE_CHECK_SIGMA_N
  printf("WARNING: Not checking that `σ_n(rot(B^(k))) < ℓ_b`.\n");
#endif

  for(long k=0; k<self->pk->kappa; k++) {
    if (!(rerand_mask & (1ULL)<<k))
      continue;


    fmpz_poly_init(self->b[k][0]);
    fmpz_poly_init(self->b[k][1]);

    while(1) {
      printf("\rk: %2ld :: !i: %4ld, !s: %4ld, !n: %4ld",k+1, fail[0], fail[1], fail[2]);
      fflush(0);
      fmpz_poly_sample_D(self->b[k][0], D, randstate);
      fmpz_poly_sample_D(self->b[k][1], D, randstate);

      if (fmpz_poly_ideal_subset(self->g, self->b[k][0], self->b[k][1]) != 0) {
        fail[0]++;
        continue;
      }

      _fmpz_vec_set(B->coeffs+0, self->b[k][0]->coeffs, n);
      _fmpz_vec_set(B->coeffs+n, self->b[k][1]->coeffs, n);


#ifdef GGHLITE_CHECK_SIGMA_N
      fmpz_poly_rot_basis_sigma_n(sigma_n, B);
      if (mpfr_cmp(sigma_n, self->pk->ell_b) < 0)
        fail[1]++;
        continue;
#endif

      _fmpz_vec_2norm_mpfr(norm, B->coeffs, 2*n);
      if (mpfr_cmp(norm, sqrtn_sigma_p) > 0) {
        fail[2]++;
        continue;
      }

      break;
    }
    printf("\n");
  }

  mpfr_clear(sqrtn_sigma_p);
  mpfr_clear(sigma_n);
  fmpz_poly_clear(B);
  gpv_mp_clear(D);
}

void _gghlite_set_pzt(gghlite_t self) {
  assert(self->pk);
  assert(self->pk->n);
  assert(fmpz_cmp_ui(self->pk->q,0)>0);
  assert(!fmpz_mod_poly_is_zero(self->z));
  assert(!fmpz_poly_is_zero(self->h));
  assert(fmpz_poly_degree(self->pk->cyclotomic_polynomial) == self->pk->n);

  fmpz_mod_poly_t modulus;
  fmpz_mod_poly_init(modulus, self->pk->q);
  fmpz_mod_poly_set_fmpz_poly(modulus, self->pk->cyclotomic_polynomial);

  fmpz_mod_poly_t z_kappa;
  fmpz_mod_poly_init(z_kappa, self->pk->q);
  fmpz_mod_poly_set(z_kappa, self->z);
  fmpz_mod_poly_powmod_ui_binexp(z_kappa, z_kappa, self->pk->kappa, modulus);

  fmpz_mod_poly_t g_inv;
  fmpz_mod_poly_init(g_inv, self->pk->q);
  fmpz_mod_poly_set_fmpq_poly(g_inv, self->g_inv);

  fmpz_mod_poly_t pzt;
  fmpz_mod_poly_init(pzt, self->pk->q);
  fmpz_mod_poly_mulmod(pzt, z_kappa, g_inv, modulus);


  fmpz_mod_poly_t h;
  fmpz_mod_poly_init(h, self->pk->q);
  fmpz_mod_poly_set_fmpz_poly(h, self->h);

  fmpz_mod_poly_mulmod(pzt, pzt, h, modulus);
  fmpz_mod_poly_init(self->pk->pzt, self->pk->q);
  fmpz_mod_poly_set(self->pk->pzt, pzt);

  fmpz_mod_poly_clear(h);
  fmpz_mod_poly_clear(pzt);
  fmpz_mod_poly_clear(z_kappa);
  fmpz_mod_poly_clear(g_inv);
  fmpz_mod_poly_clear(modulus);
}

void _gghlite_sample_g(gghlite_t self, flint_rand_t randstate) {
  assert(self->pk);
  assert(self->pk->n);

  fmpz_poly_init(self->g);
  fmpq_poly_init(self->g_inv);


  mpfr_t g_inv_norm;
  mpfr_init2(g_inv_norm, fmpz_sizeinbase(self->pk->q,2));

  mpfr_t norm;
  mpfr_init2(norm, mpfr_get_prec(self->pk->sigma));

  mpfr_t sqrtn_sigma;
  mpfr_init2(sqrtn_sigma, mpfr_get_prec(self->pk->sigma));
  mpfr_set_si(sqrtn_sigma, self->pk->n, MPFR_RNDN);
  mpfr_sqrt(sqrtn_sigma, sqrtn_sigma, MPFR_RNDN);
  mpfr_mul(sqrtn_sigma, sqrtn_sigma, self->pk->sigma, MPFR_RNDN);

  gpv_mp_t *D = _gghlite_gpv_from_n(self->pk->n, self->pk->sigma);

  long fail[3] = {0,0,0};

#ifndef GGHLITE_CHECK_PRIMALITY
  printf("WARNING: Not checking that `<g>` is prime.\n");
#endif

  while(1) {
    printf("\r    g :: !n: %4ld, !p: %4ld, !i: %4ld",fail[0], fail[1], fail[2]);
    fflush(0);

    fmpz_poly_sample_D(self->g, D, randstate);

    _fmpz_vec_2norm_mpfr(norm, self->g->coeffs, self->pk->n);
    if(mpfr_cmp(norm, sqrtn_sigma)>0) {
      fail[0]++;
      continue;
    }
#ifdef GGHLITE_CHECK_PRIMALITY
    /* 1. check if prime */
    if (!fmpz_poly_ideal_is_probaprime(self->g)) {
      fail[1]++;
      continue;
    }
#endif
    /* 2. check norm of inverse */
    fmpz_poly_invert_mod_fmpq(self->g_inv, self->g, self->pk->cyclotomic_polynomial);
    if (!_gghlite_g_inv_check(self->pk, self->g_inv)) {
      fail[2]++;
      continue;
    }
    break;
  }
  printf("\n");
  mpfr_clear(norm);
  mpfr_clear(sqrtn_sigma);
  mpfr_clear(g_inv_norm);
  gpv_mp_clear(D);
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

void _gghlite_sample_z(gghlite_t self, flint_rand_t randstate) {
  assert(self->pk);
  assert(self->pk->n);
  assert(fmpz_cmp_ui(self->pk->q,0)>0);

  fmpz_mod_poly_init(self->z, self->pk->q);
  fmpz_mod_poly_randtest(self->z, randstate, self->pk->n);
}

void _gghlite_set_n_q(gghlite_pk_t self) {
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


void _gghlite_set_sigma(gghlite_pk_t self) {
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

void _gghlite_set_sigma_p(gghlite_pk_t self) {
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


/**
   CONSTRAINTS:

   #. `σ^* ≥ n^{1.5}·ℓ_g·σ'·\sqrt{2·log(4nε_ρ^{-1})/π}`, cf. [LSS14]_, p.17, Eq. (8)
   #. `σ^* ≥ n^{1.5}·(σ')²\sqrt{8πε_d^{-1}}/ℓ_b`, cf. [LSS14]_, p.17, Eq. (9) with
   `εₑ^{-1} = O(log λ/κ)`.
*/

void _gghlite_set_sigma_s(gghlite_pk_t self) {
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
