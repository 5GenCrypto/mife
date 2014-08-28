#include <string.h>

#include "gghlite-internals.h"
#include "gghlite.h"

void _gghlite_zero(gghlite_t self) {
  memset(self, 0, sizeof(struct _gghlite_struct));
}

void _gghlite_set_y(gghlite_t self) {
  assert(!fmpz_poly_is_zero(self->a));
  assert(!fmpz_is_zero(self->pk->q));

  fmpz_mod_poly_t z_inv;
  fmpz_mod_poly_init(z_inv, self->pk->q);
  fmpz_mod_poly_invert_mod(z_inv, self->z, self->pk->modulus);

  fmpz_mod_poly_t tmp;
  fmpz_mod_poly_init(tmp, self->pk->q);
  fmpz_mod_poly_set_fmpz_poly(tmp, self->a);
  fmpz_mod_poly_mulmod(tmp, tmp, z_inv, self->pk->modulus);
  fmpz_mod_poly_truncate(tmp, self->pk->n);

  fmpz_mod_poly_init(self->pk->y, self->pk->q);
  fmpz_mod_poly_set(self->pk->y, tmp);

  fmpz_mod_poly_clear(tmp);
  fmpz_mod_poly_clear(z_inv);
}

void _gghlite_sample_a(gghlite_t self, flint_rand_t randstate) {
  assert(self->pk);
  assert(self->pk->n);
  assert(!fmpz_poly_is_zero(self->g));

  const long n = self->pk->n;
  mpfr_t *c = _mpfr_vec_init(n, _gghlite_prec(self->pk));
  mpfr_set_ui(c[0], 1, MPFR_RNDN);
  gpv_mp_t *D = _gghlite_gpv_from_poly(self->g, self->pk->sigma_p, c, GPV_INLATTICE);
  _mpfr_vec_clear(c, n);

  fmpz_poly_init(self->a);
  fmpz_poly_sample_D(self->a, D, randstate);
  gpv_mp_clear(D);
}

void _gghlite_set_x(gghlite_t self) {
  fmpz_mod_poly_t z_inv;
  fmpz_mod_poly_init(z_inv, self->pk->q);
  fmpz_mod_poly_invert_mod(z_inv, self->z, self->pk->modulus);

  fmpz_mod_poly_t tmp;
  fmpz_mod_poly_init(tmp, self->pk->q);

  fmpz_mod_poly_t acc;
  fmpz_mod_poly_init(acc, self->pk->q);

  for(long k=0; k<self->pk->kappa; k++) {
    if(fmpz_poly_is_zero(self->b[k][0]) || fmpz_poly_is_zero(self->b[k][1]))
      continue;

    fmpz_mod_poly_powmod_ui_binexp(acc, z_inv, k+1, self->pk->modulus);

    fmpz_mod_poly_init(self->pk->x[k][0], self->pk->q);
    fmpz_mod_poly_set_fmpz_poly(tmp, self->b[k][0]);
    fmpz_mod_poly_mulmod(self->pk->x[k][0], tmp, acc, self->pk->modulus);
    fmpz_mod_poly_truncate(self->pk->x[k][0], self->pk->n);

    fmpz_mod_poly_init(self->pk->x[k][1], self->pk->q);
    fmpz_mod_poly_set_fmpz_poly(tmp, self->b[k][1]);
    fmpz_mod_poly_mulmod(self->pk->x[k][1], tmp, acc, self->pk->modulus);
    fmpz_mod_poly_truncate(self->pk->x[k][0], self->pk->n);

  }
  fmpz_mod_poly_clear(tmp);
  fmpz_mod_poly_clear(acc);
  fmpz_mod_poly_clear(z_inv);
}

void _gghlite_sample_b(gghlite_t self, flint_rand_t randstate) {
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
    fmpz_poly_init(self->b[k][0]);
    fmpz_poly_init(self->b[k][1]);

    if (!(self->pk->rerand_mask & (1ULL)<<k))
      continue;

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
  assert(fmpz_mod_poly_degree(self->pk->modulus) == self->pk->n);

  fmpz_mod_poly_t z_kappa;
  fmpz_mod_poly_init(z_kappa, self->pk->q);
  fmpz_mod_poly_set(z_kappa, self->z);
  fmpz_mod_poly_powmod_ui_binexp(z_kappa, z_kappa, self->pk->kappa, self->pk->modulus);

  fmpz_mod_poly_t g_inv;
  fmpz_mod_poly_init(g_inv, self->pk->q);
  fmpz_mod_poly_set_fmpq_poly(g_inv, self->g_inv);

  fmpz_mod_poly_t pzt;
  fmpz_mod_poly_init(pzt, self->pk->q);
  fmpz_mod_poly_mulmod(pzt, z_kappa, g_inv, self->pk->modulus);

  fmpz_mod_poly_t h;
  fmpz_mod_poly_init(h, self->pk->q);
  fmpz_mod_poly_set_fmpz_poly(h, self->h);

  fmpz_mod_poly_mulmod(pzt, pzt, h, self->pk->modulus);
  fmpz_mod_poly_truncate(pzt, self->pk->n);
  fmpz_mod_poly_init(self->pk->pzt, self->pk->q);
  fmpz_mod_poly_set(self->pk->pzt, pzt);

  fmpz_mod_poly_clear(h);
  fmpz_mod_poly_clear(pzt);
  fmpz_mod_poly_clear(z_kappa);
  fmpz_mod_poly_clear(g_inv);
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

  fmpz_poly_t modulus;
  fmpz_poly_init2(modulus, self->pk->n);
  fmpz_poly_set_fmpz_mod_poly(modulus, self->pk->modulus);

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
    fmpz_poly_invert_mod_fmpq(self->g_inv, self->g, modulus);
    if (!_gghlite_g_inv_check(self->pk, self->g_inv)) {
      fail[2]++;
      continue;
    }
    break;
  }

  fmpz_poly_clear(modulus);
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

void gghlite_init_instance(gghlite_t self, flint_rand_t randstate) {
  assert(self->pk->lambda);
  assert(self->pk->kappa);

  _gghlite_sample_g(self, randstate);
  _gghlite_sample_z(self, randstate);
  _gghlite_sample_h(self, randstate);
  _gghlite_sample_b(self, randstate);
  _gghlite_sample_a(self, randstate);

  _gghlite_set_y(self);
  _gghlite_set_x(self);
  _gghlite_set_pzt(self);
  _gghlite_pk_set_D_sigma_p(self->pk);
  _gghlite_pk_set_D_sigma_s(self->pk);
}

void gghlite_init(gghlite_t self, const size_t lambda, const size_t kappa,
                  const uint64_t rerand_mask, flint_rand_t randstate) {
  _gghlite_zero(self);
  gghlite_pk_init_params(self->pk, lambda, kappa, rerand_mask);
  gghlite_init_instance(self, randstate);
}

void gghlite_clear(gghlite_t self, int clear_pk) {
  fmpz_poly_clear(self->a);
  for(long k=0; k<self->pk->kappa; k++) {
    if(self->pk->rerand_mask && (1ULL)<<k) {
      fmpz_poly_clear(self->b[k][0]);
      fmpz_poly_clear(self->b[k][1]);
    }
  }
  fmpz_poly_clear(self->h);
  fmpz_mod_poly_clear(self->z);
  fmpz_poly_clear(self->g);

  if (clear_pk)
    gghlite_pk_clear(self->pk);
}
