#include <stdarg.h>
#include <stdio.h>
#include "gghlite.h"
#include "gghlite-internals.h"

void gghlite_enc_init(gghlite_enc_t op, const gghlite_params_t self) {
  assert(!fmpz_is_zero(self->q));
  fmpz_mod_poly_init2(op, self->q, self->n);
}

void gghlite_enc_rerand(gghlite_enc_t rop, const gghlite_params_t self, const gghlite_enc_t op,
                        size_t k, size_t i, flint_rand_t randstate) {
  assert(self->D_sigma_s);

  if (k!=1)
    ggh_die("Re-randomisation for k!=1 not implemented");

  gghlite_enc_t tmp;
  gghlite_enc_init(tmp, self);
  fmpz_mod_poly_set(rop, op);

  assert(!fmpz_mod_poly_is_zero(self->x[i][k-1][0]) && !fmpz_mod_poly_is_zero(self->x[i][k-1][1]));

  for(long j=0; j<2; j++) {
    fmpz_mod_poly_sample_D(tmp, self->D_sigma_s, randstate);
    fmpz_mod_poly_oz_ntt_enc(tmp, tmp, self->ntt);
    fmpz_mod_poly_oz_ntt_mul(tmp, tmp, self->x[i][k-1][j], self->n);
    fmpz_mod_poly_oz_ntt_add(rop, rop, tmp);
  }

  fmpz_mod_poly_clear(tmp);
}

void gghlite_enc_raise(gghlite_enc_t rop, const gghlite_params_t self, const gghlite_enc_t op,
                       size_t l, size_t k, size_t i, int rerand, flint_rand_t randstate) {
  assert(k <= l);
  assert(l <= self->kappa);

  if(!gghlite_params_is_symmetric(self) && (l>1))
    ggh_die("Raising to higher levels than 1 not supported. Instead, multiply by the right combination of y_i.");

  assert(!fmpz_mod_poly_is_zero(self->y[i]));

  if (l>k) {
    gghlite_enc_t yk;
    gghlite_enc_init(yk, self);
    fmpz_mod_poly_oz_ntt_pow_ui(yk, self->y[i], l-k, self->n);
    fmpz_mod_poly_oz_ntt_mul(rop, op, yk, self->n);
    fmpz_mod_poly_clear(yk);
  }

  if(rerand && l>0) {
    gghlite_enc_rerand(rop, self, rop, l, i, randstate);
  }
}


void gghlite_enc_sample(gghlite_enc_t rop, gghlite_params_t self, size_t k, size_t i, flint_rand_t randstate) {
  assert(self->kappa);
  assert(k <= self->kappa);

  fmpz_mod_poly_sample_D(rop, self->D_sigma_p, randstate);
  fmpz_mod_poly_oz_ntt_enc(rop, rop, self->ntt);
  if (k == 0)
    return;
  gghlite_enc_raise(rop, self, rop, k, 0, i, 1, randstate);
}

void gghlite_enc_set_gghlite_clr(gghlite_enc_t rop, const gghlite_sk_t self, const gghlite_clr_t f,
                                 const size_t k, const size_t i, const int rerand,
                                 flint_rand_t randstate) {

  fmpz_poly_t t_i;  fmpz_poly_init(t_i);
  fmpz_poly_t t_o;  fmpz_poly_init(t_o);
  mpfr_t norm_i; mpfr_init2(norm_i, _gghlite_prec(self->params));
  mpfr_t norm_o; mpfr_init2(norm_o, _gghlite_prec(self->params));

  fmpz_poly_set(t_i, f);

  if (fmpz_poly_degree(f) == 0) {
    /* if f is a constant, we write it as f == ∑ 2^(i⋅b) f_i and deal with each chunk individually */
    const mp_bitcnt_t b = self->params->lambda * self->params->kappa * self->params->lambda * self->params->kappa;
    _fmpz_poly_oz_rem_small_fmpz_split(t_o, f->coeffs, self->g, self->params->n, self->g_inv, b);
  } else {
    fmpz_poly_set(t_o, t_i);
  }

  if(self->params->flags & GGHLITE_FLAGS_VERBOSE)
    printf("|f|: %10.1f, |g|: %10.1f, |f%%g|: %10.1f\n",
           fmpz_poly_eucl_norm_log2(t_i),
           fmpz_poly_eucl_norm_log2(self->g),
           fmpz_poly_eucl_norm_log2(t_o));

  /* the precision of g_inv might not be sufficient to do this in one step, hence, we repeat until
     the result does not improve any more*/

  fmpq_poly_t g_inv; fmpq_poly_init(g_inv);
  fmpq_poly_set(g_inv, self->g_inv);

  do {
    fmpz_poly_set(t_i, t_o);
    fmpz_poly_eucl_norm_mpfr(norm_i, t_i, MPFR_RNDN);
    fmpq_poly_truncate_prec(g_inv, fmpz_poly_eucl_norm_log2(t_i)/2);
    _fmpz_poly_oz_rem_small(t_o, t_i, self->g, self->params->n, g_inv);
    fmpz_poly_eucl_norm_mpfr(norm_o, t_o, MPFR_RNDN);

    if(self->params->flags & GGHLITE_FLAGS_VERBOSE)
      printf("|f|: %10.1f, |g|: %10.1f, |f%%g|: %10.1f\n",
             fmpz_poly_eucl_norm_log2(t_i),
             fmpz_poly_eucl_norm_log2(self->g),
             fmpz_poly_eucl_norm_log2(t_o));

  } while (mpfr_cmp(norm_o, norm_i) < 0);
  fmpz_poly_set(t_o, t_i);

  mpfr_clear(norm_i);
  mpfr_clear(norm_o);

  if (rerand)
    dgsl_rot_mp_call_plus_fmpz_poly(t_o, self->D_g, t_o, randstate->gmp_state);

  // encode at level zero
  fmpz_mod_poly_oz_ntt_enc_fmpz_poly(rop, t_o, self->params->ntt);

  fmpz_poly_clear(t_i);
  fmpz_poly_clear(t_o);

  // encode at level k
  if(k > 0) {
    if(!gghlite_sk_is_symmetric(self) && (k>1))
      ggh_die("Raising to higher levels than 1 not supported. Instead, multiply by the right combination of y_i.");

    for(size_t j=0; j<k; j++) // divide by z_i^k
      fmpz_mod_poly_oz_ntt_mul(rop, rop, self->z_inv[i], self->params->n);
  }
}


void gghlite_enc_extract(gghlite_clr_t rop, const gghlite_params_t self, const gghlite_enc_t op) {
  _gghlite_enc_extract_raw(rop, self, op);
  long logq = fmpz_sizeinbase(self->q,2) - 1;
  for(long i=0; i<fmpz_poly_length(rop); i++) {
    fmpz_tdiv_q_2exp(rop->coeffs+i,rop->coeffs+i, logq-self->ell+1);
  }
  _fmpz_poly_normalise(rop);
}

int gghlite_enc_is_zero(const gghlite_params_t self, const fmpz_mod_poly_t op) {
  gghlite_clr_t t;
  gghlite_clr_init(t);
  _gghlite_enc_extract_raw(t, self, op);

  mpfr_t norm;
  mpfr_init2(norm, _gghlite_prec(self));
  fmpz_poly_eucl_norm_mpfr(norm, t, MPFR_RNDN);

  mpfr_t bound;
  mpfr_init2(bound, _gghlite_prec(self));
  mpz_t q_z;
  mpz_init(q_z);
  fmpz_get_mpz(q_z, self->q);
  mpfr_set_z(bound, q_z, MPFR_RNDN);
  mpz_clear(q_z);

  mpfr_t ex;
  mpfr_init2(ex, _gghlite_prec(self));
  mpfr_set_ui(ex, 1, MPFR_RNDN);
  mpfr_sub(ex, ex, self->xi, MPFR_RNDN);
  mpfr_pow(bound, bound, ex, MPFR_RNDN);
  mpfr_clear(ex);

  int r = mpfr_cmp(norm, bound);

  mpfr_clear(bound);
  mpfr_clear(norm);
  gghlite_clr_clear(t);

  if (r<=0)
    return 1;
  else
    return 0;
}


void gghlite_params_print(const gghlite_params_t self) {
  assert(self->lambda);
  assert(self->kappa);
  assert(self->n);
  assert(!fmpz_is_zero(self->q));

  const long lambda = self->lambda;
  const long kappa = self->kappa;
  const long n = self->n;
  printf("symmetric: %9d\n", gghlite_params_is_symmetric(self));
  printf("        λ: %9ld,          k: %9ld\n",lambda, kappa);
  printf("        n: %9ld,        δ_0: %9.6f\n",n, gghlite_params_get_delta_0(self));
  printf("log(t_en): %9.2f,  log(t_sv): %9.2f\n", gghlite_params_cost_bkz_enum(self), gghlite_params_cost_bkz_sieve(self));
  printf("   log(q): %9ld,          ξ: %9.4f\n", fmpz_sizeinbase(self->q, 2), mpfr_get_d(self->xi, MPFR_RNDN));
  printf("   log(σ): %9.2f,   log(ℓ_g): %9.2f\n", log2(mpfr_get_d(self->sigma,   MPFR_RNDN)), log2(mpfr_get_d(self->ell_g,   MPFR_RNDN)));
  printf("  log(σ'): %9.2f,   log(ℓ_b): %9.2f\n", log2(mpfr_get_d(self->sigma_p, MPFR_RNDN)), log2(mpfr_get_d(self->ell_b,   MPFR_RNDN)));
  printf(" log(σ^*): %9.2f,   \n", log2(mpfr_get_d(self->sigma_s, MPFR_RNDN)));


  mpfr_t enc;
  mpfr_init2(enc, _gghlite_prec(self));

  _gghlite_params_get_q_mpfr(enc, self, MPFR_RNDN);
  mpfr_log2(enc, enc, MPFR_RNDN);
  mpfr_mul_ui(enc, enc, n, MPFR_RNDN);

  const char *units[3] = {"KB","MB","GB"};

  double sd = mpfr_get_d(enc, MPFR_RNDN)/8.0;
  int i;
  for(i=0; i<3; i++) {
    if (sd < 1024.0)
      break;
    sd = sd/1024;
  }
  printf("    |enc|: %6.1f %s,",sd, units[i-1]);

  mpfr_t par;
  mpfr_init2(par, _gghlite_prec(self));
  mpfr_set(par, enc, MPFR_RNDN);

  int count = 0;
  for(size_t k=0; k<self->kappa; k++)
    if (gghlite_params_have_rerand(self, k))
      count++;

  if (gghlite_params_is_symmetric(self))
    mpfr_mul_ui(par, par, count*2 + 1 + 1, MPFR_RNDN);
  else
    mpfr_mul_ui(par, par, count*3 + 1 + 1, MPFR_RNDN);

  mpfr_get_d(par, MPFR_RNDN);
  sd = mpfr_get_d(par, MPFR_RNDN)/8.0;
  for(i=0; i<3; i++) {
    if (sd < 1024.0)
      break;
    sd = sd/1024;
  }
  printf("      |par|: %6.1f %s\n", sd, units[i-1]);

  mpfr_clear(enc);
  mpfr_clear(par);
}

