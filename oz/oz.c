#include "flint-addons.h"
#include "util.h"
#include "oz.h"

void fmpq_poly_oz_ideal_norm(fmpq_t norm, const fmpq_poly_t f, const long n, const mpfr_prec_t prec) {
  if (prec == 0) {
    fmpq_poly_t modulus;
    fmpq_poly_oz_init_modulus(modulus, n);
    // TODO: patch FLINT to accept a bound on the norm
    fmpq_poly_resultant(norm, f, modulus);
    fmpq_poly_clear(modulus);
  } else {
    mpfr_t norm_f;
    mpfr_init2(norm_f, prec);
    fmpq_poly_eucl_norm_mpfr(norm_f, f, MPFR_RNDN);
    mpfr_pow_ui(norm_f, norm_f, n, MPFR_RNDN);
    fmpq_set_mpfr(norm, norm_f, MPFR_RNDN);
    mpfr_clear(norm_f);
  }
}

void fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const long n, const mpfr_prec_t prec) {
  if (prec == 0) {    
    fmpz_poly_t modulus;
    fmpz_poly_oz_init_modulus(modulus, n);
    fmpz_poly_ideal_norm(norm, f, modulus);
    fmpz_poly_clear(modulus);
  } else {
    mpfr_t norm_f;
    mpfr_init2(norm_f, prec);
    fmpz_poly_eucl_norm_mpfr(norm_f, f, MPFR_RNDN);
    mpfr_pow_ui(norm_f, norm_f, n, MPFR_RNDN);
    mpz_t norm_z;
    mpz_init(norm_z);
    mpfr_get_z(norm_z, norm_f, MPFR_RNDN);
    fmpz_set_mpz(norm, norm_z);
    mpfr_clear(norm_f);
    mpz_clear(norm_z);
  }
}

void _fmpq_poly_oz_invert_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, const int n, const mpfr_prec_t prec) {
  if(f_inv == f)
    oz_die("_fmpq_poly_invert_mod_cnf2pow_approx does not support parameter aliasing");

  fmpq_poly_zero(f_inv);

  fmpq_poly_t V;
  fmpq_poly_init(V);
  fmpq_poly_set(V, f);

  fmpq_poly_t U;
  fmpq_poly_init(U);

  /* F(x) = x^n+1 */
  fmpq_poly_t V2;
  fmpq_poly_init(V2);

  fmpq_poly_t Se;
  fmpq_poly_init(Se);
  fmpq_poly_t So;
  fmpq_poly_init(So);

  fmpq_t tmp;
  fmpq_init(tmp);

  int deg;
  if (n > 1) {
    /* set V2(x) = V(-x) */
    fmpq_poly_set(V2, V);
    deg = fmpq_poly_degree(V2);
    for (int i = 1 ; i <= deg ; i += 2) {
      /* negate odd coefficients */
      fmpq_poly_get_coeff_fmpq(tmp, V2, i);
      fmpq_neg(tmp,tmp);
      fmpq_poly_set_coeff_fmpq(V2, i, tmp);
    }

    /* compute Se, So */
    /* Set Se = (V(x) + V(-x))/2 */
    for (int i = 0; i <= deg/2 ; i++) {
      fmpq_poly_get_coeff_fmpq(tmp, V2, 2*i);
      fmpq_poly_set_coeff_fmpq(Se, i, tmp);
    }
    fmpq_poly_truncate(Se,deg/2+1);

    /* Set So to the compressed (V(-x) - V(x))/2 */
    for (int i = 0 ; i < (deg+1)/2 ; i++) {
      fmpq_poly_get_coeff_fmpq(tmp,V2,2*i+1);
      fmpq_poly_set_coeff_fmpq(So,i,tmp);
    }
    fmpq_poly_truncate(So,deg/2+1);

    /* V = V(x) * V(-x) mod f(x) */
    fmpq_poly_oz_mul(V, V, V2, n);

    deg = fmpq_poly_degree(V);
    /* Compress the non-zero coeffs of V */
    for (int i = 0; i <= deg/2 ; i ++) {
      fmpq_poly_get_coeff_fmpq(tmp,V,2*i);
      fmpq_poly_set_coeff_fmpq(V,i,tmp);
    }
    fmpq_poly_truncate(V,deg/2+1);

    _fmpq_poly_oz_invert_approx(V2, V, n/2, prec);

    deg = fmpq_poly_degree(V2);
    /* Te=G*Se, To = G*So */
    fmpq_poly_mul(V,V2,Se);
    fmpq_poly_mul(U,V2,So);

    deg = fmpq_poly_degree(V);

    for (int i = 0 ; i <= deg ; i ++) {
      fmpq_poly_get_coeff_fmpq(tmp,V,i);
      fmpq_poly_set_coeff_fmpq(f_inv,2*i,tmp);
    }
    deg = fmpq_poly_degree(U);
    for (int i = 0 ; i <= deg ; i ++) {
      fmpq_poly_get_coeff_fmpq(tmp,U,i);
      fmpq_poly_set_coeff_fmpq(f_inv,2*i+1,tmp);
    }
    if (deg>= 0)
      fmpq_poly_truncate(f_inv,2*deg+2);
    fmpq_poly_oz_rem(f_inv,f_inv,n);

    /* truncate results on the precision required by algorithm */
    if (prec)
      fmpq_poly_truncate_prec(f_inv, prec);
  } else {
    if (fmpz_poly_is_zero(f))
      oz_die("division by zero.");
    f_inv->length = 1;
    fmpq_poly_get_coeff_fmpq(tmp, f, 0);
    fmpq_inv(tmp, tmp);
    fmpq_poly_set_coeff_fmpq(f_inv, 0, tmp);
  }

  /* free memory */
  fmpq_clear(tmp);
  fmpq_poly_clear(So);
  fmpq_poly_clear(Se);
  fmpq_poly_clear(V2);
  fmpq_poly_clear(U);
  fmpq_poly_clear(V);
}

void fmpq_poly_oz_invert_approx(fmpq_poly_t rop, const fmpq_poly_t f, const long n,
                                const mpfr_prec_t prec, const uint64_t flags) {
  fmpq_poly_t tmp;
  fmpq_poly_init(tmp);

  fmpq_poly_t f_inv;
  fmpq_poly_init(f_inv);

  const mpfr_prec_t realprec = (2*prec < 53) ? 53 : 2*prec;

  mpfr_t norm;
  mpfr_init2(norm, realprec);

  fmpq_t c;
  fmpq_init(c);

  /** we check for distance < 2^-prec with precision 2*prec **/
  mpfr_t bound;
  mpfr_init2(bound, realprec);
  mpfr_set_ui(bound, 1, MPFR_RNDN);
  mpfr_div_2exp(bound, bound, prec, MPFR_RNDN);


  mpfr_t tmp_f;
  mpfr_init2(tmp_f, realprec);
  uint64_t t = oz_walltime(0);

  for(long b=ceil(log2(prec)); ; b++) {
    _fmpq_poly_oz_invert_approx(f_inv, f, n, (1<<b));
    fmpq_poly_oz_mul(tmp, f_inv, f, n);

    fmpq_poly_get_coeff_fmpq(c, tmp, 0);
    fmpq_sub_si(c, c, 1);
    fmpq_poly_set_coeff_fmpq(tmp, 0, c);

    fmpq_poly_eucl_norm_mpfr(norm, tmp, MPFR_RNDN);

    if (flags & OZ_VERBOSE) {
      mpfr_log2(tmp_f, norm, MPFR_RNDN);
      mpfr_fprintf(stderr, "   Computing f^-1::  b: %4d,     Δ=|f^-1·f-1|: %7.2Rf", b, tmp_f);
      mpfr_log2(tmp_f, bound, MPFR_RNDN);
      mpfr_fprintf(stderr, " <? %6.2Rf, ", tmp_f);
      fprintf(stderr, "t: %8.2fs\n",oz_walltime(t)/1000000.0);
      fflush(0);
    }

    if(mpfr_cmp(norm, bound) <= 0)
      break;
  }
  if(flags & OZ_VERBOSE)
    fprintf(stderr, "\n");
  mpfr_clear(tmp_f);

  fmpq_poly_set(rop, f_inv);
  fmpq_poly_clear(f_inv);

  mpfr_clear(bound);
  fmpq_clear(c);
  mpfr_clear(norm);
  fmpq_poly_clear(tmp);
}

int fmpq_poly_oz_sqrt_approx(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init) {
  fmpq_poly_t y;
  fmpq_poly_t y_next;

  fmpq_poly_init(y);
  fmpq_poly_init(y_next);
  if (init) {
    fmpq_poly_set(y, init);
  } else {
    fmpq_poly_set(y, f);
  }

  mpfr_t f_norm;
  mpfr_init2(f_norm, prec);
  fmpq_poly_eucl_norm_mpfr(f_norm, f, MPFR_RNDN);

  fmpq_poly_t z;
  fmpq_poly_t z_next;

  if (!(flags & OZ_BABYLONIAN)) {
    fmpq_poly_init(z);
    fmpq_poly_init(z_next);
    if (init) {
      // z = y/x
      _fmpq_poly_oz_invert_approx(z, f, n, prec);
      fmpq_poly_oz_mul(z, z, y, n);
    } else {
      fmpq_poly_set_coeff_si(z, 0, 1);
    }
  }

  mpfr_t norm;
  mpfr_init2(norm, prec);

  /* We scale by about |det(y) · det(z)|^(-1/(2n)) to make it converge faster */
  mpfr_t gamma;
  mpfr_init2(gamma, prec);

  fmpq_t gamma_q;
  fmpq_init(gamma_q);


  if (!init) {
    if((flags & OZ_BABYLONIAN)) {
      /* det(y)^(-1/n) */
      fmpq_poly_eucl_norm_mpfr(norm, y, MPFR_RNDN);
      mpfr_set(gamma, norm, MPFR_RNDN);
      mpfr_ui_div(gamma, 1, gamma, MPFR_RNDN);
      fmpq_set_mpfr(gamma_q, gamma, MPFR_RNDN);
      fmpq_poly_scalar_mul_fmpq(y, y, gamma_q);
    } else {
      /* det(y) */
      fmpq_poly_oz_ideal_norm(gamma_q, y, n, prec);
      fmpq_get_mpfr(gamma, gamma_q, MPFR_RNDN);

      /* det(y) · det(z) */
      fmpq_poly_oz_ideal_norm(gamma_q, z, n, prec);
      fmpq_get_mpfr(norm, gamma_q, MPFR_RNDN);
      mpfr_mul(gamma, gamma, norm, MPFR_RNDN);

      /* (det(y) · det(z))^(-1/(2n)) */
      mpfr_root(gamma, gamma, 2*n, MPFR_RNDN);
      mpfr_ui_div(gamma, 1, gamma, MPFR_RNDN);

      fmpq_set_mpfr(gamma_q, gamma, MPFR_RNDN);
      fmpq_poly_scalar_mul_fmpq(y, y, gamma_q);
      fmpq_poly_scalar_mul_fmpq(z, z, gamma_q);
    }
  }
  fmpq_clear(gamma_q);
  mpfr_clear(gamma);

  fmpq_poly_t f_approx;
  fmpq_poly_init(f_approx);

  mpfr_t prev_norm;
  mpfr_init2(prev_norm, prec);

  mpfr_t bound;
  mpfr_init2(bound, prec);
  mpfr_set_ui(bound, 1, MPFR_RNDN);
  mpfr_div_2exp(bound, bound, prec_bound, MPFR_RNDN);

  mpfr_t log_f;
  mpfr_init2(log_f, prec);
  uint64_t t = oz_walltime(0);

  int r = 0;
  for(long k=0; ; k++) {
    if(!(flags & OZ_BABYLONIAN)) {
#pragma omp parallel sections
      {
#pragma omp section
        {
          _fmpq_poly_oz_invert_approx(y_next, z, n, prec);
          fmpq_poly_add(y_next, y_next, y);
          fmpq_poly_scalar_div_si(y_next, y_next, 2);
        }
#pragma omp section
        {
          _fmpq_poly_oz_invert_approx(z_next, y, n, prec);
          fmpq_poly_add(z_next, z_next, z);
          fmpq_poly_scalar_div_si(z_next, z_next, 2);
        }
      }
      fmpq_poly_set(y, y_next);
      fmpq_poly_set(z, z_next);
    } else {
      _fmpq_poly_oz_invert_approx(y_next, y, n, prec);
      fmpq_poly_oz_mul(y_next, f, y_next, n);
      fmpq_poly_add(y_next, y_next, y);
      fmpq_poly_scalar_div_si(y_next, y_next, 2);
      fmpq_poly_set(y, y_next);
    }

    fmpq_poly_oz_mul(f_approx, y, y, n);
    fmpq_poly_sub(f_approx, f_approx, f);
    fmpq_poly_eucl_norm_mpfr(norm, f_approx, MPFR_RNDN);
    mpfr_div(norm, norm, f_norm, MPFR_RNDN);

    if(flags & OZ_VERBOSE) {
      mpfr_log2(log_f, norm, MPFR_RNDN);
      mpfr_fprintf(stderr, "Computing sqrt(Σ)::  k: %4d,  Δ=|sqrt(Σ)^2-Σ|: %7.2Rf", k, log_f);
      mpfr_log2(log_f, bound, MPFR_RNDN);
      mpfr_fprintf(stderr, " <? %6.2Rf, ", log_f);
      fprintf(stderr, "t: %8.2fs\n", oz_walltime(t)/1000000.0);
      fflush(0);
    }

    if(mpfr_cmp(norm, bound) < 0) {
      r = 0;
      break;
    }

    mpfr_div_ui(prev_norm, prev_norm, 2, MPFR_RNDN);
    if (k>0 && mpfr_cmp(norm, prev_norm) >= 0) {
      /*  we don't converge any more */
      r = -1;
      break;
    }
    mpfr_set(prev_norm, norm, MPFR_RNDN);
  }

  /* if(flags & OZ_VERBOSE) */
  /*   fprintf(stderr, "\n"); */
  mpfr_clear(log_f);

  fmpq_poly_set(f_sqrt, y);

  mpfr_clear(bound);
  mpfr_clear(f_norm);
  mpfr_clear(norm);
  mpfr_clear(prev_norm);
  fmpq_poly_clear(f_approx);
  fmpq_poly_clear(y_next);
  fmpq_poly_clear(y);
  if (!(flags & OZ_BABYLONIAN)) {
    fmpq_poly_clear(z_next);
    fmpq_poly_clear(z);
  }
  return r;
}

void fmpz_poly_oz_rem(fmpz_poly_t rem, const fmpz_poly_t f, const long n) {
  fmpz_poly_set(rem, f);
  fmpz_t lead; fmpz_init(lead);
  fmpz_t coef; fmpz_init(coef);
  for(long d=fmpz_poly_degree(rem); d>=n; d--) {
    fmpz_poly_get_coeff_fmpz(lead, rem, d);
    fmpz_poly_get_coeff_fmpz(coef, rem, d-n);
    fmpz_sub(coef, coef, lead);
    fmpz_poly_set_coeff_si(rem, d, 0);
    fmpz_poly_set_coeff_fmpz(rem, d-n, coef);
  }
  fmpz_clear(coef);
  fmpz_clear(lead);
}

void fmpq_poly_oz_rem(fmpq_poly_t rem, const fmpq_poly_t f, const long n) {
  // TODO: this is faster than the above naive method because FLINT avoids GCD
  // calls
  fmpq_poly_t modulus;
  fmpq_poly_oz_init_modulus(modulus, n);
  fmpq_poly_rem(rem, f, modulus);
  fmpq_poly_clear(modulus);
  return;
}

void fmpz_poly_oz_conjugate(fmpz_poly_t fT, const fmpz_poly_t f, const long n) {
  fmpz_poly_zero(fT);
  fmpz_t t0; fmpz_init(t0);
  fmpz_t t1; fmpz_init(t1);

  fmpz_poly_get_coeff_fmpz(t0, f, 0);
  fmpz_poly_set_coeff_fmpz(fT, 0, t0);

  for(int i=1; i<n/2; i++) {
    fmpz_poly_get_coeff_fmpz(t0, f,   i);
    fmpz_poly_get_coeff_fmpz(t1, f, n-i);
    fmpz_neg(t0, t0);
    fmpz_neg(t1, t1);
    fmpz_poly_set_coeff_fmpz(fT, n-i, t0);
    fmpz_poly_set_coeff_fmpz(fT, i, t1);
  }
  fmpz_poly_get_coeff_fmpz(t0, f, n/2);
  fmpz_neg(t0, t0);
  fmpz_poly_set_coeff_fmpz(fT, n/2, t0);

  fmpz_clear(t1);
  fmpz_clear(t0);
}

void fmpq_poly_oz_conjugate(fmpq_poly_t fT, const fmpq_poly_t f, const long n) {
  fmpq_poly_zero(fT);
  fmpq_t t0; fmpq_init(t0);
  fmpq_t t1; fmpq_init(t1);

  fmpq_poly_get_coeff_fmpq(t0, f, 0);
  fmpq_poly_set_coeff_fmpq(fT, 0, t0);

  for(int i=1; i<n/2; i++) {
    fmpq_poly_get_coeff_fmpq(t0, f,   i);
    fmpq_poly_get_coeff_fmpq(t1, f, n-i);
    fmpq_neg(t0, t0);
    fmpq_neg(t1, t1);
    fmpq_poly_set_coeff_fmpq(fT, n-i, t0);
    fmpq_poly_set_coeff_fmpq(fT, i, t1);
  }
  fmpq_poly_get_coeff_fmpq(t0, f, n/2);
  fmpq_neg(t0, t0);
  fmpq_poly_set_coeff_fmpq(fT, n/2, t0);

  fmpq_clear(t1);
  fmpq_clear(t0);
}
