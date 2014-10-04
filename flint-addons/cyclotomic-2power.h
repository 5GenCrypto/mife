#ifndef _CYCLOTOMIC_2POWER_H_
#define _CYCLOTOMIC_2POWER_H_
#include <stdio.h>
#include <mpfr.h>
#include "flint-addons/flint-addons.h"

static inline void fmpz_poly_init_cyc2pow_modulus(fmpz_poly_t f, const long n) {
  fmpz_poly_init(f);
  fmpz_poly_set_coeff_si(f, 0, 1);
  fmpz_poly_set_coeff_si(f, n, 1);
}


static inline void fmpq_poly_init_cyc2pow_modulus(fmpq_poly_t f, const long n) {
  fmpq_poly_init(f);
  fmpq_poly_set_coeff_si(f, 0, 1);
  fmpq_poly_set_coeff_si(f, n, 1);
}

static inline void _fmpq_poly_invert_mod_cnf2pow_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, const int n, mpfr_prec_t prec) {
  int i;

  if(f_inv == f)
    ggh_die("_fmpq_poly_invert_mod_cnf2pow_approx does not support parameter aliasing");

  fmpq_poly_zero(f_inv);

  fmpq_poly_t V;
  fmpq_poly_init(V);
  fmpq_poly_set(V, f);

  fmpq_poly_t U;
  fmpq_poly_init(U);

  /* F(x) = x^n+1 */
  fmpq_poly_t F;
  fmpq_poly_init_cyc2pow_modulus(F, n);

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
    for (i = 1 ; i <= deg ; i += 2) {
      /* negate odd coefficients */
      fmpq_poly_get_coeff_fmpq(tmp, V2, i);
      fmpq_neg(tmp,tmp);
      fmpq_poly_set_coeff_fmpq(V2, i, tmp);
    }

    /* compute Se, So */
    /* Set Se = (V(x) + V(-x))/2 */
    for (i = 0; i <= deg/2 ; i++) {
      fmpq_poly_get_coeff_fmpq(tmp, V2, 2*i);
      fmpq_poly_set_coeff_fmpq(Se, i, tmp);
    }
    fmpq_poly_truncate(Se,deg/2+1);

    /* Set So to the compressed (V(-x) - V(x))/2 */
    for (i = 0 ; i < (deg+1)/2 ; i++) {
      fmpq_poly_get_coeff_fmpq(tmp,V2,2*i+1);
      fmpq_poly_set_coeff_fmpq(So,i,tmp);
    }
    fmpq_poly_truncate(So,deg/2+1);

    /* V = V(x) * V(-x) mod f(x) */
    fmpq_poly_mulmod(V, V, V2, F);

#ifndef NDEBUG
    /* Sanity-check: verify that the odd coeffs in V are zero */
    for (i = 1; i < n ; i += 2) {
      fmpq_poly_get_coeff_fmpq(tmp,V,i);
      assert(fmpq_is_zero(tmp));
    }
#endif

    deg = fmpq_poly_degree(V);
    /* Compress the non-zero coeffs of V */
    for (i = 0; i <= deg/2 ; i ++) {
      fmpq_poly_get_coeff_fmpq(tmp,V,2*i);
      fmpq_poly_set_coeff_fmpq(V,i,tmp);
    }
    fmpq_poly_truncate(V,deg/2+1);

    _fmpq_poly_invert_mod_cnf2pow_approx(V2, V, n/2, prec);

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
    fmpq_poly_rem(f_inv,f_inv,F);

    /* truncate results on the precision required by algorithm */
    if (prec)
      fmpq_poly_truncate_prec(f_inv, prec);
  } else {
    if (fmpz_poly_is_zero(f))
      dgs_die("division by zero.");
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
  fmpq_poly_clear(F);
  fmpq_poly_clear(U);
  fmpq_poly_clear(V);
}

static inline void fmpq_poly_invert_mod_cnf2pow_approx(fmpq_poly_t rop, const fmpq_poly_t f, int n, mpfr_prec_t prec) {
  fmpq_poly_t tmp;
  fmpq_poly_init(tmp);
  fmpq_poly_t modulus;
  fmpq_poly_init_cyc2pow_modulus(modulus, n);

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

#ifndef GGHLITE_QUIET
  mpfr_t tmp_f;
  mpfr_init2(tmp_f, realprec);
  uint64_t t = ggh_walltime(0);
#endif
  for(long b=ceil(log2(prec)); ; b++) {
    _fmpq_poly_invert_mod_cnf2pow_approx(f_inv, f, n, (1<<b));
    fmpq_poly_mulmod(tmp, f_inv, f, modulus);

    fmpq_poly_get_coeff_fmpq(c, tmp, 0);
    fmpq_sub_si(c, c, 1);
    fmpq_poly_set_coeff_fmpq(tmp, 0, c);

    _fmpq_vec_2norm_mpfr(norm, tmp->coeffs, tmp->den, fmpq_poly_length(tmp));

#ifndef GGHLITE_QUIET
    mpfr_log2(tmp_f, norm, MPFR_RNDN);
    mpfr_fprintf(stderr, "\r   Computing f^-1::  b: %4d,     Δ=|f^-1·f-1|: %7.2Rf", b, tmp_f);
    mpfr_log2(tmp_f, bound, MPFR_RNDN);
    mpfr_fprintf(stderr, " <? %6.2Rf, ", tmp_f);
    fprintf(stderr, "t: %8.2fs",ggh_walltime(t)/1000000.0);
    fflush(0);
#endif

    if(mpfr_cmp(norm, bound) <= 0)
      break;
  }
#ifndef GGHLITE_QUIET
  fprintf(stderr, "\n");
  mpfr_clear(tmp_f);
#endif

  fmpq_poly_set(rop, f_inv);
  fmpq_poly_clear(f_inv);

  mpfr_clear(bound);
  fmpq_clear(c);
  mpfr_clear(norm);
  fmpq_poly_clear(modulus);
  fmpq_poly_clear(tmp);
}

// TODO: allow flag to choose method
static inline void fmpq_poly_sqrt_mod_cnf2pow_approx(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t boundbits) {

  fmpq_poly_t y;      fmpq_poly_init(y);
  fmpq_poly_t y_next; fmpq_poly_init(y_next);
  fmpq_poly_set(y, f);

  fmpq_poly_t z; fmpq_poly_init(z);
  fmpq_poly_t z_next; fmpq_poly_init(z_next);
  fmpq_poly_set_coeff_si(z, 0, 1);
  
  fmpq_poly_t f_approx; fmpq_poly_init(f_approx);

  fmpq_poly_t modulus;
  fmpq_poly_init_cyc2pow_modulus(modulus, n);

  mpfr_t norm;
  mpfr_init2(norm, prec);

  mpfr_t bound;
  mpfr_init2(bound, prec);
  mpfr_set_ui(bound, 1, MPFR_RNDN);
  mpfr_div_2exp(bound, bound, boundbits, MPFR_RNDN);

#ifndef GGHLITE_QUIET
  mpfr_t log_f;
  mpfr_init2(log_f, prec);
  uint64_t t = ggh_walltime(0);
#endif

  long i = 0;
  for(i=0; i<100; i++) {
#if 1
#pragma omp parallel sections
    {
#pragma omp section
      {
        _fmpq_poly_invert_mod_cnf2pow_approx(y_next, z, n, prec);
        fmpq_poly_add(y_next, y_next, y);
        fmpq_poly_scalar_div_si(y_next, y_next, 2);
      }
#pragma omp section
      {
        _fmpq_poly_invert_mod_cnf2pow_approx(z_next, y, n, 2*prec);
        fmpq_poly_add(z_next, z_next, z);
        fmpq_poly_scalar_div_si(z_next, z_next, 2);
      }
    }
    fmpq_poly_set(y, y_next);
    fmpq_poly_set(z, z_next);
#else
    /** The Babylonian Method is more forgiving about precision **/
    _fmpq_poly_invert_mod_cnf2pow_approx(y_next, y, n, prec);
    fmpq_poly_mulmod(y_next, f, y_next, modulus);
    fmpq_poly_add(y_next, y_next, y);
    fmpq_poly_scalar_div_si(y_next, y_next, 2);
    fmpq_poly_set(y, y_next);
#endif
    fmpq_poly_mulmod(f_approx, y, y, modulus);
    fmpq_poly_sub(f_approx, f_approx, f);
    _fmpq_vec_2norm_mpfr(norm, f_approx->coeffs, f_approx->den, f_approx->length);
    
#ifndef GGHLITE_QUIET
    mpfr_log2(log_f, norm, MPFR_RNDN);
    mpfr_fprintf(stderr, "\rComputing sqrt(Σ)::  i: %4d,  Δ=|sqrt(Σ)^2-Σ|: %7.2Rf", i, log_f);
    mpfr_log2(log_f, bound, MPFR_RNDN);
    mpfr_fprintf(stderr, " <? %6.2Rf, ", log_f);
    fprintf(stderr, "t: %8.2fs",ggh_walltime(t)/1000000.0);
    fflush(0);
#endif

    if(mpfr_cmp(norm, bound) < 0)
      break;
  }

  if(i==100)
    ggh_die("sqrt not converging fast enough, consider increasing precision or use Babylonian method.");
  
#ifndef GGHLITE_QUIET
  fprintf(stderr, "\n");
  mpfr_clear(log_f);
#endif
  fmpq_poly_set(f_sqrt, y);

  mpfr_clear(bound);
  mpfr_clear(norm);
  fmpq_poly_clear(modulus);
  fmpq_poly_clear(f_approx);
  fmpq_poly_clear(y_next);
  fmpq_poly_clear(y);
  fmpq_poly_clear(z_next);
  fmpq_poly_clear(z);
}

/*
  fT = f^T
*/

static inline void fmpq_poly_transpose_cnf2pow(fmpq_poly_t fT, const fmpq_poly_t f, const long n) {
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

#endif /* _CYCLOTOMIC_2POWER_H_ */
