#ifndef _CYCLOTOMIC_2POWER_H_
#define _CYCLOTOMIC_2POWER_H_

#include "flint-addons/flint-addons.h"

static inline void fmpq_poly_init_cyc2power_modulus(fmpq_poly_t f, const long n) {
  fmpq_poly_init(f);
  fmpq_poly_set_coeff_si(f, 0, 1);
  fmpq_poly_set_coeff_si(f, n, 1);
}

static inline void _fmpq_poly_invert_mod_cnf2pow_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, const int n, mpfr_prec_t prec) {
  int i;

  fmpq_poly_t V;
  fmpq_poly_init(V);
  fmpq_poly_set(V, f);

  fmpq_poly_t U;
  fmpq_poly_init(U);

  /* F(x) = x^n+1 */
  fmpq_poly_t F;
  fmpq_poly_init_cyc2power_modulus(F, n);

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

    /* V2 = Inverse3(V,q,n/2) */
    _fmpq_poly_invert_mod_cnf2pow_approx(V2,V,n/2, prec);

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
    fmpq_poly_truncate_prec(f_inv, prec);
  } else {
    f_inv->length = 1;
    fmpq_poly_get_coeff_fmpq(tmp, f, 0);
    fmpq_inv(tmp, tmp);
    fmpq_poly_set_coeff_fmpq(f_inv, 0, tmp);
  }

  /* free memory */
  fmpq_clear(tmp);
  fmpq_poly_clear(Se);
  fmpq_poly_clear(So);
  fmpq_poly_clear(V2);
  fmpq_poly_clear(F);
  fmpq_poly_clear(U);
  fmpq_poly_clear(V);
}

static inline void fmpq_poly_invert_mod_cnf2pow_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, int n, mp_bitcnt_t prec) {
  fmpq_poly_t tmp;
  fmpq_poly_init(tmp);
  fmpq_poly_t modulus;
  fmpq_poly_init_cyc2power_modulus(modulus, n);
  
  mpfr_t norm;
  mpfr_init2(norm, 2*prec);

  fmpq_t c;
  fmpq_init(c);

  /** we check for distance < 2^-prec with precision 2*prec **/
  mpfr_t bound;
  mpfr_init2(bound, 2*prec);
  mpfr_set_ui(bound, 1, MPFR_RNDN);
  mpfr_div_2exp(bound, bound, prec, MPFR_RNDN);

  for(long b=ceil(log2(prec)); ; b++) {
    _fmpq_poly_invert_mod_cnf2pow_approx(f_inv, f, n, (1<<b));
    fmpq_poly_mulmod(tmp, f_inv, f, modulus);

    fmpq_poly_get_coeff_fmpq(c, tmp, 0);
    fmpq_sub_si(c, c, 1);
    fmpq_poly_set_coeff_fmpq(tmp, 0, c);

    _fmpq_vec_2norm_mpfr(norm, tmp->coeffs, tmp->den, fmpq_poly_length(tmp));
    if(mpfr_cmp(norm, bound) <= 0)
      break;
  }
  printf("\n");
  fmpq_poly_clear(tmp);
  fmpq_poly_clear(modulus);
  mpfr_clear(norm);
  mpfr_clear(bound);
  fmpq_clear(c);
}


static inline void fmpq_poly_sqrt_mod_cnf2power_approx(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec) {

  fmpq_poly_t y; fmpq_poly_init(y); fmpq_poly_set(y, f);
  fmpq_poly_t z; fmpq_poly_init(z); fmpq_poly_set_coeff_si(z, 0, 1);

  fmpq_poly_t y_next; fmpq_poly_init(y_next);
  fmpq_poly_t z_next; fmpq_poly_init(z_next);

  fmpq_poly_t tmp; fmpq_poly_init(tmp);

  fmpq_poly_t modulus;
  fmpq_poly_init_cyc2power_modulus(modulus, n);
  
  mpfr_t norm;
  mpfr_init2(norm, 2*prec);

  mpfr_t tmp_f;
  mpfr_init2(tmp_f, prec);
  
  mpfr_t bound;
  mpfr_init2(bound, 2*prec);
  mpfr_set_ui(bound, 1, MPFR_RNDN);
  mpfr_div_2exp(bound, bound, prec, MPFR_RNDN);
  
  for(long i=0; ; i++) {
    uint64_t t = ggh_walltime(0);
    fmpq_poly_set(y_next, z);
    _fmpq_poly_invert_mod_cnf2pow_approx(y_next, y_next, n, 2*prec);
    fmpq_poly_add(y_next, y_next, y);
    fmpq_poly_scalar_div_si(y_next, y_next, 2);

    fmpq_poly_set(z_next, y);
    _fmpq_poly_invert_mod_cnf2pow_approx(z_next, z_next, n, 2*prec);
    fmpq_poly_add(z_next, z_next, z);
    fmpq_poly_scalar_div_si(z_next, z_next, 2);
    
    fmpq_poly_set(y, y_next);
    fmpq_poly_set(z, z_next);

    fmpq_poly_mulmod(tmp, y, y, modulus);
    fmpq_poly_sub(tmp, tmp, f);
    _fmpq_vec_2norm_mpfr(norm, tmp->coeffs, tmp->den, n);
    mpfr_log2(tmp_f, norm, MPFR_RNDN);
    mpfr_printf("i: %3d |Δ|: %10.4Rf", i, tmp_f);
    mpfr_log2(tmp_f, bound, MPFR_RNDN);
    mpfr_printf(" <? %.Rf ", tmp_f);
    printf(" t: %8.2fs\n",ggh_walltime(t)/1000000.0);
    fflush(0);
    if(mpfr_cmp(norm, bound) < 0)
      break;
  }

  mpfr_clear(tmp_f);
  fmpq_poly_set(f_sqrt, y);
  mpfr_clear(bound);
  mpfr_clear(norm);
  fmpq_poly_clear(modulus);
  fmpq_poly_clear(tmp);
  fmpq_poly_clear(y); fmpq_poly_clear(y_next);
  fmpq_poly_clear(z); fmpq_poly_clear(z_next);  
}

/*
  rop = f · f^T
*/

void fmpq_poly_mul_mod_cnf2pow_transpose(fmpq_poly_t rop, const fmpq_poly_t f, const long n) {
  fmpq_poly_t fT;
  fmpq_poly_init(fT);

  fmpq_t t0; fmpq_init(t0);

  fmpq_poly_t modulus;
  fmpq_poly_init_cyc2power_modulus(modulus, n);
  
  fmpq_poly_get_coeff_fmpq(t0, f, 0);
  fmpq_poly_set_coeff_fmpq(fT, 0, t0);

  for(int i=1; i<n/2; i++) {
    fmpq_poly_get_coeff_fmpq(t0, f,   i);
    fmpq_neg(t0, t0);
    fmpq_poly_set_coeff_fmpq(fT, n-i, t0);
    
    fmpq_poly_get_coeff_fmpq(t0, f, n-i);
    fmpq_neg(t0, t0);
    fmpq_poly_set_coeff_fmpq(fT, i, t0);
  }
  fmpq_poly_mulmod(rop, f, fT, modulus);

  fmpq_poly_clear(fT);
  fmpq_poly_clear(modulus);
  fmpq_clear(t0);
}

#endif /* _CYCLOTOMIC_2POWER_H_ */
