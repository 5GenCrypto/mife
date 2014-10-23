#include "flint-addons.h"
#include "util.h"
#include "oz.h"

static inline mp_bitcnt_t _fmpq_poly_oz_ideal_norm_bound(const fmpq_poly_t f, const long n) {
  mp_bitcnt_t bits1 = FLINT_ABS(_fmpz_vec_max_bits(f->coeffs, f->length)); 
  mp_bitcnt_t bound = 2*n*FLINT_BIT_COUNT((20*n + 26)/27) + 3;
  bound += (n - 1) + n*bits1;
  return bound;
}

void fmpq_poly_oz_ideal_norm(fmpq_t norm, const fmpq_poly_t f, const long n, const mpfr_prec_t prec) {
  fmpq_poly_t modulus;
  
  if (prec == 0) {
    mp_bitcnt_t bound = _fmpq_poly_oz_ideal_norm_bound(f, n);
    fmpq_poly_oz_init_modulus(modulus, n);
    fmpq_poly_resultant_modular_bound(norm, f, modulus, bound);
    fmpq_poly_clear(modulus);

  } else if  (prec < 0) {

    mpfr_t norm_f;
    mpfr_init2(norm_f, -prec);
    fmpq_poly_eucl_norm_mpfr(norm_f, f, MPFR_RNDN);
    mpfr_pow_ui(norm_f, norm_f, n, MPFR_RNDN);
    fmpq_set_mpfr(norm, norm_f, MPFR_RNDN);
    mpfr_clear(norm_f);

  } else {

    fmpq_poly_oz_init_modulus(modulus, n);
    fmpq_poly_t f_trunc;
    fmpq_poly_init(f_trunc);
    fmpq_poly_set(f_trunc, f);
    fmpq_poly_truncate(f_trunc, prec);

    mp_bitcnt_t bound = _fmpq_poly_oz_ideal_norm_bound(f_trunc, n);
    fmpq_poly_resultant_modular_bound(norm, f_trunc, modulus, bound);
    
    fmpq_poly_clear(modulus);
    fmpq_poly_clear(f_trunc);
  }
}

void fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const long n, const mpfr_prec_t prec) {
  mp_bitcnt_t bits = FLINT_ABS(_fmpz_vec_max_bits(f->coeffs, f->length));
  mp_bitcnt_t bound = f->length * (bits + n_clog(f->length, 2));  
  if (prec == 0) {
    fmpz_poly_t modulus;
    fmpz_poly_oz_init_modulus(modulus, n);
    fmpz_poly_resultant_modular_bound(norm, f, modulus, bound);
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
      fprintf(stderr, " <? %4ld, ", -prec);
      fprintf(stderr, "t: %8.2fs\n",oz_walltime(t)/1000000.0);
      fflush(0);
    }

    if(mpfr_cmp(norm, bound) <= 0)
      break;
  }
  /* if(flags & OZ_VERBOSE) */
  /*   fprintf(stderr, "\n"); */
  mpfr_clear(tmp_f);

  fmpq_poly_set(rop, f_inv);
  fmpq_poly_clear(f_inv);

  mpfr_clear(bound);
  fmpq_clear(c);
  mpfr_clear(norm);
  fmpq_poly_clear(tmp);
}

int _fmpq_poly_oz_sqrt_approx_break(mpfr_t norm, const fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t bound, const mpfr_prec_t prec) {
  fmpq_poly_t f_approx;
  fmpq_poly_init(f_approx);
  fmpq_poly_oz_mul(f_approx, f_sqrt, f_sqrt, n);
  fmpq_poly_sub(f_approx, f_approx, f);
  fmpq_poly_eucl_norm_mpfr(norm, f_approx, MPFR_RNDN);

  mpfr_t f_norm;
  mpfr_init2(f_norm, prec);
  fmpq_poly_eucl_norm_mpfr(f_norm, f, MPFR_RNDN);
  mpfr_div(norm, norm, f_norm, MPFR_RNDN);

  int r = 0;
  if(mpfr_cmp_si_2exp(norm, 1, -bound) < 0)
    r = 1;
  mpfr_clear(f_norm);
  fmpq_poly_clear(f_approx);
  return r;
}

void _fmpq_poly_oz_sqrt_approx_scale(fmpq_poly_t y, fmpq_poly_t z, const long n, const mpfr_prec_t prec) {    
  /* We scale by about |det(y) · det(z)|^(-1/(2n)) to make it converge faster */
  mpfr_t gamma;
  mpfr_init2(gamma, prec);

  mpfr_t tmp;
  mpfr_init2(tmp, prec);
  
  fmpq_t gamma_q;
  fmpq_init(gamma_q);

  /* det(y) */
  fmpq_poly_oz_ideal_norm(gamma_q, y, n, 10);
  fmpq_get_mpfr(gamma, gamma_q, MPFR_RNDN);

  /* det(y) · det(z) */
  fmpq_poly_oz_ideal_norm(gamma_q, z, n, 10);
  fmpq_get_mpfr(tmp, gamma_q, MPFR_RNDN);
  mpfr_mul(gamma, gamma, tmp, MPFR_RNDN);

  /* (det(y) · det(z))^(-1/(2n)) */
  mpfr_root(gamma, gamma, 2*n, MPFR_RNDN);
  mpfr_ui_div(gamma, 1, gamma, MPFR_RNDN);
  mpfr_abs(gamma, gamma, MPFR_RNDN);

  
  fmpq_set_mpfr(gamma_q, gamma, MPFR_RNDN);
  fmpq_poly_scalar_mul_fmpq(y, y, gamma_q);
  fmpq_poly_scalar_mul_fmpq(z, z, gamma_q);

  fmpq_clear(gamma_q);
  mpfr_clear(gamma);
  mpfr_clear(tmp);
}

int fmpq_poly_oz_sqrt_approx_babylonian(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init) {
  fmpq_poly_t y;      fmpq_poly_init(y);
  fmpq_poly_t y_next; fmpq_poly_init(y_next);

  mpfr_t norm;      mpfr_init2(norm, prec);
  mpfr_t prev_norm; mpfr_init2(prev_norm, prec);
  
  if (init) {
    fmpq_poly_set(y, init);
  } else {
    fmpq_poly_set(y, f);
  }

  mpfr_t log_f;
  mpfr_init2(log_f, prec);

  uint64_t t = oz_walltime(0);
  int r = 0;

  for(long k=0; ; k++) {
    _fmpq_poly_oz_invert_approx(y_next, y, n, prec);
    fmpq_poly_oz_mul(y_next, f, y_next, n);
    fmpq_poly_add(y_next, y_next, y);
    fmpq_poly_scalar_div_si(y_next, y_next, 2);
    fmpq_poly_set(y, y_next);

    r = _fmpq_poly_oz_sqrt_approx_break(norm, y, f, n, prec_bound, prec);
    
    if(flags & OZ_VERBOSE) {
      mpfr_log2(log_f, norm, MPFR_RNDN);
      mpfr_fprintf(stderr, "Computing sqrt(Σ)::  k: %4d,  Δ=|sqrt(Σ)^2-Σ|: %7.2Rf", k, log_f);
      fprintf(stderr, " <? %4ld, ", -prec_bound);
      fprintf(stderr, "t: %8.2fs\n", oz_walltime(t)/1000000.0);
      fflush(0);
    }

    if(r) {
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
  mpfr_clear(log_f);
  fmpq_poly_set(f_sqrt, y);
  mpfr_clear(norm);
  mpfr_clear(prev_norm);
  fmpq_poly_clear(y_next);
  fmpq_poly_clear(y);
  return r;
}

int fmpq_poly_oz_sqrt_approx_db(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init) {
  fmpq_poly_t y;       fmpq_poly_init(y);
  fmpq_poly_t y_next;  fmpq_poly_init(y_next);
  fmpq_poly_t z;       fmpq_poly_init(z);
  fmpq_poly_t z_next;  fmpq_poly_init(z_next);

  mpfr_t norm;       mpfr_init2(norm, prec);
  mpfr_t prev_norm;  mpfr_init2(prev_norm, prec);
  mpfr_t log_f;      mpfr_init2(log_f, prec);

  uint64_t t = oz_walltime(0);

  if (init) {
    // z = y/x
    fmpq_poly_set(y, init);
    _fmpq_poly_oz_invert_approx(z, f, n, prec);
    fmpq_poly_oz_mul(z, z, y, n);
  } else {
    fmpq_poly_set(y, f);
    fmpq_poly_set_coeff_si(z, 0, 1);
  }


  int r = 0;
  for(long k=0; ; k++) {
    if (k == 0)
      _fmpq_poly_oz_sqrt_approx_scale(y, z, n, prec);

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

    r = _fmpq_poly_oz_sqrt_approx_break(norm, y, f, n, prec_bound, prec);

    if(flags & OZ_VERBOSE) {
      mpfr_log2(log_f, norm, MPFR_RNDN);
      mpfr_fprintf(stderr, "Computing sqrt(Σ)::  k: %4d,  Δ=|sqrt(Σ)^2-Σ|: %7.2Rf", k, log_f);
      fprintf(stderr, " <? %4ld, ", -prec_bound);
      fprintf(stderr, "t: %8.2fs\n", oz_walltime(t)/1000000.0);
      fflush(0);
    }

    if(r) {
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

  mpfr_clear(log_f);
  fmpq_poly_set(f_sqrt, y);
  mpfr_clear(norm);
  mpfr_clear(prev_norm);
  fmpq_poly_clear(y_next);
  fmpq_poly_clear(y);
  fmpq_poly_clear(z_next);
  fmpq_poly_clear(z);
  return r;
}

int fmpq_poly_oz_sqrt_approx_pade(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const int p, const mpfr_prec_t prec, const mpfr_prec_t bound, uint64_t flags, const fmpq_poly_t init) {
  fmpq_poly_t y;       fmpq_poly_init(y);
  fmpq_poly_t y_next;  fmpq_poly_init(y_next);
  fmpq_poly_t z;       fmpq_poly_init(z);
  fmpq_poly_t z_next;  fmpq_poly_init(z_next);

  mpfr_t norm;      mpfr_init2(norm, prec);
  mpfr_t prev_norm; mpfr_init2(prev_norm, prec);
  mpfr_t log_f;     mpfr_init2(log_f, prec);

  if (init) {
    // z = y/x
    fmpq_poly_set(y, init);
    _fmpq_poly_oz_invert_approx(z, f, n, prec);
    fmpq_poly_oz_mul(z, z, y, n);
  } else {
    fmpq_poly_set(y, f);
    fmpq_poly_set_coeff_si(z, 0, 1);
  }

  fmpq_t *xi = (fmpq_t*)calloc(p, sizeof(fmpq_t));
  fmpq_t *a2 = (fmpq_t*)calloc(p, sizeof(fmpq_t));
  fmpq_t *c  = (fmpq_t*)calloc(p, sizeof(fmpq_t));
  fmpq_poly_t *t_ = (fmpq_poly_t*)calloc(p, sizeof(fmpq_poly_t));
  fmpq_poly_t *s_ = (fmpq_poly_t*)calloc(p, sizeof(fmpq_poly_t));

  mpfr_t pi;  mpfr_init2(pi, 4*prec);
  mpfr_const_pi(pi, MPFR_RNDN);
  
#pragma omp parallel for
  for(int i=0; i<p; i++) {
    mpfr_t xi_r; mpfr_init2(xi_r, 4*prec);
    mpfr_t a2_r; mpfr_init2(a2_r, 4*prec);

    /*  ζ_i = 1/2 * (1 + cos( (2·i -1)·π/(2·p) )) */
    mpfr_set_si(xi_r, 2*i+1, MPFR_RNDN);
    mpfr_mul(xi_r, xi_r, pi, MPFR_RNDN);
    mpfr_div_si(xi_r, xi_r, 2*p, MPFR_RNDN);
    mpfr_cos(xi_r, xi_r, MPFR_RNDN);
    mpfr_add_si(xi_r, xi_r, 1, MPFR_RNDN);
    mpfr_div_si(xi_r, xi_r, 2, MPFR_RNDN);

    /* α_i^2 = 1/ζ_i -1 */
    mpfr_set_si(a2_r, 1, MPFR_RNDN);
    mpfr_div(a2_r, a2_r, xi_r, MPFR_RNDN);
    mpfr_sub_si(a2_r, a2_r, 1, MPFR_RNDN);

    fmpq_init(xi[i]);
    fmpq_init(a2[i]);
    fmpq_set_mpfr(xi[i], xi_r, MPFR_RNDN);
    fmpq_set_mpfr(a2[i], a2_r, MPFR_RNDN);

    fmpq_init(c[i]);
    fmpq_poly_init(t_[i]);
    fmpq_poly_init(s_[i]);

    mpfr_clear(xi_r);
    mpfr_clear(a2_r);
  }

  mpfr_clear(pi);
  
  uint64_t t = oz_walltime(0);

  int r = 0;
  int cont = 1;
  for(long  k=0; cont; k++) {
    if (k == 0)
      _fmpq_poly_oz_sqrt_approx_scale(y, z, n, prec);

    /*   T = sum([1/xi[i] * ~(Z*Y + a2[i]) for i in range(p)]) */
#pragma omp parallel for
  for(int i=0; i<p; i++) {
    fmpq_poly_oz_mul(t_[i], z, y, n);
    fmpq_poly_get_coeff_fmpq(c[i], t_[i], 0);
    fmpq_add(c[i], c[i], a2[i]);
    fmpq_poly_set_coeff_fmpq(t_[i], 0, c[i]);
    fmpq_poly_scalar_mul_fmpq(t_[i], t_[i], xi[i]);
    _fmpq_poly_oz_invert_approx(s_[i], t_[i], n, prec);
  }
    
    for(int i=1; i<p; i++)
      fmpq_poly_add(s_[0], s_[0], s_[i]);

#pragma omp parallel sections
    {
#pragma omp section
      {
        fmpq_poly_oz_mul(y_next, y, s_[0], n);
        fmpq_poly_scalar_div_si(y_next, y_next, p);
        fmpq_poly_set(y, y_next);
      }
#pragma omp section
      {
        fmpq_poly_oz_mul(z_next, z, s_[0], n);
        fmpq_poly_scalar_div_si(z_next, z_next, p);
        fmpq_poly_set(z, z_next);
      }
    }
    cont = !_fmpq_poly_oz_sqrt_approx_break(norm, y, f, n, bound, prec);

    if(flags & OZ_VERBOSE) {
      mpfr_log2(log_f, norm, MPFR_RNDN);
      mpfr_fprintf(stderr, "Computing sqrt(Σ)::  k: %4d,  Δ=|sqrt(Σ)^2-Σ|: %7.2Rf", k, log_f);
      fprintf(stderr, " <? %4ld, ", -bound);
      fprintf(stderr, "t: %8.2fs\n", oz_walltime(t)/1000000.0);
      fflush(0);
    }

    if (cont) {
      if (k>0 && mpfr_cmp(norm, prev_norm) >= 0) {
        /*  we don't converge any more */
        r = -1;
        break;
      }
      mpfr_set(prev_norm, norm, MPFR_RNDN);
    }
  }

  for(int i=0; i<p; i++) {
    fmpq_clear(xi[i]);
    fmpq_clear(a2[i]);
    fmpq_clear(c[i]);
    fmpq_poly_clear(t_[i]);
    fmpq_poly_clear(s_[i]);
  }
  free(xi);
  free(a2);
  free(c);
  free(t_);
  free(s_);
  
  mpfr_clear(log_f);
  fmpq_poly_set(f_sqrt, y);
  mpfr_clear(norm);
  mpfr_clear(prev_norm);
  fmpq_poly_clear(y_next);
  fmpq_poly_clear(y);
  fmpq_poly_clear(z_next);
  fmpq_poly_clear(z);
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

mp_limb_t *_fmpz_poly_oz_ideal_is_probaprime_small_primes(const long n, const int k) {
  mp_limb_t q = 1;
  mp_limb_t *small_primes = (mp_limb_t*)calloc(sizeof(mp_limb_t), k);
  small_primes[0] = 2;
  
  for(int i=1; i<k; ) {
    q += n;
    if (n_is_probabprime(q)) {
      small_primes[i] = q;
      i++;
    }
  }
  return small_primes;
}

int fmpz_poly_oz_ideal_is_probaprime(fmpz_poly_t f, fmpz_poly_t g, int sloppy, const int k, const mp_limb_t *small_primes) {
  int r = 1;

  nmod_poly_t f_mod;
  nmod_poly_t g_mod;

  for(int i=0; i<k; i++) {
    mp_limb_t p = small_primes[i];

    // TODO: use _nmod_vec and check conditions manually
    nmod_poly_init(f_mod, p); fmpz_poly_get_nmod_poly(f_mod, f);
    nmod_poly_init(g_mod, p); fmpz_poly_get_nmod_poly(g_mod, g);

    if (nmod_poly_resultant(f_mod, g_mod) == 0) {
      r = 0;
      nmod_poly_clear(f_mod);
      nmod_poly_clear(g_mod);
      break;
    }
    nmod_poly_clear(f_mod);
    nmod_poly_clear(g_mod);
  }

  if (sloppy) {
    return r;
  }

  if (r) {
    fmpz_t norm;
    fmpz_init(norm);
    fmpz_poly_oz_ideal_norm(norm, f, fmpz_poly_degree(g), 0);
    r = fmpz_is_probabprime(norm);
    fmpz_clear(norm);
  }
  return r;
}
