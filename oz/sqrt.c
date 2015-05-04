#include "sqrt.h"
#include "oz.h"
#include "util.h"
#include "flint-addons.h"

int _fmpq_poly_oz_sqrt_approx_break(mpfr_t norm, const fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t bound, const mpfr_prec_t prec) {
  fmpq_poly_t f_approx;
  fmpq_poly_init(f_approx);
  fmpq_poly_oz_mul(f_approx, f_sqrt, f_sqrt, n);
  fmpq_poly_sub(f_approx, f_approx, f);
  fmpq_poly_2norm_mpfr(norm, f_approx, MPFR_RNDN);

  mpfr_t f_norm;
  mpfr_init2(f_norm, prec);
  fmpq_poly_2norm_mpfr(f_norm, f, MPFR_RNDN);
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
  fmpq_poly_oz_ideal_norm(gamma_q, y, n, 1);
  fmpq_get_mpfr(gamma, gamma_q, MPFR_RNDN);

  /* det(y) · det(z) */
  fmpq_poly_oz_ideal_norm(gamma_q, z, n, 1);
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

int fmpq_poly_oz_sqrt_approx_babylonian(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t bound, oz_flag_t flags, const fmpq_poly_t init) {
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

    r = _fmpq_poly_oz_sqrt_approx_break(norm, y, f, n, bound, prec);

    if(flags & OZ_VERBOSE) {
      mpfr_log2(log_f, norm, MPFR_RNDN);
      mpfr_fprintf(stderr, "Computing sqrt(Σ)::  k: %4d,  Δ=|sqrt(Σ)^2-Σ|: %7.2Rf", k, log_f);
      fprintf(stderr, " <? %4ld, ", -bound);
      fprintf(stderr, "t: %8.2fs\n", oz_seconds(oz_walltime(t)));
      fflush(0);
    }

    if(r) {
      r = 0;
      break;
    }

    if (k>0 && mpfr_cmp_ui_2exp(norm, 1, bound) >= 0) {
      /* something went really wrong */
      r = -1;
      break;
    }

    mpfr_div_ui(prev_norm, prev_norm, 2, MPFR_RNDN);
    if (k>0 && mpfr_cmp(norm, prev_norm) >= 0) {
      /*  we don't converge any more */
      r = 1;
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

int fmpq_poly_oz_sqrt_approx_db(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t bound, oz_flag_t flags, const fmpq_poly_t init) {
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
    if (k == 0 || mpfr_cmp_ui(prev_norm, 1) > 0)
      _fmpq_poly_oz_sqrt_approx_scale(y, z, n, prec);

#pragma omp parallel sections
    {
#pragma omp section
      {
        _fmpq_poly_oz_invert_approx(y_next, z, n, prec);
        fmpq_poly_add(y_next, y_next, y);
        fmpq_poly_scalar_div_si(y_next, y_next, 2);
        flint_cleanup();
      }
#pragma omp section
      {
        _fmpq_poly_oz_invert_approx(z_next, y, n, prec);
        fmpq_poly_add(z_next, z_next, z);
        fmpq_poly_scalar_div_si(z_next, z_next, 2);
        flint_cleanup();
      }
    }
    fmpq_poly_set(y, y_next);
    fmpq_poly_set(z, z_next);

    r = _fmpq_poly_oz_sqrt_approx_break(norm, y, f, n, bound, prec);

    if(flags & OZ_VERBOSE) {
      mpfr_log2(log_f, norm, MPFR_RNDN);
      mpfr_fprintf(stderr, "Computing sqrt(Σ)::  k: %4d,  Δ=|sqrt(Σ)^2-Σ|: %7.2Rf", k, log_f);
      fprintf(stderr, " <? %4ld, ", -bound);
      fprintf(stderr, "t: %8.2fs\n", oz_seconds(oz_walltime(t)));
      fflush(0);
    }

    if(r) {
      r = 0;
      break;
    }

    if (k>0 && mpfr_cmp_ui_2exp(norm, 1, bound) >= 0) {
      /* something went really wrong */
      r = -1;
      break;
    }

    mpfr_div_ui(prev_norm, prev_norm, 2, MPFR_RNDN);
    if (k>0 && mpfr_cmp(norm, prev_norm) >= 0) {
      /*  we don't converge any more */
      r = 1;
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

int fmpq_poly_oz_sqrt_approx_pade(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const int p, const mpfr_prec_t prec, const mpfr_prec_t bound, oz_flag_t flags, const fmpq_poly_t init) {
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
    if (k == 0 || mpfr_cmp_ui(prev_norm, 1) > 0)
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
    fmpq_poly_add(s_[0],   s_[0], s_[i]);

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
      fprintf(stderr, "t: %8.2fs\n", oz_seconds(oz_walltime(t)));
      fflush(0);
    }

    if (cont) {
      if (k>0 && mpfr_cmp_ui_2exp(norm, 1, bound) >= 0) {
        /* something went really wrong */
        r = -1;
        break;
      }
      if (k>0 && mpfr_cmp(norm, prev_norm) >= 0) {
        /*  we don't converge any more */
        r = 1;
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
