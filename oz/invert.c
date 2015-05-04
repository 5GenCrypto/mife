#include "invert.h"
#include "util.h"
#include "oz.h"
#include "flint-addons.h"

void fmpz_mod_poly_oz_invert(fmpz_mod_poly_t f_inv, const fmpz_mod_poly_t f, const long n) {
  assert(1<<n_clog(n,2) == n);
  if(f_inv == f)
    oz_die("fmpz_mod_poly_oz_invert does not support parameter aliasing");

  fmpz_mod_poly_t V;  fmpz_mod_poly_init(V, fmpz_mod_poly_modulus(f));
  fmpz_mod_poly_set(V, f);

  fmpz_mod_poly_t U;  fmpz_mod_poly_init(U, fmpz_mod_poly_modulus(f));

  /* F(x) = x^n+1 */
  fmpz_mod_poly_t V2;  fmpz_mod_poly_init(V2, fmpz_mod_poly_modulus(f));

  fmpz_mod_poly_t Se;  fmpz_mod_poly_init(Se, fmpz_mod_poly_modulus(f));
  fmpz_mod_poly_t So;  fmpz_mod_poly_init(So, fmpz_mod_poly_modulus(f));

  fmpz_t tmp;
  fmpz_init(tmp);

  if (n > 1) {
    /* set V2(x) = V(-x) */
    fmpz_mod_poly_set(V2, V);
    int deg = fmpz_mod_poly_degree(V2);
    for (int i = 1 ; i <= deg ; i += 2) {
      /* negate odd coefficients */
      fmpz_mod_poly_get_coeff_fmpz(tmp,V2,i);
      fmpz_neg(tmp,tmp);
      fmpz_mod_poly_set_coeff_fmpz(V2,i,tmp);
    }

    /* compute Se, So */
    /* Set Se = (V(x) + V(-x))/2 */
    for (int i = 0; i <= deg/2 ; i++) {
      fmpz_mod_poly_get_coeff_fmpz(tmp, V2, 2*i);
      fmpz_mod_poly_set_coeff_fmpz(Se, i, tmp);
    }
    fmpz_mod_poly_truncate(Se,deg/2+1);

    /* Set So to the compressed (V(-x) - V(x))/2 */
    for (int i = 0 ; i < (deg+1)/2 ; i++) {
      fmpz_mod_poly_get_coeff_fmpz(tmp, V2, 2*i+1);
      fmpz_mod_poly_set_coeff_fmpz(So,i,tmp);
    }
    fmpz_mod_poly_truncate(So,deg/2+1);

    /* V = V(x) * V(-x) mod f(x) */
    fmpz_mod_poly_oz_mul(V, V, V2, n);

    deg = fmpz_mod_poly_degree(V);
    /* Compress the non-zero coeffs of V */
    for (int i = 0; i <= deg/2 ; i ++) {
      fmpz_mod_poly_get_coeff_fmpz(tmp, V, 2*i);
      fmpz_mod_poly_set_coeff_fmpz(V, i, tmp);
    }
    fmpz_mod_poly_truncate(V,deg/2+1);
    fmpz_mod_poly_oz_invert(V2, V, n/2);

    /* Te=G*Se, To = G*So */
    fmpz_mod_poly_mul(V,V2,Se);
    fmpz_mod_poly_mul(U,V2,So);

    deg = fmpz_mod_poly_degree(V);
    fmpz_mod_poly_zero(f_inv);
    for (int i = 0 ; i <= deg ; i ++) {
      fmpz_mod_poly_get_coeff_fmpz(tmp, V, i);
      fmpz_mod_poly_set_coeff_fmpz(f_inv, 2*i, tmp);
    }
    deg = fmpz_mod_poly_degree(U);
    for (int i = 0 ; i <= deg ; i ++) {
      fmpz_mod_poly_get_coeff_fmpz(tmp, U, i);
      fmpz_mod_poly_set_coeff_fmpz(f_inv, 2*i+1, tmp);
    }
    fmpz_mod_poly_oz_rem(f_inv,f_inv,n);

    /* truncate results on the precision required by algorithm */
  } else {
    fmpz_mod_poly_zero(f_inv);
    fmpz_mod_poly_get_coeff_fmpz(tmp, f, 0);
    fmpz_invmod(tmp, tmp, fmpz_mod_poly_modulus(f));
    fmpz_mod_poly_set_coeff_fmpz(f_inv, 0, tmp);
    fmpz_mod_poly_truncate(f_inv,1);
  }

  /* free memory */
  fmpz_clear(tmp);
  fmpz_mod_poly_clear(So);
  fmpz_mod_poly_clear(Se);
  fmpz_mod_poly_clear(V2);
  fmpz_mod_poly_clear(U);
  fmpz_mod_poly_clear(V);
}

void _fmpq_poly_oz_invert_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, const long n, const mpfr_prec_t prec) {
  if(f_inv == f)
    oz_die("_fmpq_poly_oz_invert_approx does not support parameter aliasing");

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

  mpq_t *tmp_q = (mpq_t*)calloc(2*n, sizeof(mpq_t));
  for(int i=0; i<2*n; i++)
    mpq_init(tmp_q[i]);

  if (n > 1) {
    /* set V2(x) = V(-x) */
    fmpq_poly_set(V2, V);
    int deg = fmpq_poly_degree(V2);
    for (int i = 1 ; i <= deg ; i += 2) {
      /* negate odd coefficients */
      fmpz_neg(fmpq_poly_numref(V2) + i, fmpq_poly_numref(V2) + i);
    }

    /* compute Se, So */
    /* Set Se = (V(x) + V(-x))/2 */
    for (int i = 0; i <= deg/2 ; i++) {
      fmpq_poly_get_coeff_mpq(tmp_q[i], V2, 2*i);
      /* fmpq_poly_set_coeff_fmpq(Se, i, tmp); */
    }
    fmpq_poly_set_array_mpq(Se, (const mpq_t*)tmp_q, deg/2+1);
    fmpq_poly_truncate(Se,deg/2+1);

    /* Set So to the compressed (V(-x) - V(x))/2 */
    mpq_set_si(tmp_q[(deg+1)/2], 0, 1);
    for (int i = 0 ; i < (deg+1)/2 ; i++) {
      fmpq_poly_get_coeff_mpq(tmp_q[i], V2, 2*i+1);
      /* fmpq_poly_set_coeff_fmpq(So,i,tmp); */
    }
    fmpq_poly_set_array_mpq(So,(const mpq_t*)tmp_q, deg/2+1);
    fmpq_poly_truncate(So,deg/2+1);

    /* V = V(x) * V(-x) mod f(x) */
    fmpq_poly_oz_mul(V, V, V2, n);

    deg = fmpq_poly_degree(V);
    /* Compress the non-zero coeffs of V */
    for (int i = 0; i <= deg/2 ; i ++) {
      fmpq_poly_get_coeff_mpq(tmp_q[i],V,2*i);
      //fmpq_poly_set_coeff_fmpq(V,i,tmp);
    }
    fmpq_poly_set_array_mpq(V, (const mpq_t*)tmp_q, deg/2+1);
    fmpq_poly_truncate(V,deg/2+1);

    if (prec)
      fmpq_poly_truncate_prec(V, prec);
    _fmpq_poly_oz_invert_approx(V2, V, n/2, prec);

    deg = fmpq_poly_degree(V2);
    /* Te=G*Se, To = G*So */
    fmpq_poly_mul(V,V2,Se);
    fmpq_poly_mul(U,V2,So);

    fmpq_poly_zero(f_inv);
    for(int i=0; i<2*n; i++)
      mpq_set_si(tmp_q[i], 0, 1);

    deg = fmpq_poly_degree(V);
    for (int i = 0 ; i <= deg ; i ++) {
      fmpq_poly_get_coeff_mpq(tmp_q[2*i], V, i);
      /* fmpq_poly_set_coeff_fmpq(f_inv,2*i,tmp); */
    }
    deg = fmpq_poly_degree(U);
    for (int i = 0 ; i <= deg ; i ++) {
      fmpq_poly_get_coeff_mpq(tmp_q[2*i+1], U, i);
      /* fmpq_poly_set_coeff_fmpq(f_inv,2*i+1,tmp); */
    }
    fmpq_poly_set_array_mpq(f_inv, (const mpq_t*)tmp_q, 2*deg+3);
    fmpq_poly_oz_rem(f_inv,f_inv,n);

    /* truncate results on the precision required by algorithm */
    if (prec)
      fmpq_poly_truncate_prec(f_inv, prec);
  } else {
    fmpq_poly_zero(f_inv);
    if (fmpz_poly_is_zero(f))
      oz_die("division by zero.");
    f_inv->length = 1;
    fmpq_poly_get_coeff_fmpq(tmp, f, 0);
    fmpq_inv(tmp, tmp);
    fmpq_poly_set_coeff_fmpq(f_inv, 0, tmp);
  }

  /* free memory */
  for(int i=0; i<2*n; i++)
    mpq_clear(tmp_q[i]);
  free(tmp_q);

  fmpq_clear(tmp);
  fmpq_poly_clear(So);
  fmpq_poly_clear(Se);
  fmpq_poly_clear(V2);
  fmpq_poly_clear(U);
  fmpq_poly_clear(V);
}

void fmpq_poly_oz_invert_approx(fmpq_poly_t rop, const fmpq_poly_t f, const long n,
                                const mpfr_prec_t prec, const oz_flag_t flags) {

  if (prec == 0) {
    _fmpq_poly_oz_invert_approx(rop, f, n, prec);
    return;
  }
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

    fmpq_poly_2norm_mpfr(norm, tmp, MPFR_RNDN);

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
  mpfr_clear(tmp_f);

  fmpq_poly_set(rop, f_inv);
  fmpq_poly_clear(f_inv);

  mpfr_clear(bound);
  fmpq_clear(c);
  mpfr_clear(norm);
  fmpq_poly_clear(tmp);
}
