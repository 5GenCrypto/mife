#include "invert.h"
#include "util.h"
#include "oz.h"
#include "flint-addons.h"


void _fmpq_poly_oz_invert_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, const int n, const mpfr_prec_t prec) {
  if(f_inv == f)
    oz_die("_fmpq_poly_oz_invert_approx does not support parameter aliasing");

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

  mpq_t *tmp_q = (mpq_t*)calloc(2*n, sizeof(mpq_t));
  for(int i=0; i<2*n; i++)
    mpq_init(tmp_q[i]);

  int deg;
  if (n > 1) {
    /* set V2(x) = V(-x) */
    fmpq_poly_set(V2, V);
    deg = fmpq_poly_degree(V2);
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
  mpfr_clear(tmp_f);

  fmpq_poly_set(rop, f_inv);
  fmpq_poly_clear(f_inv);

  mpfr_clear(bound);
  fmpq_clear(c);
  mpfr_clear(norm);
  fmpq_poly_clear(tmp);
}
