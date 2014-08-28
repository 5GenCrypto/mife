#ifndef MISC__H
#define MISC__H

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include <gpv/gpv.h>
#include <mpfr.h>

#define S_TO_SIGMA 0.398942280401433 /* 1/sqrt(2*pi) */

static inline void ggh_die(const char *msg, ...) {
  va_list lst;
  va_start(lst, msg);
  vfprintf(stderr, msg, lst);
  va_end(lst);
  abort();
}

static inline void *ggh_malloc(size_t size) {
  void *ret = malloc(size);
  if (!ret)
    ggh_die("Out of memory");
  return ret;
}

static inline void ggh_free(void *condemned) {
  free(condemned);
}

static inline void *ggh_calloc(size_t nmemb, size_t size) {
  void *ret = calloc(nmemb, size);
  if (!ret)
    ggh_die("Out of memory");
  return ret;
}

static inline gpv_mp_t *_gghlite_gpv_from_poly(fmpz_poly_t g, mpfr_t sigma, mpfr_t *c, gpv_alg_t algorithm) {
  fmpz_mat_t B;
  const long n = fmpz_poly_length(g);
  fmpz_mat_init(B, 1, n);
  _fmpz_vec_set(B->rows[0], g->coeffs, n);

  /* GPV samples proportionally to `\exp(-(x-c)²/(2σ²))` but GGHLite is
     specifiied with respect to `\exp(-π(x-c)²/σ²)`. So we divide by \sqrt{2π}
  */

  mpfr_t sigma_;
  mpfr_init2(sigma_, mpfr_get_prec(sigma));
  mpfr_set(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, S_TO_SIGMA, MPFR_RNDN);
  
  gpv_mp_t *D = gpv_mp_init(B, sigma_, c, algorithm);
  fmpz_mat_clear(B);
  mpfr_clear(sigma_);
  return D;
}

static inline gpv_mp_t *_gghlite_gpv_from_n(const long n, mpfr_t sigma) {
  fmpz_mat_t I;
  fmpz_mat_init(I, 1, n);
  fmpz_mat_one(I);

  /* GPV samples proportionally to `\exp(-(x-c)²/(2σ²))` but GGHLite is
     specifiied with respect to `\exp(-π(x-c)²/σ²)`. So we divide by \sqrt{2π}
  */

  mpfr_t sigma_;
  mpfr_init2(sigma_, mpfr_get_prec(sigma));
  mpfr_set(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, S_TO_SIGMA, MPFR_RNDN);
 
  gpv_mp_t *D = gpv_mp_init(I, sigma_, NULL, GPV_IDENTITY);
  fmpz_mat_clear(I);
  mpfr_clear(sigma_);
  return D;
}

#endif //MISC__H
