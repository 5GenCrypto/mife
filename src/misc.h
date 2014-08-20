#ifndef MISC__H
#define MISC__H

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include <gpv/gpv.h>
#include <mpfr.h>

static inline void ggh_die(const char *msg, ...) {
  va_list lst;
  va_start(lst, msg);
  vfprintf(stderr, msg, lst);
  va_end(lst);
  abort();
};

static inline void *ggh_malloc(size_t size) {
  void *ret = malloc(size);
  if (!ret)
    ggh_die("Out of memory");
  return ret;
};

static inline void ggh_free(void *condemned) {
  free(condemned);
};

static inline void *ggh_calloc(size_t nmemb, size_t size) {
  void *ret = calloc(nmemb, size);
  if (!ret)
    ggh_die("Out of memory");
  return ret;
};

static inline gpv_mp_t *_gghlite_gpv_from_poly(fmpz_poly_t g, mpfr_t sigma, const fmpz *c, gpv_alg_t algorithm) {
  fmpz_mat_t B;
  const long n = fmpz_poly_length(g);
  fmpz_mat_init(B, 1, n);
  _fmpz_vec_set(B->rows[0], g->coeffs, n);
  gpv_mp_t *D = gpv_mp_init(B, sigma, c, algorithm);
  fmpz_mat_clear(B);
  return D;
}

static inline gpv_mp_t *_gghlite_gpv_from_n(const long n, mpfr_t sigma) {
  fmpz_mat_t I;
  fmpz_mat_init(I, 1, n);
  fmpz_mat_one(I);
  gpv_mp_t *D = gpv_mp_init(I, sigma, NULL, GPV_IDENTITY);
  fmpz_mat_clear(I);
  return D;
}

#endif //MISC__H
