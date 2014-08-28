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

gpv_mp_t *_gghlite_gpv_from_poly(fmpz_poly_t g, mpfr_t sigma, mpfr_t *c, gpv_alg_t algorithm);

gpv_mp_t *_gghlite_gpv_from_n(const long n, mpfr_t sigma);

#endif //MISC__H
