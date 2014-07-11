#ifndef MISC__H
#define MISC__H

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

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

#endif //MISC__H
