#ifndef MISC__H
#define MISC__H

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>

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


#include <sys/time.h>
#include <unistd.h>

static inline uint64_t ggh_walltime(uint64_t t0) {
  static time_t base_sec;
  struct timeval tp;
  gettimeofday(&tp, NULL);
  if (base_sec == 0)
    base_sec = tp.tv_sec;
  return ((uint64_t)(tp.tv_sec - base_sec)) * 1000000 + (uint64_t)tp.tv_usec - t0;
}

#endif //MISC__H
