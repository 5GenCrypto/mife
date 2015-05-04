#ifndef _UTIL_H_
#define _UTIL_H_

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>

static inline void oz_die(const char *msg, ...) {
  va_list lst;
  va_start(lst, msg);
  vfprintf(stderr, msg, lst);
  va_end(lst);
  abort();
}

#include <sys/time.h>
#include <unistd.h>

static inline uint64_t oz_walltime(uint64_t t0) {
  static time_t base_sec;
  struct timeval tp;
  gettimeofday(&tp, NULL);
  if (base_sec == 0)
    base_sec = tp.tv_sec;
  return ((uint64_t)(tp.tv_sec - base_sec)) * 1000000 + (uint64_t)tp.tv_usec - t0;
 }

static inline double oz_seconds(uint64_t t) {
  return t/1000000.0;
}

#endif /* _UTIL_H_ */
