/**
   @file misc.h
   @brief Utility Functions
*/

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

static inline double ggh_seconds(uint64_t t) {
  return t/1000000.0;
}

#include <gghlite/gghlite-defs.h>


/**
   @brief Print if GGHLITE_FLAGS_VERBOSE is set

   @param self   GGHITE `params`
   @param msg    format string for printf()
   @param ...    arguments for format string `msg`
*/

void ggh_printf_v(const gghlite_params_t self, const char *msg, ...);

/**
   @brief Print if GGHLITE_FLAGS_QUIET is not set

   @param self   GGHITE `params`
   @param msg    format string for printf()
   @param ...    arguments for format string `msg`
*/

void ggh_printf(const gghlite_params_t self, const char *msg, ...);

#endif //MISC__H
