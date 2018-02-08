#ifndef _MIFE_UTILS_H
#define _MIFE_UTILS_H

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_poly.h>

#define ALLOC_FAILS(path, len) (NULL == ((path) = malloc((len) * sizeof(*(path)))))
#define AES_SEED_BYTE_SIZE 32

extern int PRINT_TIMERS;
extern uint64_t T;

typedef struct {
	char *path;
	bool stack_allocated;
} location;

void location_free(location loc);
location location_append(const location loc, const char *const path);

typedef enum {
	PARSE_SUCCESS,
	PARSE_INVALID,
	PARSE_OUT_OF_MEMORY,
	PARSE_IO_ERROR
} parse_result;

void check_parse_result(parse_result result, void usage(int), int problem);
bool create_directory_if_missing(char *dir);

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

static inline void start_timer (void) {
    T = ggh_walltime(0);
}

static inline void timer_printf(const char *msg, ...) {
    if(PRINT_TIMERS) {
        va_list lst;
        va_start(lst, msg);
        vfprintf(stdout, msg, lst);
        va_end(lst);
        fflush(0);
    }
}

static inline void print_timer(void) {
    timer_printf("%8.2fs", ggh_seconds(ggh_walltime(T)));
}

int fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);

#endif /* ifndef _MIFE_UTILS_H */
