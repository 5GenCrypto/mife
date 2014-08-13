#ifndef GSO__H
#define GSO__H

#include <flint/fmpz_mat.h>
#include <mpfr.h>

typedef struct {
    mpfr_t * entries;
    long r;
    long c;
    mpfr_t ** rows;
} mpfr_mat_struct;

typedef mpfr_mat_struct mpfr_mat_t[1];

void mpfr_mat_init(mpfr_mat_t mat, long rows, long cols, mpfr_prec_t prec);
void mpfr_mat_set_fmpz_mat(mpfr_mat_t rop, const fmpz_mat_t op);
void mpfr_mat_clear(mpfr_mat_t mat);
mpfr_prec_t mpfr_mat_get_prec(mpfr_mat_t mat);
void mpfr_mat_gso(mpfr_mat_t mat);

#endif
