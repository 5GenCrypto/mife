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

static inline void _mpfr_vec_dot_product(mpfr_t rop, const mpfr_t *op1, const mpfr_t *op2, const long n) {
  mpfr_set_si(rop, 0, MPFR_RNDN);
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(rop));
  for(long i=0; i<n; i++) {
    mpfr_mul(tmp, op1[i], op2[i], MPFR_RNDN);
    mpfr_add(rop, rop, tmp, MPFR_RNDN);
  }
  mpfr_clear(tmp);
}

static inline void _mpfr_vec_submul(mpfr_t *rop, mpfr_t scalar, mpfr_t *op, const long n) {
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(rop[0]));
  for(long i=0; i<n; i++) {
    mpfr_mul(tmp, scalar, op[i], MPFR_RNDN);
    mpfr_sub(rop[i], rop[i], tmp, MPFR_RNDN);
  }
  mpfr_clear(tmp);
}

static inline void _mpfr_vec_2norm(mpfr_t rop, const mpfr_t *vec, long n) {
  _mpfr_vec_dot_product(rop, vec, vec, n);
  mpfr_sqrt(rop, rop, MPFR_RNDN);
}

#endif
