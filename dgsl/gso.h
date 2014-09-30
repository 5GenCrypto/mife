#ifndef GSO__H
#define GSO__H

#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include "flint-addons/flint-addons.h"

typedef struct {
    mpfr_t * entries;
    long r;
    long c;
    mpfr_t ** rows;
} mpfr_mat_struct;

typedef mpfr_mat_struct mpfr_mat_t[1];

void mpfr_mat_init(mpfr_mat_t mat, long rows, long cols, mpfr_prec_t prec);
void mpfr_mat_set_fmpz_mat(mpfr_mat_t rop, const fmpz_mat_t op);
void mpfr_mat_set_fmpz_mat_rot(mpfr_mat_t rop, const fmpz_mat_t op);
void mpfr_mat_clear(mpfr_mat_t mat);

static inline int mpfr_mat_is_empty(const mpfr_mat_t mat) {
  return (mat->r==0) || (mat->c==0);
}

mpfr_prec_t mpfr_mat_get_prec(mpfr_mat_t mat);
void mpfr_mat_gso(mpfr_mat_t mat, mpfr_rnd_t rnd);

static inline mpfr_t * _mpfr_vec_init(const long n, mpfr_prec_t prec) {
  mpfr_t *ret = (mpfr_t*)calloc(n, sizeof(mpfr_t));
  if (!ret)
    dgs_die("out of memory");
  for(long i=0; i<n; i++) {
    mpfr_init2(ret[i], prec);
  }
  return ret;
}

static inline void _mpfr_vec_clear(mpfr_t *op, const long n) {
  for(long i=0; i<n; i++) {
    mpfr_clear(op[i]);
  }
  free(op);
}

static inline void _mpfr_vec_set(mpfr_t *rop, mpfr_t *op, const long n, const mpfr_rnd_t rnd) {
  for(long i=0; i<n; i++) {
    mpfr_set(rop[i], op[i], MPFR_RNDN);
  }
}

static inline void _mpfr_vec_set_fmpz_vec(mpfr_t *rop, fmpz *op, const long n, const mpfr_rnd_t rnd) {
  mpz_t tmp;
  mpz_init(tmp);
  for(long i=0; i<n; i++) {
    fmpz_get_mpz(tmp, op + i);
    mpfr_set_z(rop[i], tmp, MPFR_RNDN);
  }
  mpz_clear(tmp);
}

static inline int _mpfr_vec_is_zero(mpfr_t *op, const long n) {
  for(long i=0; i<n; i++) {
    if (!mpfr_zero_p(op[i]))
      return 0;    
  }
  return 1;
}

static inline void _mpfr_vec_dot_product(mpfr_t rop, mpfr_t *op1, mpfr_t *op2, const long n, const mpfr_rnd_t rnd) {
  mpfr_set_si(rop, 0, MPFR_RNDN);
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(rop));
  for(long i=0; i<n; i++) {
    mpfr_mul(tmp, op1[i], op2[i], rnd);
    mpfr_add(rop, rop, tmp, rnd);
  }
  mpfr_clear(tmp);
}

static inline void _mpfr_vec_scalar_addmul_mpfr(mpfr_t *rop, mpfr_t *op, const long n, mpfr_t scalar, const mpfr_rnd_t rnd) {
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(op[0]));
  for(long i=0; i<n; i++) {
    mpfr_mul(tmp, scalar, op[i], rnd);
    mpfr_add(rop[i], rop[i], tmp, rnd);
  }
  mpfr_clear(tmp);
}

static inline void _mpfr_vec_2norm(mpfr_t rop, mpfr_t *vec, long n, const mpfr_rnd_t rnd) {
  _mpfr_vec_dot_product(rop, vec, vec, n, rnd);
  mpfr_sqrt(rop, rop, rnd);
}

static inline void _mpfr_vec_add(mpfr_t *rop, mpfr_t *op1, mpfr_t *op2, const long n, mpfr_rnd_t rnd) {
  for(long i=0; i<n; i++) {
    mpfr_add(rop[i], op1[i], op2[i], rnd);
  }
}

#endif
