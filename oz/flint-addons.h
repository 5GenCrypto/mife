#ifndef FLINT_ADDONS__H
#define FLINT_ADDONS__H

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <mpfr.h>
#include <assert.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_poly.h>
#include <math.h>

/* functions dealing with fmpz types and matrix multiplications mod fmpz_t */
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);
void fmpz_mat_scalar_mul_modp(fmpz_mat_t m, fmpz_t scalar, fmpz_t modp);
void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n,
    fmpz_t p);
void fmpz_init_exp(fmpz_t exp, int base, int n);

/* functions for computing matrix inverse mod fmpz_t */
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);
void fmpz_mat_det_modp(fmpz_t det, fmpz_mat_t a, int n, fmpz_t p);
void fmpz_mat_cofactor_modp(fmpz_mat_t b, fmpz_mat_t a, int n, fmpz_t p);
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);

static inline mp_limb_t n_prevprime(mp_limb_t n, int proved) {
  n--;
  if (n <= 1)
    abort();
  if (n <= 3)
    return n;
  if (n%2 == 0)
    n -= 1;
  while ((!proved && !n_is_probabprime(n)) || (proved && !n_is_prime(n)))
    n -= 2;
  return n;
}

static inline void fmpq_set_mpfr(fmpq_t rop, const mpfr_t op, mpfr_rnd_t rnd) {
  mpf_t  op_f;
  mpf_init2(op_f, mpfr_get_prec(op));
  mpfr_get_f(op_f, op, rnd);
  mpq_t  op_q;
  mpq_init(op_q);
  mpq_set_f(op_q, op_f);
  mpf_clear(op_f);
  fmpq_set_mpq(rop, op_q);
  mpq_clear(op_q);
}

/**
   fmpz_mat_t
*/

static inline void _fmpz_vec_mul(fmpz_t rop, fmpz *a, fmpz *b, long len) {
  fmpz_zero(rop);
  for (long i = 0; i < len; i++)
    fmpz_addmul(rop, a + i, b + i);
}

static inline void _fmpz_vec_rot_left(fmpz *rop, fmpz *op, long len, long shift) {
  fmpz* tmp = _fmpz_vec_init(shift);
  for(long i=0; i<shift; i++) {
    fmpz_set(tmp + i, op + len - shift + i);
  }
  for(long i=len-shift-1; i>=0; i--)
    fmpz_set(rop + shift + i, op + i);
  for(long i=0; i<shift; i++) {
    fmpz_set(rop + i, tmp + i);
  }
  _fmpz_vec_clear(tmp, shift);
}

static inline void _fmpz_vec_rot_left_neg(fmpz *rop, fmpz *op, long len) {
  _fmpz_vec_rot_left(rop, op, len, 1);
  fmpz_neg(rop + 0, rop + 0);
}

// TODO: This function should accept a rounding mode.
static inline void _fmpz_vec_set_mpfr_vec(fmpz *rop, mpfr_t *op, const long len) {
  mpz_t tmp;
  mpz_init(tmp);
  for(long i=0; i<len; i++) {
    mpfr_get_z(tmp, op[i], MPFR_RNDN);
    fmpz_set_mpz(rop + i, tmp);
  }
  mpz_clear(tmp);
}

void _fmpz_vec_eucl_norm_mpfr(mpfr_t rop, const fmpz *vec, const long len, const mpfr_rnd_t rnd);

static inline void fmpz_poly_2norm_mpfr(mpfr_t rop, const fmpz_poly_t poly, const mpfr_rnd_t rnd) {
  _fmpz_vec_eucl_norm_mpfr(rop, poly->coeffs, poly->length, rnd);
}

static inline double fmpz_poly_2norm_log2(const fmpz_poly_t poly) {
  if (fmpz_poly_is_zero(poly))
    return -1;
  mpfr_t tmp;
  mpfr_init2(tmp, labs(fmpz_poly_max_bits(poly)));
  _fmpz_vec_eucl_norm_mpfr(tmp, poly->coeffs, poly->length, MPFR_RNDN);
  mpfr_log2(tmp, tmp, MPFR_RNDN);
  double r = mpfr_get_d(tmp, MPFR_RNDN);
  mpfr_clear(tmp);
  return r;
}


static inline void fmpz_mod_poly_eucl_norm_mpfr(mpfr_t rop, const fmpz_mod_poly_t poly, const mpfr_rnd_t rnd) {
  _fmpz_vec_eucl_norm_mpfr(rop, poly->coeffs, poly->length, rnd);
}

void _fmpq_vec_eucl_norm_mpfr(mpfr_t rop, const fmpz *num, const fmpz_t den, const long len, const mpfr_rnd_t rnd);

static inline void fmpq_poly_2norm_mpfr(mpfr_t rop, const fmpq_poly_t poly, const mpfr_rnd_t rnd) {
  _fmpq_vec_eucl_norm_mpfr(rop, poly->coeffs, poly->den, poly->length, rnd);
}

static inline void fmpz_poly_mulmod(fmpz_poly_t rop, const fmpz_poly_t op1, const fmpz_poly_t op2, const fmpz_poly_t modulus) {
  fmpz_poly_mul(rop, op1, op2);
  fmpz_poly_rem(rop, rop, modulus);
}

static inline void fmpq_poly_mulmod(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2, const fmpq_poly_t modulus) {
  fmpq_poly_mul(rop, op1, op2);
  fmpq_poly_rem(rop, rop, modulus);
}

void fmpq_poly_truncate_prec(fmpq_poly_t op, const mp_bitcnt_t prec);

static inline void fmpq_poly_invert_mod(fmpq_poly_t f_inv, const fmpq_poly_t f, const fmpq_poly_t g) {
  fmpq_poly_t r; fmpq_poly_init(r);
  fmpq_poly_t t; fmpq_poly_init(t);

  fmpq_poly_xgcd(r, f_inv, t, f, g);
  assert(fmpq_poly_is_one(r));

  fmpq_poly_clear(t);
  fmpq_poly_clear(r);
}

static inline void fmpz_poly_invert_mod_fmpq(fmpq_poly_t f_inv, fmpz_poly_t f, const fmpz_poly_t g) {
  fmpq_poly_t r; fmpq_poly_init(r);
  fmpq_poly_t t; fmpq_poly_init(t);

  fmpq_poly_t f_; fmpq_poly_init(f_); fmpq_poly_set_fmpz_poly(f_, f);
  fmpq_poly_t g_; fmpq_poly_init(g_); fmpq_poly_set_fmpz_poly(g_, g);

  fmpq_poly_xgcd(r, f_inv, t, f_, g_);
  assert(fmpq_poly_is_one(r));

  fmpq_poly_clear(g_);
  fmpq_poly_clear(f_);
  fmpq_poly_clear(t);
  fmpq_poly_clear(r);
}

static inline void fmpz_poly_ideal_rot_basis(fmpz_mat_t rop, fmpz_poly_t op) {
  assert(fmpz_mat_nrows(rop) == fmpz_poly_length(op));
  assert(fmpz_mat_ncols(rop) == fmpz_poly_length(op));

  printf("WARNING: called fmpz_poly_ideal_rot_basis, O(n^ω) algorithms ahead.\n");

  const long n = fmpz_poly_length(op);

  fmpz* v = _fmpz_vec_init(n);
  _fmpz_vec_set(v, op->coeffs, n);
  for(long i=0; i<n; i++) {
    for (long j=0; j<n; j++)
      fmpz_set(fmpz_mat_entry(rop, i, j), v+j);
    _fmpz_vec_rot_left_neg(v, v, n);
  }
  _fmpz_vec_clear(v, n);
}

void fmpz_poly_resultant_modular_bound(fmpz_t res, const fmpz_poly_t poly1,
                                       const fmpz_poly_t poly2, const mp_bitcnt_t bound);

void fmpq_poly_resultant_modular_bound(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g,
                                       const mp_bitcnt_t bound);

void _fmpz_poly_resultant_modular_bound(fmpz_t res, const fmpz * poly1, slong len1,
                                        const fmpz * poly2, slong len2, mp_bitcnt_t bound);


void _fmpq_poly_resultant_modular(fmpz_t rnum, fmpz_t rden, 
                                  const fmpz *poly1, const fmpz_t den1, slong len1, 
                                  const fmpz *poly2, const fmpz_t den2, slong len2,
                                  const mp_bitcnt_t bound);

static inline void fmpz_poly_set_fmpz_mod_poly(fmpz_poly_t rop, fmpz_mod_poly_t op) {
  long n = fmpz_mod_poly_length(op);

  fmpz *q = fmpz_mod_poly_modulus(op);

  fmpz_t q2;
  fmpz_init(q2);
  fmpz_fdiv_q_2exp(q2,q, 1);

  fmpz_t tmp;
  fmpz_init(tmp);

  for(long i=0; i<n; i++) {
    fmpz_mod_poly_get_coeff_fmpz(tmp, op, i);
    if (fmpz_cmp(tmp, q2)>=0)
      fmpz_sub(tmp, tmp, q);
    fmpz_poly_set_coeff_fmpz(rop, i, tmp);
  }
  fmpz_clear(tmp);
  fmpz_clear(q2);
}

/**
   fmpz_mod_poly_t
 */

static inline void fmpz_mod_poly_invert_mod(fmpz_mod_poly_t f_inv, fmpz_mod_poly_t f, fmpz_mod_poly_t modulus) {
  fmpz_mod_poly_t r; fmpz_mod_poly_init(r, fmpz_mod_poly_modulus(modulus));
  fmpz_mod_poly_t t; fmpz_mod_poly_init(t, fmpz_mod_poly_modulus(modulus));

  fmpz_mod_poly_xgcd(r, f_inv, t, f, modulus);
  if (fmpz_mod_poly_degree(r) != 0 || !fmpz_is_one(r->coeffs + 0)) {
      fmpz_mod_poly_zero(f_inv);
  }

  fmpz_mod_poly_clear(t);
  fmpz_mod_poly_clear(r);
}

static inline void fmpz_mod_poly_set_fmpq_poly(fmpz_mod_poly_t rop, fmpq_poly_t op) {
  fmpz_t tmp;
  fmpz_t den_inv;

  fmpz *q = fmpz_mod_poly_modulus(rop);

  fmpz_init(tmp);
  fmpz_init(den_inv);

  fmpz_set(den_inv, fmpq_poly_denref(op));
  fmpz_invmod(den_inv, den_inv, q);

  for(long i=0; i<fmpq_poly_length(op); i++) {
    fmpz_set(tmp, fmpq_poly_numref(op)+i);
    fmpz_mul(tmp, tmp, den_inv);
    fmpz_mod(tmp, tmp, q);
    fmpz_mod_poly_set_coeff_fmpz(rop, i, tmp);
  }
  fmpz_clear(tmp);
  fmpz_clear(den_inv);
}

/**
   flint_rand_t
*/

static inline void flint_randinit_seed(flint_rand_t randstate, mp_limb_t seed, int gmp) {
  flint_randinit(randstate);
  randstate->__randval  = seed;
  randstate->__randval2 = seed ^ 0x5555555555555555ULL;
  if (gmp) {
    _flint_rand_init_gmp(randstate);
    gmp_randseed_ui(randstate->gmp_state, seed ^ 0xDEADBEEFDEADBEEFULL);
  }
}

#endif //FLINT_ADDONS__H
