#ifndef FLINT_ADDONS__H
#define FLINT_ADDONS__H

#include <assert.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>

#include <gghlite/misc.h>
/**
   fmpz_mat_t
*/

static inline int fmpz_mat_rot_is_one(const fmpz_mat_t mat) {
  assert(mat->c);
  if (mat->r != 1)
    return 0;
  if(!fmpz_is_one(&mat->rows[0][0]))
    return 0;
  for(long i=1; i<mat->c; i++)
    if(!fmpz_is_zero(&mat->rows[0][i]))
      return 0;
  return 1;
}

static inline int fmpz_mat_is_one(const fmpz_mat_t mat) {
  if (!fmpz_mat_is_square(mat))
    return 0;
  for(long i=0; i<mat->r; i++) {
    for(long j=0; j<i; j++)
      if(!fmpz_is_zero(&mat->rows[i][j]))
        return 0;
    if(!fmpz_is_one(&mat->rows[i][i]))
      return 0;
    for(long j=i+1; j<mat->c; j++)
      if(!fmpz_is_zero(&mat->rows[i][j]))
        return 0;
  }
  return 1;
}

/**
   fmpz*
*/

// TODO: This function should accept a rounding mode.
static inline void _fmpz_vec_2norm_mpfr(mpfr_t rop, fmpz *vec, long len) {
  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_zero(tmp);
  for (long i = 0; i < len; i++)
    fmpz_addmul(tmp, vec + i, vec + i);
  mpz_t tmp_g;
  mpz_init(tmp_g);
  fmpz_get_mpz(tmp_g, tmp);
  fmpz_clear(tmp);

  mpfr_set_z(rop, tmp_g, MPFR_RNDN);
  mpfr_sqrt(rop, rop, MPFR_RNDN);
  mpz_clear(tmp_g);
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
/**
   fmpq*

*/

// TODO: This function should accept a rounding mode.
static inline void _fmpq_vec_2norm_mpfr(mpfr_t rop, fmpz *num, fmpz_t den, long len) {
  fmpz_t acc_num;
  fmpz_init(acc_num);
  fmpz_zero(acc_num);
  for (long i = 0; i < len; i++)
    fmpz_addmul(acc_num, num + i, num + i);

  fmpz_t acc_den;
  fmpz_init(acc_den);
  fmpz_mul(acc_den, den, den);
  fmpz_mul_ui(acc_den, acc_den, len);

  fmpq_t acc;
  fmpq_init(acc);

  fmpq_set_fmpz_frac(acc, acc_num, acc_den);
  fmpq_get_mpfr(rop, acc, MPFR_RNDN);
  mpfr_sqrt(rop, rop, MPFR_RNDN);

  fmpq_clear(acc);
  fmpz_clear(acc_den);
  fmpz_clear(acc_num);
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

/**
   fmpz_poly_t
*/

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

static inline int fmpz_poly_ideal_norm(fmpz_t norm, fmpz_poly_t f, fmpz_poly_t modulus) {
  fmpz_poly_resultant(norm, f, modulus);
}

/**
   Decide if <b_0,b_1> = <g>
 */

static inline int fmpz_poly_ideal_subset(fmpz_poly_t g, fmpz_poly_t b0, fmpz_poly_t b1, fmpz_poly_t modulus) {
  const long n = fmpz_poly_degree(modulus);
  fmpz_t det;
  fmpz_init(det);
  fmpz_poly_ideal_norm(det, g, modulus);

  fmpz_t det_b0, det_b1;
  fmpz_init(det_b0);
  fmpz_init(det_b1);

  fmpz_poly_ideal_norm(det_b0, b0, modulus);
  fmpz_poly_ideal_norm(det_b1, b1, modulus);

  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_gcd(tmp, det_b0, det_b1);
  fmpz_clear(det_b0);
  fmpz_clear(det_b1);

  int r = fmpz_cmp(det, tmp);

  fmpz_clear(det);
  fmpz_clear(tmp);
  return r;
}

/**

   
*/

static inline int fmpz_poly_ideal_is_probaprime(fmpz_poly_t op, fmpz_poly_t modulus) {
  fmpz_t det;
  fmpz_init(det);

  // TODO: Implement multi-modular resultant
  fmpz_poly_ideal_norm(det, op, modulus);

  int r = fmpz_is_probabprime(det);
  fmpz_clear(det);
  return r;
}

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
  assert(fmpz_mod_poly_degree(r) == 0 && fmpz_is_one(r->coeffs + 0));

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
