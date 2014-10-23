#ifndef _OZ_H_
#define _OZ_H_

#include <stdio.h>
#include <stdint.h>
#include <mpfr.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>

typedef enum {
  OZ_VERBOSE    = 0x1, //< print debug messages
  OZ_BABYLONIAN = 0x2, //< use Babylonian method for sqrt
} oz_flag_t;

static inline void fmpz_poly_oz_init_modulus(fmpz_poly_t f, const long n) {
  fmpz_poly_init(f);
  fmpz_poly_set_coeff_si(f, 0, 1);
  fmpz_poly_set_coeff_si(f, n, 1);
}


static inline void fmpq_poly_oz_init_modulus(fmpq_poly_t f, const long n) {
  fmpq_poly_init(f);
  fmpq_poly_set_coeff_si(f, 0, 1);
  fmpq_poly_set_coeff_si(f, n, 1);
}

void fmpz_poly_oz_rem(fmpz_poly_t rem, const fmpz_poly_t f, const long n);

void fmpq_poly_oz_rem(fmpq_poly_t rem, const fmpq_poly_t f, const long n);

static inline void fmpz_poly_oz_mul(fmpz_poly_t rop, const fmpz_poly_t op1, const fmpz_poly_t op2, const long n) {
  fmpz_poly_mul(rop, op1, op2);
  fmpz_poly_oz_rem(rop, rop, n);
}

static inline void fmpq_poly_oz_mul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2, const long n) {
  fmpq_poly_mul(rop, op1, op2);
  fmpq_poly_oz_rem(rop, rop, n);
}

void _fmpq_poly_oz_invert_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, const int n, const mpfr_prec_t prec);

void fmpq_poly_oz_invert_approx(fmpq_poly_t rop, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const uint64_t flags);

int fmpq_poly_oz_sqrt_approx_db(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init);
int fmpq_poly_oz_sqrt_approx_babylonian(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init);
int fmpq_poly_oz_sqrt_approx_pade(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const int p, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init);

/*
  Set `fT` to the conjugate of `f`.

  This is equivalent to taking the transpose of the matrix associated with `f`.
*/

void fmpz_poly_oz_conjugate(fmpz_poly_t fT, const fmpz_poly_t f, const long n);

/*
  Set `fT` to the conjugate of `f`.

  This is equivalent to taking the transpose of the matrix associated with `f`.
*/

void fmpq_poly_oz_conjugate(fmpq_poly_t fT, const fmpq_poly_t f, const long n);

void fmpq_poly_oz_ideal_norm(fmpq_t norm, const fmpq_poly_t f, const long n, const mpfr_prec_t prec);

void fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const long n, const mpfr_prec_t prec);

mp_limb_t *_fmpz_poly_oz_ideal_is_probaprime_small_primes(const long n, const int k);

int fmpz_poly_oz_ideal_is_probaprime(fmpz_poly_t f, fmpz_poly_t g, int sloppy, const int k, const mp_limb_t *small_primes);

/**
   Decide if <b_0,b_1> = <g>
 */

static inline int fmpz_poly_oz_ideal_subset(fmpz_poly_t g, fmpz_poly_t b0, fmpz_poly_t b1, const long n) {
  fmpz_t det;
  fmpz_init(det);
  fmpz_poly_oz_ideal_norm(det, g, n, 0);

  fmpz_t det_b0, det_b1;
  fmpz_init(det_b0);
  fmpz_init(det_b1);

  fmpz_poly_oz_ideal_norm(det_b0, b0, n, 0);
  fmpz_poly_oz_ideal_norm(det_b1, b1, n, 0);

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


#endif /* _OZ_H_ */
