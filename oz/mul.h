#ifndef _MUL_H_
#define _MUL_H_

#include <stdint.h>
#include <stdio.h>
#include <mpfr.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpq_poly.h>

/**
   Set `r` to `f` modulo `x^n + 1`

   :param r: return polynomial in coefficient representation
   :param f: polynomial in coefficient representation
   :param n: power of two
*/

void fmpz_poly_oz_rem(fmpz_poly_t r, const fmpz_poly_t f, const long n);
void fmpq_poly_oz_rem(fmpq_poly_t r, const fmpq_poly_t f, const long n);
void fmpz_mod_poly_oz_rem(fmpz_mod_poly_t rem, const fmpz_mod_poly_t f, const long n);

/**
   Set `r` to `f · g` modulo `x^n + 1`

   :param r: return polynomial in coefficient representation
   :param f: multiplicant in coefficient representation
   :param g: multiplicant in coefficient representation
   :param n: power of two
*/

static inline void fmpz_poly_oz_mul(fmpz_poly_t r, const fmpz_poly_t f, const fmpz_poly_t g, const long n) {
  fmpz_poly_mul(r, f, g);
  fmpz_poly_oz_rem(r, r, n);
}

static inline void fmpq_poly_oz_mul(fmpq_poly_t r, const fmpq_poly_t f, const fmpq_poly_t g, const long n) {
  fmpq_poly_mul(r, f, g);
  fmpq_poly_oz_rem(r, r, n);
}

static inline void fmpz_mod_poly_oz_mul(fmpz_mod_poly_t r, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const long n) {
  fmpz_mod_poly_mul(r, f, g);
  fmpz_mod_poly_oz_rem(r, r, n);
}

#endif /* _MUL_H_ */
