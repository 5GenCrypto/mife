#ifndef _MUL_H_
#define _MUL_H_

/******************************************************************************
*
*              OZ: Operations in 2^k-th Cyclotomic Number Fields
*
*    Copyright (C) 2014 Martin Albrecht <martinralbrecht+oz@googlemail.com>
*
*  Distributed under the terms of the GNU General Public License (GPL)
*  version 2 or higher.
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/

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

/**
   Precomputed data for number-theoretic transform

   FIELDS:

   - `n` - dimension, must be a  power of two
   - `w` - a vector holding `ω_n^i` at index `i` where `ω_n` as a `n`-th root of
     unity
   - `w` - a vector holding `ω_n^-i` at index `i` where `ω_n` as a `n`-th root
     of unity
   - `phi` - a vector holding `φ^i` at index `i` where `φ = sqrt(ω_n) mod q`
   - `phi_inv` - a vector holding `φ^-i` at index `i` where `φ = sqrt(ω_n) mod q`
*/

struct fmpz_mod_poly_oz_ntt_precomp_struct {
  size_t n;
  fmpz_mod_poly_t w;
  fmpz_mod_poly_t w_inv;
  fmpz_mod_poly_t phi;
  fmpz_mod_poly_t phi_inv;
};

typedef struct fmpz_mod_poly_oz_ntt_precomp_struct fmpz_mod_poly_oz_ntt_precomp_t[1];

void fmpz_mod_poly_oz_ntt_precomp_init(fmpz_mod_poly_oz_ntt_precomp_t op, const size_t n, const fmpz_t q);
void fmpz_mod_poly_oz_ntt_precomp_clear(fmpz_mod_poly_oz_ntt_precomp_t op);

void fmpz_mod_poly_oz_ntt (fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const size_t n);
void fmpz_mod_poly_oz_intt(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const size_t n);

void _fmpz_mod_poly_oz_ntt_dec(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp);
void _fmpz_mod_poly_oz_ntt_set_ui(fmpz_mod_poly_t op, const unsigned long c, const size_t n);
void _fmpz_mod_poly_oz_ntt_mul(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n);
#define _fmpz_mod_poly_oz_ntt_add fmpz_mod_poly_add
void _fmpz_mod_poly_oz_ntt_pow_ui(fmpz_mod_poly_t rop, const fmpz_mod_poly_t f, unsigned long e, const size_t n);
void _fmpz_mod_poly_oz_ntt_enc(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp);

void _fmpz_mod_poly_oz_ntt(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_t w, const size_t n);

void _fmpz_mod_poly_oz_mul_nttnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_poly_oz_ntt_precomp_t precomp);
void fmpz_mod_poly_oz_mul_nttnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n);

#endif /* _MUL_H_ */
