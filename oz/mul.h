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

void fmpz_poly_oz_rem(fmpz_poly_t rem, const fmpz_poly_t f, const long n);
void fmpq_poly_oz_rem(fmpq_poly_t rem, const fmpq_poly_t f, const long n);
void fmpz_mod_poly_oz_rem(fmpz_mod_poly_t rem, const fmpz_mod_poly_t f, const long n);

static inline void fmpz_poly_oz_mul(fmpz_poly_t rop, const fmpz_poly_t op1, const fmpz_poly_t op2, const long n) {
  fmpz_poly_mul(rop, op1, op2);
  fmpz_poly_oz_rem(rop, rop, n);
}

static inline void fmpq_poly_oz_mul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2, const long n) {
  fmpq_poly_mul(rop, op1, op2);
  fmpq_poly_oz_rem(rop, rop, n);
}

static inline void fmpz_mod_poly_oz_mul(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op1, const fmpz_mod_poly_t op2, const long n) {
  fmpz_mod_poly_mul(rop, op1, op2);
  fmpz_mod_poly_oz_rem(rop, rop, n);
}

struct fmpz_mod_poly_oz_fft_precomp_struct {
  size_t n;
  fmpz_mod_poly_t w;
  fmpz_mod_poly_t w_inv;
  fmpz_mod_poly_t phi;
  fmpz_mod_poly_t phi_inv;
};

typedef struct fmpz_mod_poly_oz_fft_precomp_struct fmpz_mod_poly_oz_fft_precomp_t[1];

void fmpz_mod_poly_oz_fft_precomp_init(fmpz_mod_poly_oz_fft_precomp_t op, const size_t n, const fmpz_t q);
void fmpz_mod_poly_oz_fft_precomp_clear(fmpz_mod_poly_oz_fft_precomp_t op);


void fmpz_mod_poly_oz_fft (fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const size_t n);
void fmpz_mod_poly_oz_ifft(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const size_t n);

void _fmpz_mod_poly_oz_fft    (fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_t w, const size_t n);

void _fmpz_mod_poly_oz_mul_fftnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_poly_oz_fft_precomp_t precomp);
void fmpz_mod_poly_oz_mul_fftnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n);

#endif /* _MUL_H_ */
