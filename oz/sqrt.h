#ifndef _SQRT_H_
#define _SQRT_H_

/******************************************************************************
*
*              OZ: Operations in 2^k-th Cyclotomic Number Fields
*
*    Copyright (C) 2014 Martin Albrecht <martinralbrecht+oz@googlemail.com>
*    Copyright (C) 2014 Catalin Cocis <catalincocis@yahoo.com >
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
#include <flint/fmpq_poly.h>

int fmpq_poly_oz_sqrt_approx_db(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init);
int fmpq_poly_oz_sqrt_approx_babylonian(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init);
int fmpq_poly_oz_sqrt_approx_pade(fmpq_poly_t f_sqrt, const fmpq_poly_t f, const long n, const int p, const mpfr_prec_t prec, const mpfr_prec_t prec_bound, uint64_t flags, const fmpq_poly_t init);

#endif /* _SQRT_H_ */
