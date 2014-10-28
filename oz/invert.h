#ifndef _INVERT_H_
#define _INVERT_H_

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

#include <stdio.h>
#include <stdint.h>
#include <mpfr.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mod_poly.h>

void _fmpq_poly_oz_invert_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, const long n, const mpfr_prec_t prec);

void fmpq_poly_oz_invert_approx(fmpq_poly_t rop, const fmpq_poly_t f, const long n, const mpfr_prec_t prec, const uint64_t flags);

void fmpz_mod_poly_oz_invert(fmpz_mod_poly_t rop, const fmpz_mod_poly_t f, const long n);

#endif /* _INVERT_H_ */
