/**
    @file rem.h
    @brief Computing remainders modulo $g \\in \\ZZ[x]/(x^n+1)$.
*/

#ifndef _REM_H
#define _REM_H

#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <oz/flags.h>

/**
   @brief Return a small representative of $f \\mod \\ideal{g}$.

   @param rem           return value, a small representative of $f \bmod \ideal{g}$.
   @param f             an element $f$ in $\\R$
   @param g             an element $g$ in $\\R$
   @param n             degree of cyclotomic polynomial, must be power of two
   @param ginv          pre-computed approximate inverse of $g$ in $\\R$.
 */

void _fmpz_poly_oz_rem_small(fmpz_poly_t rem, const fmpz_poly_t f, const fmpz_poly_t g, const long n, const fmpq_poly_t ginv);

/**
   @brief Return a small representative of $f \\mod \\ideal{g}$ with $f \\in \\Z$.

   @param rem           return value, a small representative of $f \bmod \ideal{g}$.
   @param f             an element $f$ in $\\Z$
   @param g             an element $g$ in $\\R$
   @param n             degree of cyclotomic polynomial, must be power of two
   @param ginv          pre-computed approximate inverse of $g$ in $\\R$.
   @param bound         log_2 of bound on $|rem|_∞$ (set to zero to disable)
 */

void _fmpz_poly_oz_rem_small_fmpz(fmpz_poly_t rem, const fmpz_t f, const fmpz_poly_t g, const long n,
                                  const fmpq_poly_t g_inv, const mp_bitcnt_t bound);

/**
   @brief Return a small representative of $f \\mod \\ideal{g}$ with $f \\in \\Z$.

   @param rem           return value, a small representative of $f \bmod \ideal{g}$.
   @param f             an element $f$ in $\\Z$
   @param g             an element $g$ in $\\R$
   @param n             degree of cyclotomic polynomial, must be power of two
   @param ginv          pre-computed approximate inverse of $g$ in $\\R$.
   @param b             process $f$ in chunks of size $b$ bits.
 */

void _fmpz_poly_oz_rem_small_fmpz_split(fmpz_poly_t rem, const fmpz_t f, const fmpz_poly_t g,
                                        const long n, const fmpq_poly_t g_inv, const mp_bitcnt_t b);

/**
   @brief Return a small representative of $f \\mod \\ideal{g}$.

   @param rem           return value, a small representative of $f \bmod \ideal{g}$.
   @param f             an element $f$ in \\R$
   @param g             an element $g$ in \\R$
   @param n             degree of cyclotomic polynomial, must be power of two
 */

void fmpz_poly_oz_rem_small(fmpz_poly_t rem, const fmpz_poly_t f, const fmpz_poly_t g, const long n);

/**
   @brief Return a small representative of $f \\mod \\ideal{g}$.

   @param rem           return value, a small representative of $f \bmod \ideal{g}$.
   @param f             an element $f$ in $\\R$
   @param g             an element $g$ in $\\R$
   @param n             degree of cyclotomic polynomial, must be power of two
   @param ginv          pre-computed approximate inverse of $g$ in $\\R$.
   @param prec          process $f$ in chunks of size $prec$ bits.
   @param flags         flags controlling verbosity et al.
 */

void _fmpz_poly_oz_rem_small_iter(fmpz_poly_t rem,
                                  const fmpz_poly_t f, const fmpz_poly_t g, const long n, const fmpq_poly_t ginv,
                                  const mp_bitcnt_t b, const oz_flag_t flags);

#endif /* _REM_H */
