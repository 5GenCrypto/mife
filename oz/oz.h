/** @file oz.h
*/
#ifndef _OZ_H_
#define _OZ_H_

#include <stdio.h>
#include <stdint.h>
#include <mpfr.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mod_poly.h>

typedef enum {
  OZ_VERBOSE    = 0x1, //< print debug messages
} oz_flag_t;

static inline void fmpz_poly_init_oz_modulus(fmpz_poly_t f, const long n) {
  fmpz_poly_init(f);
  fmpz_poly_set_coeff_si(f, 0, 1);
  fmpz_poly_set_coeff_si(f, n, 1);
}

static inline void fmpq_poly_init_oz_modulus(fmpq_poly_t f, const long n) {
  fmpq_poly_init(f);
  fmpq_poly_set_coeff_si(f, 0, 1);
  fmpq_poly_set_coeff_si(f, n, 1);
}

static inline void fmpz_mod_poly_init_oz_modulus(fmpz_mod_poly_t f, const fmpz_t q, const long n) {
  fmpz_mod_poly_init(f, q);
  fmpz_mod_poly_set_coeff_ui(f, 0, 1);
  fmpz_mod_poly_set_coeff_ui(f, n, 1);
}

/**
   Set $f^T$ to $f$'s conjugate.

   @param fT  pre-allocated output for \f$f^T \in \R\f$
   @param f   input \f$f \in \R\f$
   @param n   degree of cyclotomic polynomial, must be power of two
*/

void fmpz_poly_oz_conjugate(fmpz_poly_t fT, const fmpz_poly_t f, const long n);

void fmpq_poly_oz_conjugate(fmpq_poly_t fT, const fmpq_poly_t f, const long n);

void fmpq_poly_oz_ideal_norm(fmpq_t norm, const fmpq_poly_t f, const long n, const mpfr_prec_t prec);

void fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const long n, const mpfr_prec_t prec);

/**
   Return list of likely prime factors of ideal norms.
*/


mp_limb_t *_fmpz_poly_oz_ideal_small_primes(const long n, const int k);

/**
   \brief Return false if \f$\ideal{f}\f$ is not a prime ideal. If `sloppy=0` return true if it is
   probably a prime ideal.

   @param f            \f$f \in \R\f$
   @param n            degree of cyclotomic polynomial, must be power of two
   @param sloppy       only check modulo primes in `small_primes`
   @param k            number of primes in `small_primes`, must be > 1
   @param small_primes a list of small primes which are checked first

  1. This function first computes \f$\res{f,x^n+1} \bmod p_i\f$ where $p_i$ with $0 ≤ i < k$ is a
     small prime in the list `small_primes`. If any of those computations produces $0$ we return
     false. If `sloppy = 1` then true is returned after all small primes were exhausted and no zero
     was found.

  2. If `sloppy = 1` then we compute \f$\res{f,x^n+1}\f$ and check if it is probably prime.

  @see _fmpz_poly_oz_ideal_small_primes
*/

int fmpz_poly_oz_ideal_is_probaprime(const fmpz_poly_t f, const long n, int sloppy, const int k, const mp_limb_t *small_primes);

/**
   \brief Return true if @f$\ideal{b_0, b_1} = \ideal{g}@f$.

   @param g             an element $g$ in $\\R$
   @param b0            an element in $\\langle g \\rangle$, otherwise behaviour is undefined
   @param b1            an element in $\\ideal{ g }$, otherwise behaviour is undefined
   @param n             degree of cyclotomic polynomial, must be power of two
   @param sloppy        only check modulo primes in `small_primes`
   @param k             number of primes in `small_primes`, must be > 1
   @param small_primes  a list of small primes which are checked first

   1. This function checks if the resultants of $b_0$ resp. $b_1$ with $x^n+1$ is zero mod $p_i$
      where $p_i$ with $0 ≤ i < k$ is a small prime in the list `small_primes`. We assume that $g$
      does not have small prime factors in `small_primes` and return false in this case.

   2. If this test passes and `sloppy = 0` we compute the the resultants \f$r_0 = \res{b_0,
      x^n+1}\f$, \f$r_1 = \res{b_1,x^n+1}\f$ and \f$r_g = \mbox{res}(g,x^n+1)\f$ and return true if
      \f$\gcd{r_0, r_1} = r_g = \N{g}\f$.

   @warning This function returns false negatives, i.e. it might return `0` even though
   \f$\ideal{b_0,b_1} = \ideal{g}\f$ holds, even when `sloppy=0`. This function checks if
   \f$\gcd{\N{b_0}, \N{b_1}} = \N{g}\f$, but \f$\ideal{b_0,b_1} = \ideal{g}\f$ may hold even if the
   gcd is a strict multiple of $\\N{g}$.
 */

int fmpz_poly_oz_ideal_span(const fmpz_poly_t g, const fmpz_poly_t b0, const fmpz_poly_t b1, const long n,
                            const int sloppy, const int k, const mp_limb_t *small_primes);


/**
   \brief Return true if $b_0$ and $b_1$ are coprime, false otherwise.

   @param b0  a polynomial in \f$\ZZ[x]\f$
   @param b1  a polynomial in \f$\ZZ[x]\f$

   This function computes \f$\gcd{b_0,b_1}\f$ and compares the result with 1.
*/

static inline int fmpz_poly_oz_coprime(fmpz_poly_t b0, fmpz_poly_t b1) {
  fmpz_poly_t t;  fmpz_poly_init(t);
  fmpz_poly_gcd(t, b0, b1);
  int r = (fmpz_poly_degree(t) == 0);
  if (r)
    r = (fmpz_cmp_ui(t->coeffs, 1) == 0);
  fmpz_poly_clear(t);
  return r;
}


#include "mul.h"
#include "ntt.h"
#include "invert.h"
#include "sqrt.h"

#endif /* _OZ_H_ */
