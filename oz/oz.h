/** @file oz.h
    @brief Computing in \f$\R\f$ and \f$\QQ[x]/\ideal{x^n+1}\f$, main include file.
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

/**
   @brief Initialise $f$ to $x^n+1 \\in \\ZZ[x]$.
*/

static inline void fmpz_poly_init_oz_modulus(fmpz_poly_t f, const long n) {
  fmpz_poly_init(f);
  fmpz_poly_set_coeff_si(f, 0, 1);
  fmpz_poly_set_coeff_si(f, n, 1);
}

/**
   @brief Initialise $f$ to $x^n+1 \\in \\QQ[x]$.
*/


static inline void fmpq_poly_init_oz_modulus(fmpq_poly_t f, const long n) {
  fmpq_poly_init(f);
  fmpq_poly_set_coeff_si(f, 0, 1);
  fmpq_poly_set_coeff_si(f, n, 1);
}

/**
   @brief Initialise $f$ to $x^n+1 \\in \\ZZ_q[x]$.
*/

static inline void fmpz_mod_poly_init_oz_modulus(fmpz_mod_poly_t f, const fmpz_t q, const long n) {
  fmpz_mod_poly_init(f, q);
  fmpz_mod_poly_set_coeff_ui(f, 0, 1);
  fmpz_mod_poly_set_coeff_ui(f, n, 1);
}

/**
   @brief Set $f^T$ to $f$'s conjugate in $\\R$.
*/

void fmpz_poly_oz_conjugate(fmpz_poly_t fT, const fmpz_poly_t f, const long n);

/**
   @brief Set $f^T$ to $f$'s conjugate in $\\QQ[x]/(x^n+1)$.
*/

void fmpq_poly_oz_conjugate(fmpq_poly_t fT, const fmpq_poly_t f, const long n);


/**
   @brief Return array with $k$ likely prime factors of ideal norms in $R$.
*/


mp_limb_t *_fmpz_poly_oz_ideal_probable_prime_factors(const long n, const size_t k);

/**
   @brief Return array with all primes smaller than `bound` with probable prime factors of ideals in $\\R$ first.
*/

mp_limb_t *_fmpz_poly_oz_ideal_small_prime_factors(const long n, const mp_limb_t bound);

/**
   \brief Return false if \f$\ideal{f}\f$ is not a prime ideal. If `sloppy=0` return true if it is
   probably a prime ideal.

   @param f            \f$f \in \R\f$
   @param n            degree of cyclotomic polynomial, must be power of two
   @param sloppy       only check modulo primes in `small_primes`
   @param primes       an array of probable prime factors which are checked first

  1. This function first computes \f$\res{f,x^n+1} \bmod p_i\f$ where $p_i$ with $0 ≤ i <
  |\\mbox{primes}|$ is a small prime in the list `primes`. If any of those computations produces $0$
  we return false. If `sloppy = 1` then true is returned after all small primes were exhausted and
  no zero was found.

  2. If `sloppy = 0` then we compute \f$\res{f,x^n+1}\f$ and check if it is probably prime.

  @see _fmpz_poly_oz_ideal_probable_prime_factors
*/

int fmpz_poly_oz_ideal_is_probaprime(const fmpz_poly_t f, const long n, int sloppy, const mp_limb_t *primes);

/**
   @brief Return true if \f$\N{f}\f$ has none of the elements of `primes` as a prime factor.

   @param f            \f$f \in \R\f$
   @param n            degree of cyclotomic polynomial, must be power of two
   @param primes       an array of possible prime factors

  This function computes \f$\res{f,x^n+1} \bmod p_i\f$ where $p_i$ with $0 ≤ i < |\\mbox{primes}|$
  is a prime in the list `primes`. If any of those computations produces $0$ we return
  false. Otherwise, return true.

  @see _fmpz_poly_oz_ideal_small_prime_factors
*/

int fmpz_poly_oz_ideal_not_prime_factors(const fmpz_poly_t f, const long n, const mp_limb_t *primes);

/**
   \brief Return true if @f$\ideal{b_0, b_1} = \ideal{g}@f$.

   @param g             an element $g$ in $\\R$
   @param b0            an element in $\\langle g \\rangle$, otherwise behaviour is undefined
   @param b1            an element in $\\ideal{ g }$, otherwise behaviour is undefined
   @param n             degree of cyclotomic polynomial, must be power of two
   @param sloppy        only check modulo primes in `small_primes`
   @param primes        a list of small primes which are checked first

   1. This function first checks if the resultants of $b_0$ resp. $b_1$ with $x^n+1$ is zero mod
      $p_i$ where $p_i$ with $0 ≤ i < |\\mbox{primes}|$ is a prime in the list `primes`. If `sloppy
      = 1` we assume that $g$ does not have small prime factors and return false.

   2. If this test passes and `sloppy = 0` we compute the the resultants \f$r_0 = \res{b_0,
      x^n+1}\f$, \f$r_1 = \res{b_1,x^n+1}\f$ and \f$r_g = \mbox{res}(g,x^n+1)\f$ and return true if
      \f$\gcd{r_0, r_1} = r_g = \N{g}\f$.

   @warning This function returns false negatives, i.e. it might return `0` even though
   \f$\ideal{b_0,b_1} = \ideal{g}\f$ holds, even when `sloppy=0`. This function checks if
   \f$\gcd{\N{b_0}, \N{b_1}} = \N{g}\f$, but \f$\ideal{b_0,b_1} = \ideal{g}\f$ may hold even if the
   gcd is a strict multiple of $\\N{g}$.
 */

int fmpz_poly_oz_ideal_span(const fmpz_poly_t g, const fmpz_poly_t b0, const fmpz_poly_t b1, const long n,
                            const int sloppy, const mp_limb_t *primes);


/**
   \brief Return true if @f$\ideal{b_0} + \ideal{b_1} = R@f$.

   @param b0            an element
   @param b1            an element
   @param n             degree of cyclotomic polynomial, must be power of two
   @param sloppy        only check modulo primes in `primes`
   @param primes        an array list of probable primes which are checked first

   1. This function checks if \f$\res{b_0,x^n+1} = 0 \bmod p_i\f$ and \f$\res{b_1,x^n+1} = 0 \bmod
      p_i\f$ for $0 ≤ i < k$ where $p_i$ is a small prime in the list `small_primes`. We return
      false in this case.

   2. If this test passes and `sloppy = 0` we compute the the resultants \f$r_0 = \res{b_0,
      x^n+1}\f$ and \f$r_1 = \res{b_1,x^n+1}\f$ and return true if \f$\gcd{r_0, r_1} = 1\f$ and
      false otherwise.

   @see _fmpz_poly_oz_ideal_probable_prime_factors
 */

int fmpz_poly_oz_coprime(const fmpz_poly_t b0, const fmpz_poly_t b1, const long n,
                         const int sloppy, const mp_limb_t *small_primes);


/**
   \brief Return true if the norm of @f$\ideal{b_0}@f$ and @f$det_{b_1}@f$ are co-prime

   @param b0            an element
   @param det_b1        the norm of another element
   @param n             degree of cyclotomic polynomial, must be power of two
   @param sloppy        only check modulo primes in `primes`
   @param primes        an array list of probable primes which are checked first

   1. This function checks if \f$\res{b_0,x^n+1} = 0 \bmod p_i\f$ for $0 ≤ i < k$ where $p_i$ is a
      small prime in the list `small_primes`. We return false in this case.

   2. If this test passes and `sloppy = 0` we compute the the resultants \f$r_0 = \res{b_0,
      x^n+1}\f$ and return true if \f$\gcd{r_0, det_{b_1}} = 1\f$ and false otherwise.

   @see _fmpz_poly_oz_ideal_probable_prime_factors
 */


int fmpz_poly_oz_coprime_det(const fmpz_poly_t b0, const fmpz_t det_b1, const long n,
                             const int sloppy, const mp_limb_t *primes);

#include "mul.h"
#include "ntt.h"
#include "invert.h"
#include "sqrt.h"
#include "norm.h"

#endif /* _OZ_H_ */
