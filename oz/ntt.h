/**
   @file ntt.h
   @brief Computing with the number-theoretic transform.

   Let @f$ω_n@f$ be an $n$th root of unity in @f$\ZZ_q@f$ and @f$φ^2 = ω_n@f$.  Let @f$a =
   \sum_{i=0}^{n-1} a_i x^i@f$, @f$b = \sum_{i=0}^{n-1} b_i x^i@f$ and @f$c = a · b \in
   \ZZ_q[x]/\ideal{x^n+1}@f$. Let @f$\overline{a} = (a_0, φa_1, \dots, φ^{n-1}a_{n-1})@f$ and define
   @f$\overline{b}@f$ and @f$\overline{c}@f$ analogously. Then @f$\overline{c} = 1/n ·
   \mbox{NTT}_{ω_n}^{-1}(\mbox{NTT}_{ω_n}(\overline{a}) \odot \mbox{NTT}_{ω_n}(\overline{b}))@f$
   where @f$\mbox{NTT}_{ω_n}(·)@f$ is the number-theoretic transform and
   @f$\mbox{NTT}_{ω_n}^{-1}(·)@f$ is its inverse.
 */

#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include <stdio.h>
#include <mpfr.h>
#include <flint/fmpz_mod_poly.h>

/**
   @brief Pre-computed data for number-theoretic transform
*/

struct fmpz_mod_poly_oz_ntt_precomp_struct {
  size_t n;                   //!< dimension, must be a  power of two
  fmpz_mod_poly_t w;          //!< a vector holding $ω_n^i$ at index $i$ where $ω_n$ as an $n$-th root of unity.
  fmpz_mod_poly_t w_inv;      //!< a vector holding $ω_n^{-i}$ at index $i$ where $ω_n$ as an $n$-th root of unity.
  fmpz_mod_poly_t phi;        //!< a vector holding $φ^i$ at index $i$ where @f$φ = \sqrt{ω_n} \bmod q@f$.
  fmpz_mod_poly_t phi_inv;    //!< a vector holding $φ^{-i}$ at index $i$ where @f$φ = \sqrt{ω_n} \bmod q@f$.
};

typedef struct fmpz_mod_poly_oz_ntt_precomp_struct fmpz_mod_poly_oz_ntt_precomp_t[1];

/**
   @brief Pre-compute NTT data for $\\ZZ_q[x]/\\ideal{x^n+1}$.
*/

void fmpz_mod_poly_oz_ntt_precomp_init(fmpz_mod_poly_oz_ntt_precomp_t op, const size_t n, const fmpz_t q);

/**
   @brief Clear pre-computed data.
*/

void fmpz_mod_poly_oz_ntt_precomp_clear(fmpz_mod_poly_oz_ntt_precomp_t op);

/**
   @brief Compute @f$\mbox{rop} = \NTT{\mbox{op}}@f$.
*/

void fmpz_mod_poly_oz_ntt (fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const size_t n);

/**
   @brief Compute @f$\mbox{rop} = \INTT{\mbox{op}}@f$.
*/

void fmpz_mod_poly_oz_intt(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const size_t n);

/**
   @brief Compute @f$\mbox{rop} = \NTT{\mbox{op}}@f$ using `precomp`.
*/

void fmpz_mod_poly_oz_ntt_enc(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp);

/**
   @brief Compute @f$\mbox{rop} = \NTT{\mbox{op}}@f$ using `precomp` for @f$\mbox{op} \in \ZZ[x]/\ideal{x^n+1}@f$.
*/

void fmpz_mod_poly_oz_ntt_enc_fmpz_poly(fmpz_mod_poly_t rop, const fmpz_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp);

/**
   @brief Compute @f$\mbox{rop} = \INTT{\mbox{op}}@f$ using `precomp`.
*/

void fmpz_mod_poly_oz_ntt_dec(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp);

/**
   @brief Compute @f$\mbox{rop} = \NTT{c}@f$ using `precomp`.
*/

void fmpz_mod_poly_oz_ntt_set_ui(fmpz_mod_poly_t rop, const unsigned long c, const size_t n);

/**
   @brief Compute $h = \\NTT{f'+g'}$ from $f = \\NTT{f'}$ and $g = \\NTT{g'}$.
*/

#define fmpz_mod_poly_oz_ntt_add fmpz_mod_poly_add

/**
   @brief Compute $h = \\NTT{f' · g'}$  where $f',g' \\in \\ZZ_q[x]/\\ideal{x^n+1}$ from $f = \\NTT{f'}$ and $g = \\NTT{g'}$
*/

void fmpz_mod_poly_oz_ntt_mul(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n);

/**
   @brief Compute $h = \\NTT{f'^{-1}}$  where $f' \\in \\ZZ_q[x]/\\ideal{x^n+1}$ from $f = \\NTT{f'}$.
*/

void fmpz_mod_poly_oz_ntt_inv(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const size_t n);

/**
   @brief Compute $h = \\NTT{f'^e}$  where $f' \\in \\ZZ_q[x]/\\ideal{x^n+1}$ from $f = \\NTT{f'}$.
*/

void fmpz_mod_poly_oz_ntt_pow_ui(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, unsigned long e, const size_t n);

/**
   @brief Perform $\\mbox{rop} = \\NTT{\\mbox{op}}$ given $w = (1,ω,ω^2,…,ω^{n-1})$.
*/

void _fmpz_mod_poly_oz_ntt(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_t w, const size_t n);

/**
   @brief Compute $h = f · g$ using the number-theoretic transform using `precomp`.
*/

void _fmpz_mod_poly_oz_mul_nttnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_poly_oz_ntt_precomp_t precomp);

/**
   @brief Compute $h = f · g$ using the number-theoretic transform.
*/

void fmpz_mod_poly_oz_mul_nttnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n);

#endif /* NTT_H */
