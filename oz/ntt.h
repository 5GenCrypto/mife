#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include <stdio.h>
#include <mpfr.h>
#include <flint/fmpz_mod_poly.h>

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

void fmpz_mod_poly_oz_ntt_enc(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp);
void fmpz_mod_poly_oz_ntt_dec(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp);
void fmpz_mod_poly_oz_ntt_enc_fmpz_poly(fmpz_mod_poly_t rop, const fmpz_poly_t op, const fmpz_mod_poly_oz_ntt_precomp_t precomp);

#define fmpz_mod_poly_oz_ntt_add fmpz_mod_poly_add

void fmpz_mod_poly_oz_ntt_set_ui(fmpz_mod_poly_t op, const unsigned long c, const size_t n);
void fmpz_mod_poly_oz_ntt_mul(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n);
void fmpz_mod_poly_oz_ntt_inv(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const size_t n);
void fmpz_mod_poly_oz_ntt_pow_ui(fmpz_mod_poly_t rop, const fmpz_mod_poly_t f, unsigned long e, const size_t n);

void _fmpz_mod_poly_oz_ntt(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, const fmpz_mod_poly_t w, const size_t n);

void _fmpz_mod_poly_oz_mul_nttnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_poly_oz_ntt_precomp_t precomp);
void fmpz_mod_poly_oz_mul_nttnwc(fmpz_mod_poly_t h, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const size_t n);

#endif /* NTT_H */
