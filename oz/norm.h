#ifndef NORM_H
#define NORM_H

#include <flint/ulong_extras.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>

static inline mp_limb_t _n_next_oz_good_probaprime(mp_limb_t p, const mp_limb_t n) {
  p += n;
  while (!n_is_probabprime(p)) {
    p += n;
  }
  return p;
}

static inline mp_limb_t _n_prev_oz_good_probaprime(mp_limb_t p, const mp_limb_t n) {
  p -= n;
  while (!n_is_probabprime(p)) {
    if (p>n)
      p -= n;
    else
      return 0;
  }
  return p;
}


static inline mp_limb_t _nmod_nth_root(const long n, mp_limb_t p) {
  mp_limb_t a = n_primitive_root_prime(p);
  mp_limb_t e = (p-1)/n;
  a = n_powmod(a, e, p);
  return a;
}

void nmod_poly_oz_set_powers(nmod_poly_t op, const size_t n, const mp_limb_t w);
void _nmod_poly_oz_ntt(nmod_poly_t rop, const nmod_poly_t op, const nmod_poly_t w, const size_t n);
mp_limb_t nmod_poly_oz_resultant(const nmod_poly_t a, const long n);

void fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const long n, mpfr_prec_t prec);
void _fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const long n);
void fmpq_poly_oz_ideal_norm(fmpq_t norm, const fmpq_poly_t f, const long n, const mpfr_prec_t prec);

#endif /* NORM_H */
