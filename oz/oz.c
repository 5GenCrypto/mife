#include <omp.h>
#include "flint-addons.h"
#include "util.h"
#include "oz.h"

static inline mp_bitcnt_t _fmpq_poly_oz_ideal_norm_bound(const fmpq_poly_t f, const long n) {
  mp_bitcnt_t bits1 = FLINT_ABS(_fmpz_vec_max_bits(f->coeffs, f->length));
  mp_bitcnt_t bound = 2*n*FLINT_BIT_COUNT((20*n + 26)/27) + 3;
  bound += (n - 1) + n*bits1;
  return bound;
}

void fmpq_poly_oz_ideal_norm(fmpq_t norm, const fmpq_poly_t f, const long n, const mpfr_prec_t prec) {
  fmpq_poly_t modulus;

  if (prec == 0) {
    mp_bitcnt_t bound = _fmpq_poly_oz_ideal_norm_bound(f, n);
    fmpq_poly_init_oz_modulus(modulus, n);
    fmpq_poly_resultant_modular_bound(norm, f, modulus, bound);
    fmpq_poly_clear(modulus);

  } else if  (prec < 0) {

    mpfr_t norm_f;
    mpfr_init2(norm_f, -prec);
    fmpq_poly_eucl_norm_mpfr(norm_f, f, MPFR_RNDN);
    mpfr_pow_ui(norm_f, norm_f, n, MPFR_RNDN);
    fmpq_set_mpfr(norm, norm_f, MPFR_RNDN);
    mpfr_clear(norm_f);

  } else {

    fmpq_poly_init_oz_modulus(modulus, n);
    fmpq_poly_t f_trunc;
    fmpq_poly_init(f_trunc);
    fmpq_poly_set(f_trunc, f);
    fmpq_poly_truncate(f_trunc, prec);

    mp_bitcnt_t bound = _fmpq_poly_oz_ideal_norm_bound(f_trunc, n);
    fmpq_poly_resultant_modular_bound(norm, f_trunc, modulus, bound);

    fmpq_poly_clear(modulus);
    fmpq_poly_clear(f_trunc);
  }
}

void fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const long n, const mpfr_prec_t prec) {
  mp_bitcnt_t bits = FLINT_ABS(_fmpz_vec_max_bits(f->coeffs, f->length));
  mp_bitcnt_t bound = f->length * (bits + n_clog(f->length, 2));
  if (prec == 0) {
    fmpz_poly_t modulus;
    fmpz_poly_init_oz_modulus(modulus, n);
    fmpz_poly_resultant_modular_bound(norm, f, modulus, bound);
    fmpz_poly_clear(modulus);
  } else {
    mpfr_t norm_f;
    mpfr_init2(norm_f, prec);
    fmpz_poly_eucl_norm_mpfr(norm_f, f, MPFR_RNDN);
    mpfr_pow_ui(norm_f, norm_f, n, MPFR_RNDN);
    mpz_t norm_z;
    mpz_init(norm_z);
    mpfr_get_z(norm_z, norm_f, MPFR_RNDN);
    fmpz_set_mpz(norm, norm_z);
    mpfr_clear(norm_f);
    mpz_clear(norm_z);
  }
}

void fmpz_poly_oz_rem(fmpz_poly_t rem, const fmpz_poly_t f, const long n) {
  fmpz_poly_set(rem, f);
  fmpz_t lead; fmpz_init(lead);
  fmpz_t coef; fmpz_init(coef);
  for(long d=fmpz_poly_degree(rem); d>=n; d--) {
    fmpz_poly_get_coeff_fmpz(lead, rem, d);
    fmpz_poly_get_coeff_fmpz(coef, rem, d-n);
    fmpz_sub(coef, coef, lead);
    fmpz_poly_set_coeff_si(rem, d, 0);
    fmpz_poly_set_coeff_fmpz(rem, d-n, coef);
  }
  fmpz_clear(coef);
  fmpz_clear(lead);
}

void fmpz_mod_poly_oz_rem(fmpz_mod_poly_t rem, const fmpz_mod_poly_t f, const long n) {
  fmpz_mod_poly_set(rem, f);
  fmpz_t lead; fmpz_init(lead);
  fmpz_t coef; fmpz_init(coef);
  for(long d=fmpz_mod_poly_degree(rem); d>=n; d--) {
    fmpz_mod_poly_get_coeff_fmpz(lead, rem, d);
    fmpz_mod_poly_get_coeff_fmpz(coef, rem, d-n);
    fmpz_sub(coef, coef, lead);
    fmpz_mod_poly_set_coeff_ui(rem, d, 0);
    fmpz_mod_poly_set_coeff_fmpz(rem, d-n, coef);
  }
  fmpz_clear(coef);
  fmpz_clear(lead);
}

void fmpq_poly_oz_rem(fmpq_poly_t rem, const fmpq_poly_t f, const long n) {

  const long d = fmpq_poly_degree(f);
  if(d > 2*n-1) {
    // TODO: Don't be so lazy and take care of this case as well
    fmpq_poly_t modulus;
    fmpq_poly_init_oz_modulus(modulus, n);
    fmpq_poly_rem(rem, f, modulus);
    fmpq_poly_clear(modulus);
  } else if (d>=n) {
    assert(d-n+1 <= n);
    fmpq_poly_t r;
    fmpq_poly_init(r);
    fmpq_poly_realloc(r, n);    
    mpq_t *z = (mpq_t*)calloc(n, sizeof(mpq_t));
    for (int i = 0; i < n; i++) {
      mpq_init(z[i]);
      fmpq_poly_get_coeff_mpq(z[i], f, i);
    }
    fmpq_poly_set_array_mpq(r, (const mpq_t*)z, n);
              
    for (int i=0; i<d-n+1; i++)
      fmpq_poly_get_coeff_mpq(z[i], f, n+i);
    fmpq_poly_t t;
    fmpq_poly_init(t);
    fmpq_poly_set_array_mpq(t, (const mpq_t*)z, d-n+1);
    fmpq_poly_sub(r, r, t);
    fmpq_poly_set(rem, r);

    fmpq_poly_clear(t);    
    fmpq_poly_clear(r);
    for (int i=0; i<n; i++)
      mpq_clear(z[i]);
    free(z);
  }
  return;
}

void fmpz_poly_oz_conjugate(fmpz_poly_t fT, const fmpz_poly_t f, const long n) {
  fmpz_poly_zero(fT);
  fmpz_t t0; fmpz_init(t0);
  fmpz_t t1; fmpz_init(t1);

  fmpz_poly_get_coeff_fmpz(t0, f, 0);
  fmpz_poly_set_coeff_fmpz(fT, 0, t0);

  for(int i=1; i<n/2; i++) {
    fmpz_poly_get_coeff_fmpz(t0, f,   i);
    fmpz_poly_get_coeff_fmpz(t1, f, n-i);
    fmpz_neg(t0, t0);
    fmpz_neg(t1, t1);
    fmpz_poly_set_coeff_fmpz(fT, n-i, t0);
    fmpz_poly_set_coeff_fmpz(fT, i, t1);
  }
  fmpz_poly_get_coeff_fmpz(t0, f, n/2);
  fmpz_neg(t0, t0);
  fmpz_poly_set_coeff_fmpz(fT, n/2, t0);

  fmpz_clear(t1);
  fmpz_clear(t0);
}

void fmpq_poly_oz_conjugate(fmpq_poly_t fT, const fmpq_poly_t f, const long n) {
  fmpq_poly_zero(fT);
  fmpq_t t0; fmpq_init(t0);
  fmpq_t t1; fmpq_init(t1);

  fmpq_poly_get_coeff_fmpq(t0, f, 0);
  fmpq_poly_set_coeff_fmpq(fT, 0, t0);

  for(int i=1; i<n/2; i++) {
    fmpq_poly_get_coeff_fmpq(t0, f,   i);
    fmpq_poly_get_coeff_fmpq(t1, f, n-i);
    fmpq_neg(t0, t0);
    fmpq_neg(t1, t1);
    fmpq_poly_set_coeff_fmpq(fT, n-i, t0);
    fmpq_poly_set_coeff_fmpq(fT, i, t1);
  }
  fmpq_poly_get_coeff_fmpq(t0, f, n/2);
  fmpq_neg(t0, t0);
  fmpq_poly_set_coeff_fmpq(fT, n/2, t0);

  fmpq_clear(t1);
  fmpq_clear(t0);
}

mp_limb_t *_fmpz_poly_oz_ideal_is_probaprime_small_primes(const long n, const int k) {
  mp_limb_t q = 1;
  mp_limb_t *small_primes = (mp_limb_t*)calloc(sizeof(mp_limb_t), k);
  small_primes[0] = 2;

  for(int i=1; i<k; ) {
    q += n;
    if (n_is_probabprime(q)) {
      small_primes[i] = q;
      i++;
    }
  }
  return small_primes;
}

int fmpz_poly_oz_ideal_is_probaprime(fmpz_poly_t f, const long n, int sloppy, const int k, const mp_limb_t *small_primes) {
  fmpz_poly_t g; fmpz_poly_init_oz_modulus(g, n);
  int num_threads = omp_get_max_threads();

  nmod_poly_t a[num_threads];
  nmod_poly_t b[num_threads];
  int r[num_threads];

  for(int j=0; j<num_threads; j++)
    r[j] = 1;

  for(int i=0; i<k; i+=num_threads) {
    if (k-i < num_threads)
      num_threads = k-i;
#pragma omp parallel for   
    for (int j=0; j<num_threads; j++) {
      mp_limb_t p = small_primes[i+j];
      nmod_poly_init(a[j], p); fmpz_poly_get_nmod_poly(a[j], f);
      nmod_poly_init(b[j], p); fmpz_poly_get_nmod_poly(b[j], g);
      r[j] = nmod_poly_resultant(a[j], b[j]);
      nmod_poly_clear(a[j]);
      nmod_poly_clear(b[j]);
    }
    for(int j=0; j<num_threads; j++)
      if (r[j]==0) {
        r[0] = 0;
        break;
      }
    if (r[0] == 0)
      break;
  }
  fmpz_poly_clear(g);
  if (sloppy)
    return r[0];

  if (r[0]) {
    fmpz_t norm;
    fmpz_init(norm);
    fmpz_poly_oz_ideal_norm(norm, f, n, 0);
    r[0] = fmpz_is_probabprime(norm);
    fmpz_clear(norm);
  }
  return r[0];
}
