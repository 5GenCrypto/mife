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

mp_limb_t *_fmpz_poly_oz_ideal_small_primes(const long n, const int k) {
  assert(k>=1);
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

int fmpz_poly_oz_ideal_is_probaprime(const fmpz_poly_t f, const long n, int sloppy, const int k, const mp_limb_t *small_primes) {
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


int fmpz_poly_oz_ideal_span(const fmpz_poly_t g, const fmpz_poly_t b0, const fmpz_poly_t b1, const long n,
                            const int sloppy, const int k, const mp_limb_t *small_primes) {

  fmpz_poly_t mod; fmpz_poly_init_oz_modulus(mod, n);
  int num_threads = omp_get_max_threads();

  nmod_poly_t n0[num_threads];
  nmod_poly_t n1[num_threads];
  nmod_poly_t nm[num_threads];

  int r0[num_threads];
  int r1[num_threads];

  for(int j=0; j<num_threads; j++) {
    r0[j] = 1;
    r1[j] = 1;
  }
  for(int i=0; i<k; i+=num_threads) {
    if (k-i < num_threads)
      num_threads = k-i;
#pragma omp parallel for
    for (int j=0; j<num_threads; j++) {
      mp_limb_t p = small_primes[i+j];
      nmod_poly_init(n0[j], p); fmpz_poly_get_nmod_poly(n0[j], b0);
      nmod_poly_init(n1[j], p); fmpz_poly_get_nmod_poly(n1[j], b1);
      nmod_poly_init(nm[j], p); fmpz_poly_get_nmod_poly(nm[j], mod);
      r0[j] = nmod_poly_resultant(n0[j], nm[j]);
      r1[j] = nmod_poly_resultant(n1[j], nm[j]);
      nmod_poly_clear(n0[j]);
      nmod_poly_clear(n1[j]);
    }
    for(int j=0; j<num_threads; j++) {
      /* if both resultants are zero we're in a sub-ideal as g is expected to not to be divisible by
         any small prime */
      if (r0[j] == 0 && r1[j] == 0) {
        r0[0] = 0;
        break;
      } else
        r0[0] = 1;
    }
    if (r0[0] == 0)
      break;
  }
  fmpz_poly_clear(mod);

  if (sloppy || r0[0] == 0)
    return (r0[0] != 0);

  fmpz_t det;
  fmpz_init(det);
  fmpz_poly_oz_ideal_norm(det, g, n, 0);

  fmpz_t det_b0, det_b1;
  fmpz_init(det_b0);
  fmpz_init(det_b1);

  fmpz_poly_oz_ideal_norm(det_b0, b0, n, 0);
  fmpz_poly_oz_ideal_norm(det_b1, b1, n, 0);

  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_gcd(tmp, det_b0, det_b1);
  fmpz_clear(det_b0);
  fmpz_clear(det_b1);

  r0[0] = fmpz_equal(det, tmp);

  fmpz_clear(det);
  fmpz_clear(tmp);
  return r0[0];
}
