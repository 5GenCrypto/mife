#include <assert.h>
#include <omp.h>
#include "norm.h"
#include "util.h"
#include "oz.h"
#include "flint-addons.h"

void _nmod_vec_oz_set_powers(mp_ptr op, const size_t n, const mp_limb_t w, const nmod_t q) {
  mp_limb_t acc = 1;
  op[0] = 1;
  for(size_t i=1; i<n; i++) {
    acc = n_mulmod2_preinv(acc, w, q.n, q.ninv);
    op[i] = acc;
  }
}

static const unsigned char bit_reverse_table_256[] =  {
  0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
  0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
  0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
  0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
  0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
  0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
  0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
  0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
  0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
  0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
  0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
  0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
  0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
  0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
  0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
  0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};


void _nmod_vec_oz_ntt(mp_ptr rop, const mp_ptr op, const mp_ptr w, const size_t n, const nmod_t q) {
  const size_t k = n_flog(n,2);

  mp_ptr a = _nmod_vec_init(n);

#if 0
  for (size_t i = 0; i < n; i++) {
    size_t ii = i, r=0;
    for (size_t h = 0; h < k; h++) {
      r = (r << 1) | (ii & 1);
      ii >>= 1;
    }
    a[r] = op[i];
  }
#else
  assert(k <= 24);
  assert(sizeof(int) >= 4);
  for (unsigned int i = 0; i < n; i++) {
    unsigned int r;
    r = (bit_reverse_table_256[i & 0xff] << 16) | \
      (bit_reverse_table_256[(i >>  8) & 0xff] << 8) | \
      (bit_reverse_table_256[(i >> 16) & 0xff]);
    r >>= (24 - k);
    a[r] = op[i];
  }
#endif

  mp_ptr b = _nmod_vec_init(n);

  const double ninv = n_precompute_inverse(q.n);
  for(size_t i=0; i<k; i++) {
    const mp_limb_t tkm = ~(((1UL)<<(k-1-i)) - 1);
    for(size_t j=0; j<n/2; j++) {
      const size_t pij = j & tkm;
      mp_limb_t tmp = n_mulmod_precomp(a[2*j+1], w[pij], q.n, ninv);
      b[j]      = n_addmod(a[2*j], tmp, q.n);
      b[j+n/2]  = n_submod(a[2*j], tmp, q.n);
    }
    if(i!=k-1)
      _nmod_vec_set(a, b, n);
  }
  _nmod_vec_set(rop, b, n);
  _nmod_vec_clear(b);
  _nmod_vec_clear(a);
}


mp_limb_t _nmod_vec_oz_resultant(const mp_ptr a, const long n, nmod_t q) {
  const mp_limb_t w_ = _nmod_nth_root(2*n, q.n);
  mp_ptr w = _nmod_vec_init(2*n);
  mp_ptr t = _nmod_vec_init(2*n);

  _nmod_vec_oz_set_powers(w, 2*n, w_, q);
  _nmod_vec_oz_ntt(t, a, w, 2*n, q);

  mp_limb_t acc = 1;
  for(int i=1; i<2*n; i+=2)
    acc = n_mulmod2_preinv(acc, t[i], q.n, q.ninv);

  _nmod_vec_clear(w);
  _nmod_vec_clear(t);
  return acc;
}


mp_limb_t nmod_poly_oz_resultant(const nmod_poly_t a, const long n) {
  nmod_t q;
  nmod_init(&q, nmod_poly_modulus(a));
  mp_ptr t = _nmod_vec_init(2*n);
  _nmod_vec_set(t, a->coeffs, a->length);
  for(long i=a->length; i<2*n; i++)
    t[i] = 0;
  mp_limb_t res = _nmod_vec_oz_resultant(t, n, q);
  _nmod_vec_clear(t);
  return res;
}

void _fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const long n) {
  mp_bitcnt_t bits = FLINT_ABS(_fmpz_vec_max_bits(f->coeffs, f->length));
  mp_bitcnt_t bound = f->length * (bits + n_clog(f->length, 2));

  mp_bitcnt_t pbits;
  slong i, num_primes;
  fmpz_comb_t comb;
  fmpz_comb_temp_t comb_temp;

  /* compute content of f */
  fmpz_t fc;  fmpz_init(fc);
  _fmpz_vec_content(fc, f->coeffs, n);

  /* divide f by content */
  fmpz *F = _fmpz_vec_init(n);
  _fmpz_vec_scalar_divexact_fmpz(F, f->coeffs, n, fc);

  /* get product of leading coefficients */
  fmpz_t l;  fmpz_init(l);
  fmpz_set(l, f->coeffs + n-1);

  /* set size of first prime */
  pbits = FLINT_D_BITS - 1;

  num_primes = (bound + pbits - 1)/pbits;
  mp_ptr parr = _nmod_vec_init(num_primes);
  mp_ptr rarr = _nmod_vec_init(num_primes);

  fmpz_zero(norm);

  /* make space for polynomials mod p */
  const int num_threads = omp_get_max_threads();

  mp_ptr a[num_threads];

  for(i=0; i<num_threads; i++) {
    a[i] = _nmod_vec_init(2*n);
    for(size_t j=0; j<2*n; j++)
      a[i][j] = 0;
  }

  mp_limb_t p = (UWORD(1)<<pbits) + 1;
  for(i=0; i<num_primes;) {
    p = _n_next_oz_good_probaprime(p, 2*n);
    if (fmpz_fdiv_ui(l, p) == 0)
      continue;
    parr[i++] = p;
  }

#pragma omp parallel for
  for (i = 0; i<num_primes; i++) {
    nmod_t mod;
    nmod_init(&mod, parr[i]);

    const int id = omp_get_thread_num();
    /* reduce polynomials modulo p */
    _fmpz_vec_get_nmod_vec(a[id], F, n, mod);
    /* compute resultant over Z/pZ */
    rarr[i] = _nmod_vec_oz_resultant(a[id], n, mod);
    flint_cleanup();
  }

  fmpz_comb_init(comb, parr, num_primes);
  fmpz_comb_temp_init(comb_temp, comb);

  fmpz_multi_CRT_ui(norm, rarr, comb, comb_temp, 1);

  fmpz_comb_temp_clear(comb_temp);
  fmpz_comb_clear(comb);

  for(i=0; i<num_threads; i++) {
    _nmod_vec_clear(a[i]);
  }

  _nmod_vec_clear(parr);
  _nmod_vec_clear(rarr);

  /* finally multiply by powers of content */
  if (!fmpz_is_one(fc)) {
    fmpz_pow_ui(l, fc, n - 1);
    fmpz_mul(norm, norm, l);
  }

  fmpz_clear(l);

  _fmpz_vec_clear(F, n);
  fmpz_clear(fc);
}


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
  if (prec == 0) {
    _fmpz_poly_oz_ideal_norm(norm, f, n);
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
