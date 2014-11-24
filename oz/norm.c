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

void _nmod_vec_oz_ntt(mp_ptr rop, const mp_ptr op, const mp_ptr w, const size_t n, const nmod_t q) {
  const size_t k = n_flog(n,2);

  mp_ptr a = _nmod_vec_init(n);
  for(size_t i=0; i<n; i++)
    a[i] = 0;

  for (size_t i = 0; i < n; i++) {
    size_t ii = i, r=0;
    for (size_t h = 0; h < k; h++) {
      r = (r << 1) | (ii & 1);
      ii >>= 1;
    }
    a[r] = op[i];
  }

  mp_ptr b = _nmod_vec_init(n);
  for(size_t i=0; i<n; i++)
    b[i] = 0;

  for(size_t i=0; i<k; i++) {
    for(size_t j=0; j<n/2; j++) {
      const size_t tk  = (1UL<<(k-1-i));
      const size_t pij = (j/tk) * tk;
      mp_limb_t tmp = n_mulmod2_preinv(a[2*j+1], w[pij], q.n, q.ninv);
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

void _fmpz_poly_oz_ideal_norm(fmpz_t norm, const fmpz_poly_t f, const size_t n) {
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
#if 0
    mp_bitcnt_t bits = FLINT_ABS(_fmpz_vec_max_bits(f->coeffs, f->length));
    mp_bitcnt_t bound = f->length * (bits + n_clog(f->length, 2));
    fmpz_poly_t modulus;
    fmpz_poly_init_oz_modulus(modulus, n);
    fmpz_poly_resultant_modular_bound(norm, f, modulus, bound);
    fmpz_poly_clear(modulus);
#else
    _fmpz_poly_oz_ideal_norm(norm, f, n);
#endif
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
