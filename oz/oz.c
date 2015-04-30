#include <omp.h>
#include "flint-addons.h"
#include "util.h"
#include "oz.h"


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


mp_limb_t *_fmpz_poly_oz_ideal_probable_prime_factors(const long n, const size_t k) {
  assert(k>=1);
  mp_limb_t q = 1;
  mp_limb_t *primes = (mp_limb_t*)calloc(sizeof(mp_limb_t), k+1);
  primes[0] = k;
  primes[1] = 2;

  for(size_t i=2; i<k+1; ) {
    q += n;
    if (n_is_probabprime(q)) {
      primes[i] = q;
      i++;
    }
  }
  return primes;
}


mp_limb_t *_fmpz_poly_oz_ideal_small_prime_factors(const long n, const mp_limb_t bound) {
  size_t kmax = ceil(bound/log((double)bound));
  size_t i;

  mp_limb_t *primes = (mp_limb_t*)malloc(sizeof(mp_limb_t) * kmax+1);
  mp_limb_t *t = NULL;
  primes[1] = 2;

  mp_limb_t q = 1;
  for(i=2; i<kmax+1 && q < bound; ) {
    q += n;
    if (n_is_probabprime(q)) {
      primes[i] = q;
      i++;
      if (i == kmax) {
        kmax = 2*kmax;
        t = realloc(primes, sizeof(mp_limb_t) * kmax);
        if (t == NULL)
          oz_die("Not enough memory");
        else
          primes = t;
      }
    }
  }
  q = 2;
  for( ; q < bound; ) {
    q = n_nextprime(q, 0);
    if (q%n != 1) {
      primes[i] = q;
      i++;
      if (i == kmax) {
        kmax = 2*kmax;
        t = realloc(primes, sizeof(mp_limb_t) * kmax);
        if (t == NULL)
          oz_die("Not enough memory");
        else
          primes = t;
      }
    }
  }
  t = realloc(primes, sizeof(mp_limb_t) * i);
  if (t == NULL)
    oz_die("memory allocation error");
  else
    primes = t;
  primes[0] = i-1;
  return primes;
}

int fmpz_poly_oz_ideal_is_probaprime(const fmpz_poly_t f, const long n, int sloppy, const mp_limb_t *primes) {
  int r = fmpz_poly_oz_ideal_not_prime_factors(f, n, primes);
  if (r) {
    fmpz_t norm;
    fmpz_init(norm);
    fmpz_poly_oz_ideal_norm(norm, f, n, 0);
    r = fmpz_is_probabprime(norm);
    fmpz_clear(norm);
  }
  return r;
}

int fmpz_poly_oz_ideal_not_prime_factors(const fmpz_poly_t f, const long n, const mp_limb_t *primes) {
  fmpz_poly_t g; fmpz_poly_init_oz_modulus(g, n);
  int num_threads = omp_get_max_threads();

  const size_t k = primes[0];

  nmod_poly_t a[num_threads];
  nmod_poly_t b[num_threads];
  int r[num_threads];

  for(int j=0; j<num_threads; j++)
    r[j] = 1;

  for(size_t i=0; i<k; i+=num_threads) {
    if (k-i < (unsigned long)num_threads)
      num_threads = k-i;
#pragma omp parallel for
    for (int j=0; j<num_threads; j++) {
      mp_limb_t p = primes[1+i+j];
      nmod_poly_init(a[j], p); fmpz_poly_get_nmod_poly(a[j], f);
      nmod_poly_init(b[j], p); fmpz_poly_get_nmod_poly(b[j], g);
      if(p%(2*n) == 1)
        r[j] = nmod_poly_oz_resultant(a[j], n);
      else
        r[j] = nmod_poly_resultant(a[j], b[j]);
      nmod_poly_clear(a[j]);
      nmod_poly_clear(b[j]);
      flint_cleanup();
    }
    for(int j=0; j<num_threads; j++)
      if (r[j] == 0) {
        r[0] = 0;
        break;
      }
    if (r[0] == 0)
      break;
  }
  fmpz_poly_clear(g);
  return r[0];
}

int fmpz_poly_oz_ideal_span(const fmpz_poly_t g, const fmpz_poly_t b0, const fmpz_poly_t b1, const long n,
                            const int sloppy, const mp_limb_t *primes) {

  fmpz_poly_t mod; fmpz_poly_init_oz_modulus(mod, n);
  int num_threads = omp_get_max_threads();

  const long k = primes[0];

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
      mp_limb_t p = primes[1+i+j];
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

void _fmpz_poly_oz_rem_small(fmpz_poly_t rem, const fmpz_poly_t f, const fmpz_poly_t g, const long n, const fmpq_poly_t g_inv) {
  fmpz_poly_t fc; fmpz_poly_init(fc); fmpz_poly_set(fc, f);
  fmpq_poly_t fq; fmpq_poly_init(fq);
  fmpq_poly_set_fmpz_poly(fq, f);
  fmpq_poly_oz_rem(fq, fq, n);

  fmpq_poly_oz_mul(fq, g_inv, fq, n);

  fmpz_t t; fmpz_init(t);
  for(int i=0; i<fmpq_poly_length(fq); i++) {
    fmpz_tdiv_q(t, fq->coeffs + i, fq->den);
    fmpz_poly_set_coeff_fmpz(rem, i, t);
  }
  fmpz_clear(t);

  fmpz_poly_oz_mul(rem, rem, g, n);
  fmpz_poly_sub(rem, fc, rem);

  fmpz_poly_clear(fc);
  fmpq_poly_clear(fq);
}


void _fmpz_poly_oz_rem_small_fmpz(fmpz_poly_t rem, const fmpz_t f, const fmpz_poly_t g, const long n,
                                  const fmpq_poly_t g_inv, const mp_bitcnt_t bound) {
  if (fmpz_is_zero(f)) {
    fmpz_poly_set_ui(rem, 0);
    return;
  }
  fmpz_t fc; fmpz_init_set(fc, f);
  fmpq_poly_t fq; fmpq_poly_init(fq);
  /* fmpq_poly_scalar_mul_fmpz(fq, g_inv, f); but without all the gcd computations */
  fmpq_poly_realloc(fq, fmpq_poly_length(g_inv));
  _fmpz_vec_scalar_mul_fmpz(fq->coeffs, g_inv->coeffs, fmpq_poly_length(g_inv), f);

  fmpz_t t; fmpz_init(t);
  for(int i=0; i<fmpq_poly_length(g_inv); i++) {
    fmpz_tdiv_q(t, fq->coeffs + i, g_inv->den);
    fmpz_poly_set_coeff_fmpz(rem, i, t);
  }
  fmpz_clear(t);

  if (bound) {
    /* we know the result is small, so compute modulo the bound on the size */
    fmpz_poly_t tmp; fmpz_poly_init(tmp);
    fmpz_t q;
    fmpz_init(q);
    fmpz_setbit(q, bound);
    fmpz_mod_poly_t rem_q;
    while (1) {
      fmpz_mod_poly_init(rem_q, q);
      fmpz_mod_poly_set_fmpz_poly(rem_q, rem);

      fmpz_mod_poly_t g_q;
      fmpz_mod_poly_init(g_q, q);
      fmpz_mod_poly_set_fmpz_poly(g_q, g);

      fmpz_mod_poly_oz_mul(rem_q, rem_q, g_q, n);
      fmpz_mod_poly_neg(rem_q, rem_q);
      if(fmpz_mod_poly_degree(rem_q)>-1) {
        fmpz_add(rem_q->coeffs, fc, rem_q->coeffs);
        fmpz_mod(rem_q->coeffs, rem_q->coeffs, q);
      } else {
        fmpz_mod_poly_set_coeff_fmpz(rem_q, 0,fc);
        fmpz_mod(rem_q->coeffs, rem_q->coeffs, q);
      }

      fmpz_poly_set_fmpz_mod_poly(tmp, rem_q);
      fmpz_mod_poly_clear(rem_q);
      fmpz_mod_poly_clear(g_q);

      /* if there was no overflow, break*/
      if(labs(fmpz_poly_max_bits(tmp)) < (long)fmpz_sizeinbase(q, 2)-2)
        break;
      fmpz_mul_2exp(q, q, bound);
    }
    fmpz_poly_set(rem, tmp);
    fmpz_poly_clear(tmp);
    fmpz_clear(q);
  } else {
    fmpz_poly_oz_mul(rem, rem, g, n);
    fmpz_poly_neg(rem, rem);
    if(fmpz_poly_degree(rem)>-1)
      fmpz_add(rem->coeffs, fc, rem->coeffs);
    else
      fmpz_poly_set_coeff_fmpz(rem, 0,fc);
  }
  fmpq_poly_clear(fq);
  fmpz_clear(fc);
}

void _fmpz_poly_oz_rem_small_fmpz_split(fmpz_poly_t rem, const fmpz_t f, const fmpz_poly_t g,
                                        const long n, const fmpq_poly_t g_inv, const mp_bitcnt_t b) {

  const size_t num_threads = omp_get_max_threads();

  fmpz_t F; fmpz_init_set(F, f);
  fmpz_t H; fmpz_init(H);
  fmpz_poly_t t; fmpz_poly_init(t);
  fmpz_poly_set_ui(t, 1);
  fmpz_poly_t acc; fmpz_poly_init(acc);

  fmpz_t H_[num_threads];
  fmpz_poly_t f_[num_threads];
  fmpz_poly_t t_[num_threads];

  for(size_t j=0; j<num_threads; j++) {
    fmpz_init(H_[j]);
    fmpz_poly_init(f_[j]);
    fmpz_poly_init(t_[j]);
  }

  const mp_bitcnt_t B = num_threads*b;
  const mp_bitcnt_t rem_bound = log2(n)/2 * labs(fmpz_poly_max_bits(g));

  fmpz_poly_t powb; fmpz_poly_init(powb);
  fmpz_poly_set_coeff_ui(powb, 0, 2); // powb ~= 2^b
  fmpz_pow_ui(powb->coeffs, powb->coeffs, b);
  _fmpz_poly_oz_rem_small_fmpz(powb, powb->coeffs, g, n, g_inv, rem_bound);

  fmpz_poly_t powB; fmpz_poly_init(powB);
  fmpz_poly_set_ui(powB, 1);
  for(size_t j=0; j<num_threads; j++)
    fmpz_poly_oz_mul(powB, powB, powb, n);
  _fmpz_poly_oz_rem_small(powB, powB, g, n, g_inv);

  const size_t nparts = (fmpz_sizeinbase(f, 2)/B) + ((fmpz_sizeinbase(f, 2)%B) ? 1 : 0);

  for(size_t i=0; i<nparts; i++) {
    fmpz_set(H, F);
    fmpz_fdiv_r_2exp(H, H, B); // H = H % 2^B

    fmpz_poly_set(t_[0], t);
    for(size_t j=1; j<num_threads; j++)
      fmpz_poly_oz_mul(t_[j], t_[j-1], powb, n);

#pragma omp parallel for
    for(size_t j=0; j<num_threads; j++) {
      fmpz_set(H_[j], H);
      fmpz_fdiv_q_2exp(H_[j], H_[j], j*b);
      fmpz_fdiv_r_2exp(H_[j], H_[j], b); // H_j = (H >> j*b) % 2^b

      _fmpz_poly_oz_rem_small_fmpz(f_[j], H_[j], g, n, g_inv, rem_bound); // f_j ~= H_j
      fmpz_poly_oz_mul(f_[j], t_[j], f_[j], n); // f_j ~= 2^(b*j) * H_j
      flint_cleanup();
    }

    for(size_t j=0; j<num_threads; j++)
      fmpz_poly_add(acc, acc, f_[j]);

    fmpz_poly_oz_mul(t, t, powB, n);
    if (labs(fmpz_poly_max_bits(t)) > (long)b/2)
      _fmpz_poly_oz_rem_small(t, t, g, n, g_inv);

    fmpz_fdiv_q_2exp(F, F, B); // F >> B
  }

  fmpz_poly_set(rem, acc);

  fmpz_clear(F);
  fmpz_clear(H);
  fmpz_poly_clear(acc);
  fmpz_poly_clear(t);
  fmpz_poly_clear(powb);
  fmpz_poly_clear(powB);

  for(size_t j=0; j<num_threads; j++) {
    fmpz_clear(H_[j]);
    fmpz_poly_clear(f_[j]);
    fmpz_poly_clear(t_[j]);
  }

}

void fmpz_poly_oz_rem_small(fmpz_poly_t rem, const fmpz_poly_t f, const fmpz_poly_t g, const long n) {
  fmpq_poly_t gq; fmpq_poly_init(gq);
  fmpq_poly_set_fmpz_poly(gq, g);
  const mpfr_prec_t prec = labs(_fmpz_vec_max_bits(f->coeffs, fmpz_poly_length(f)));
  fmpq_poly_t g_inv; fmpq_poly_init(g_inv);
  fmpq_poly_oz_invert_approx(g_inv, gq, n, prec, 0);
  fmpq_poly_clear(gq);

  _fmpz_poly_oz_rem_small(rem, f, g, n, g_inv);

  fmpq_poly_clear(g_inv);
}

int fmpz_poly_oz_coprime(const fmpz_poly_t b0, const fmpz_poly_t b1, const long n,
                         const int sloppy, const mp_limb_t *primes) {

  fmpz_poly_t mod; fmpz_poly_init_oz_modulus(mod, n);
  int num_threads = omp_get_max_threads();

  /* If one operand is much larger than the other consider it mod the other */
  const mp_bitcnt_t s0 = labs(fmpz_poly_max_bits(b0));
  const mp_bitcnt_t s1 = labs(fmpz_poly_max_bits(b1));
  fmpz_poly_t v0; fmpz_poly_init(v0);
  fmpz_poly_t v1; fmpz_poly_init(v1);

  if (s0 > 1.1*s1) {
    fmpz_poly_oz_rem_small(v0, b0, b1, n);
    fmpz_poly_set(v1, b1);
  } else if (s1 > 1.1*s0) {
    fmpz_poly_set(v0, b0);
    fmpz_poly_oz_rem_small(v1, b1, b0, n);
  } else {
    fmpz_poly_set(v0, b0);
    fmpz_poly_set(v1, b1);
  }

  const size_t k = primes[0];

  nmod_poly_t n0[num_threads];
  nmod_poly_t n1[num_threads];
  nmod_poly_t nm[num_threads];

  int r0[num_threads];
  int r1[num_threads];

  for(int j=0; j<num_threads; j++) {
    r0[j] = 1;
    r1[j] = 1;
  }
  for(size_t i=0; i<k; i+=num_threads) {
    if (k-i < (unsigned long)num_threads)
      num_threads = k-i;
#pragma omp parallel for
    for (int j=0; j<num_threads; j++) {
      mp_limb_t p = primes[1+i+j];

      nmod_poly_init(n0[j], p); fmpz_poly_get_nmod_poly(n0[j], v0);
      nmod_poly_init(n1[j], p); fmpz_poly_get_nmod_poly(n1[j], v1);
      if(p%(2*n) == 1) {
        r0[j] = nmod_poly_oz_resultant(n0[j], n);
        r1[j] = nmod_poly_oz_resultant(n1[j], n);
      } else {
        nmod_poly_init(nm[j], p); fmpz_poly_get_nmod_poly(nm[j], mod);
        r0[j] = nmod_poly_resultant(n0[j], nm[j]);
        r1[j] = nmod_poly_resultant(n1[j], nm[j]);
        nmod_poly_clear(nm[j]);
      }
      nmod_poly_clear(n0[j]);
      nmod_poly_clear(n1[j]);
      flint_cleanup();
    }
    for(int j=0; j<num_threads; j++) {
      /* if both resultants are zero they share the prime factor p */
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

  /* run expensive test if we're not sloppy and we haven't ruled out co-primality yet */
  if (!sloppy && r0[0] == 1) {
    fmpz_t det_v0, det_v1;
    fmpz_init(det_v0);
    fmpz_init(det_v1);

    fmpz_poly_oz_ideal_norm(det_v0, v0, n, 0);
    fmpz_poly_oz_ideal_norm(det_v1, v1, n, 0);

    fmpz_t tmp;
    fmpz_init(tmp);
    fmpz_gcd(tmp, det_v0, det_v1);

    r0[0] = fmpz_equal_si(tmp, 1);
    fmpz_clear(tmp);

    fmpz_clear(det_v0);
    fmpz_clear(det_v1);
  }
  fmpz_poly_clear(v0);
  fmpz_poly_clear(v1);

  return r0[0];
}

int fmpz_poly_oz_coprime_det(const fmpz_poly_t b0, const fmpz_t det_b1, const long n,
                             const int sloppy, const mp_limb_t *primes) {

  fmpz_poly_t mod; fmpz_poly_init_oz_modulus(mod, n);
  int num_threads = omp_get_max_threads();

  const size_t k = primes[0];

  nmod_poly_t n0[num_threads];
  nmod_poly_t nm[num_threads];

  int r0[num_threads];

  for(int j=0; j<num_threads; j++) {
    r0[j] = 1;
  }
  for(size_t i=0; i<k; i+=num_threads) {
    if (k-i < (unsigned long)num_threads)
      num_threads = k-i;
#pragma omp parallel for
    for (int j=0; j<num_threads; j++) {
      mp_limb_t p = primes[1+i+j];
      nmod_poly_init(n0[j], p); fmpz_poly_get_nmod_poly(n0[j], b0);
      nmod_poly_init(nm[j], p); fmpz_poly_get_nmod_poly(nm[j], mod);
      r0[j] = nmod_poly_resultant(n0[j], nm[j]);
      nmod_poly_clear(n0[j]);
      nmod_poly_clear(nm[j]);
      flint_cleanup();
    }
    for(int j=0; j<num_threads; j++) {
      /* if both resultants are zero they share the prime factor p */
      if (r0[j] == 0) {
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

  fmpz_t det_b0;
  fmpz_init(det_b0);

  fmpz_poly_oz_ideal_norm(det_b0, b0, n, 0);

  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_gcd(tmp, det_b0, det_b1);
  fmpz_clear(det_b0);
  r0[0] = fmpz_equal_si(tmp, 1);
  fmpz_clear(tmp);
  return r0[0];
}
