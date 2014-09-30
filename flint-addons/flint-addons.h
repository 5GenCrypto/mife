#ifndef FLINT_ADDONS__H
#define FLINT_ADDONS__H

#include <assert.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>

#include <gghlite/misc.h>
#include <math.h>

static inline mp_limb_t n_prevprime(mp_limb_t n, int proved) {
  n--;
  if (n <= 1)
    abort();
  if (n <= 3)
    return n;
  if (n%2 == 0)
    n -= 1;
  while ((!proved && !n_is_probabprime(n)) || (proved && !n_is_prime(n)))
    n -= 2;
  return n;
}

/**
   fmpz_mat_t
*/

static inline int fmpz_mat_rot_is_one(const fmpz_mat_t mat) {
  assert(mat->c);
  if (mat->r != 1)
    return 0;
  if(!fmpz_is_one(&mat->rows[0][0]))
    return 0;
  for(long i=1; i<mat->c; i++)
    if(!fmpz_is_zero(&mat->rows[0][i]))
      return 0;
  return 1;
}

static inline int fmpz_mat_is_one(const fmpz_mat_t mat) {
  if (!fmpz_mat_is_square(mat))
    return 0;
  for(long i=0; i<mat->r; i++) {
    for(long j=0; j<i; j++)
      if(!fmpz_is_zero(&mat->rows[i][j]))
        return 0;
    if(!fmpz_is_one(&mat->rows[i][i]))
      return 0;
    for(long j=i+1; j<mat->c; j++)
      if(!fmpz_is_zero(&mat->rows[i][j]))
        return 0;
  }
  return 1;
}

/**
   fmpz*
*/

// TODO: This function should accept a rounding mode.
static inline void _fmpz_vec_2norm_mpfr(mpfr_t rop, fmpz *vec, long len) {
  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_zero(tmp);
  for (long i = 0; i < len; i++)
    fmpz_addmul(tmp, vec + i, vec + i);
  mpz_t tmp_g;
  mpz_init(tmp_g);
  fmpz_get_mpz(tmp_g, tmp);
  fmpz_clear(tmp);

  mpfr_set_z(rop, tmp_g, MPFR_RNDN);
  mpfr_sqrt(rop, rop, MPFR_RNDN);
  mpz_clear(tmp_g);
}

static inline void _fmpz_vec_mul(fmpz_t rop, fmpz *a, fmpz *b, long len) {
  fmpz_zero(rop);
  for (long i = 0; i < len; i++)
    fmpz_addmul(rop, a + i, b + i);
}

static inline void _fmpz_vec_rot_left(fmpz *rop, fmpz *op, long len, long shift) {
  fmpz* tmp = _fmpz_vec_init(shift);
  for(long i=0; i<shift; i++) {
    fmpz_set(tmp + i, op + len - shift + i);
  }
  for(long i=len-shift-1; i>=0; i--)
    fmpz_set(rop + shift + i, op + i);
  for(long i=0; i<shift; i++) {
    fmpz_set(rop + i, tmp + i);
  }
  _fmpz_vec_clear(tmp, shift);
}

static inline void _fmpz_vec_rot_left_neg(fmpz *rop, fmpz *op, long len) {
  _fmpz_vec_rot_left(rop, op, len, 1);
  fmpz_neg(rop + 0, rop + 0);
}

/*
   fmpq*
*/

// TODO: This function should accept a rounding mode.
static inline void _fmpq_vec_2norm_mpfr(mpfr_t rop, fmpz *num, fmpz_t den, long len) {
  fmpz_t acc_num;
  fmpz_init(acc_num);
  fmpz_zero(acc_num);
  for (long i = 0; i < len; i++)
    fmpz_addmul(acc_num, num + i, num + i);

  fmpz_t acc_den;
  fmpz_init(acc_den);
  fmpz_mul(acc_den, den, den);
  fmpz_mul_ui(acc_den, acc_den, len);

  fmpq_t acc;
  fmpq_init(acc);

  fmpq_set_fmpz_frac(acc, acc_num, acc_den);
  fmpq_get_mpfr(rop, acc, MPFR_RNDN);
  mpfr_sqrt(rop, rop, MPFR_RNDN);

  fmpq_clear(acc);
  fmpz_clear(acc_den);
  fmpz_clear(acc_num);
}

// TODO: This function should accept a rounding mode.
static inline void _fmpz_vec_set_mpfr_vec(fmpz *rop, mpfr_t *op, const long len) {
  mpz_t tmp;
  mpz_init(tmp);
  for(long i=0; i<len; i++) {
    mpfr_get_z(tmp, op[i], MPFR_RNDN);
    fmpz_set_mpz(rop + i, tmp);
  }
  mpz_clear(tmp);
}

static inline void fmpq_poly_truncate_prec(fmpq_poly_t op, mp_bitcnt_t prec) {
  mpq_t tmp_q; mpq_init(tmp_q);
  mpf_t tmp_f; mpf_init2(tmp_f, prec);

  for (int i=0; i<fmpq_poly_length(op); i ++) {
    fmpq_poly_get_coeff_mpq(tmp_q, op, i);
    mpf_set_q(tmp_f, tmp_q);
    mpq_set_f(tmp_q, tmp_f);
    fmpq_poly_set_coeff_mpq(op, i, tmp_q);
  }

  mpf_clear(tmp_f);
  mpq_clear(tmp_q);
}

static inline void fmpq_poly_mulmod(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2, fmpq_poly_t modulus) {
  fmpq_poly_mul(rop, op1, op2);
  fmpq_poly_rem(rop, rop, modulus);
}

static inline void _fmpq_poly_invert_mod_cnf2pow_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, const int n, mpfr_prec_t prec) {
  int i;

  fmpq_poly_t V;
  fmpq_poly_init(V);
  fmpq_poly_set(V, f);

  fmpq_poly_t U;
  fmpq_poly_init(U);

  /* F(x) = x^n+1 */
  fmpq_poly_t F;
  fmpq_poly_init(F);
  fmpq_poly_set_coeff_si(F, 0, 1);
  fmpq_poly_set_coeff_si(F, n, 1);

  fmpq_poly_t V2;
  fmpq_poly_init(V2);

  fmpq_poly_t Se;
  fmpq_poly_init(Se);
  fmpq_poly_t So;
  fmpq_poly_init(So);

  fmpq_t tmp;
  fmpq_init(tmp);

  int deg;

  if (n > 1) {
    /* set V2(x) = V(-x) */
    fmpq_poly_set(V2, V);
    deg = fmpq_poly_degree(V2);
    for (i = 1 ; i <= deg ; i += 2) {
      /* negate odd coefficients */
      fmpq_poly_get_coeff_fmpq(tmp, V2, i);
      fmpq_neg(tmp,tmp);
      fmpq_poly_set_coeff_fmpq(V2, i, tmp);
    }

    /* compute Se, So */
    /* Set Se = (V(x) + V(-x))/2 */
    for (i = 0; i <= deg/2 ; i++) {
      fmpq_poly_get_coeff_fmpq(tmp, V2, 2*i);
      fmpq_poly_set_coeff_fmpq(Se, i, tmp);
    }
    fmpq_poly_truncate(Se,deg/2+1);

    /* Set So to the compressed (V(-x) - V(x))/2 */
    for (i = 0 ; i < (deg+1)/2 ; i++) {
      fmpq_poly_get_coeff_fmpq(tmp,V2,2*i+1);
      fmpq_poly_set_coeff_fmpq(So,i,tmp);
    }
    fmpq_poly_truncate(So,deg/2+1);

    /* V = V(x) * V(-x) mod f(x) */
    fmpq_poly_mulmod(V, V, V2, F);

#ifndef NDEBUG
    /* Sanity-check: verify that the odd coeffs in V are zero */
    for (i = 1; i < n ; i += 2) {
      fmpq_poly_get_coeff_fmpq(tmp,V,i);
      assert(fmpq_is_zero(tmp));
    }
#endif

    deg = fmpq_poly_degree(V);
    /* Compress the non-zero coeffs of V */
    for (i = 0; i <= deg/2 ; i ++) {
      fmpq_poly_get_coeff_fmpq(tmp,V,2*i);
      fmpq_poly_set_coeff_fmpq(V,i,tmp);
    }
    fmpq_poly_truncate(V,deg/2+1);

    /* V2 = Inverse3(V,q,n/2) */
    _fmpq_poly_invert_mod_cnf2pow_approx(V2,V,n/2, prec);

    deg = fmpq_poly_degree(V2);
    /* Te=G*Se, To = G*So */
    fmpq_poly_mul(V,V2,Se);
    fmpq_poly_mul(U,V2,So);

    deg = fmpq_poly_degree(V);

    for (int i = 0 ; i <= deg ; i ++) {
      fmpq_poly_get_coeff_fmpq(tmp,V,i);
      fmpq_poly_set_coeff_fmpq(f_inv,2*i,tmp);
    }
    deg = fmpq_poly_degree(U);
    for (int i = 0 ; i <= deg ; i ++) {
      fmpq_poly_get_coeff_fmpq(tmp,U,i);
      fmpq_poly_set_coeff_fmpq(f_inv,2*i+1,tmp);
    }
    fmpq_poly_truncate(f_inv,2*deg+2);
    fmpq_poly_rem(f_inv,f_inv,F);

    /* truncate results on the precision required by algorithm */
    fmpq_poly_truncate_prec(f_inv, prec);
  } else {
    f_inv->length = 1;
    fmpz_set(f_inv->coeffs, f->den);
    fmpz_set(f_inv->den, f->coeffs);
  }

  /* free memory */
  fmpq_clear(tmp);
  fmpq_poly_clear(Se);
  fmpq_poly_clear(So);
  fmpq_poly_clear(V2);
  fmpq_poly_clear(F);
  fmpq_poly_clear(U);
  fmpq_poly_clear(V);
}

static inline void fmpq_poly_invert_mod_cnf2pow_approx(fmpq_poly_t f_inv, const fmpq_poly_t f, int n, mp_bitcnt_t prec) {
  fmpq_poly_t tmp;
  fmpq_poly_init(tmp);
  fmpq_poly_t modulus;
  fmpq_poly_init(modulus);
  fmpq_poly_set_coeff_si(modulus, 0, 1);
  fmpq_poly_set_coeff_si(modulus, n, 1);

  mpfr_t norm;
  mpfr_init2(norm, 2*prec);

  fmpq_t c;
  fmpq_init(c);

  /** we check for distance < 2^-prec with precision 2*prec **/
  mpfr_t bound;
  mpfr_init2(bound, 2*prec);
  mpfr_set_ui(bound, 1, MPFR_RNDN);
  mpfr_div_2exp(bound, bound, prec, MPFR_RNDN);

  for(long b=ceil(log2(prec)); ; b++) {
    _fmpq_poly_invert_mod_cnf2pow_approx(f_inv, f, n, (1<<b));
    fmpq_poly_mulmod(tmp, f_inv, f, modulus);

    fmpq_poly_get_coeff_fmpq(c, tmp, 0);
    fmpq_sub_si(c, c, 1);
    fmpq_poly_set_coeff_fmpq(tmp, 0, c);

    _fmpq_vec_2norm_mpfr(norm, tmp->coeffs, tmp->den, fmpq_poly_length(tmp));
    if(mpfr_cmp(norm, bound) <= 0)
      break;
  }
  fmpq_poly_clear(tmp);
  fmpq_poly_clear(modulus);
  mpfr_clear(norm);
  mpfr_clear(bound);
  fmpq_clear(c);
}


static inline void fmpq_polq_invert_mod(fmpq_poly_t f_inv, fmpq_poly_t f, const fmpz_poly_t g) {
  fmpq_poly_t r; fmpq_poly_init(r);
  fmpq_poly_t t; fmpq_poly_init(t);

  fmpq_poly_t g_; fmpq_poly_init(g_); fmpq_poly_set_fmpz_poly(g_, g);

  fmpq_poly_xgcd(r, f_inv, t, f, g_);
  assert(fmpq_poly_is_one(r));

  fmpq_poly_clear(g_);
  fmpq_poly_clear(t);
  fmpq_poly_clear(r);
}


/**
   fmpz_poly_t
*/

static inline void fmpz_poly_invert_mod_fmpq(fmpq_poly_t f_inv, fmpz_poly_t f, const fmpz_poly_t g) {
  fmpq_poly_t r; fmpq_poly_init(r);
  fmpq_poly_t t; fmpq_poly_init(t);

  fmpq_poly_t f_; fmpq_poly_init(f_); fmpq_poly_set_fmpz_poly(f_, f);
  fmpq_poly_t g_; fmpq_poly_init(g_); fmpq_poly_set_fmpz_poly(g_, g);

  fmpq_poly_xgcd(r, f_inv, t, f_, g_);
  assert(fmpq_poly_is_one(r));

  fmpq_poly_clear(g_);
  fmpq_poly_clear(f_);
  fmpq_poly_clear(t);
  fmpq_poly_clear(r);
}

static inline void fmpz_poly_ideal_rot_basis(fmpz_mat_t rop, fmpz_poly_t op) {
  assert(fmpz_mat_nrows(rop) == fmpz_poly_length(op));
  assert(fmpz_mat_ncols(rop) == fmpz_poly_length(op));

  printf("WARNING: called fmpz_poly_ideal_rot_basis, O(n^ω) algorithms ahead.\n");

  const long n = fmpz_poly_length(op);

  fmpz* v = _fmpz_vec_init(n);
  _fmpz_vec_set(v, op->coeffs, n);
  for(long i=0; i<n; i++) {
    for (long j=0; j<n; j++)
      fmpz_set(fmpz_mat_entry(rop, i, j), v+j);
    _fmpz_vec_rot_left_neg(v, v, n);
  }
  _fmpz_vec_clear(v, n);
}

static inline void _fmpz_poly_resultant_modular_bound(fmpz_t res, const fmpz * poly1, slong len1,
                                                      const fmpz * poly2, slong len2, mp_bitcnt_t bound) {
  mp_bitcnt_t pbits, curr_bits = 0;
  slong i, num_primes;
  fmpz_comb_t comb;
  fmpz_comb_temp_t comb_temp;
  fmpz_t ac, bc, l, modulus;
  fmpz * A, * B, * lead_A, * lead_B;
  mp_ptr a, b, rarr, parr;
  mp_limb_t p;
  nmod_t mod;

  /* special case, one of the polys is a constant */
  if (len2 == 1) /* if len1 == 1 then so does len2 */ {
    fmpz_pow_ui(res, poly2, len1 - 1);
    return;
  }

  fmpz_init(ac);
  fmpz_init(bc);

  /* compute content of poly1 and poly2 */
  _fmpz_vec_content(ac, poly1, len1);
  _fmpz_vec_content(bc, poly2, len2);

  /* divide poly1 and poly2 by their content */
  A = _fmpz_vec_init(len1);
  B = _fmpz_vec_init(len2);
  _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, ac);
  _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, bc);

  /* get product of leading coefficients */
  fmpz_init(l);

  lead_A = A + len1 - 1;
  lead_B = B + len2 - 1;
  fmpz_mul(l, lead_A, lead_B);

  /* set size of first prime */
  pbits = FLINT_BITS -1;
  p = (UWORD(1)<<pbits);

  num_primes = (bound + pbits - 1)/pbits;
  parr = _nmod_vec_init(num_primes);
  rarr = _nmod_vec_init(num_primes);

  fmpz_init(modulus);
  fmpz_set_ui(modulus, 1);
  fmpz_zero(res);

  /* make space for polynomials mod p */
  a = _nmod_vec_init(len1);
  b = _nmod_vec_init(len2);

  for (i = 0; curr_bits < bound; ) {
    /* get new prime and initialise modulus */
    p = n_prevprime(p, 0);
    if (fmpz_fdiv_ui(l, p) == 0)
      continue;

    curr_bits += pbits;

    nmod_init(&mod, p);

    /* reduce polynomials modulo p */
    _fmpz_vec_get_nmod_vec(a, A, len1, mod);
    _fmpz_vec_get_nmod_vec(b, B, len2, mod);

    /* compute resultant over Z/pZ */
    parr[i] = p;
    rarr[i++] = _nmod_poly_resultant(a, len1, b, len2, mod);
  }

  fmpz_comb_init(comb, parr, num_primes);
  fmpz_comb_temp_init(comb_temp, comb);

  fmpz_multi_CRT_ui(res, rarr, comb, comb_temp, 1);

  fmpz_clear(modulus);
  fmpz_comb_temp_clear(comb_temp);
  fmpz_comb_clear(comb);

  _nmod_vec_clear(a);
  _nmod_vec_clear(b);

  _nmod_vec_clear(parr);
  _nmod_vec_clear(rarr);

  /* finally multiply by powers of content */
  if (!fmpz_is_one(ac)) {
    fmpz_pow_ui(l, ac, len2 - 1);
    fmpz_mul(res, res, l);
  }

  if (!fmpz_is_one(bc)) {
    fmpz_pow_ui(l, bc, len1 - 1);
    fmpz_mul(res, res, l);
  }

  fmpz_clear(l);

  _fmpz_vec_clear(A, len1);
  _fmpz_vec_clear(B, len2);

  fmpz_clear(ac);
  fmpz_clear(bc);
}

static inline void fmpz_poly_resultant_modular_bound(fmpz_t res, const fmpz_poly_t poly1,
                                                     const fmpz_poly_t poly2, const mp_bitcnt_t bound) {
  slong len1 = poly1->length;
  slong len2 = poly2->length;

  if (len1 == 0 || len2 == 0)
    fmpz_zero(res);
  else if (len1 >= len2)
    _fmpz_poly_resultant_modular_bound(res, poly1->coeffs, len1, poly2->coeffs, len2, bound);
  else {
    _fmpz_poly_resultant_modular_bound(res, poly2->coeffs, len2, poly1->coeffs, len1, bound);
    if ((len1 > 1) && (!(len1 & WORD(1)) & !(len2 & WORD(1))))
      fmpz_neg(res, res);
  }
}

static inline int fmpz_poly_ideal_norm(fmpz_t norm, fmpz_poly_t f, fmpz_poly_t modulus) {
  mp_bitcnt_t bits = FLINT_ABS(_fmpz_vec_max_bits(f->coeffs, f->length));
  mp_bitcnt_t bound = f->length * (bits + n_clog(f->length, 2));
  fmpz_poly_resultant_modular_bound(norm, f, modulus, bound);
}

/**
   Decide if <b_0,b_1> = <g>
 */

static inline int fmpz_poly_ideal_subset(fmpz_poly_t g, fmpz_poly_t b0, fmpz_poly_t b1, fmpz_poly_t modulus) {
  const long n = fmpz_poly_degree(modulus);
  fmpz_t det;
  fmpz_init(det);
  fmpz_poly_ideal_norm(det, g, modulus);

  fmpz_t det_b0, det_b1;
  fmpz_init(det_b0);
  fmpz_init(det_b1);

  fmpz_poly_ideal_norm(det_b0, b0, modulus);
  fmpz_poly_ideal_norm(det_b1, b1, modulus);

  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_gcd(tmp, det_b0, det_b1);
  fmpz_clear(det_b0);
  fmpz_clear(det_b1);

  int r = fmpz_cmp(det, tmp);

  fmpz_clear(det);
  fmpz_clear(tmp);
  return r;
}

/**


*/

static inline int fmpz_poly_ideal_is_probaprime(fmpz_poly_t op, fmpz_poly_t modulus) {
  fmpz_t det;
  fmpz_init(det);

  fmpz_poly_ideal_norm(det, op, modulus);

  int r = fmpz_is_probabprime(det);
  fmpz_clear(det);
  return r;
}

static inline void fmpz_poly_set_fmpz_mod_poly(fmpz_poly_t rop, fmpz_mod_poly_t op) {
  long n = fmpz_mod_poly_length(op);

  fmpz *q = fmpz_mod_poly_modulus(op);

  fmpz_t q2;
  fmpz_init(q2);
  fmpz_fdiv_q_2exp(q2,q, 1);

  fmpz_t tmp;
  fmpz_init(tmp);

  for(long i=0; i<n; i++) {
    fmpz_mod_poly_get_coeff_fmpz(tmp, op, i);
    if (fmpz_cmp(tmp, q2)>=0)
      fmpz_sub(tmp, tmp, q);
    fmpz_poly_set_coeff_fmpz(rop, i, tmp);
  }
  fmpz_clear(tmp);
  fmpz_clear(q2);
}

/**
   fmpz_mod_poly_t
 */

static inline void fmpz_mod_poly_invert_mod(fmpz_mod_poly_t f_inv, fmpz_mod_poly_t f, fmpz_mod_poly_t modulus) {
  fmpz_mod_poly_t r; fmpz_mod_poly_init(r, fmpz_mod_poly_modulus(modulus));
  fmpz_mod_poly_t t; fmpz_mod_poly_init(t, fmpz_mod_poly_modulus(modulus));

  fmpz_mod_poly_xgcd(r, f_inv, t, f, modulus);
  assert(fmpz_mod_poly_degree(r) == 0 && fmpz_is_one(r->coeffs + 0));

  fmpz_mod_poly_clear(t);
  fmpz_mod_poly_clear(r);
}

static inline void fmpz_mod_poly_set_fmpq_poly(fmpz_mod_poly_t rop, fmpq_poly_t op) {
  fmpz_t tmp;
  fmpz_t den_inv;

  fmpz *q = fmpz_mod_poly_modulus(rop);

  fmpz_init(tmp);
  fmpz_init(den_inv);

  fmpz_set(den_inv, fmpq_poly_denref(op));
  fmpz_invmod(den_inv, den_inv, q);

  for(long i=0; i<fmpq_poly_length(op); i++) {
    fmpz_set(tmp, fmpq_poly_numref(op)+i);
    fmpz_mul(tmp, tmp, den_inv);
    fmpz_mod(tmp, tmp, q);
    fmpz_mod_poly_set_coeff_fmpz(rop, i, tmp);
  }
  fmpz_clear(tmp);
  fmpz_clear(den_inv);
}

/**
   flint_rand_t
*/

static inline void flint_randinit_seed(flint_rand_t randstate, mp_limb_t seed, int gmp) {
  flint_randinit(randstate);
  randstate->__randval  = seed;
  randstate->__randval2 = seed ^ 0x5555555555555555ULL;
  if (gmp) {
    _flint_rand_init_gmp(randstate);
    gmp_randseed_ui(randstate->gmp_state, seed ^ 0xDEADBEEFDEADBEEFULL);
  }
}

#endif //FLINT_ADDONS__H
