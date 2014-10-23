#include <omp.h>
#include "flint-addons.h"

void _fmpz_vec_eucl_norm_mpfr(mpfr_t rop, const fmpz *vec, const long len, const mpfr_rnd_t rnd) {
  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_zero(tmp);
  for (long i = 0; i < len; i++)
    fmpz_addmul(tmp, vec + i, vec + i);
  mpz_t tmp_g;
  mpz_init(tmp_g);
  fmpz_get_mpz(tmp_g, tmp);
  fmpz_clear(tmp);

  mpfr_set_z(rop, tmp_g, rnd);
  mpfr_sqrt(rop, rop, rnd);
  mpz_clear(tmp_g);
}

void _fmpq_vec_eucl_norm_mpfr(mpfr_t rop, const fmpz *num, const fmpz_t den, const long len, const mpfr_rnd_t rnd) {
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
  fmpq_get_mpfr(rop, acc, rnd);
  mpfr_sqrt(rop, rop, rnd);

  fmpq_clear(acc);
  fmpz_clear(acc_den);
  fmpz_clear(acc_num);
}

void fmpq_poly_truncate_prec(fmpq_poly_t op, const mp_bitcnt_t prec) {
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

void _fmpz_poly_resultant_modular_bound(fmpz_t res, const fmpz * poly1, const slong len1,
                                        const fmpz * poly2, const slong len2, const mp_bitcnt_t bound) {
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

