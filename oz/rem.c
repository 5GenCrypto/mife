#include <omp.h>
#include "flint-addons.h"
#include "util.h"
#include "rem.h"
#include "invert.h"

void _fmpz_poly_oz_rem_small_fmpz(fmpz_poly_t rem, const fmpz_t f, const fmpz_poly_t g, const long n,
                                  const fmpq_poly_t g_inv, const mp_bitcnt_t bound) {
  if (fmpz_is_zero(f)) {
    fmpz_poly_set_ui(rem, 0);
    return;
  }
  fmpz_t fc; fmpz_init_set(fc, f);

  mp_bitcnt_t den_log_approx = fmpz_sizeinbase(g_inv->den, 2)-1;
  fmpz_t t; fmpz_init(t);
  for(int i=0; i<fmpq_poly_length(g_inv); i++) {
    fmpz_mul_tdiv_q_2exp(t, fc, g_inv->coeffs + i, den_log_approx);
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
      fprintf(stderr, "overflow |t|: %5ld >= |q|: %5ld\n", labs(fmpz_poly_max_bits(tmp)), (long)fmpz_sizeinbase(q, 2)-2);
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
  const mp_bitcnt_t rem_bound = log2(n)/2*labs(fmpz_poly_max_bits(g));

  fmpz_poly_t powb; fmpz_poly_init(powb);
  fmpz_poly_set_coeff_ui(powb, 0, 2); // powb ~= 2^b
  fmpz_pow_ui(powb->coeffs, powb->coeffs, b);
  _fmpz_poly_oz_rem_small_fmpz(powb, powb->coeffs, g, n, g_inv, rem_bound);

  fmpz_poly_t powB; fmpz_poly_init(powB);
  fmpz_poly_set_ui(powB, 1);
  for(size_t j=0; j<num_threads; j++)
    fmpz_poly_oz_mul(powB, powB, powb, n);
  /* invest a bit more here as it keeps everything below small */
  _fmpz_poly_oz_rem_small_iter(powB, powB, g, n, g_inv, 0, 0);

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
  assert(fmpz_is_zero(F));
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


void _fmpz_poly_oz_rem_small_iter(fmpz_poly_t rem, const fmpz_poly_t f, const fmpz_poly_t g,
                                  const long n, const fmpq_poly_t ginv, const mp_bitcnt_t b, const oz_flag_t flags) {

  mp_bitcnt_t prec = (b) ? b : labs(_fmpz_vec_max_bits(ginv->coeffs, fmpq_poly_length(ginv)))/2;
  fmpz_poly_t t_i;  fmpz_poly_init(t_i);
  fmpz_poly_t t_o;  fmpz_poly_init(t_o);
  mpfr_t norm_i; mpfr_init2(norm_i, prec);
  mpfr_t norm_o; mpfr_init2(norm_o, prec);

  fmpz_poly_set(t_i, f);
  fmpq_poly_t g_inv; fmpq_poly_init(g_inv);
  fmpq_poly_set(g_inv, ginv);

  if (fmpz_poly_degree(f) == 0) {
    uint64_t t = oz_walltime(0);
    _fmpz_poly_oz_rem_small_fmpz_split(t_o, f->coeffs, g, n, g_inv, prec);
    t = oz_walltime(t);

    if (flags & OZ_VERBOSE) {
      fprintf(stderr, "|f|: %10.1f, |g|: %10.1f, |f%%g|: %10.1f, t: %10.6f\n",
             fmpz_poly_2norm_log2(t_i), fmpz_poly_2norm_log2(g), fmpz_poly_2norm_log2(t_o),
             oz_seconds(t));
      fflush(stderr);
    }
  } else {
    fmpz_poly_set(t_o, t_i);
  }

  /* the precision of g_inv might not be sufficient to do this in one step, hence, we repeat until
     the result does not improve any more*/

  do {
    uint64_t t = oz_walltime(0);
    fmpz_poly_set(t_i, t_o);
    fmpz_poly_2norm_mpfr(norm_i, t_i, MPFR_RNDN);
    fmpq_poly_truncate_prec(g_inv, fmpz_poly_2norm_log2(t_i)/2);
    _fmpz_poly_oz_rem_small(t_o, t_i, g, n, g_inv);
    t = oz_walltime(t);
    fmpz_poly_2norm_mpfr(norm_o, t_o, MPFR_RNDN);

    if (flags & OZ_VERBOSE) {
      fprintf(stderr, "|f|: %10.1f, |g|: %10.1f, |f%%g|: %10.1f, t: %10.6f\n",
             fmpz_poly_2norm_log2(t_i), fmpz_poly_2norm_log2(g), fmpz_poly_2norm_log2(t_o),
             oz_seconds(t));
      fflush(stderr);
    }
  } while (mpfr_cmp(norm_o, norm_i) < 0);

  fmpz_poly_set(rem, t_i);
  mpfr_clear(norm_i);
  mpfr_clear(norm_o);
  fmpq_poly_clear(g_inv);
  fmpz_poly_clear(t_i);
  fmpz_poly_clear(t_o);
}
