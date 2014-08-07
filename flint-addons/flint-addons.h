#ifndef FLINT_ADDONS__H
#define FLINT_ADDONS__H

#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <gpv/gpv.h>

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


static inline void fmpz_poly_sample_D(fmpz_poly_t f, dgs_disc_gauss_lattice_mp_t *D, flint_rand_t randstate) {
  assert(D);
  assert(randstate->gmp_init);

  const long n = fmpz_mat_ncols(D->B);
  fmpz_poly_realloc(f, n);
  D->call(f->coeffs, D, randstate->gmp_state);
  assert(!fmpz_is_zero(f->coeffs +(n-1)));
  _fmpz_poly_set_length(f, n);
}

static inline void fmpz_poly_sample_sigma(fmpz_poly_t f, long len, mpfr_t sigma, flint_rand_t randstate) {
  fmpz_mat_t I;
  fmpz_mat_init(I, 1, len);
  fmpz_mat_one(I);
  dgs_disc_gauss_lattice_mp_t *D = dgs_disc_gauss_lattice_mp_init(I, sigma, NULL, DGS_LATTICE_IDENTITY);

  fmpz_poly_sample_D(f, D, randstate);

  dgs_disc_gauss_lattice_mp_clear(D);
  fmpz_mat_clear(I);
}

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


#endif //FLINT_ADDONS__H
