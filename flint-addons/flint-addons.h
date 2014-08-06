#ifndef FLINT_ADDONS__H
#define FLINT_ADDONS__H

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
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

static inline void fmpz_poly_sample(fmpz_poly_t f, dgs_disc_gauss_lattice_mp_t *D, flint_rand_t randstate) {
  assert(D);
  assert(randstate->gmp_init);

  const long n = fmpz_mat_ncols(D->B);
  fmpz_poly_realloc(f, n);
  D->call(f->coeffs, D, randstate->gmp_state);
  assert(!fmpz_is_zero(f->coeffs +(n-1)));
  _fmpz_poly_set_length(f, n);
}

#endif //FLINT_ADDONS__H
