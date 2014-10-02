#include "gso.h"
#include <gmp.h>
#include <assert.h>

void mpfr_mat_init(mpfr_mat_t mat, long rows, long cols, mpfr_prec_t prec) {
  if ((rows) && (cols)) {
    mat->entries = (mpfr_t *) calloc(rows * cols, sizeof(mpfr_t));
    if (!mat->entries)
      dgs_die("out of memory");
    mat->rows = (mpfr_t **) calloc(rows, sizeof(mpfr_t *));
    if (!mat->rows)
      dgs_die("out of memory");

    for (long i = 0; i < rows; i++) {
      mat->rows[i] = mat->entries + i * cols;
      for(long j=0; j< cols; j++)
        mpfr_init2(mat->rows[i][j], prec);
    }
  } else {
    mat->entries = NULL;
  }

  mat->r = rows;
  mat->c = cols;
}

void mpfr_mat_set_fmpz_poly(mpfr_mat_t rop, const fmpz_poly_t op) {
  assert(rop->r == rop->c);
  assert(fmpz_poly_length(op) <= rop->r);
  const long n = rop->r;

  mpz_t t_g;
  mpz_init(t_g);

  fmpz_poly_t tmp;
  fmpz_poly_init(tmp);
  fmpz_poly_set(tmp, op);
  fmpz_poly_realloc(tmp, rop->r);

  for(long i=0; i<n; i++) {
    for(long j=0; j<n; j++) {
      fmpz_get_mpz(t_g, tmp->coeffs + j);
      mpfr_set_z(rop->rows[i][j], t_g, MPFR_RNDN);
    }
    _fmpz_vec_rot_left_neg(tmp->coeffs, tmp->coeffs, n);
  }
  fmpz_poly_clear(tmp);
  mpz_clear(t_g);
}

void mpfr_mat_set_fmpz_mat(mpfr_mat_t rop, const fmpz_mat_t op) {
  assert((rop->r == op->r) && (rop->c == op->c));
  mpz_t t_g;
  mpz_init(t_g);
  for(long i=0; i<op->r; i++) {
    for(long j=0; j<op->c; j++) {
      fmpz_get_mpz(t_g, &op->rows[i][j]);
      mpfr_set_z(rop->rows[i][j], t_g, MPFR_RNDN);
    }
  }
  mpz_clear(t_g);
}

void mpfr_mat_clear(mpfr_mat_t mat) {
  if (mat->entries) {
    for(long i = 0; i < mat->r * mat->c; i++) {
      mpfr_clear(mat->entries[i]);   /* Clear all coefficients */
    }
    free(mat->entries);     /* Clean up array of entries */
    flint_free(mat->rows);        /* Clean up row array */
  }
}

mpfr_prec_t mpfr_mat_get_prec(mpfr_mat_t mat) {
  if(mat->entries)
    return mpfr_get_prec(mat->entries[0]);
  else
    return 53;
}

void mpfr_mat_gso(mpfr_mat_t mat, mpfr_rnd_t rnd)  {
  const long m = mat->r;
  const long n = mat->c;
  const mpfr_prec_t prec = mpfr_mat_get_prec(mat);

  mpfr_mat_t dot_products;
  mpfr_mat_init(dot_products, 1, m, prec);

  mpfr_t tmp;
  mpfr_init2(tmp, prec);

  for(long k = 0; k < m; k++) {
    for(long i = 0; i < k; i++) {
      _mpfr_vec_dot_product(tmp, mat->rows[i], mat->rows[k], n, rnd);
      mpfr_div(tmp, tmp, dot_products->rows[0][i], rnd);
      mpfr_neg(tmp, tmp, rnd);
      _mpfr_vec_scalar_addmul_mpfr(mat->rows[k], mat->rows[i], n, tmp, rnd);
    }
    _mpfr_vec_dot_product(dot_products->rows[0][k], mat->rows[k], mat->rows[k], n, rnd);
  }
  mpfr_clear(tmp);
  mpfr_mat_clear(dot_products);
}
