#include "gso.h"
#include <gmp.h>
#include <assert.h>

void mpfr_mat_init(mpfr_mat_t mat, long rows, long cols, mpfr_prec_t prec) {
  if ((rows) && (cols)) {
    mat->entries = (mpfr_t *) calloc(rows * cols, sizeof(mpfr_t));
    mat->rows = (mpfr_t **) flint_malloc(rows * sizeof(mpfr_t *));    /* Initialise rows */

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

void mpfr_mat_set_fmpz_mat_rot(mpfr_mat_t rop, const fmpz_mat_t op) {
  assert((rop->r == op->c) && (rop->c == op->c));
  const long n = op->c;  

  mpz_t t_g;
  mpz_init(t_g);

  fmpz *gen = _fmpz_vec_init(n);
  _fmpz_vec_set(gen, op->rows[0], n);
  
  for(long i=0; i<n; i++) {
    for(long j=0; j<n; j++) {
      fmpz_get_mpz(t_g, gen +j );
      mpfr_set_z(rop->rows[i][j], t_g, MPFR_RNDN);
    }
    _fmpz_vec_rot_left_neg(gen, gen, n);
  }
  _fmpz_vec_clear(gen, n);
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

void mpfr_mat_gso(mpfr_mat_t mat)  {
  const long m = mat->r;
  const long n = mat->c;
  const mpfr_prec_t prec = mpfr_mat_get_prec(mat);

  mpfr_mat_t dot_products;
  mpfr_mat_init(dot_products, 1, m, prec);

  mpfr_t tmp;
  mpfr_init2(tmp, prec);
  
  for(long k = 0; k < m; k++) {
    for(long i = 0; i < k; i++) {
      _mpfr_vec_dot_product(tmp, (const mpfr_t *)mat->rows[i], (const mpfr_t *)mat->rows[k], n);
      mpfr_div(tmp, tmp, dot_products->rows[0][i], MPFR_RNDN);
      _mpfr_vec_submul(mat->rows[k], tmp, mat->rows[i], n);      
    }
    _mpfr_vec_dot_product(dot_products->rows[0][k], (const mpfr_t *)mat->rows[k], (const mpfr_t *)mat->rows[k], n);
  }
  mpfr_clear(tmp);
  mpfr_mat_clear(dot_products);
}
