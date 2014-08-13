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

static inline void _mpfr_vec_dot_product(mpfr_t rop, const mpfr_t *op1, const mpfr_t *op2, const long n) {
  mpfr_set_si(rop, 0, MPFR_RNDN);
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(rop));
  for(long i=0; i<n; i++) {
    mpfr_mul(tmp, op1[i], op2[i], MPFR_RNDN);
    mpfr_add(rop, rop, tmp, MPFR_RNDN); 
  }
  mpfr_clear(tmp);
}

static inline void _mpfr_vec_submul(mpfr_t *rop, mpfr_t scalar, mpfr_t *op, const long n) {
  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(rop[0]));
  for(long i=0; i<n; i++) {
    mpfr_mul(tmp, scalar, op[i], MPFR_RNDN);
    mpfr_sub(rop[i], rop[i], tmp, MPFR_RNDN);
  }
  mpfr_clear(tmp);
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
