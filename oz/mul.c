#include "mul.h"
#include "ntt.h"
#include "oz.h"
#include "util.h"
#include "flint-addons.h"

void fmpz_poly_oz_rem(fmpz_poly_t rem, const fmpz_poly_t f, const long n) {
  fmpz_poly_set(rem, f);
  fmpz_t lead; fmpz_init(lead);
  fmpz_t coef; fmpz_init(coef);
  for(long d=fmpz_poly_degree(rem); d>=n; d--) {
    fmpz_poly_get_coeff_fmpz(lead, rem, d);
    fmpz_poly_get_coeff_fmpz(coef, rem, d-n);
    fmpz_sub(coef, coef, lead);
    fmpz_poly_set_coeff_si(rem, d, 0);
    fmpz_poly_set_coeff_fmpz(rem, d-n, coef);
  }
  fmpz_clear(coef);
  fmpz_clear(lead);
}

void fmpz_mod_poly_oz_rem(fmpz_mod_poly_t rem, const fmpz_mod_poly_t f, const long n) {
  fmpz_mod_poly_set(rem, f);
  fmpz_t lead; fmpz_init(lead);
  fmpz_t coef; fmpz_init(coef);
  for(long d=fmpz_mod_poly_degree(rem); d>=n; d--) {
    fmpz_mod_poly_get_coeff_fmpz(lead, rem, d);
    fmpz_mod_poly_get_coeff_fmpz(coef, rem, d-n);
    fmpz_sub(coef, coef, lead);
    fmpz_mod_poly_set_coeff_ui(rem, d, 0);
    fmpz_mod_poly_set_coeff_fmpz(rem, d-n, coef);
  }
  fmpz_clear(coef);
  fmpz_clear(lead);
}

void fmpq_poly_oz_rem(fmpq_poly_t rem, const fmpq_poly_t f, const long n) {
  const long d = fmpq_poly_degree(f);
  if(d > 2*n-1) {
    // TODO: Don't be so lazy and take care of this case as well
    fmpq_poly_t modulus;
    fmpq_poly_init_oz_modulus(modulus, n);
    fmpq_poly_rem(rem, f, modulus);
    fmpq_poly_clear(modulus);
  } else if (d>=n) {
    assert(d-n+1 <= n);
    fmpq_poly_t r;
    fmpq_poly_init(r);
    fmpq_poly_realloc(r, n);
    mpq_t *z = (mpq_t*)calloc(n, sizeof(mpq_t));
    for (int i = 0; i < n; i++) {
      mpq_init(z[i]);
      fmpq_poly_get_coeff_mpq(z[i], f, i);
    }
    fmpq_poly_set_array_mpq(r, (const mpq_t*)z, n);

    for (int i=0; i<d-n+1; i++)
      fmpq_poly_get_coeff_mpq(z[i], f, n+i);
    fmpq_poly_t t;
    fmpq_poly_init(t);
    fmpq_poly_set_array_mpq(t, (const mpq_t*)z, d-n+1);
    fmpq_poly_sub(r, r, t);
    fmpq_poly_set(rem, r);

    fmpq_poly_clear(t);
    fmpq_poly_clear(r);
    for (int i=0; i<n; i++)
      mpq_clear(z[i]);
    free(z);
  }
  return;
}

