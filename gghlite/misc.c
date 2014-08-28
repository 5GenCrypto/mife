#include "misc.h"

gpv_mp_t *_gghlite_gpv_from_poly(fmpz_poly_t g, mpfr_t sigma, mpfr_t *c, gpv_alg_t algorithm) {
  fmpz_mat_t B;
  const long n = fmpz_poly_length(g);
  fmpz_mat_init(B, 1, n);
  _fmpz_vec_set(B->rows[0], g->coeffs, n);

  /* GPV samples proportionally to `\exp(-(x-c)²/(2σ²))` but GGHLite is
     specifiied with respect to `\exp(-π(x-c)²/σ²)`. So we divide by \sqrt{2π}
  */

  mpfr_t sigma_;
  mpfr_init2(sigma_, mpfr_get_prec(sigma));
  mpfr_set(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, S_TO_SIGMA, MPFR_RNDN);
  
  gpv_mp_t *D = gpv_mp_init(B, sigma_, c, algorithm);
  fmpz_mat_clear(B);
  mpfr_clear(sigma_);
  return D;
}

gpv_mp_t *_gghlite_gpv_from_n(const long n, mpfr_t sigma) {
  fmpz_mat_t I;
  fmpz_mat_init(I, 1, n);
  fmpz_mat_one(I);

  /* GPV samples proportionally to `\exp(-(x-c)²/(2σ²))` but GGHLite is
     specifiied with respect to `\exp(-π(x-c)²/σ²)`. So we divide by \sqrt{2π}
  */

  mpfr_t sigma_;
  mpfr_init2(sigma_, mpfr_get_prec(sigma));
  mpfr_set(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, S_TO_SIGMA, MPFR_RNDN);
 
  gpv_mp_t *D = gpv_mp_init(I, sigma_, NULL, GPV_IDENTITY);
  fmpz_mat_clear(I);
  mpfr_clear(sigma_);
  return D;
}
