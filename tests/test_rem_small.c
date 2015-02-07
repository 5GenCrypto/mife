#include <dgsl/dgsl.h>
#include <oz/oz.h>
#include <oz/util.h>
#include <mpfr.h>
#include <math.h>

int test_fmpz_poly_oz_rem_small(const long n, const mp_bitcnt_t bits, flint_rand_t state) {

  /** make it **/
  mpfr_t sigma;
  mpfr_init(sigma);

  fmpz_poly_t g; fmpz_poly_init(g);
  mpfr_set_d(sigma, (double)n, MPFR_RNDN);
  fmpz_poly_sample_sigma(g, n, sigma, state);

  fmpz_poly_t h; fmpz_poly_init(h);
  mpfr_set_si_2exp(sigma, 1, bits, MPFR_RNDN);
  fmpz_poly_sample_sigma(h, n, sigma, state);

  fmpz_poly_t small; fmpz_poly_init(small);
  fmpz_poly_oz_rem_small(small, h, g, n);

  /** check it **/
  fmpz_poly_t t; fmpz_poly_init(t);
  fmpz_poly_sub(t, h, small);

  fmpq_poly_t tq; fmpq_poly_init(tq);
  fmpq_poly_set_fmpz_poly(tq, t);

  fmpq_poly_t gq; fmpq_poly_init(gq);
  fmpq_poly_set_fmpz_poly(gq, g);

  fmpq_poly_t ginv; fmpq_poly_init(ginv);
  fmpq_poly_oz_invert_approx(ginv, gq, n, 0, 0);

  fmpq_poly_oz_mul(tq, tq, ginv, n);

  /* fmpq_poly_print_pretty(tq, "x"); */
  /* printf("\n"); */

  int r = (fmpz_is_one(tq->den)) ? 0 : 1;

  printf("n: %4ld, bits: %4ld:", n, bits);
  if (r == 0)
    printf("PASS\n");
  else
    printf("FAIL\n");

  fmpq_poly_clear(ginv);
  fmpq_poly_clear(gq);
  fmpq_poly_clear(tq);
  fmpz_poly_clear(t);
  fmpz_poly_clear(small);
  fmpz_poly_clear(h);
  fmpz_poly_clear(g);
  mpfr_clear(sigma);
  return r;
}


int main(int argc, char *argv[]) {

  flint_rand_t state;
  flint_randinit_seed(state, 0x1337, 1);

  int status = 0;

  long n[5] = {32,64,128,256,0};

  for(int i=0; n[i]; i++)
    for(mp_bitcnt_t bits=2; bits<=2*n[i]; bits=2*bits)
      status += test_fmpz_poly_oz_rem_small(n[i], bits, state);

  flint_randclear(state);
  flint_cleanup();
  return status;
}
