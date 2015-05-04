#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include <oz/oz.h>

int main(int argc, char *argv[]) {
  assert(argc>=3);
  const long n = atol(argv[1]);
  const long prec = atol(argv[2]);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, 0x1337, 1);

  mpfr_t sigma;
  mpfr_init2(sigma, prec);
  mpfr_set_d(sigma, _gghlite_sigma(n), MPFR_RNDN);

  fmpz_poly_t g;
  fmpz_poly_init(g);
  fmpz_poly_sample_sigma(g, n, sigma, randstate);

  fmpq_poly_t gq; fmpq_poly_init(gq);
  fmpq_poly_set_fmpz_poly(gq, g);
  fmpq_poly_t g_inv; fmpq_poly_init(g_inv);

  fmpq_poly_t modulus;
  fmpq_poly_init_oz_modulus(modulus, n);

  printf(" n: %4ld, log σ': %7.2f, ",n, log2(_gghlite_sigma(n)));

  uint64_t t = ggh_walltime(0);
  fmpq_poly_invert_mod(g_inv, gq, modulus);
  t = ggh_walltime(t);
  printf("xgcd: %7.3f, ", ggh_seconds(t));
  fflush(0);

  t = ggh_walltime(0);
  _fmpq_poly_oz_invert_approx(g_inv, gq, n, 160);
  t = ggh_walltime(t);
  printf("%4ld: %7.3f, ", prec, ggh_seconds(t));
  fflush(0);

  t = ggh_walltime(0);
  fmpq_poly_oz_invert_approx(g_inv, gq, n, 160, 0);
  t = ggh_walltime(t);
  printf("%4lditer: %7.3f, ", prec, ggh_seconds(t));
  fflush(0);

  t = ggh_walltime(0);
  _fmpq_poly_oz_invert_approx(g_inv, gq, n, 0);
  t = ggh_walltime(t);
  printf("∞: %7.3f.", ggh_seconds(t));


  printf("\n");

  mpfr_clear(sigma);
  fmpz_poly_clear(g);
  fmpq_poly_clear(gq);
  fmpq_poly_clear(g_inv);
  fmpq_poly_clear(modulus);
  flint_randclear(randstate);
  flint_cleanup();
  return 0;
}
