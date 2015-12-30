#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include <oz/oz.h>

int main(int argc, char *argv[]) {
  assert(argc>=3);
  const long n = atol(argv[1]);
  const long m = atol(argv[2]);

  aes_randstate_t randstate;
  aes_randinit(randstate);

  fmpz_poly_t g;
  fmpz_poly_init(g);

  mpfr_t sigma;
  mpfr_init2(sigma, 80);
  mpfr_set_d(sigma, _gghlite_sigma(n), MPFR_RNDN);

  fmpz_poly_t modulus;
  fmpz_poly_init_oz_modulus(modulus, n);

  uint64_t t = 0;

  int k = ceil((log2(_gghlite_sigma(n)) + log2(n)/2.0) * n/100.0/(FLINT_BITS -1));
  if (k < 20)
    k = 20;

  mp_limb_t *primes = _fmpz_poly_oz_ideal_probable_prime_factors(n, k);

  int r = 0;
  for(long i=0; i<m; i++) {
    fmpz_poly_sample_sigma(g, n, sigma, randstate);
    uint64_t t0 = ggh_walltime(0);
    r += fmpz_poly_oz_ideal_is_probaprime(g, n, 0, primes);
    t += ggh_walltime(t0);
    printf("\rm: %6ld, #prime: %d,",i+1, r); fflush(0);
  }
  printf(" n: %4ld, log σ: %7.2f, t: %.6f\n",n, log2(_gghlite_sigma(n)), ggh_seconds(t)/m);

  free(primes);
  mpfr_clear(sigma);
  fmpz_poly_clear(g);
  aes_randclear(randstate);
  flint_cleanup();
  return 0;
}
