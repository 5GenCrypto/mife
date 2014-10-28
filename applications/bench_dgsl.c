#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include "common.h"

int main(int argc, char *argv[]) {

  if(argc != 3)
    ggh_die("parameters n and ntrials needed");

  const long n = atol(argv[1]);
  const long m = atol(argv[2]);
  
  flint_rand_t randstate;
  flint_randinit_seed(randstate, 0x1337, 1);

  mpfr_t sigma;
  mpfr_init2(sigma, 80);
  mpfr_set_d(sigma, _gghlite_sigma(n) * 0.398942280401433, MPFR_RNDN);

  mpfr_t sigma_p;
  mpfr_init2(sigma_p, 80);
  mpfr_set_d(sigma_p, _gghlite_sigma_p(n) * 0.398942280401433, MPFR_RNDN);

  fmpz_poly_t g;
  fmpz_poly_init2(g, n);
  fmpz_poly_sample_sigma(g, n, sigma, randstate);

  dgsl_rot_mp_t *D = dgsl_rot_mp_init(n, g, sigma_p, NULL, DGSL_INLATTICE, 0);

  fmpz_poly_t f;
  fmpz_poly_init(f);
  uint64_t t = ggh_walltime(0);
  for(int i=0; i<m; i++) {
    D->call(f, D, randstate->gmp_state);
  }
  t = ggh_walltime(t);
  fmpz_poly_clear(f);

  printf("n: %4ld, σ: %.1f, σ': %.1f, time: %.2f ms, rate: %.2f\n", n, mpfr_get_d(sigma, MPFR_RNDN),
         mpfr_get_d(sigma_p, MPFR_RNDN), t/1000.0/m, (1000000.0*m)/t);

  dgsl_rot_mp_clear(D);
  fmpz_poly_clear(g);
  flint_randclear(randstate);
  mpfr_free_cache();
}
