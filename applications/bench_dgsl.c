#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include <oz/util.h>
#include "common.h"

int main(int argc, char *argv[]) {

  if(argc != 4)
    ggh_die("parameters n, ntrials and prec needed");

  const long n    = atol(argv[1]);
  const long m    = atol(argv[2]);
  const long prec = atol(argv[3]);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, 0x1337, 1);

  mpfr_t sigma;
  mpfr_init2(sigma, prec);
  mpfr_set_d(sigma, _gghlite_sigma(n) * 0.398942280401433, MPFR_RNDN);

  mpfr_t sigma_p;
  mpfr_init2(sigma_p, prec);
  mpfr_set_d(sigma_p, _gghlite_sigma_p(n) * 0.398942280401433, MPFR_RNDN);

  fmpz_poly_t g;
  fmpz_poly_init2(g, n);
  fmpz_poly_sample_sigma(g, n, sigma, randstate);

  uint64_t t_sqrt = oz_walltime(0);
  dgsl_rot_mp_t *D = dgsl_rot_mp_init(n, g, sigma_p, NULL, DGSL_INLATTICE, OZ_VERBOSE);
  t_sqrt = oz_walltime(t_sqrt);

  fmpz_poly_t f;
  fmpz_poly_init(f);
  uint64_t t_sample = oz_walltime(0);
  for(int i=0; i<m; i++) {
    D->call(f, D, randstate->gmp_state);
  }
  t_sample = oz_walltime(t_sample);
  fmpz_poly_clear(f);

  printf("prec: %4ld, n: %4ld, log σ: %.1f, log σ': %.1f, sqrt time: %.2fs, sample time: %.2f ms, rate: %.2f\n", prec, n, log2(mpfr_get_d(sigma, MPFR_RNDN)),
         log2(mpfr_get_d(sigma_p, MPFR_RNDN)), t_sqrt/1000000.0, t_sample/1000.0/m, (1000000.0*m)/t_sample);

  dgsl_rot_mp_clear(D);
  fmpz_poly_clear(g);
  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
}
