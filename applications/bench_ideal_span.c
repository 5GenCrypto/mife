#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include "common.h"

int main(int argc, char *argv[]) {
  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "Playing with primes", NULL);

  print_header("Playing with primes", params);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, params->seed, 1);

  gghlite_t self;
  gghlite_pk_init_params(self->pk, params->lambda, params->kappa, params->rerand, params->flags);
  gghlite_print_params(self->pk);

  printf("\n---\n");

  gghlite_init_instance(self, randstate);

  const int nsp = _gghlite_nsmall_primes(self->pk);
  mp_limb_t *small_primes = _fmpz_poly_oz_ideal_small_primes(self->pk->n, nsp);

  printf("number of small primes: %4d\n",nsp);

  fmpz_poly_t b0;
  fmpz_poly_init2(b0, self->pk->n);

  fmpz_poly_t b1;
  fmpz_poly_init2(b1, self->pk->n);

  uint64_t t0, t1, T0, T1;
  int r0, r1;

  T0 = 0;
  T1 = 0;
  int c = 0;
  for(int i=0; i<100; i++) {
    fmpz_poly_sample_D(b0, self->D_g, randstate);
    fmpz_poly_sample_D(b1, self->D_g, randstate);

    t0 = ggh_walltime(0);
    r0 = fmpz_poly_oz_ideal_span(self->g, b0, b1, self->pk->n, 1, nsp, small_primes);
    t0 = ggh_walltime(t0);

    t1 = ggh_walltime(0);
    r1 = fmpz_poly_oz_ideal_span(self->g, b0, b1, self->pk->n, 0, nsp, small_primes);
    t1 = ggh_walltime(t1);

    T0 += t0;
    T1 += t1;
    c += (r0<=r1);
    printf("\ri: %4d, eq: %4d, t0: %6.2fms, t1: %6.2fs", i+1, c, T0/1000.0/(i+1), T1/1000000.0/(i+1));
    fflush(0);
  }
  printf("\n");

  fmpz_poly_clear(b0);
  fmpz_poly_clear(b1);

  gghlite_clear(self, 1);

  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
}
