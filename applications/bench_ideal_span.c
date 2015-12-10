#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include "common.h"

int main(int argc, char *argv[]) {
  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "Playing with primes", NULL);

  print_header("Playing with primes", params);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, params->seed, 1);

  gghlite_sk_t self;
  gghlite_params_init(self->params, params->lambda, params->kappa, params->gamma, params->rerand, params->flags);
  gghlite_params_print(self->params);

  printf("\n---\n");

  gghlite_sk_init(self, randstate);

  const int nsp = _gghlite_nsmall_primes(self->params);
  mp_limb_t *primes = _fmpz_poly_oz_ideal_probable_prime_factors(self->params->n, nsp);

  printf("number of small primes: %4d\n",nsp);

  fmpz_poly_t b0;
  fmpz_poly_init2(b0, self->params->n);

  fmpz_poly_t b1;
  fmpz_poly_init2(b1, self->params->n);

  uint64_t t0, t1, T0, T1;
  int r0, r1;

  T0 = 0;
  T1 = 0;
  int c = 0;
  int d = 0;
  for(int i=0; i<100; i++) {
    fmpz_poly_sample_D(b0, self->D_g, randstate);
    fmpz_poly_sample_D(b1, self->D_g, randstate);

    t0 = ggh_walltime(0);
    r0 = fmpz_poly_oz_ideal_span(self->g, b0, b1, self->params->n, 1, primes);
    t0 = ggh_walltime(t0);

    t1 = ggh_walltime(0);
    r1 = fmpz_poly_oz_ideal_span(self->g, b0, b1, self->params->n, 0, primes);
    t1 = ggh_walltime(t1);

    T0 += t0;
    T1 += t1;
    c += (r0==r1);
    d += (r1!=0);
    printf("\ri: %4d, eq: %4d, #1: %4d, t0: %6.2fms, t1: %6.2fs", i+1, c, d, T0/1000.0/(i+1), ggh_seconds(T1)/(i+1));
    fflush(0);
  }
  printf("\n");

  free(primes);
  fmpz_poly_clear(b0);
  fmpz_poly_clear(b1);

  gghlite_sk_clear(self, 1);

  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
}
