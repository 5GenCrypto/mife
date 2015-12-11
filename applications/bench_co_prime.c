#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include <oz/oz.h>
#include "common.h"

int main(int argc, char *argv[]) {
  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "Bench co-primality test", NULL);

  print_header("Bench co-primality test", params);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, params->seed, 1);

  gghlite_sk_t self;
  gghlite_params_init(self->params, params->lambda, params->kappa, params->rerand, params->flags);
  gghlite_params_print(self->params);

  printf("\n---\n");

  gghlite_sk_init(self, randstate);

  printf("\n---\n");

  const int nsp = _gghlite_nsmall_primes(self->params);
  mp_limb_t *primes = _fmpz_poly_oz_ideal_probable_prime_factors(self->params->n, nsp);

  {
    uint64_t t1 = ggh_walltime(0);
    fmpz_poly_t h_mod; fmpz_poly_init(h_mod);
    fmpz_poly_oz_rem_small(h_mod, self->h, self->g, self->params->n);
    fmpz_poly_oz_coprime(self->g, h_mod, self->params->n, 0, primes);
    t1 = ggh_walltime(t1);
    printf("gcd(N(g), N(h%%g)): %.2fs (%.1f, %.1f), ", ggh_seconds(t1),
           fmpz_poly_2norm_log2(self->g),
           fmpz_poly_2norm_log2(h_mod));
    fmpz_poly_clear(h_mod);
    fflush(0);
  }

  {
    uint64_t t0 = ggh_walltime(0);
    fmpz_poly_oz_coprime(self->g, self->h, self->params->n, 0, primes);
    t0 = ggh_walltime(t0);
    printf("gcd(N(g), N(h)): %.2fs (%.1f, %.1f), ", ggh_seconds(t0),
           fmpz_poly_2norm_log2(self->g),
           fmpz_poly_2norm_log2(self->h));
    fflush(0);
  }

  free(primes);

  printf("\n");
  gghlite_sk_clear(self, 1);

  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
}
