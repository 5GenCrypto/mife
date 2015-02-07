#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include <oz/oz.h>
#include "common.h"

double fmpz_poly_eucl_norm_log2(fmpz_poly_t f) {
  mpfr_t norm;
  mpfr_init(norm);
  fmpz_poly_eucl_norm_mpfr(norm, f, MPFR_RNDN);
  mpfr_log2(norm, norm, MPFR_RNDN);
  double d = mpfr_get_d(norm, MPFR_RNDN);
  mpfr_clear(norm);
  return d;
}

int main(int argc, char *argv[]) {
  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "Bench co-primality test", NULL);

  print_header("Bench co-primality test", params);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, params->seed, 1);

  gghlite_t self;
  gghlite_pk_init_params(self->pk, params->lambda, params->kappa, params->rerand, params->flags);
  gghlite_print_params(self->pk);

  printf("\n---\n");

  gghlite_init_instance(self, randstate);

  printf("\n---\n");

  const int nsp = _gghlite_nsmall_primes(self->pk);
  mp_limb_t *primes = _fmpz_poly_oz_ideal_probable_prime_factors(self->pk->n, nsp);

  {
    uint64_t t1 = ggh_walltime(0);
    fmpz_poly_t h_mod; fmpz_poly_init(h_mod);
    fmpz_poly_oz_rem_small(h_mod, self->h, self->g, self->pk->n);
    fmpz_poly_oz_coprime(self->g, h_mod, self->pk->n, 0, primes);
    t1 = ggh_walltime(t1);
    printf("gcd(N(g), N(h%%g)): %.2fs (%.1f, %.1f), ", t1/1000000.0,
           fmpz_poly_eucl_norm_log2(self->g),
           fmpz_poly_eucl_norm_log2(h_mod));
    fmpz_poly_clear(h_mod);
    fflush(0);
  }

  {
    uint64_t t0 = ggh_walltime(0);
    fmpz_poly_oz_coprime(self->g, self->h, self->pk->n, 0, primes);
    t0 = ggh_walltime(t0);
    printf("gcd(N(g), N(h)): %.2fs (%.1f, %.1f), ", t0/1000000.0,
           fmpz_poly_eucl_norm_log2(self->g),
           fmpz_poly_eucl_norm_log2(self->h));
    fflush(0);
  }

  free(primes);

  printf("\n");
  gghlite_clear(self, 1);

  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
}
