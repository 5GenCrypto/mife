#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include "common.h"

int main(int argc, char *argv[]) {
  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "GGHLite Check Primality", NULL);

  print_header("GGHLite False Positives Check", params);

  if (params->flags & GGHLITE_FLAGS_PRIME_G)
    ggh_die("Setting GGHLITE_FLAGS_PRIME_G defeats the purpose of this application.");

  aes_randstate_t randstate;
  aes_randinit_seed(randstate, params->shaseed);

  gghlite_sk_t self;
  gghlite_params_init(self->params, params->lambda, params->kappa, params->rerand, params->flags);
  gghlite_params_print(self->params);
  printf("\n---\n");
  gghlite_sk_init(self, randstate);
  printf("---\n\n");

  mpfr_t sqrt_q;
  mpfr_init2(sqrt_q, fmpz_sizeinbase(self->params->q, 2));
  _gghlite_params_get_q_mpfr(sqrt_q, self->params, MPFR_RNDN);
  mpfr_sqrt(sqrt_q, sqrt_q, MPFR_RNDN);
  mpfr_mul_d(sqrt_q, sqrt_q, S_TO_SIGMA, MPFR_RNDN);

  mp_limb_t *primes = _fmpz_poly_oz_ideal_small_prime_factors(self->params->n, self->params->lambda);

  int count = 0;
  while(1) {
    fmpz_poly_sample_sigma(self->g, self->params->n, self->params->sigma, randstate);
    fmpz_poly_sample_sigma(self->h, self->params->n, sqrt_q, randstate);

    if(fmpz_poly_oz_coprime(self->g, self->h, self->params->n, 0, primes))
      printf("-");
    else {
      printf("+");
      break;
    }
    fflush(0);
    count++;
  }
  printf(" prob: %7.5f%%\n", 1.0 - 1.0/(double)count);
  free(primes);
  mpfr_clear(sqrt_q);

  gghlite_enc_t *u = calloc(params->kappa, sizeof(gghlite_enc_t));

  gghlite_enc_t tmp;
  gghlite_enc_init(tmp, self->params);
  for(int i=0; i<params->kappa; i++)
      gghlite_enc_init(u[i], self->params);

  uint64_t t = ggh_walltime(0);
  int m = (1<<16);
  for(int j=0; j<m; j++) {
    printf("\rtesting: %d",j);
    fflush(0);
    for(int i=0; i<params->kappa; i++) {
      gghlite_enc_sample(u[i], self->params, 1, 0, randstate);
    }
    gghlite_enc_set_ui0(tmp, 1, self->params);
    for(int i=0; i<params->kappa; i++) {
      gghlite_enc_mul(tmp, self->params, tmp, u[i]);
    }
    assert(!gghlite_enc_is_zero(self->params, tmp));
  }
  for(int i=0; i<params->kappa; i++)
      gghlite_enc_clear(u[i]);
  gghlite_enc_clear(tmp);
  t = ggh_walltime(t);
  printf(" PASS, wall time: %8.2fs, per sample: %8.2fms\n", ggh_seconds(t), t/1000.0/m);
  free(u);
  gghlite_sk_clear(self, 1);
  flint_cleanup();
  return 0;
}
