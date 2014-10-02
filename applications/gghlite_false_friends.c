#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include "common.h"

int main(int argc, char *argv[]) {

  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "GGHLite Check Primality", NULL);

  print_header("GGHLite Check Primality", params);

  if (!(params->flags & GGHLITE_FLAGS_SLOPPY))
    ggh_die("Not setting GGHLITE_FLAGS_SLOPPY defeats the purpose of this program.");

  flint_rand_t randstate;
  flint_randinit_seed(randstate, params->seed, 1);

  gghlite_t self;
  gghlite_pk_init_params(self->pk, params->lambda, params->kappa, 1<<0, params->flags);
  gghlite_print_params(self->pk);
  printf("\n---\n");
  gghlite_init_instance(self, randstate);
  printf("---\n\n");

  gghlite_enc_t *u = calloc(params->kappa, sizeof(gghlite_enc_t));

  gghlite_enc_t tmp;

  for(int j=0; j<(1<<20); j++) {
    printf("\rtesting: %d",j);
    fflush(0);
    for(int i=0; i<params->kappa; i++) {
      gghlite_enc_init(u[i], self->pk);
      gghlite_sample(u[i], self->pk, 1, randstate);
    }
    gghlite_enc_init(tmp, self->pk);
    gghlite_enc_set_ui(tmp, 1);
    for(int i=0; i<params->kappa; i++) {
      gghlite_mult(tmp, self->pk, tmp, u[i]);
    }
    assert(!gghlite_is_zero(self->pk, tmp));
    gghlite_enc_clear(tmp);
    for(int i=0; i<params->kappa; i++) {
      gghlite_enc_clear(u[i]);
    }
  }
  printf(" PASS\n");
  free(u);
  gghlite_clear(self, 1);
  return 0;
}
