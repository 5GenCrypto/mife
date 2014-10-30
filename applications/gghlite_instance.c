#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include "common.h"

int main(int argc, char *argv[]) {

  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "GGHLite Instance Generator", NULL);

  print_header("GGHLite Instance Generator", params);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, params->seed, 1);

  gghlite_t self;
  gghlite_pk_init_params(self->pk, params->lambda, params->kappa, 1<<0, params->flags);
  gghlite_print_params(self->pk);

  printf("\n---\n");

  gghlite_init_instance(self, randstate);

  printf("\n---\n");

  gghlite_print_norms(self);

  printf("\n---\n");

  gghlite_print_times(self);

  gghlite_clear(self, 1);

  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
}
