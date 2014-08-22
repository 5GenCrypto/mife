#include <gghlite/gghlite.h>
#include "common.h"

int main(int argc, char *argv[]) {

  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "GGHLite Instance Generator", NULL);

  print_header("GGHLite Instance Generator", params);
  
  flint_rand_t randstate;
  flint_randinit(randstate);
  randstate->__randval  = params->seed;
  randstate->__randval2 = params->seed ^ 0x5555555555555555ULL;
  _flint_rand_init_gmp(randstate);
  gmp_randseed_ui(randstate->gmp_state, params->seed ^ 0xDEADBEEFDEADBEEFULL);
    
  gghlite_t self;
  gghlite_init(self, params->lambda, params->kappa, 1<<0, randstate);

  printf("\n");
  gghlite_print_params(self->pk);
  gghlite_clear(self, 1);

  flint_randclear(randstate);
  mpfr_free_cache();
}
