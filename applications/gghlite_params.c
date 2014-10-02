#include "gghlite/gghlite.h"
#include "common.h"

int main(int argc, char *argv[]) {
  cmdline_params_t params;
  const char *name = "GGHLite Parameters";
  parse_cmdline(params, argc, argv, name, NULL);
  print_header(name, params);
  gghlite_t self;

  gghlite_pk_init_params(self->pk, params->lambda, params->kappa, 1<<0, params->flags);
  gghlite_print_params(self->pk);

  gghlite_clear(self, 1);
  mpfr_free_cache();
}
