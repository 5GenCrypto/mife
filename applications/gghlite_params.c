#include "gghlite/gghlite.h"
#include "common.h"

int main(int argc, char *argv[]) {
  cmdline_params_t params;
  const char *name = "GGHLite Parameters";
  parse_cmdline(params, argc, argv, name, NULL);
  print_header(name, params);
  gghlite_sk_t self;

  gghlite_params_init(self->params, params->lambda, params->kappa, params->rerand, params->flags);
  gghlite_params_print(self->params);

  gghlite_sk_clear(self, 1);
  flint_cleanup();
  mpfr_free_cache();
}
