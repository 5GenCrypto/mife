#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include "common.h"

#define NTRIALS 256

int main(int argc, char *argv[]) {

  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "GGHLite Instance Generator", NULL);

  print_header("GGHLite Noise Distribution", params);

  aes_randstate_t randstate;
  aes_randinit_seed(randstate, params->shaseed, NULL);


  gghlite_sk_t self;

  gghlite_params_init(self->params, params->lambda, params->kappa, params->rerand, params->flags);
  gghlite_params_print(self->params);

  printf("\n---\n");

  gghlite_sk_init(self, randstate);

  printf("\n---\n");

  gghlite_sk_print_norms(self);

  gghlite_enc_t tmp;
  gghlite_enc_init(tmp, self->params);
  gghlite_enc_t *u = calloc(params->kappa, sizeof(gghlite_enc_t));

  gghlite_clr_t out;
  gghlite_clr_init(out);

  mpfr_t norm;
  mpfr_init2(norm, _gghlite_prec(self->params));

  mpfr_t acc;
  mpfr_init2(acc, _gghlite_prec(self->params));
  mpfr_set_ui(acc, 0, MPFR_RNDN);

  mpfr_t max;
  mpfr_init2(max, _gghlite_prec(self->params));
  mpfr_set_ui(max, 0, MPFR_RNDN);

  for(int k=0; k<params->kappa; k++)
      gghlite_enc_init(u[k], self->params);

  for(int i=0; i<NTRIALS; i++) {
    gghlite_enc_set_ui0(tmp, 1, self->params);
    for(int k=0; k<params->kappa; k++) {
      gghlite_enc_set_ui0(u[k], 0, self->params);
      gghlite_enc_raise0(u[k], self->params, u[k], 1, randstate);
      gghlite_enc_mul(tmp, self->params, tmp, u[k]);
    }
    _gghlite_enc_extract_raw(out, self->params, tmp);
    fmpz_poly_2norm_mpfr(norm, out, MPFR_RNDN);
    mpfr_add(acc, acc, norm, MPFR_RNDN);
    if (mpfr_cmp(norm, max)>0)
      mpfr_set(max, norm, MPFR_RNDN);
  }

  mpfr_div_ui(acc, acc, NTRIALS, MPFR_RNDN);

  mpfr_log2(max, max, MPFR_RNDN);
  mpfr_log2(acc, acc, MPFR_RNDN);


  const double bound = (1.0-mpfr_get_d(self->params->xi, MPFR_RNDN)) * fmpz_sizeinbase(self->params->q,2);

  printf("\n  log(avg): %6.1f\n  log(max): %6.1f ?< %6.1f\n", mpfr_get_d(acc, MPFR_RNDN), mpfr_get_d(max, MPFR_RNDN), bound);

  for(int k=0; k<params->kappa; k++)
      gghlite_enc_clear(u[k]);

  free(u);

  gghlite_sk_clear(self, 1);

  aes_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
}
