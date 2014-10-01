#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include "common.h"

#define NTRIALS 256

int main(int argc, char *argv[]) {

  cmdline_params_t params;
  parse_cmdline(params, argc, argv, "GGHLite Instance Generator", NULL);

  print_header("GGHLite Noise Distribution", params);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, params->seed, 1);

  gghlite_t self;
  gghlite_pk_init_params(self->pk, params->lambda, params->kappa, 1<<0);
  gghlite_print_params(self->pk);

  printf("\n---\n");

  gghlite_init_instance(self, randstate);

  printf("\n---\n");

  gghlite_print_norms(self);

  gghlite_enc_t tmp;
  gghlite_enc_init(tmp, self->pk);
  gghlite_enc_t *u = calloc(params->kappa, sizeof(gghlite_enc_t));

  gghlite_clr_t out;
  gghlite_clr_init(out);

  mpfr_t norm;
  mpfr_init2(norm, _gghlite_prec(self->pk));

  mpfr_t acc;
  mpfr_init2(acc, _gghlite_prec(self->pk));
  mpfr_set_ui(acc, 0, MPFR_RNDN);

  mpfr_t max;
  mpfr_init2(max, _gghlite_prec(self->pk));
  mpfr_set_ui(max, 0, MPFR_RNDN);

  for(int k=0; k<params->kappa; k++)
      gghlite_enc_init(u[k], self->pk);

  for(int i=0; i<NTRIALS; i++) {
    gghlite_enc_set_ui(tmp, 1);
    for(int k=0; k<params->kappa; k++) {
      gghlite_enc_set_ui(u[k], 0);
      gghlite_enc(u[k], self->pk, u[k], 1, 1, randstate);
      gghlite_mult(tmp, self->pk, tmp, u[k]);
    }

    fmpz_mod_poly_mulmod(tmp, self->pk->pzt, tmp, self->pk->modulus);
    fmpz_poly_set_fmpz_mod_poly(out, tmp);
    _fmpz_vec_2norm_mpfr(norm, out->coeffs, self->pk->n);
    mpfr_add(acc, acc, norm, MPFR_RNDN);
    if (mpfr_cmp(norm, max)>0)
      mpfr_set(max, norm, MPFR_RNDN);
  }

  mpfr_div_ui(acc, acc, NTRIALS, MPFR_RNDN);

  mpfr_log2(max, max, MPFR_RNDN);
  mpfr_log2(acc, acc, MPFR_RNDN);


  const double bound = (1.0-mpfr_get_d(self->pk->ell_q, MPFR_RNDN)) * fmpz_sizeinbase(self->pk->q,2);

  printf("\n  log(avg): %6.1f\n  log(max): %6.1f ?< %6.1f\n", mpfr_get_d(acc, MPFR_RNDN), mpfr_get_d(max, MPFR_RNDN), bound);

  for(int k=0; k<params->kappa; k++)
      gghlite_enc_clear(u[k]);

  free(u);

  gghlite_clear(self, 1);

  flint_randclear(randstate);
  mpfr_free_cache();
}
