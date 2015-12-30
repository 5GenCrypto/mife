#include <gghlite/gghlite.h>
#include "common.h"

int main(int argc, char *argv[]) {
  cmdline_params_t cmdline_params;

  const char *name =  "Jigsaw Puzzles";
  parse_cmdline(cmdline_params, argc, argv, name, NULL);

  print_header(name, cmdline_params);

  aes_randstate_t randstate;
  aes_randinit_seed(randstate, cmdline_params->shaseed);

  uint64_t t = ggh_walltime(0);
  uint64_t t_total = ggh_walltime(0);

  uint64_t t_gen = 0;

  gghlite_sk_t self;


  gghlite_jigsaw_init(self,
                      cmdline_params->lambda,
                      cmdline_params->kappa,
                      cmdline_params->flags,
                      randstate);

  printf("\n");
  gghlite_params_print(self->params);
  printf("\n---\n\n");

  t_gen = ggh_walltime(t);
  printf("1. GGH InstGen wall time:                 %8.2f s\n", ggh_seconds(t_gen));

  t = ggh_walltime(0);

  fmpz_t p; fmpz_init(p);
  fmpz_poly_oz_ideal_norm(p, self->g, self->params->n, 0);

  fmpz_t a[5];
  for(long k=0; k<5; k++)
    fmpz_init(a[k]);

  fmpz_t acc;  fmpz_init(acc);
  fmpz_set_ui(acc, 1);

  // reducing these elements to something small mod g is expensive, for benchmarketing purposes we
  // hence avoid this costly step most of the time by doing it only twice and by then computing the
  // remaining elements from those two elements using cheap operations. The reported times still
  // refer to the expensive step for fairness.
  fmpz_randm_aes(a[0], randstate, p);
  fmpz_mul(acc, acc, a[0]);
  fmpz_mod(acc, acc, p);

  fmpz_randm_aes(a[1], randstate, p);
  fmpz_mul(acc, acc, a[1]);
  fmpz_mod(acc, acc, p);

  fmpz_set(a[3], a[0]);
  fmpz_set(a[4], a[1]);

  for(long k=2; k<cmdline_params->kappa; k++) {
    fmpz_mul(a[2], a[0], a[1]);
    fmpz_add_ui(a[2], a[2], k);
    fmpz_mod(a[2], a[2], p);
    fmpz_mul(acc, acc, a[2]);
    fmpz_mod(acc, acc, p);
    fmpz_set(a[1], a[0]);
    fmpz_set(a[0], a[2]);
  }

  printf("2. Sampling from U(Z_p) wall time:        %8.2f s\n", ggh_seconds(ggh_walltime(t)));

  gghlite_clr_t e[3];
  gghlite_clr_init(e[0]);
  gghlite_clr_init(e[1]);
  gghlite_clr_init(e[2]);

  gghlite_enc_t u_k;
  gghlite_enc_init(u_k, self->params);

  gghlite_enc_t left;
  gghlite_enc_init(left, self->params);
  gghlite_enc_set_ui0(left, 1, self->params);

  fmpz_poly_set_coeff_fmpz(e[0], 0, a[3]);

  int GAMMA = 20;
  uint64_t t_enc = ggh_walltime(0);
	int group0[GAMMA];
	memset(group0, 0, GAMMA * sizeof(int));
	group0[0] = 1;
  gghlite_enc_set_gghlite_clr(u_k, self, e[0], 1, group0, 1, randstate);
  t_enc = ggh_walltime(t_enc);

  printf("3. Encoding wall time:                    %8.2f s (per elem)\n", ggh_seconds(t_enc));
  fflush(0);

  uint64_t t_mul = 0;
  t = ggh_walltime(0);
  gghlite_enc_mul(left, self->params, left, u_k);
  t_mul += ggh_walltime(t);

  fmpz_poly_set_coeff_fmpz(e[1], 0, a[4]);
	int group1[GAMMA];
	memset(group1, 0, GAMMA * sizeof(int));
	group1[1] = 1;
  gghlite_enc_set_gghlite_clr(u_k, self, e[1], 1, group1, 1, randstate);
  t = ggh_walltime(0);
  gghlite_enc_mul(left, self->params, left, u_k);
  t_mul += ggh_walltime(t);

  const mp_bitcnt_t prec = (self->params->n/4 < 8192) ? 8192 : self->params->n/4;
  const oz_flag_t flags = (self->params->flags & GGHLITE_FLAGS_VERBOSE) ? OZ_VERBOSE : 0;
  _fmpz_poly_oz_rem_small_iter(e[0], e[0], self->g, self->params->n, self->g_inv, prec, flags);
  _fmpz_poly_oz_rem_small_iter(e[1], e[1], self->g, self->params->n, self->g_inv, prec, flags);

  for(long k=2; k<cmdline_params->kappa; k++) {
    fmpz_poly_oz_mul(e[2], e[0], e[1], self->params->n);
    assert(fmpz_poly_degree(e[2])>=0);
    fmpz_add_ui(e[2]->coeffs, e[2]->coeffs, k);
    _fmpz_poly_oz_rem_small_iter(e[2], e[2], self->g, self->params->n, self->g_inv, prec, flags);
		int groupk[GAMMA];
		memset(groupk, 0, GAMMA * sizeof(int));
		groupk[k] = 1;
		gghlite_enc_set_gghlite_clr(u_k, self, e[2], 1, groupk, 1, randstate);
    t = ggh_walltime(0);
    gghlite_enc_mul(left, self->params, left, u_k);
    t_mul += ggh_walltime(t);
    fmpz_poly_set(e[1], e[0]);
    fmpz_poly_set(e[0], e[2]);
  }

  printf("4. Multiplication wall time:              %8.4f s\n", ggh_seconds(t_mul));
  fflush(0);

  t = ggh_walltime(0);
  gghlite_enc_t rght;
  gghlite_enc_init(rght, self->params);
  gghlite_enc_set_ui0(rght, 1, self->params);

  fmpz_poly_t tmp; fmpz_poly_init(tmp);
  fmpz_poly_set_coeff_fmpz(tmp, 0, acc);
  gghlite_enc_set_gghlite_clr0(rght, self, tmp, randstate);

  for(long k=0; k<cmdline_params->kappa; k++) {
    gghlite_enc_mul(rght, self->params, rght, self->z_inv[k]);
  }

  printf("5. RHS generation wall time:              %8.2f s\n", ggh_seconds(ggh_walltime(t)));
  t = ggh_walltime(0);

  gghlite_enc_sub(rght, self->params, rght, left);
  int status = 1 - gghlite_enc_is_zero(self->params, rght);

  gghlite_clr_t clr; gghlite_clr_init(clr);
  gghlite_enc_extract(clr, self->params, rght);
  double size = fmpz_poly_2norm_log2(clr);
  gghlite_clr_clear(clr);

  printf("6. Checking correctness wall time:        %8.2f s\n", ggh_seconds(ggh_walltime(t)));
  printf("   Correct:                               %8s (%8.2f)\n\n", (status == 0) ? "TRUE" : "FALSE", size);

  for(long i=0; i<5; i++) {
    fmpz_clear(a[i]);
  }

  gghlite_clr_clear(e[0]);
  gghlite_clr_clear(e[1]);
  gghlite_clr_clear(e[2]);
  gghlite_enc_clear(u_k);

  fmpz_clear(acc);
  gghlite_enc_clear(left);
  gghlite_enc_clear(rght);
  gghlite_clr_clear(tmp);
  fmpz_clear(p);
  gghlite_sk_clear(self, 1);

  t_total = ggh_walltime(t_total);
  printf("λ: %3ld, κ: %2ld, n: %6ld, seed: 0x%08lx, success: %d, gen: %10.2fs, enc: %8.2fs, mul: %8.4fs, time: %10.2fs\n",
         cmdline_params->lambda, cmdline_params->kappa, self->params->n, cmdline_params->seed,
         status==0,
         ggh_seconds(t_gen), ggh_seconds(t_enc), ggh_seconds(t_mul), ggh_seconds(t_total));

  aes_randclear(randstate);
  mpfr_free_cache();
  flint_cleanup();
  return status;
}
