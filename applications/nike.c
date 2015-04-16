#include <gghlite/gghlite.h>
#include "nike.h"

int main(int argc, char *argv[]) {
  cmdline_params_t cmdline_params;

  parse_cmdline(cmdline_params, argc, argv);
  int verbose = (cmdline_params->flags & GGHLITE_FLAGS_VERBOSE);
  int ret = 0;

  printf("####################################################################\n");
  if (verbose)
    printf("%s\n", name);
  else
    printf("NIKE\n");
  printf(" λ: %3ld, N: %2ld                              seed: 0x%016lx\n",
         cmdline_params->lambda,
         cmdline_params->N,
         cmdline_params->seed);
  printf("#############################################all logs are base two##\n\n");

  flint_rand_t randstate;
  flint_randinit_seed(randstate, cmdline_params->seed, 1);


  if (verbose)
    print_intro();

  printf("-----------------------------------------------------\n");
  printf("Step 1: GCHQ runs Setup:\n");
  printf("-----------------------------------------------------\n");

  uint64_t t = ggh_walltime(0);
  uint64_t t_total = ggh_walltime(0);

  gghlite_sk_t self;

  long mlm_lambda;
  if (cmdline_params->flags & GGHLITE_FLAGS_GDDH_HARD)
    mlm_lambda = cmdline_params->lambda;
  else
    mlm_lambda = 2*cmdline_params->lambda+1;

  gghlite_params_init(self->params, mlm_lambda, cmdline_params->N-1, 1<<0,
                         cmdline_params->flags);
  gghlite_params_print(self->params);
  printf("\n---\n");
  gghlite_sk_init(self, randstate);

  gghlite_params_t params;
  gghlite_params_ref(params, self);
  gghlite_sk_clear(self, 0);

  printf("\n");
  printf("wall time: %.2f s\n\n", ggh_walltime(t)/1000000.0);

  printf("-----------------------------------------------------\n");
  printf("Step 2: Publish");
  if (verbose)
    printf("\n-----------------------------------------------------\n");

  t = ggh_walltime(0);

  gghlite_enc_t *e = calloc(2, sizeof(gghlite_enc_t));
  gghlite_enc_t *u = calloc(2, sizeof(gghlite_enc_t));
  gghlite_enc_t *b = calloc(2, sizeof(gghlite_enc_t));
  gghlite_clr_t *s = calloc(2, sizeof(gghlite_clr_t));

  for(int i=0; i<2; i++) {
    if (verbose)
      printf("%8s samples e_%d, ",agents[i],i);
    gghlite_enc_init(e[i], params);
    gghlite_enc_sample(e[i], params, 0, 0, randstate);

    if (verbose)
      printf("computes u_%d and publishes it\n", i);
    gghlite_enc_init(u[i], params);
    gghlite_enc_raise0(u[i], params, e[i], 1, randstate);
  }

  if (verbose)
    printf("\n… and so on\n\n");
  else
    printf("  ");
  t = ggh_walltime(t);
  printf("wall time: %.2f s, per party: %.2f s\n", t/1000000.0,t/1000000.0/2);

  printf("-----------------------------------------------------\n");
  printf("Step 3: KeyGen");
  if (verbose)
    printf("\n-----------------------------------------------------\n");

  uint64_t t_sample = 0;

  gghlite_enc_init(b[0], params);
  gghlite_enc_mul(b[0], params, e[0], u[1]);

  t = ggh_walltime(0);
  gghlite_enc_init(b[1], params);
  gghlite_enc_mul(b[1], params, e[1], u[0]);

  if(verbose)
    printf("%8s computes e_%d·u_%d", agents[0], 0, 1);

  gghlite_enc_t u_i;
  gghlite_enc_init(u_i, params);
  gghlite_enc_t acc;
  gghlite_enc_init(acc, params);
  gghlite_enc_set_ui0(acc, 1, params);
  for(int j=2; j<cmdline_params->N; j++) {
    uint64_t t_s = ggh_walltime(0);
    gghlite_enc_sample(u_i, params, 0, 0, randstate);
    gghlite_enc_raise0(u_i, params, u_i, 1, randstate);
    t_s = ggh_walltime(t_s);
    t_sample += t_s;

    if(verbose)
      printf("·u_%d", j);

    gghlite_enc_mul(acc, params, acc, u_i);
  }
  gghlite_enc_clear(u_i);

  gghlite_enc_mul(b[0], params, acc, b[0]);
  gghlite_clr_init(s[0]);
  gghlite_enc_extract(s[0], params, b[0]);
  t = ggh_walltime(t) - t_sample;

  gghlite_enc_mul(b[1], params, acc, b[1]);
  gghlite_clr_init(s[1]);
  gghlite_enc_extract(s[1], params, b[1]);

  /* make sure we're rerandomising, i.e. not computing trivially the same thing */
  if (fmpz_mod_poly_equal(b[0], b[1])) ret = -1;

  /* the difference should be zero. We don't call gghlite_enc_is_zero() here because this might fail
     for small parameters: we are multiplying by e_i as well which makes our encodings of zero
     slightly bigger than they should be */
  gghlite_enc_t diff;
  gghlite_clr_t diff_c;
  gghlite_enc_init(diff, params);
  gghlite_clr_init(diff_c);
  gghlite_enc_sub(diff, params, b[0], b[1]);
  gghlite_enc_extract(diff_c, params, diff);
  if (!fmpz_poly_is_zero(diff_c))    ret = -1;
  gghlite_enc_clear(diff);
  gghlite_clr_clear(diff_c);

  if (verbose)
    printf("\n\n… and so on\n\n");
  else
    printf("  ");

  printf("wall time: %.2f s, per party: %.2f s\n", t/1000000.0,t/1000000.0);


  printf("-----------------------------------------------------\n");
  printf(" Check: ");
  if (verbose)
    printf("\n-----------------------------------------------------\n");

  if(verbose)
    printf("s_0 == s_1: ");
  assert(!fmpz_poly_is_zero(s[0]));

  if (!gghlite_clr_equal(s[1], s[0])) {
    ret = -1;
  }

  if(ret == 0) {
    verbose ? printf("TRUE\n") : printf ("+");
  } else {
    verbose ? printf("FALSE\n") : printf("-");
  }

  printf("\n");
  if (!verbose)
    printf("-----------------------------------------------------\n");


  printf("λ: %3ld, N: %2ld, n: %6ld, seed: 0x%08lx, prime: %u, success: %d, time: %10.2fs\n",
         cmdline_params->lambda, cmdline_params->N, params->n, cmdline_params->seed, cmdline_params->flags & GGHLITE_FLAGS_PRIME_G,
         ret==0, ggh_walltime(t_total)/1000000.0);

  for(int i=0; i<2; i++) {
    gghlite_enc_clear(e[i]);
    gghlite_enc_clear(u[i]);
    gghlite_clr_clear(s[i]);
    gghlite_enc_clear(b[i]);
  }
  free(e);
  free(u);
  free(s);
  free(b);
  gghlite_params_clear(params);

  flint_randclear(randstate);
  mpfr_free_cache();
  flint_cleanup();
  return ret;
}
