#include <gghlite/gghlite.h>

#include "nike.h"

int main(int argc, char *argv[]) {  
  cmdline_params_t cmdline_params;

  parse_cmdline(cmdline_params, argc, argv);
  
  printf("####################################################################\n");
  printf("%s\n", name);
  printf(" λ: %3d, N: %2d                              seed: 0x%016lx\n",
         cmdline_params->lambda,
         cmdline_params->N,
         cmdline_params->seed);
  printf("#############################################all logs are base two##\n\n");
  
  flint_rand_t randstate;
  flint_randinit_seed(randstate, cmdline_params->seed, 1);
  
  print_intro();

  printf("------------------------------------\n");
  printf("Step 1: GCHQ runs Setup:\n");
  printf("------------------------------------\n");

  unsigned long long t = walltime(0);
  
  gghlite_t self;
  gghlite_pk_init_params(self->pk, cmdline_params->lambda, cmdline_params->N-1, 1<<0);
  gghlite_print_params(self->pk);
  printf("\n---\n");  
  gghlite_init_instance(self, randstate);

  gghlite_pk_t params;
  gghlite_pk_ref(params, self);
  gghlite_clear(self, 0);

  printf("\n");
  printf("wall time: %.2f s\n\n", walltime(t)/1000000.0); 

  printf("------------------------------------\n");
  printf("Step 2: \n");
  printf("------------------------------------\n");

  t = walltime(0);
  
  gghlite_enc_t *e = calloc(cmdline_params->N, sizeof(gghlite_enc_t));
  gghlite_enc_t *u = calloc(cmdline_params->N, sizeof(gghlite_enc_t));
  gghlite_clr_t *s = calloc(cmdline_params->N, sizeof(gghlite_clr_t));

  for(int i=0; i<cmdline_params->N; i++) {
    printf("%8s samples e_%d,",agents[i],i);
    gghlite_enc_init(e[i], params);
    gghlite_sample(e[i], params, 0, randstate);

    printf("computes u_%d and publishes it\n", i);
    gghlite_enc_init(u[i], params);
    gghlite_elevate(u[i], params, e[i], 1, 0, 1, randstate);
  }

  printf("\n");
  printf("wall time: %.2f s\n\n", walltime(t)/1000000.0); 

  printf("------------------------------------\n");
  printf("Step 3: \n");
  printf("------------------------------------\n");

  t = walltime(0);
  
  gghlite_enc_t tmp;
  for(int i=0; i<cmdline_params->N; i++) {
    printf("%8s computes: s_%d = ", agents[i], i);
    gghlite_enc_init(tmp, params);
    gghlite_enc_set_ui(tmp, 1);
    for(int j=0; j<cmdline_params->N; j++) {
      if (i==j) {
        gghlite_mult(tmp, params, tmp, e[j]); printf("e_%d",j);
      } else {
        gghlite_mult(tmp, params, tmp, u[j]); printf("u_%d",j);
      }
      if(j<cmdline_params->N-1)
        printf("·");
    }
    printf("\n");
    gghlite_clr_init(s[i]);
    gghlite_extract(s[i], params, tmp);
    gghlite_enc_clear(tmp);
  }

  printf("\n");
  printf("wall time: %.2f s\n\n", walltime(t)/1000000.0); 
  printf("------------------------------------\n");
  printf("Check: \n");
  printf("------------------------------------\n");

int ret = 0;
  for(int i=1; i<cmdline_params->N; i++) {
    printf("s_%d == s_%d: ",0, i);
    if (gghlite_clr_equal(s[i], s[0])) {
      printf("TRUE\n");
    } else {
      printf("FALSE\n");
      ret = -1;
    }
  }
  printf("\n");

  for(int i=0; i<cmdline_params->N; i++) {
    gghlite_enc_clear(e[i]);
    gghlite_enc_clear(u[i]);
    gghlite_clr_clear(s[i]);
  }
  gghlite_pk_clear(params);
  
  flint_randclear(randstate);
  mpfr_free_cache();
  return ret;
}
