#include <sys/time.h>
#include <unistd.h>
#include <gghlite/gghlite.h>
#include "nike_header.h"

unsigned long long walltime(unsigned long long t0) {
  static time_t base_sec;
  struct timeval tp;
  gettimeofday(&tp, NULL);
  if (base_sec == 0)
    base_sec = tp.tv_sec;
  return (tp.tv_sec - base_sec) * 1000000 + tp.tv_usec - t0;
}

struct _cmdline_params_struct{
  long lambda;
  long N;
  uint64_t seed;
};

typedef struct _cmdline_params_struct cmdline_params_t[1];

#define DEFAULT_N       3
#define DEFAULT_LAMBDA 16
#define DEFAULT_SEED    0

static inline void print_help_and_exit(const char *name, const char *extra) {
  printf("####################################################################\n");
  printf(" %s\n",name);
  printf("####################################################################\n");
  printf("-l   security parameter λ > 0 (default: %d)\n", DEFAULT_LAMBDA);
  printf("-N   number of parties N > 2 (default: %d)\n", DEFAULT_N);
  printf("-s   seed (default: %d)\n",DEFAULT_SEED);
  if (extra)
    printf("%s\n", extra);
  abort();
}

void print_intro() {
  printf("Alice: Guys, guys! There is this thing now where we can agree on a key\n");
  printf("       by sending only one broadcast message each!\n");
  printf("  Bob: Wow, this sounds great!\n");
  printf("Alice: Yeah, and it is also really really efficient\n");
  printf("  Bob: Do tell!\n");
  printf("Alice: If κ = poly(log(λ)) then it is 'asymptotically close to optimal,\n");
  printf("       namely quasi-linear in the security parameter λ'!\n");
  printf("  Bob: Wow. I can't think of anything better, ever!\n");
  printf("Alice: It get's even better: security is defined in the\n");
  printf("       'common reference string model'!\n");
  printf("  Bob: That sounds very innocent and reassuring … hang on, what is it?\n");
  printf("Alice: It means that there is a *trusted setup*, so someone you really\n");
  printf("       really trust - like the government - gets to pick the parameters\n");
  printf("       for you! … and this trusted party can check out our shared key\n");
  printf("       to make sure we're doing it right and are not abusing\n");
  printf("       our civil liberties!\n");
  printf("  Bob: This is *so* cool! This makes me feel so much more secure already.\n");
  printf("Alice: Let's try it!\n");
  printf("\n");
}

int main(int argc, char *argv[]) {  
  cmdline_params_t cmdline_params;

  cmdline_params->N      =  DEFAULT_N;
  cmdline_params->lambda =  DEFAULT_LAMBDA;
  cmdline_params->seed   =  DEFAULT_SEED;

  int c;
  while ((c = getopt(argc, argv, "l:N:s:")) != -1) {
    switch(c) {
    case 'l':
      cmdline_params->lambda = (long)atol(optarg);
      break;
    case 'N':
      cmdline_params->N = (long)atol(optarg);
      break;
    case 's':
      cmdline_params->seed = (long)atol(optarg);
      break;
    case ':':  /* without operand */
      print_help_and_exit(name, NULL);
    case '?':
      print_help_and_exit(name, NULL);
    }
  }

  if (cmdline_params->N<3)
    print_help_and_exit(name, NULL);
  if (cmdline_params->lambda<1)
    print_help_and_exit(name, NULL);

  printf("####################################################################\n");
  printf("%s\n", name);
  printf(" λ: %3d, N: %2d                              seed: 0x%016lx\n",cmdline_params->lambda, cmdline_params->N, cmdline_params->seed);
  printf("#############################################all logs are base two##\n\n");

  
  flint_rand_t randstate;
  flint_randinit(randstate);
  randstate->__randval  = cmdline_params->seed;
  randstate->__randval2 = cmdline_params->seed ^ 0x5555555555555555ULL;
  _flint_rand_init_gmp(randstate);
  gmp_randseed_ui(randstate->gmp_state, cmdline_params->seed ^ 0xDEADBEEFDEADBEEFULL);

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
  
  fmpz_mod_poly_t  *e = calloc(cmdline_params->N, sizeof(fmpz_mod_poly_t));
  fmpz_mod_poly_t  *u = calloc(cmdline_params->N, sizeof(fmpz_mod_poly_t));
  fmpz_poly_t      *s = calloc(cmdline_params->N, sizeof(fmpz_poly_t));

  for(int i=0; i<cmdline_params->N; i++) {
    printf("%8s samples e_%d,",agents[i],i);
    fmpz_mod_poly_init_gghlite(e[i], params);
    gghlite_sample(e[i], params, 0, randstate);

    printf("computes u_%d and publishes it\n", i);
    fmpz_mod_poly_init_gghlite(u[i], params);
    gghlite_elevate(u[i], params, e[i], 1, 0, 1, randstate);
  }

  printf("\n");
  printf("wall time: %.2f s\n\n", walltime(t)/1000000.0); 
  printf("------------------------------------\n");
  printf("Step 3: \n");
  printf("------------------------------------\n");

  t = walltime(0);
  
  fmpz_mod_poly_t tmp;
  for(int i=0; i<cmdline_params->N; i++) {
    printf("%8s computes: s_%d = ", agents[i], i);
    fmpz_mod_poly_init_gghlite(tmp, params);
    fmpz_mod_poly_set_coeff_ui(tmp, 0, 1);
    for(int j=0; j<cmdline_params->N; j++) {
      if (i==j) {
        gghlite_mult(tmp, params, tmp, e[j]);
        printf("e_%d",j);
      } else {
        gghlite_mult(tmp, params, tmp, u[j]);
        printf("u_%d",j);
      }
      if(j<cmdline_params->N-1)
        printf("·");
    }
    printf("\n");
    fmpz_poly_init(s[i]);
    gghlite_extract(s[i], params, tmp, 1);
    fmpz_mod_poly_clear(tmp);
  }

  printf("\n");
  printf("wall time: %.2f s\n\n", walltime(t)/1000000.0); 
  printf("------------------------------------\n");
  printf("Check: \n");
  printf("------------------------------------\n");

int ret = 0;
  for(int i=1; i<cmdline_params->N; i++) {
    printf("s_%d == s_%d: ",0, i);
    if (fmpz_poly_equal(s[i], s[0])) {
      printf("TRUE\n");
    } else {
      printf("FALSE\n");
      ret = -1;
    }
  }
  printf("\n");
  
  flint_randclear(randstate);
  mpfr_free_cache();
  return ret;
}
