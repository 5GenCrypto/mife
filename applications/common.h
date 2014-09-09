#ifndef _COMMON_H_
#define _COMMON_H_

#include <unistd.h>

#define DEFAULT_KAPPA   2
#define DEFAULT_LAMBDA 16
#define DEFAULT_SEED    0

static inline void print_help_and_exit(const char *name, const char *extra) {
  printf("####################################################################\n");
  printf(" %s\n",name);
  printf("####################################################################\n");
  printf("-l   security parameter λ > 0 (default: %d)\n", DEFAULT_LAMBDA);
  printf("-k   multi-linearity parameter k > 1 (default: %d)\n", DEFAULT_KAPPA);
  printf("-s   seed (default: %d)\n",DEFAULT_SEED);
  if (extra)
    printf("%s\n", extra);
  abort();
}

struct _cmdline_params_struct{
  long lambda;
  long kappa;
  mp_limb_t seed;
};

typedef struct _cmdline_params_struct cmdline_params_t[1];

static inline void print_header(const char *name, cmdline_params_t params) {
#ifdef GGHLITE_HEURISTICS
  int heuristics = 1;
#else
  int heuristics = 0;
#endif

  printf("####################################################################\n");
  printf("%s\n", name);
  printf(" λ: %3d, κ: %2d, heutistics: %d               seed: 0x%016lx\n",params->lambda, params->kappa, heuristics,  params->seed);
  printf("#############################################all logs are base two##\n\n");
}

static inline int parse_cmdline(cmdline_params_t params, int argc, char *argv[], const char *name, const char *extra) {
  params->kappa  =  DEFAULT_KAPPA;
  params->lambda =  DEFAULT_LAMBDA;
  params->seed   =  DEFAULT_SEED;

  int c;
  while ((c = getopt(argc, argv, "l:k:s:")) != -1) {
    switch(c) {
    case 'l':
      params->lambda = (long)atol(optarg);
      break;
    case 'k':
      params->kappa = (long)atol(optarg);
      break;
    case 's':
      params->seed = (long)atol(optarg);
      break;
    case ':':  /* without operand */
      print_help_and_exit(name, extra);
    case '?':
      print_help_and_exit(name, extra);
    }
  }

  if (params->kappa<2)
    print_help_and_exit(name, extra);
  if (params->lambda<1)
    print_help_and_exit(name, extra);
}


#endif /* _COMMON_H_ */
