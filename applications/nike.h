#ifndef _NIKE_H_
#define _NIKE_H_

#include <sys/time.h>
#include <unistd.h>


static inline unsigned long long walltime(unsigned long long t0) {
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
  mp_limb_t seed;
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

const char *name = "                                                                           \n                                                               ._a(        \n              _                                             _ayZ\"`         \n              r                                         _swmQT^            \n             j(                                     ._wmQQBT`              \n            ]W                                   _wwQQQQ@\"                 \n           _QD                               _awQQQQQQT\"                   \n           dQk                           ._wmQQWWQQWY`                     \n          ]QQh                        _wyQQQQQQQQ@\"`                       \n         .QQQQ                    _awQQQWQQQQQQT\"                          \n         jWQQQ(               ._wQQWWWQQQQQQWP'                            \n         QWQQQQ,           _wmQQQWQQQQQQQQ@!                               \n        ]QQQQQWQ%      _awQQQWQQQQQQQQQQP\"                                 \n        dQQQQQQWQQwaawWQQQQQQQQQQQQQQWP'                                   \n       .QQQQQQQQQQWWWWQQQQQQQQQQQQQW!                                      \n       .QQQQQQQQQQQQQQQQQQQQQQQQQP\"                                        \n       =QQQQQQQQQQQQQQQQQQQQQQWP'                                          \n       .QQQQQQQQQQQQQQQQQQQQW!                                             \n        WQQQQQQQQQQQQQQQQQP\"                                               \n        ]WQQQQQQQQQQQQQWT'                                                 \n         )QQQQQQQQQQQ@!                                                    \n          -4QQQQQQ@?^                                                      \n             \"?!\"`                                                         \n                                                                           \n";

static inline int parse_cmdline(cmdline_params_t cmdline_params, int argc, char *argv[]) {
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
}


const char *agents[26] = {"Alice",
                          "Bob",
                          "Charlie",
                          "Dan",
                          "Eve",
                          "Frank",
                          "Georgia",
                          "Hannah",
                          "Ilena",
                          "Jacob",
                          "Karl",
                          "Ludwig",
                          "Mallory",
                          "Nick",
                          "Oscar",
                          "Peggy",
                          "Q",
                          "Robert",
                          "Sybil",
                          "Trent",
                          "Uwe",
                          "Viet",
                          "Wendy",
                          "Xavier",
                          "Y",
                          "Zorro"};




static inline void print_intro() {
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

#endif /* _NIKE_H_ */
