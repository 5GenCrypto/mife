#ifndef _NIKE_H_
#define _NIKE_H_

#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <gghlite/gghlite-defs.h>

struct _cmdline_params_struct{
  long lambda;
  long N;
  uint64_t flags;
  mp_limb_t seed;
};

typedef struct _cmdline_params_struct cmdline_params_t[1];

#define DEFAULT_N       3
#define DEFAULT_LAMBDA 20
#define DEFAULT_SEED    0

static inline void print_help_and_exit(const char *name, const char *extra) {
  printf("####################################################################\n");
  printf(" %s\n",name);
  printf("####################################################################\n");
  printf("-l   security parameter λ > 0 (default: %d)\n", DEFAULT_LAMBDA);
  printf("-N   number of parties N > 2 (default: %d)\n", DEFAULT_N);
  printf("-s   seed (default: %d)\n",DEFAULT_SEED);
  printf("-p   enforce that g generates a prime ideal (default: False)\n");
  printf("-d   pick parameters to make GDDH hard (default: False)\n");
  printf("-v   be more verbose (default: False)\n");
  if (extra)
    printf("%s\n", extra);
  abort();
}

//const char *name = "                                                                           \n                                                               ._a(        \n              _                                             _ayZ\"`         \n              r                                         _swmQT^            \n             j(                                     ._wmQQBT`              \n            ]W                                   _wwQQQQ@\"                 \n           _QD                               _awQQQQQQT\"                   \n           dQk                           ._wmQQWWQQWY`                     \n          ]QQh                        _wyQQQQQQQQ@\"`                       \n         .QQQQ                    _awQQQWQQQQQQT\"                          \n         jWQQQ(               ._wQQWWWQQQQQQWP'                            \n         QWQQQQ,           _wmQQQWQQQQQQQQ@!                               \n        ]QQQQQWQ%      _awQQQWQQQQQQQQQQP\"                                 \n        dQQQQQQWQQwaawWQQQQQQQQQQQQQQWP'                                   \n       .QQQQQQQQQQWWWWQQQQQQQQQQQQQW!                                      \n       .QQQQQQQQQQQQQQQQQQQQQQQQQP\"                                        \n       =QQQQQQQQQQQQQQQQQQQQQQWP'                                          \n       .QQQQQQQQQQQQQQQQQQQQW!                                             \n        WQQQQQQQQQQQQQQQQQP\"                                               \n        ]WQQQQQQQQQQQQQWT'                                                 \n         )QQQQQQQQQQQ@!                                                    \n          -4QQQQQQ@?^                                                      \n             \"?!\"`                                                         \n                                                                           \n";

//const char *name = "\n[8m               [0;10m[1mygggr[0;10m[8m [0;10mj[1mgggzjggg[0;10m6[1myggg[0;10m'_[1mwgggngggggggr[0;10mj[1mu[0;10mp[8m           [0;10m\n[8m              [0;10m[1mjQQQQ[[0;10m[8m [0;10m[1mmQQQ\\QQQPjQQQ[[0;10mq[1mQQQD[0;10m5[1mQQQWBUBB`[0;10m[8m-[0;10m-[8m:           [0;10m\n[8m             [0;10mj[1mQQQQQ[]QQQfQQQQ\\QQQDwQQQ![0;10m[8m [0;10m[1mjQQQ'[0;10m[8m                   [0;10m\n[8m            :[0;10m[1mmQQQQQ[0;10mq[1mQQQ@jQQQ[mQQQWWQP[0;10m'[8m [0;10mj[1mQQQP[0;10m[8mI                   [0;10m\n[8m            [0;10m[1mjQQQQQQmQQQ\\QQQPjWQQQQQf[0;10m[8m  [0;10m_[1mQQQQQQQ@[0;10m'[8m                [0;10m\n[8m           [0;10mj[1mWQQWQQQQQQFmQQQ\\QQQQQQQL[0;10m[8m  [0;10m[1mjQQQBBHB'[0;10m[8mn                [0;10m\n[8m        [0;10m,[8m: [0;10m[1myQQQ4WQQQQ@jQQQfdQQQ5WQQk[0;10m[8m [0;10mj[1mQQQP[0;10m[8m.               [0;10m[1m_[0;10mg[1ma,!`[0;10m\n[8m      [0;10m_[1m%[0;10m[8m|.[0;10m[1m]QQQF]QQQQQ\\QQQ@]WQQF[0;10mj[1mQQQk[0;10m[8m [0;10m[1mmQQQ'[0;10m[8m:        [0;10m_[1m_aawmD?^[0;10m`[8m   [0;10m\n[8m     [0;10mj[1mm'[0;10m[8m+[0;10m[1m_QQQ@[0;10m`[1mjWQQQPyQQQ\\WWQ@[0;10m'j[1mQQQ#jQQQQgmgg[0;10mLg[1maawQQQP?^[0;10m[8m        [0;10m\n[8m    [0;10mq[1mQ@[0;10m[8m  [0;10m[1mjQQQ'[0;10m[8mI[0;10m[1mjQQQW]QQQfjQQQ'[0;10m[8m [0;10mj[1mQQQQQQQQQQQQQQQQWD?\"[0;10m`[8m           [0;10m\n[8m   [0;10mq[1mQQm[0;10m[8m  [0;10m=:[8mIIIi[0;10m-[8mXI+ [0;10m-[8mnvv=[0;10m;[8m#ii  =Q[0;10mgg[1mawWQQQQQQQV?^[0;10m[8m                [0;10m\n[8m  [0;10mQ[1mQQQQ,[0;10m[8m  --             .=[0;10mq[1mawwQQQQQQQQQQ@?\"[0;10m[8m                    [0;10m\n[8m [0;10mj[1mQQQQQQ[0;10mp[8m [0;10m[1mI}[0;10m[8m        [0;10m[1m_[0;10mg[1mawmQQQQWWWQQQQQU?^[0;10m`[8m                       [0;10m\n_[1mQQQQQQQQmaa[0;10mgg[1maaawmQQQQWQQQQQQQQQD?\"[0;10m[8mv                           [0;10m\n[1m]QQQQQQQQQQQQQQWWWQQQQQQQQQQQ8?^[0;10m'[8m                               [0;10m\n[1mjQQQQQQQQQQQQQQQQQQQQQQQQD?\"[0;10m[8m                                    [0;10m\n[1m]QQQQQQQQQQQQQQQQQQQQB?\"[0;10m'[8m                                       [0;10m\n[8m=[0;10m[1m4QQQQQQQQQQQQQQQD?\"[0;10m[8m-  -                                        [0;10m\n[8m  [0;10m[1m\"4WQQQQQQQWP?\"[0;10m'[8m                                               [0;10m\n[8m      [0;10m[1m\"!!^[0;10m`[8m                                                     [0;10m\n";

const char *name = "\n               ygggr _gggzjggg(wggg`.wgggngggggggr_o/           \n              jQQQQ[ mQQQ\\QQQPjWQQ[_QQQP\\WQQWBUBB`  .           \n             _QQQQQ[]QQQfQQQQ\\WQQDwQQQ! jQQQ'                   \n             mQQQQQ<QQQ@jQQQ[mQQQQWQP` <QQQP                    \n            jQQQQQQmQQQ\\QQQPjWQQQQQf   QQQQQQQ@                 \n           _QQQQQQQQQQfmQQQ\\QQQQQQQL  jQQQBHBB'                 \n        .  yQQQ4WQQQQ@jQQQfdQQQ5WQQk _QQQP                _aa,!`\n       J  ]QQQF]WQQQQ\\QQQ@]QQQF)WQQk mQQQ'         ._aawmD?^    \n     _m' _QQQ@ jQQQQPyQQQ\\WWQ@ )WQQEjQQQQgmgg=saawQQQP?^        \n    <Q@  jQQQ'.jQQQW]QQQfjQQQ' =QQQQQQQQQQQQQQWWWD?\"            \n   <QQm  ..  .  .    . . ..      =sawWQQQQQQQV?^                \n  <QQQQ,                   _awwQQQQWWWQQQ@?\"                    \n _QQQQQQ, I}        __awmQQQQQWWQQQQQU?^                        \n.QQQQQQQQmaaa_aaawQQQQQQQQQQQQQQQD?\"                            \n]QQQQQQQQQQWWQQQWWWQQQQQQQQQQ8?^                                \njQQQQQQQQQQQQQQQQQQQQQQQQD?\"                                    \n]QQQQQQQQQQQQQQQQQQQQB?\"                                        \n 4QQQQQQQQQQQQQQQD?\"                                            \n  \"4WQQQQQQQWP?\"                                                \n      \"!!^                                                      \n";



static inline int parse_cmdline(cmdline_params_t cmdline_params, int argc, char *argv[]) {
  cmdline_params->N      =  DEFAULT_N;
  cmdline_params->lambda =  DEFAULT_LAMBDA;
  cmdline_params->seed   =  DEFAULT_SEED;
  cmdline_params->flags  =  GGHLITE_FLAGS_DEFAULT;

  int c;
  while ((c = getopt(argc, argv, "l:N:s:pvd")) != -1) {
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
    case 'p':
      cmdline_params->flags |= GGHLITE_FLAGS_PRIME_G;
      break;
    case 'v':
      cmdline_params->flags |= GGHLITE_FLAGS_VERBOSE;
      break;
    case 'd':
      cmdline_params->flags |= GGHLITE_FLAGS_GDDH_HARD;
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
  return 0;
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
                          "Silke",
                          "Trent",
                          "Uwe",
                          "Viet",
                          "Wendy",
                          "Xavier",
                          "Y",
                          "Ziggy"};




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
