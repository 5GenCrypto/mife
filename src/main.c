#include <gghlite.h>
#include <unistd.h>

void print_help_and_exit() {
  printf("GGHLite\n");
  printf("-l   security parameter λ > 0 (default: 32)\n");
  printf("-k   multi-linearity parameter k > 1 (default: 2)\n");
  abort();
};



int main(int argc, char *argv[]) {

  long kappa = 2;
  long lamba = 32;

  int c;
  while ((c = getopt(argc, argv, "l:k:s:")) != -1) {
    switch(c) {
    case 'l':
      lamba = (long)atol(optarg);
      break;
    case 'k':
      kappa = (long)atol(optarg);
      break;
    case ':':  /* -l or -k without operand */
      print_help_and_exit();
    case '?':
      print_help_and_exit();
    }
  }

  if (kappa<2)
    print_help_and_exit();
  if (lamba<1)
    print_help_and_exit();


  gghlite_t self;
  flint_rand_t state;
  flint_randinit(state);
  _flint_rand_init_gmp(state);
  gghlite_init_step1(self, lamba, kappa);
  printf("GGHLite\n");
  gghlite_print_params(self->pk);
  printf("---\n");
  
  gghlite_init_step2(self, 1ULL<<0, state);
  gghlite_clear(self, 1);
}
