#include <gghlite.h>

int main(int argc, char *argv[]) {
  gghlite_t self;
  flint_rand_t state;
  flint_randinit(state);
  _flint_rand_init_gmp(state);
  gghlite_init_step1(self, 8, 2);
  gghlite_print(self);
  gghlite_init_step2(self, state);
  gghlite_clear(self, 1);
}
