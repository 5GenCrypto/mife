#include <gghlite.h>

int main(int argc, char *argv[]) {
  gghlite_t self;
  gghlite_init(self, 32, 2);
  gghlite_print(self);
  gghlite_clear(self);
}
