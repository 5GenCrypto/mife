#include <mife/mife.h>
#include <mmap/mmap_gghlite.h>

int main(int argc, char **argv)
  { return mife_keygen_main(&gghlite_vtable, argc, argv); }
