#include "mife.h"
#include "mife_glue.h"

int main(int argc, char **argv)
  { return mife_eval_main(&gghlite_vtable, argc, argv); }
