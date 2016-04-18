#include <mife/mife.h>
#include "mife_glue_gghlite.h"

int main(int argc, char **argv)
  { return mife_eval_main(&gghlite_vtable, argc, argv); }
