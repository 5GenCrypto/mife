#ifndef _ORE_CMDLINE_H
#define _ORE_CMDLINE_H

#include <gghlite/gghlite.h>
#include "matrix.h"
#include "mife.h"
#include "mife_glue.h"
#include "parse.h"
#include "util.h"

bool template_to_mife_pp(mife_pp_t pp, const template *const template, template_stats *const stats);
parse_result load_seed(location private_location, char *context, aes_randstate_t seed);

#endif /* ifndef _ORE_CMDLINE_H */
