#ifndef _ORE_CMDLINE_H
#define _ORE_CMDLINE_H

#include <gghlite/gghlite.h>
#include "mbp_types.h"
#include "mife.h"
#include "mbp_glue.h"
#include "util.h"

bool mbp_template_to_mife_pp(mife_pp_t pp, const mbp_template *const template, mbp_template_stats *const stats);
parse_result load_seed(location private_location, char *context, aes_randstate_t seed);

#endif /* ifndef _ORE_CMDLINE_H */
