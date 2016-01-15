#include "cmdline.h"

const mife_flag_t mife_flags = MIFE_DEFAULT;

bool template_to_mife_pp(mife_pp_t pp, const template *const template, template_stats *const stats) {
	if(!template_to_template_stats(template, stats)) return false;
	mife_init_params(pp, mife_flags);
	mife_mbp_set(stats, pp, stats->positions_len,
		template_stats_to_params,
		template_stats_to_dimensions,
		template_stats_to_position,
		template_stats_to_cleartext,
		template_stats_to_result);
	return true;
}
