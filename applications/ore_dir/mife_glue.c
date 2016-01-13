#include "mife_glue.h"
#include "util.h"

void template_stats_free(template_stats stats) {
	int i, j;
	const int len = stats.positions_len;
	free(stats.positions);
	if(stats.indexes != NULL)
		for(i = 0; i < len; i++)
			free(stats.indexes[i]);
	free(stats.indexes);
	free(stats.indexes_lens);
}

/* TODO: this function is pretty wasteful in terms of memory: it allocates
 * lists that are as long as the number of steps in the template, which is
 * always safe but probably also always much longer than necessary */
bool template_to_template_stats(const template *const template, template_stats *stats) {
	int i, j;
	const int len = template->steps_len;
	*stats = (template_stats) { 0, NULL, NULL, NULL };

	if(ALLOC_FAILS(stats->positions   , len) ||
	   ALLOC_FAILS(stats->indexes     , len) ||
	   ALLOC_FAILS(stats->indexes_lens, len))
		goto fail;

	for(i = 0; i < template->steps_len; i++) {
		for(j = 0; j < stats->positions_len; j++)
			if(!strcmp(template->steps[i].position, stats->positions[j]))
				break;

		if(j < stats->positions_len)
			/* found a match: update it */
			stats->indexes[j][stats->indexes_lens[j]++] = i;
		else {
			/* no match: allocate a new position and initialize it */
			if(ALLOC_FAILS(stats->indexes[j], len))
				goto fail;
			stats->positions[j] = template->steps[i].position;
			stats->indexes[j][0] = i;
			stats->indexes_lens[j] = 1;
			stats->positions_len++;
		}
	}

	return true;

fail:
	template_stats_free(*stats);
	return false;
}

int template_stats_to_params(mife_pp_t pp, int i) {
	const template_stats *const stats = pp->mbp_params;
	assert(i < stats->positions_len);
	return stats->indexes_lens[i];
}

void template_stats_to_dimensions(mife_pp_t pp, int *out) {
	const template_stats *const stats    = pp->mbp_params;
	const template       *const template = stats->template;
	int i;
	for(i = 0; i < template->steps_len-1; i++)
		out[i] = template->steps[i].matrix[0].num_cols;
}

/* TODO: efficiency */
void template_stats_to_position(mife_pp_t pp, int global_index, int *out_position, int *out_local_index) {
	const template_stats *const stats = pp->mbp_params;
	int i, j;
	for(i = 0; i < stats->positions_len; i++)
		for(j = 0; j < stats->indexes_lens[i]; j++)
			if(stats->indexes[i][j] == global_index)
				break;
	*out_position = i;
	*out_local_index = j;
}

void template_stats_to_cleartext(mife_pp_t pp, mife_mat_clr_t cleartext, fmpz_t bits) {
	const template_stats *const stats = pp->mbp_params;
	/* TODO */
}

/* TODO: hm! maybe the outputs should be a matrix rather than an array */
/* TODO: are we sure that templates handed to us will produce just one output
 *       every time? */
int template_stats_to_result(mife_pp_t pp, char **raw_result) {
	const template_stats *const stats = pp->mbp_params;
	int i, result, num_nonzeros = 0;
	for(i = 0; i < stats->template->outputs_len; i++)
		if(raw_result[0][i]) {
			result = i;
			num_nonzeros++;
		}
	return (1 == num_nonzeros) ? result : -1;
}
