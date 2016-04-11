#include <flint/fmpz.h>

#include "mbp_glue.h"
#include "util.h"

void mbp_template_stats_free(mbp_template_stats stats) {
	free(stats.position_index);
	free(stats.local_index);
	free(stats.step_lens);
	free(stats.positions);
}

bool mbp_template_to_mbp_template_stats(const mbp_template *const template, mbp_template_stats *stats) {
	int i, j;
	const int len = template->steps_len;
	*stats = (mbp_template_stats) { 0, NULL, NULL, NULL, NULL, NULL };

	if(ALLOC_FAILS(stats->positions     , len) ||
	   ALLOC_FAILS(stats->position_index, len) ||
	   ALLOC_FAILS(stats->local_index   , len) ||
	   ALLOC_FAILS(stats->step_lens     , len)) {
		mbp_template_stats_free(*stats);
		return false;
	}

	for(i = 0; i < len; i++) {
		for(j = 0; j < stats->positions_len; j++)
			if(!strcmp(template->steps[i].position, stats->positions[j]))
				break;

		stats->step_lens[i] = 0;
		stats->position_index[i] = j;
		stats->local_index[i] = stats->step_lens[j]++;

		if(j == stats->positions_len)
			/* never seen this position before, record it */
			stats->positions[stats->positions_len++] = template->steps[i].position;
	}

	stats->template = template;
	return true;
}

int mbp_template_stats_to_params(mife_pp_t pp, int i) {
	const mbp_template_stats *const stats = pp->mbp_params;
	assert(i < stats->positions_len);
	return stats->step_lens[i];
}

void mbp_template_stats_to_dimensions(mife_pp_t pp, int *out) {
	const mbp_template_stats *const stats    = pp->mbp_params;
	const mbp_template       *const template = stats->template;
	int i;
	for(i = 0; i < template->steps_len-1; i++)
		out[i] = template->steps[i].matrix[0].num_cols;
}

void mbp_template_stats_to_position(mife_pp_t pp, int global_index, int *out_position, int *out_local) {
	const mbp_template_stats *const stats = pp->mbp_params;
	*out_position = stats->position_index[global_index];
	*out_local    = stats->local_index[global_index];
}

/* TODO: it would be good to not call assert */
void mbp_template_stats_to_cleartext(mife_pp_t pp, mife_mat_clr_t cleartext, void *cleartext_raw_untyped) {
	const mbp_template_stats *const stats = pp->mbp_params;
	const mbp_plaintext      *const cleartext_raw = cleartext_raw_untyped;
	f2_mbp mbp;
	int i;

	bool tmp;
	tmp = mbp_template_instantiate(stats->template, cleartext_raw, &mbp);
	assert(tmp);

	/* except these asserts; the asserts in this block are okay */
	/* (but should at least print some diagnostic information) */
	assert(pp->num_inputs == stats->positions_len);
	assert(mbp.matrices_len == stats->template->steps_len);

	if(ALLOC_FAILS(cleartext->clr, pp->num_inputs)) assert(false);
	for(i = 0; i < stats->positions_len; i++)
		if(ALLOC_FAILS(cleartext->clr[i], pp->n[i]))
			assert(false);

	for(i = 0; i < mbp.matrices_len; i++) {
		int j, k;

		/* abbreviations */
		const int position = stats->position_index[i], local = stats->local_index[i];
		const f2_matrix f2_m = mbp.matrices[i];
		fmpz_mat_struct *fmpz_m = cleartext->clr[position][local];

		fmpz_mat_init(fmpz_m, f2_m.num_rows, f2_m.num_cols);
		for(j = 0; j < f2_m.num_rows; j++)
			for(k = 0; k < f2_m.num_cols; k++)
				fmpz_set_ui(fmpz_mat_entry(fmpz_m, j, k), f2_m.elems[j][k]);
	}

	f2_mbp_free(mbp);
}

/* TODO: are we sure that templates handed to us will produce just one output
 *       every time? */
int mbp_template_stats_to_result(mife_pp_t pp, f2_matrix raw_result) {
	int i, j, result, num_nonzeros = 0;
	fprintf(stderr, "The impossible happened! Deprecated mbp_glue function mbp_template_stats_to_result with insufficiently expressive type was called. We'll do our best, but something is wrong.\n");
	for(i = 0; i < raw_result.num_rows; i++)
		for(j = 0; j < raw_result.num_cols; j++)
			if(raw_result.elems[i][j]) {
				result = i*raw_result.num_cols + j;
				num_nonzeros++;
			}
	return (1 == num_nonzeros) ? result : -1;
}
