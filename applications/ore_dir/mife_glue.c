#include <fmpz.h>

#include "mife_glue.h"
#include "util.h"

void template_stats_free(template_stats stats) {
	free(stats.position_index);
	free(stats.local_index);
	free(stats.step_lens);
}

bool template_to_template_stats(const template *const template, template_stats *stats) {
	int i, j;
	const int len = template->steps_len;
	char **positions = NULL;
	*stats = (template_stats) { 0, NULL, NULL, NULL, NULL };

	if(ALLOC_FAILS(positions            , len) ||
	   ALLOC_FAILS(stats->position_index, len) ||
	   ALLOC_FAILS(stats->local_index   , len) ||
	   ALLOC_FAILS(stats->step_lens     , len))
		goto fail;

	for(i = 0; i < len; i++) {
		for(j = 0; j < stats->positions_len; j++)
			if(!strcmp(template->steps[i].position, positions[j]))
				break;

		stats->step_lens[i] = 0;
		stats->position_index[i] = j;
		stats->local_index[i] = stats->step_lens[j]++;

		if(j == stats->positions_len)
			/* never seen this position before, record it */
			positions[stats->positions_len++] = template->steps[i].position;
	}

	stats->template = template;
	free(positions);
	return true;

fail:
	free(positions);
	template_stats_free(*stats);
	return false;
}

int template_stats_to_params(mife_pp_t pp, int i) {
	const template_stats *const stats = pp->mbp_params;
	assert(i < stats->positions_len);
	return stats->step_lens[i];
}

void template_stats_to_dimensions(mife_pp_t pp, int *out) {
	const template_stats *const stats    = pp->mbp_params;
	const template       *const template = stats->template;
	int i;
	for(i = 0; i < template->steps_len-1; i++)
		out[i] = template->steps[i].matrix[0].num_cols;
}

void template_stats_to_position(mife_pp_t pp, int global_index, int *out_position, int *out_local) {
	const template_stats *const stats = pp->mbp_params;
	*out_position = stats->position_index[global_index];
	*out_local    = stats->local_index[global_index];
}

/* TODO: it would be good to not call assert */
void template_stats_to_cleartext(mife_pp_t pp, mife_mat_clr_t cleartext, void *cleartext_raw_untyped) {
	const template_stats *const stats = pp->mbp_params;
	const plaintext      *const cleartext_raw = cleartext_raw_untyped;
	f2_mbp mbp;
	int i;

	assert(template_instantiate(stats->template, cleartext_raw, &mbp));

	/* except these asserts; the asserts in this block are okay */
	assert(pp->num_inputs == stats->positions_len);
	assert(mbp.matrices_len == stats->template->steps_len);

	assert(!ALLOC_FAILS(cleartext->clr, pp->num_inputs));
	for(i = 0; i < stats->positions_len; i++)
		assert(!ALLOC_FAILS(cleartext->clr[i], pp->n[i]));

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
