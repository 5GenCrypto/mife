#ifndef _MIFE_GLUE_H
#include "matrix.h"
#include "mife.h"

/* template_stats structures should not be used after the template they were
 * constructed from are free'd! */
typedef struct {
	int positions_len;
	char **positions;   /* each element references the position string from a step in the template */
	int **indexes;      /* which steps in the template have the given position? */
	int *indexes_lens;  /* how many steps in the template have the given position? */
	template *template; /* a reference to the template in question */
} template_stats;
void template_stats_free(template_stats stats);
bool template_to_template_stats(const template *const template, template_stats *stats);

/* These glue functions assume that the mife_pp_t structures have a pointer to
 * a template_stats structure in their mbp_params. */
int  template_stats_to_params    (mife_pp_t pp, int i);
void template_stats_to_dimensions(mife_pp_t pp, int *out);
void template_stats_to_position  (mife_pp_t pp, int global_index, int *out_position, int *out_local_index);
void template_stats_to_cleartext (mife_pp_t pp, mife_mat_clr_t cleartext, void *cleartext_raw_untyped);
int  template_stats_to_result    (mife_pp_t pp, char **raw_result);
#endif
