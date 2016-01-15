#ifndef _MIFE_GLUE_H
#include "matrix.h"
#include "mife.h"

/* the template field is not owned by the template_stats structure, so don't
 * use it after the template it was constructed from is free'd! */
typedef struct {
	int positions_len;   /* what is the arity of the function being encoded? */
	int *position_index; /* for each step, map its position to a number from 0 to positions_len-1 */
	int *local_index;    /* for each step, how many previous steps had the same position? */
	int *step_lens;      /* for each function position, how many steps use that position? */
	const template *template; /* a reference to the template in question */
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
