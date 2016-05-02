#ifndef _MBP_GLUE_H
#define _MBP_GLUE_H

#include "mbp_types.h"
#include "mife.h"

/* the template field and elements of the positions field are not owned by the
 * mbp_template_stats structure, so don't use them after the template it was
 * constructed from is free'd! */
typedef struct {
	unsigned int positions_len;         /* what is the arity of the function being encoded? */
	int *position_index;                /* for each step, map its position to a number from 0 to positions_len-1 */
	int *local_index;                   /* for each step, how many previous steps had the same position? */
	int *step_lens;                     /* for each function position, how many steps use that position? */
	const char **positions;             /* a name for each position */
	const mbp_template *template;       /* a reference to the template in question */
} mbp_template_stats;
void mbp_template_stats_free(mbp_template_stats stats);
bool mbp_template_to_mbp_template_stats(const mbp_template *const template, mbp_template_stats *stats);

/* These glue functions assume that the mife_pp_t structures have a pointer to
 * a mbp_template_stats structure in their mbp_params. */
int  mbp_template_stats_to_params    (mife_pp_t pp, int i);
void mbp_template_stats_to_dimensions(mife_pp_t pp, int *out);
void mbp_template_stats_to_position  (mife_pp_t pp, int global_index, int *out_position, int *out_local_index);
void mbp_template_stats_to_cleartext (mife_pp_t pp, mife_mat_clr_t cleartext, void *cleartext_raw_untyped);
int  mbp_template_stats_to_result    (mife_pp_t pp, f2_matrix raw_result);
#endif /* ifndef _MBP_GLUE_H */
