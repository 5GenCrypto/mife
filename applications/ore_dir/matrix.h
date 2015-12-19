#ifndef _MATRIX_H
#define _MATRIX_H

#include <stdbool.h>

/* A brief note about this comment: there are two meanings of "field"
 * appropriate here; one is the mathematical term for a structured set, and the
 * other is the C term for a part of a struct. We will use "math field" for the
 * former, and "C field" for the latter.
 *
 * Naming convention: matrix branching programs over a given math field are
 * called <math field>_mbp, matrices are called <math field>_matrix, and the C
 * fields are named <math field>_foo. The math fields include (at least) f2 and
 * ggh. The C fields are:
 *
 * * unsigned int num_rows, num_cols: the number of rows and columns in a
 *   matrix
 * * <math field> **elems: a two-dimensional array of the appropriate
 *   dimensions (row index first, then column index)
 * * unsigned int len: how many matrices there are in the program
 * * <math field>_matrix *matrices: the matrices in the program; an invariant
 *   is that `matrices[i].num_cols == matrices[i+1].num_rows`
 */
typedef struct {
	unsigned int f2_num_rows, f2_num_cols;
	bool **f2_elems;
} f2_matrix;

typedef struct {
	unsigned int f2_len;
	f2_matrix *f2_matrices;
} f2_mbp;

void f2_matrix_free(f2_matrix m);
void f2_mbp_free(f2_mbp mbp);

#endif
