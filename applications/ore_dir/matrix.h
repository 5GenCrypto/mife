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
 * * unsigned int matrices_len: how many matrices there are in the program
 * * <math field>_matrix *matrices: the matrices in the program; an invariant
 *   is that `matrices[i].num_cols == matrices[i+1].num_rows`
 */
typedef struct {
	unsigned int f2_num_rows, f2_num_cols;
	bool **f2_elems;
} f2_matrix;

typedef struct {
	unsigned int f2_matrices_len;
	f2_matrix *f2_matrices;
} f2_mbp;

void f2_matrix_free(f2_matrix m);
void f2_mbp_free(f2_mbp mbp);

/* A template describes how to choose an MBP from a plaintext. A template has:
 * * a sequence of steps
 * * a description of the result of the computation
 *
 * Each step of the computation can be thought of as taking a step in a finite
 * state machine. In the language of finite state machines, a step has:
 * * a collection of input symbols that can be seen at that step
 * * a state transition for each symbol, represented as the adjacency matrix
 *
 * For our purposes, we will also want to track which "position" each symbol of
 * the input comes from -- that is, when encoding a function f(x,y,z), each
 * symbol of the input will have one position out of "x", "y", or "z".
 */
typedef struct {
	unsigned int symbols_len;
	char *position;
	char **symbols;
	f2_matrix *matrix;
} step;

typedef struct {
	unsigned int steps_len, outputs_len;
	step *steps;
	char **outputs;
} template;

void step_free(step s);
void template_free(template t);

#endif
