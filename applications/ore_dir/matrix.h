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
	unsigned int num_rows, num_cols;
	bool **elems;
} f2_matrix;

typedef struct {
	unsigned int matrices_len;
	f2_matrix *matrices;
} f2_mbp;

/* returns true iff dest is successfully initialized to the contents of src */
bool f2_matrix_copy(f2_matrix *const dest, const f2_matrix src);
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

/* For our purposes, a plaintext is a sequence of symbols and nothing more.
 * Particular applications may wish to provide an application-specific
 * interface for producing these things.
 *
 * For example, consider the case of base-2 ORE, using the compression
 * technique where we draw bits from the left and right arguments in the order
 * "left right right left left right right left ...". In this case, the
 * "semantic plaintext" number 13 would be represented by the "syntactic
 * plaintext":
 *
 *       .-------.--------- 1101 fills in the "right bits"
 *      v       v
 * [1, 11, 10, 01, 1]
 *  ^       ^      ^
 *   `-------`------`------ 1101 fills in the "left bits"
 *
 * One might certainly want a convenient way for the user to specify the
 * semantic plaintext and get the syntactic plaintext out the other end. We do
 * not try to model these "semantic plaintexts" in this type, as that is an
 * application-specific operation; instead, we operate directly on symbol
 * sequences.
 */
typedef struct {
	unsigned int symbols_len;
	char **symbols;
} plaintext;

void plaintext_free(plaintext pt);

/* Choose the matrices from a template that correspond to a particular plaintext.
 * Inputs a template t and plaintext pt.
 * Outputs a matrix branching program mbp.
 * Returns true if all the plaintext symbols were found in the template and the
 * 	mbp is initialized; otherwise leaves the mbp uninitialized and returns
 * 	false.
 */
bool template_instantiate(const template *const t, const plaintext *const pt, f2_mbp *const mbp);

/* A ciphertext mapping assigns a uid to each position, which can be used to
 * locate the particular matrices necessary for evaluation.
 */
typedef struct {
	unsigned int positions_len;
	char **positions;
	char **uids;
} ciphertext_mapping;

void ciphertext_mapping_free(ciphertext_mapping m);

/* returns the uid associated with a given position, or NULL if the position
 * isn't assigned to any uid */
char *uid_from_position(const ciphertext_mapping m, const char *const position);

#endif
