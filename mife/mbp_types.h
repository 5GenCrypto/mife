#ifndef _MBP_TYPES_H
#define _MBP_TYPES_H

#include <stdbool.h>

typedef struct {
  unsigned int num_rows, num_cols;
  bool **elems;
} f2_matrix;

/* returns true iff dest is successfully initialized to the contents of src */
bool f2_matrix_copy(f2_matrix *const dest, const f2_matrix src);
bool f2_matrix_zero(f2_matrix *const dest, const unsigned int num_rows, const unsigned int num_cols);
void f2_matrix_free(f2_matrix m);

/* an invariant is that `matrices[i].num_cols == matrices[i+1].num_rows` */
typedef struct {
	unsigned int matrices_len;
	f2_matrix *matrices;
} f2_mbp;

void f2_mbp_free(f2_mbp mbp);

typedef struct {
	unsigned int num_rows, num_cols;
	char ***elems;
} string_matrix;

void string_matrix_free(string_matrix m);

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
} mbp_step;

typedef struct {
	unsigned int steps_len;
	mbp_step *steps;
	string_matrix outputs;
} mbp_template;

void mbp_step_free(mbp_step s);
void mbp_template_free(mbp_template t);

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
} mbp_plaintext;

void mbp_plaintext_free(mbp_plaintext pt);

/* Choose the matrices from a template that correspond to a particular plaintext.
 * Inputs a template t and plaintext pt.
 * Outputs a matrix branching program mbp.
 * Returns true if all the plaintext symbols were found in the template and the
 * 	mbp is initialized; otherwise leaves the mbp uninitialized and returns
 * 	false.
 */
bool mbp_template_instantiate(const mbp_template *const t, const mbp_plaintext *const pt, f2_mbp *const mbp);

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

#endif /* ifndef _MBP_TYPES_H */
