#include <stdlib.h>
#include "matrix.h"

void f2_matrix_free(f2_matrix m) {
	int i;
	if(NULL != m.f2_elems) {
		for(i = 0; i < m.f2_num_rows; i++)
			free(m.f2_elems[i]);
		free(m.f2_elems);
	}
}

void f2_mbp_free(f2_mbp mbp) {
	int i;
	if(NULL != mbp.f2_matrices) {
		for(i = 0; i < mbp.f2_matrices_len; i++)
			f2_matrix_free(mbp.f2_matrices[i]);
		free(mbp.f2_matrices);
	}
}

void step_free(step s) {
	int i;
	if(NULL != s.symbols) {
		for(i = 0; i < s.symbols_len; i++)
			if(NULL != s.symbols[i])
				free(s.symbols[i]);
		free(s.symbols);
	}
	if(NULL != s.matrix) {
		for(i = 0; i < s.symbols_len; i++)
			f2_matrix_free(s.matrix[i]);
		free(s.matrix);
	}
	if(NULL != s.position) free(s.position);
}

void template_free(template t) {
	int i;
	if(NULL != t.steps) {
		for(i = 0; i < t.steps_len; i++)
			step_free(t.steps[i]);
		free(t.steps);
	}
	if(NULL != t.outputs) {
		for(i = 0; i < t.outputs_len; i++)
			if(NULL != t.outputs[i])
				free(t.outputs[i]);
		free(t.outputs);
	}
}
