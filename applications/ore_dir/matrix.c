#include <stdlib.h>
#include "matrix.h"

void f2_matrix_free(f2_matrix m) {
	int i;
	for(i = 0; i < m.f2_num_rows; i++) free(m.f2_elems[i]);
	free(m.f2_elems);
}

void f2_mbp_free(f2_mbp mbp) {
	int i;
	for(i = 0; i < mbp.f2_len; i++) f2_matrix_free(mbp.f2_matrices[i]);
	free(mbp.f2_matrices);
}
