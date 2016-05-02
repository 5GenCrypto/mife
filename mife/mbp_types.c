#include <stdlib.h>
#include <string.h>
#include "mbp_types.h"
#include "util.h"

bool f2_matrix_copy(f2_matrix *const dest, const f2_matrix src) {
  if(ALLOC_FAILS(dest->elems, src.num_rows))
    return false;
  int *i = &dest->num_rows; /* to shorten some lines */
  int j;
  for(*i = 0; *i < src.num_rows; (*i)++) {
    if(ALLOC_FAILS(dest->elems[*i], src.num_cols)) {
      f2_matrix_free(*dest);
      return false;
    }

    for(j = 0; j < src.num_cols; j++)
      dest->elems[*i][j] = src.elems[*i][j];
  }
  dest->num_cols = src.num_cols;
  return true;
}

bool f2_matrix_zero(f2_matrix *const dest, const unsigned int num_rows, const unsigned int num_cols) {
  if(ALLOC_FAILS(dest->elems, num_rows))
    return false;
  int j, *i = &dest->num_rows;
  for(*i = 0; *i < num_rows; (*i)++) {
    if(ALLOC_FAILS(dest->elems[*i], num_cols)) {
      f2_matrix_free(*dest);
      return false;
    }

    for(j = 0; j < num_cols; j++)
      dest->elems[*i][j] = false;
  }
  dest->num_cols = num_cols;
  return true;
}

void f2_matrix_free(f2_matrix m) {
  int i;
  if(NULL != m.elems) {
    for(i = 0; i < m.num_rows; i++)
      free(m.elems[i]);
    free(m.elems);
  }
}

void f2_mbp_free(f2_mbp mbp) {
	int i;
	if(NULL != mbp.matrices) {
		for(i = 0; i < mbp.matrices_len; i++)
			f2_matrix_free(mbp.matrices[i]);
		free(mbp.matrices);
	}
}

void string_matrix_free(string_matrix m) {
	int i, j;
	if(NULL != m.elems) {
		for(i = 0; i < m.num_rows; i++) {
			if(NULL != m.elems[i]) {
				for(j = 0; j < m.num_cols; j++)
					free(m.elems[i][j]);
				free(m.elems[i]);
			}
		}
	free(m.elems);
	}
}

void mbp_step_free(mbp_step s) {
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

void mbp_template_free(mbp_template t) {
	int i;
	if(NULL != t.steps) {
		for(i = 0; i < t.steps_len; i++)
			mbp_step_free(t.steps[i]);
		free(t.steps);
	}
	string_matrix_free(t.outputs);
}

void mbp_plaintext_free(mbp_plaintext pt) {
	int i;
	if(NULL != pt.symbols) {
		for(i = 0; i < pt.symbols_len; i++)
			free(pt.symbols[i]);
		free(pt.symbols);
	}
}

bool mbp_template_instantiate(const mbp_template *const t, const mbp_plaintext *const pt, f2_mbp *const mbp) {
	if(t->steps_len != pt->symbols_len) return false;
	if(ALLOC_FAILS(mbp->matrices, t->steps_len))
		return false;

	int *i = &mbp->matrices_len;
	for(*i = 0; *i < t->steps_len; (*i)++) {
		int j;
		bool found = false;
		for(j = 0; j < t->steps[*i].symbols_len; j++) {
			if(!strcmp(pt->symbols[*i], t->steps[*i].symbols[j])) {
				if(found) /* found a second hit! bail out */ {
					f2_matrix_free(mbp->matrices[*i]);
					f2_mbp_free(*mbp);
					return false;
				}
				found = true;

				if(!f2_matrix_copy(mbp->matrices + *i, t->steps[*i].matrix[j])) {
					f2_mbp_free(*mbp);
					return false;
				}
			}
		}
		if(!found) {
			f2_mbp_free(*mbp);
			return false;
		}
	}
	return true;
}

void ciphertext_mapping_free(ciphertext_mapping m) {
	unsigned int i;

	if(NULL != m.positions)
		for(i = 0; i < m.positions_len; i++)
			free(m.positions[i]);
	free(m.positions);

	if(NULL != m.uids)
		for(i = 0; i < m.positions_len; i++)
			free(m.uids[i]);
	free(m.uids);
}

char *uid_from_position(const ciphertext_mapping m, const char *const position) {
	unsigned int i;
	for(i = 0; i < m.positions_len; i++)
		if(!strcmp(m.positions[i], position))
			return m.uids[i];
	return NULL;
}
