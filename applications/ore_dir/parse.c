#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

#include "parse.h"

void location_free(location loc) { if(!loc.stack_allocated && NULL != loc.path) free(loc.path); }

location location_append(const location loc, const char *const path) {
	location result = { NULL, false };
	const unsigned int loc_len = strlen(loc.path), path_len = strlen(path);
	const unsigned int result_len = loc_len+1+path_len, result_size = result_len+1;
	if(NULL != (result.path = malloc(result_size * sizeof(*result.path))))
		if(snprintf(result.path, result_size, "%s/%s", loc.path, path) != result_len) {
			free(result.path);
			result.path = NULL;
		}
	return result;
}

bool jsmn_parse_f2_elem(const char *const json_string, const jsmntok_t **const json_tokens, bool *const elem) {
	if((*json_tokens)->type != JSMN_PRIMITIVE) return false;

	const int pos = (*json_tokens)->start;
	const char c = json_string[pos];

	/* for numbers, die if there are multiple digits */
	if((c == '0' || c == '1') && (*json_tokens)->end - pos > 1) {
		fprintf(stderr, "at position %d\nexpecting 0 or 1, but found a longer number\n", (*json_tokens)->start);
		return false;
	}

	switch(c) {
		case '0':
		case 'f':
			*elem = false;
			break;
		case '1':
		case 't':
			*elem = true;
			break;
		default:
			fprintf(stderr, "at position %d\nexpecting 0, 1, true, or false, found token starting with '%c'\n", (*json_tokens)->start, json_string[(*json_tokens)->start]);
			return false;
	}
	return true;
}

/* The next three functions are almost identical. I miss polymorphism. =( */
int jsmn_parse_f2_row(const char *const json_string, const jsmntok_t **const json_tokens, bool **const row, int expected_num_cols) {
	/* demand an array of elements of the appropriate length */
	if((*json_tokens)->type != JSMN_ARRAY) {
		fprintf(stderr, "at position %d\nexpecting matrix row, found non-array\n", (*json_tokens)->start);
		return -1;
	}
	if(expected_num_cols < 0) expected_num_cols = (*json_tokens)->size;
	if(expected_num_cols != (*json_tokens)->size) {
		fprintf(stderr, "at position %d\nexpecting row of length %d, found row of length %d\n",
		        (*json_tokens)->start, expected_num_cols, (*json_tokens)->size);
		return -1;
	}

	/* reserve some memory */
	if(NULL == (*row = malloc(expected_num_cols * sizeof(bool)))) {
		fprintf(stderr, "out of memory in jsmn_parse_f2_row\n");
		return -1;
	}

	/* parse each element */
	int i;
	for(i = 0; i < expected_num_cols; i++) {
		++(*json_tokens);
		if(!jsmn_parse_f2_elem(json_string, json_tokens, *row+i)) {
			free(*row);
			return -1;
		}
	}

	return expected_num_cols;
}

bool jsmn_parse_f2_matrix(const char *const json_string, const jsmntok_t **const json_tokens, f2_matrix *const matrix) {
	/* demand an array of rows */
	matrix->f2_num_rows = (*json_tokens)->size;
	if((*json_tokens)->type != JSMN_ARRAY) {
		fprintf(stderr, "at position %d\nexpecting matrix, found non-array\n", (*json_tokens)->start);
		return false;
	}

	/* reserve some memory */
	if(NULL == (matrix->f2_elems = malloc(matrix->f2_num_rows * sizeof(*matrix->f2_elems)))) {
		fprintf(stderr, "out of memory in jsmn_parse_f2_matrix\n");
		return false;
	}

	/* parse each row */
	int num_cols = -1;
	int i;
	for(i = 0; i < matrix->f2_num_rows; i++) {
		++(*json_tokens);
		if(0 > (num_cols = jsmn_parse_f2_row(json_string, json_tokens, matrix->f2_elems+i, num_cols))) {
			int j;
			for(j = 0; j < i; j++)
				free(matrix->f2_elems[j]);
			free(matrix->f2_elems);
			return false;
		}
	}
	matrix->f2_num_cols = num_cols;

	return true;
}

bool jsmn_parse_f2_mbp(const char *const json_string, const jsmntok_t **const json_tokens, f2_mbp *const mbp) {
	/* demand an array of at least one matrix */
	mbp->f2_matrices_len = (*json_tokens)->size;
	if((*json_tokens)->type != JSMN_ARRAY) {
		fprintf(stderr, "at position %d\nexpecting matrix branching program, found non-array\n", (*json_tokens)->start);
		return false;
	}
	if(mbp->f2_matrices_len < 1) return false;

	/* reserve some memory */
	if(NULL == (mbp->f2_matrices = malloc(mbp->f2_matrices_len * sizeof(*mbp->f2_matrices)))) {
		fprintf(stderr, "out of memory in jsmn_parse_f2_mbp\n");
		return false;
	}

	/* parse each matrix */
	int i;
	for(i = 0; i < mbp->f2_matrices_len; i++) {
		++(*json_tokens);
		if(!jsmn_parse_f2_matrix(json_string, json_tokens, mbp->f2_matrices+i)) {
			int j;
			for(j = 0; j < i; j++)
				f2_matrix_free(mbp->f2_matrices[j]);
			free(mbp->f2_matrices);
			return false;
		}
	}

	/* check that the dimensions match up */
	for(i = 0; i < mbp->f2_matrices_len-1; i++) {
		if(mbp->f2_matrices[i].f2_num_cols != mbp->f2_matrices[i+1].f2_num_rows) {
			fprintf(stderr, "matrices %d and %d have incompatible shapes\n", i, i+1);
			f2_mbp_free(*mbp);
			return false;
		}
	}

	return true;
}

bool jsmn_parse_step(const char *const json_string, const jsmntok_t **const json_tokens, step *const step) {
	int i;

	/* demand an object with at least two parts */
	step->symbols_len = (*json_tokens)->size - 1;
	if((*json_tokens)->type != JSMN_OBJECT || step->symbols_len < 1) {
		fprintf(stderr, "at position %d\nexpecting step (JSON object with key \"position\" and one key per input symbol), found non-object\n", (*json_tokens)->start);
		return false;
	}

	/* reserve some memory; reserve one extra slot in case the user forgot to
	 * specify a position */
	step->position = NULL;
	step->symbols  = NULL;
	step->matrix   = NULL;
	if(NULL == (step->symbols = malloc((step->symbols_len+1) * sizeof(*step->symbols)))) {
		fprintf(stderr, "out of memory in jsmn_parse_step\n");
		return false;
	}
	for(i = 0; i < step->symbols_len; i++) step->symbols[i] = NULL;
	if(NULL == (step->matrix = malloc((step->symbols_len+1) * sizeof(*step->matrix)))) {
		fprintf(stderr, "out of memory in jsmn_parse_step\n");
		step_free(*step);
		return false;
	}

	/* parse each key-value pair */
	int next_symbol = 0;
	for(i = 0; i < step->symbols_len+1; i++) {
		++(*json_tokens);
		char *key;
		if(!jsmn_parse_string(json_string, json_tokens, &key)) {
			step_free(*step);
			return false;
		}

		++(*json_tokens);
		if(strcmp(key, "position") == 0) {
			free(key);
			if(step->position != NULL) {
				fprintf(stderr, "at position %d\nduplicate position information\n", (*json_tokens)->start);
				step_free(*step);
				return false;
			}
			if(!jsmn_parse_string(json_string, json_tokens, &step->position)) {
				step_free(*step);
				return false;
			}
		}
		else {
			step->symbols[next_symbol] = key;
			if(!jsmn_parse_f2_matrix(json_string, json_tokens, step->matrix + next_symbol)) {
				free(key);
				step->symbols[next_symbol] = NULL;
				step_free(*step);
				return false;
			}
			++next_symbol;
		}
	}

	/* make sure we saw a "position" key */
	if(next_symbol == step->symbols_len+1) {
		fprintf(stderr, "at position %d\nnever saw a \"position\" key in this object\n", (*json_tokens)->start);
		++step->symbols_len;
		step_free(*step);
		return false;
	}

	/* check that all the matrices have the same size */
	for(i = 1; i < step->symbols_len; i++) {
		if(step->matrix[i].f2_num_rows != step->matrix[0].f2_num_rows ||
		   step->matrix[i].f2_num_cols != step->matrix[0].f2_num_cols) {
			fprintf(stderr, "before position %d\ndimension mismatch in matrices 0 and %d\n", (*json_tokens)->end, i);
			step_free(*step);
			return false;
		}
	}

	return true;
}

bool jsmn_parse_steps(const char *const json_string, const jsmntok_t **const json_tokens, template *const template) {
	/* demand an array with at least one step in it */
	template->steps_len = (*json_tokens)->size;
	if((*json_tokens)->type != JSMN_ARRAY || template->steps_len <= 1) {
		fprintf(stderr, "at position %d\nexpecting array of steps\n", (*json_tokens)->start);
		return false;
	}

	/* reserve some space */
	if(NULL == (template->steps = malloc(template->steps_len * sizeof(*template->steps)))) {
		fprintf(stderr, "out of memory in jsmn_parse_steps\n");
		return false;
	}

	/* parse each step */
	int i;
	for(i = 0; i < template->steps_len; i++) {
		++(*json_tokens);
		if(!jsmn_parse_step(json_string, json_tokens, template->steps+i)) {
			/* caller is responsible for cleaning up; template_free will
			 * inspect steps_len to decide how many steps to free */
			template->steps_len = i;
			return false;
		}
	}

	/* check that the steps' dimensions mesh */
	for(i = 1; i < template->steps_len; i++) {
		if(template->steps[i-1].matrix[0].f2_num_cols != template->steps[i].matrix[0].f2_num_rows) {
			fprintf(stderr, "before position %d\ndimension mismatch in steps %d and %d\n", (*json_tokens)->end, i-1, i);
			return false;
		}
	}

	return true;
}

bool jsmn_parse_string(const char *const json_string, const jsmntok_t **const json_tokens, char **const string) {
	/* demand a string */
	if((*json_tokens)->type != JSMN_STRING) {
		fprintf(stderr, "at position %d\nexpecting string\n", (*json_tokens)->start);
		return false;
	}

	/* reserve some space */
	const unsigned int string_len = (*json_tokens)->end - (*json_tokens)->start;
	if(NULL == (*string = malloc((string_len+1) * sizeof(**string)))) {
		fprintf(stderr, "out of memory in jsmn_parse_string\n");
		return false;
	}

	/* copy */
	memcpy(*string, json_string+(*json_tokens)->start, string_len);
	(*string)[string_len] = '\0';
	return true;
}

bool jsmn_parse_outputs(const char *const json_string, const jsmntok_t **const json_tokens, template *const template) {
	/* demand an array */
	template->outputs_len = (*json_tokens)->size;
	if((*json_tokens)->type != JSMN_ARRAY) {
		fprintf(stderr, "at position %d\nexpecting array of output strings\n", (*json_tokens)->start);
		return false;
	}

	/* reserve some space */
	if(NULL == (template->outputs = malloc(template->outputs_len * sizeof(*template->outputs)))) {
		fprintf(stderr, "out of memory in jsmn_parse_outputs\n");
		return false;
	}

	/* parse each output */
	int i;
	for(i = 0; i < template->outputs_len; i++) {
		++(*json_tokens);
		if(!jsmn_parse_string(json_string, json_tokens, template->outputs+i)) {
			/* caller is responsible for cleaning up; template_free will
			 * inspect outputs_len to decide how many outputs to free */
			template->outputs_len = i;
			return false;
		}
	}

	return true;
}

bool jsmn_parse_template(const char *const json_string, const jsmntok_t **const json_tokens, template *const template) {
	int i;

	/* demand an object with exactly two parts */
	if((*json_tokens)->type != JSMN_OBJECT || (*json_tokens)->size != 2) {
		fprintf(stderr, "at position %d\nexpecting template (JSON object with keys \"steps\" and \"outputs\"), found non-object\n", (*json_tokens)->start);
		return false;
	}

	template->steps   = NULL;
	template->outputs = NULL;

	for(i = 0; i < 2; i++) {
		char *key;
		++(*json_tokens);
		switch(json_string[(*json_tokens)->start]) {
			case 's': /* steps */
				++(*json_tokens);
				if(!jsmn_parse_steps(json_string, json_tokens, template)) {
					template_free(*template);
					return false;
				}
				break;

			case 'o': /* outputs */
				++(*json_tokens);
				if(!jsmn_parse_outputs(json_string, json_tokens, template)) {
					template_free(*template);
					return false;
				}
				break;

			default:
				if(jsmn_parse_string(json_string, json_tokens, &key)) {
					fprintf(stderr, "at position %d\nexpecting \"steps\" or \"outputs\", found %s\n", (*json_tokens)->start, key);
					free(key);
				}
				template_free(*template);
				return false;
		}
	}

	if(NULL == template->steps) {
		fprintf(stderr, "before position %d\nmissing key \"steps\"\n", (*json_tokens)->end);
		template_free(*template);
		return false;
	}
	if(NULL == template->outputs) {
		fprintf(stderr, "before position %d\nmissing key \"outputs\"\n", (*json_tokens)->end);
		template_free(*template);
		return false;
	}
	if(template->steps[template->steps_len-1].matrix[0].f2_num_cols != template->outputs_len) {
		fprintf(stderr, "before position %d\ndimension mismatch between final step and outputs\n", (*json_tokens)->end);
		template_free(*template);
		return false;
	}

	return true;
}

bool jsmn_parse_template_location(const location loc, template *const template) {
	int fd;
	struct stat fd_stat;
	char *json_string;
	jsmn_parser parser;
	jsmntok_t *json_tokens;
	bool status = false;

	/* read plaintext */
	if(-1 == (fd = open(loc.path, O_RDONLY, 0))) {
		fprintf(stderr, "could not open matrix branching program '%s'\n", loc.path);
		goto fail_none;
	}
	if(-1 == fstat(fd, &fd_stat)) {
		fprintf(stderr, "could not stat matrix branching program '%s'\n", loc.path);
		goto fail_close;
	}
	if(NULL == (json_string = mmap(NULL, fd_stat.st_size, PROT_READ, MAP_PRIVATE, fd, 0))) {
		fprintf(stderr, "mmap failed, out of address space?\n");
		goto fail_close;
	}

	/* first pass: figure out how much stuff is in the string */
	jsmn_init(&parser);
	int json_tokens_len = jsmn_parse(&parser, json_string, fd_stat.st_size, NULL, 0);
	if(json_tokens_len < 1) {
		fprintf(stderr, "found no JSON tokens\n");
		goto fail_unmap;
	}
	if(NULL == (json_tokens = malloc(json_tokens_len * sizeof(jsmntok_t)))) {
		fprintf(stderr, "out of memory when allocation tokens in jsmn_parse_template_location\n");
		goto fail_unmap;
	}

	/* second pass: record the results of parsing */
	jsmn_init(&parser);
	int tokens_parsed = jsmn_parse(&parser, json_string, fd_stat.st_size, json_tokens, json_tokens_len);
	if(tokens_parsed != json_tokens_len) {
		fprintf(stderr, "The impossible happened: parsed the same string twice and got two\ndifferent token counts (%d first time, %d second).\n", json_tokens_len, tokens_parsed);
		goto fail_free_tokens;
	}

	/* convert the parsed JSON to our custom type */
	const jsmntok_t **const state = malloc(sizeof(*state));
	if(NULL == state) {
		fprintf(stderr, "out of memory when allocating state in jsmn_parse_template_location\n");
		goto fail_free_tokens;
	}
	*state = json_tokens;
	status = jsmn_parse_template(json_string, state, template);

fail_free_state:
	free(state);
fail_free_tokens:
	free(json_tokens);
fail_unmap:
	munmap(json_string, fd_stat.st_size);
fail_close:
	close(fd);
fail_none:
	return status;
}
