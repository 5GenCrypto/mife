#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

#include "parse.h"

void location_free(location loc) { if(!loc.stack_allocated) free(loc.path); }

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
	mbp->f2_len = (*json_tokens)->size;
	if((*json_tokens)->type != JSMN_ARRAY) {
		fprintf(stderr, "at position %d\nexpecting matrix branching program, found non-array\n", (*json_tokens)->start);
		return false;
	}
	if(mbp->f2_len < 1) return false;

	/* reserve some memory */
	if(NULL == (mbp->f2_matrices = malloc(mbp->f2_len * sizeof(*mbp->f2_matrices)))) {
		fprintf(stderr, "out of memory in jsmn_parse_f2_mbp\n");
		return false;
	}

	/* parse each matrix */
	int i;
	for(i = 0; i < mbp->f2_len; i++) {
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
	for(i = 0; i < mbp->f2_len-1; i++) {
		if(mbp->f2_matrices[i].f2_num_cols != mbp->f2_matrices[i+1].f2_num_rows) {
			fprintf(stderr, "matrices %d and %d have incompatible shapes\n", i, i+1);
			f2_mbp_free(*mbp);
			return false;
		}
	}

	return true;
}

bool jsmn_parse_f2_mbp_location(const location loc, f2_mbp *const mbp) {
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
		fprintf(stderr, "out of memory when allocation tokens in jsmn_parse_f2_mbp_location\n");
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
		fprintf(stderr, "out of memory when allocating state in jsmn_parse_f2_mbp_location\n");
		goto fail_free_tokens;
	}
	*state = json_tokens;
	status = jsmn_parse_f2_mbp(json_string, state, mbp);

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
