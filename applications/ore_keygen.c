#include <fcntl.h>
#include <getopt.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
/* TODO: work out correct way to include jsmn in this project */
#include "jsmn.h"

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

typedef struct {
	int sec_param, log_db_size;
	f2_mbp plaintext;
} keygen_inputs;

typedef struct {
	char *path;
	bool stack_allocated;
} location;

typedef struct {
	location public, private;
} keygen_output_locations;

/* TODO: is `foo * const x` the right thing? */
void parse_cmdline(int argc, char **argv, keygen_inputs *ins, keygen_output_locations *outs);
void cleanup(keygen_inputs * const ins, keygen_output_locations * const outs);

int main(int argc, char **argv) {
	keygen_inputs ins;
	keygen_output_locations outs;
	parse_cmdline(argc, argv, &ins, &outs);

	/* TODO: do something interesting with this data */
	printf("security parameter %d\ndb size %d\npublic output to %s\nprivate output to %s\n"
	      , ins.sec_param, ins.log_db_size, outs.public.path, outs.private.path
	      );
	printf("matrix branching program\n");
	int i, j, k;
	for(i = 0; i < ins.plaintext.f2_len; i++) {
		const unsigned int r = ins.plaintext.f2_matrices[i].f2_num_rows;
		const unsigned int c = ins.plaintext.f2_matrices[i].f2_num_cols;
		printf("[\n");
		for(j = 0; j < r; j++) {
			printf("[");
			for(k = 0; k < c-1; k++)
				printf("%d, ", ins.plaintext.f2_matrices[i].f2_elems[j][k]);
			if(c > 0)
				printf("%d"  , ins.plaintext.f2_matrices[i].f2_elems[j][c-1]);
			if(j < r-1) printf("],\n");
			else        printf("]\n");
		}
		if(i < ins.plaintext.f2_len-1) printf("],\n");
		else                           printf("]\n");
	}

	/* TODO: write to the output locations */

	cleanup(&ins, &outs);
	return 0;
}

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

bool jsmn_parse_f2_elem(char * const json_string, jsmntok_t **json_tokens, bool *elem) {
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
int jsmn_parse_f2_row(char * const json_string, jsmntok_t **json_tokens, int expected_num_cols, bool **row) {
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

bool jsmn_parse_f2_matrix(char * const json_string, jsmntok_t **json_tokens, f2_matrix *matrix) {
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
		if(0 > (num_cols = jsmn_parse_f2_row(json_string, json_tokens, num_cols, matrix->f2_elems+i))) {
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

bool jsmn_parse_f2_mbp(char * const json_string, jsmntok_t **json_tokens, f2_mbp *mbp) {
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

bool jsmn_parse_f2_mbp_string(char * const json_string, const int json_string_len, f2_mbp *mbp) {
	jsmn_parser parser;
	jsmntok_t *json_tokens;

	/* first pass: figure out how much stuff is in the string */
	jsmn_init(&parser);
	int json_tokens_len = jsmn_parse(&parser, json_string, json_string_len, NULL, 0);
	if(json_tokens_len < 1) return false;
	if(NULL == (json_tokens = malloc(json_tokens_len * sizeof(jsmntok_t)))) return false;

	/* second pass: record the results of parsing */
	jsmn_init(&parser);
	int tokens_parsed = jsmn_parse(&parser, json_string, json_string_len, json_tokens, json_tokens_len);
	if(tokens_parsed != json_tokens_len) {
		fprintf(stderr, "The impossible happened: parsed the same string twice and got two\ndifferent token counts (%d first time, %d second).\n", json_tokens_len, tokens_parsed);
		return false;
	}

	/* convert the parsed JSON to our custom type */
	jsmntok_t **state;
	if(NULL == (state = malloc(sizeof(*state)))) {
		free(json_tokens);
		return false;
	}
	*state = json_tokens;

	bool result = jsmn_parse_f2_mbp(json_string, state, mbp);
	free(json_tokens);
	return result;
}

void location_free(location loc) { if(!loc.stack_allocated) free(loc.path); }

void location_init(keygen_inputs * const ins, location *loc, char * const prefix, int code) {
	/* three characters per byte of number should be a safe
	 * overapproximation of how much space is needed; note we use
	 * sizeof(int) instead of sizeof(ins->sec_param) because of the implicit
	 * cast that may happen during the call to snprintf if sec_param's type
	 * changes during a later refactor
	 */
	const int prefix_len = strlen(prefix), sec_param_len = 3*sizeof(int);
	loc->path = malloc((prefix_len+sec_param_len+1)*sizeof(char));
	if(NULL == loc->path) {
		fprintf(stderr, "Out of memory\n");
		exit(code);
	}
	memcpy(loc->path, prefix, prefix_len);
	if(snprintf(loc->path+prefix_len, sec_param_len+1, "%d", ins->sec_param) > sec_param_len+1) {
		fprintf(stderr, "The impossible happened: %d bytes was not enough to store %d.\n", sec_param_len+1, ins->sec_param);
		exit(code);
	}
	loc->stack_allocated = false;
}

void usage(int code) {
	printf(
		"The key generation phase initializes the parameters that are shared across\n"
		"an entire database. Brackets indicate default values for each argument.\n"
		"\n"
		"Common options:\n"
		"  -h, --help         Display this usage information\n"
		"  -l, --plaintext    An input file containing a sample plaintext [plaintext.json]\n"
		"  -r, --private      An output directory for private parameters [private-$secparam]\n"
		"  -u, --public       An output directory for public parameters [public-$secparam]\n"
		"\n"
		"Keygen-specific options:\n"
		"  -s, --secparam     Security parameter [80]\n"
		"  -n, --dbsize       Allow up to 2^n records [80]\n"
		);
	exit(code);
}

void parse_cmdline(int argc, char **argv, keygen_inputs *ins, keygen_output_locations *outs) {
	int opt_index = 0, error = 0;
	bool done = false;

	/* set defaults; the NULLs in outs will be overwritten */
	location plaintext_location = (location) { "plaintext.json", true };
	ins->sec_param = 80;
	ins->log_db_size = 80;
	*outs = (keygen_output_locations) { { NULL, false}, { NULL, false } };

	struct option long_opts[] =
		{ {"help"     ,       no_argument, NULL, 'h'}
		, {"plaintext", required_argument, NULL, 'l'}
		, {"dbsize"   , required_argument, NULL, 'n'}
		, {"private"  , required_argument, NULL, 'r'}
		, {"secparam" , required_argument, NULL, 's'}
		, {"public"   , required_argument, NULL, 'u'}
		, {}
		};

	while(!done) {
		int c = getopt_long(argc, argv, "hl:n:r:s:u:", long_opts, &opt_index);
		switch(c) {
			case  -1: done = true; break;
			case   0: break;
			case '?': usage(1); break; /* braking is good defensive driving */
			case 'h': usage(0); break;
			case 'l':
				plaintext_location.path = optarg;
				break;
			case 'n':
				if((ins->log_db_size = atoi(optarg)) < 1) {
					fprintf(stderr, "%s: unparseable database size '%s', should be positive number\n", *argv, optarg);
					usage(2);
				}
				break;
			case 'r':
				outs->private.path = optarg;
				outs->private.stack_allocated = true;
				break;
			case 's':
				if((ins->sec_param = atoi(optarg)) < 1) {
					fprintf(stderr, "%s: unparseable security parameter '%s', should be positive number\n", *argv, optarg);
					usage(3);
				}
				break;
			case 'u':
				outs->public.path = optarg;
				outs->public.stack_allocated = true;
				break;
			default:
				fprintf(stderr, "The impossible happened! getopt returned %d (%c)\n", c, c);
				exit(-1);
				break;
		}
	}

	if(optind < argc) {
		fprintf(stderr, "%s: unexpected non-option argument %s\n", *argv, argv[optind]);
		exit(4);
	}

	if(NULL == outs-> public.path) location_init(ins, &outs-> public,  "public-", 5);
	if(NULL == outs->private.path) location_init(ins, &outs->private, "private-", 6);

	/* read plaintext */
	int fd = open(plaintext_location.path, O_RDONLY, 0);
	if(-1 == fd) {
		fprintf(stderr, "%s: couldn't open plaintext file '%s'\n", *argv, plaintext_location.path);
		usage(7);
	}
	struct stat fd_stat;
	if(-1 == fstat(fd, &fd_stat)) {
		fprintf(stderr, "%s: couldn't stat plaintext file '%s'\n", *argv, plaintext_location.path);
		exit(8);
	}
	char *plaintext_contents = mmap(NULL, fd_stat.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if(NULL == plaintext_contents) {
		fprintf(stderr, "mmap failed, out of address space?\n");
		exit(9);
	}
	if(!jsmn_parse_f2_mbp_string(plaintext_contents, fd_stat.st_size, &ins->plaintext)) {
		fprintf(stderr, "%s: could not parse '%s' as a\nJSON representation of a matrix branching program over the field F_2\n", *argv, plaintext_location.path);
		usage(10);
	}
	location_free(plaintext_location);
}

void cleanup(keygen_inputs * const ins, keygen_output_locations * const outs) {
	location_free(outs-> public);
	location_free(outs->private);
	f2_mbp_free(ins->plaintext);
}
