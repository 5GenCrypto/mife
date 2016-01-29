#include <getopt.h>

#include "cmdline.h"
#include "types.h"
#include "mife.h"
#include "mife_glue.h"
#include "parse.h"
#include "util.h"

typedef struct {
	mife_pp_t pp;
	location database_location;
	ciphertext_mapping mapping;
} eval_inputs;

const template_stats *template_stats_from_eval_inputs(const eval_inputs ins) { return ins.pp->mbp_params; }
const template       *template_from_eval_inputs      (const eval_inputs ins) { return template_stats_from_eval_inputs(ins)->template; }

void parse_cmdline(int argc, char **argv, eval_inputs *const ins);
f2_matrix evaluate(const eval_inputs ins);
void print_outputs(const template t, const f2_matrix m);
void cleanup(eval_inputs ins, f2_matrix m);

int main(int argc, char **argv) {
	eval_inputs ins;
	f2_matrix m;
	bool success;

	parse_cmdline(argc, argv, &ins);
	m = evaluate(ins);
	success = NULL != m.elems;
	if(success) print_outputs(*template_from_eval_inputs(ins), m);
	cleanup(ins, m);

	return success ? 0 : -1;
}

void usage(const int code) {
	/* separate the diagnostic information from the usage information a little bit */
	if(0 != code) printf("\n\n");
	printf(
		"USAGE: eval [OPTIONS] MAPPING\n"
		"The evaluation operation plugs individual database records into each position\n"
		"of the function being computed. The mapping should be a JSON object whose\n"
		"field names should be positions from the public function template, and whose\n"
		"values should be strings giving uids of database records.\n"
		"\n"
		"Prints one line for each non-zero in the result matrix. The line will contain\n"
		"the appropriate string from the `outputs` field of the function template.\n"
		"\n"
		"Brackets indicate default values for each argument.\n"
		"\n"
		"Common options:\n"
		"  -h, --help               Display this usage information\n"
		"  -u, --public             A directory for public parameters [public]\n"
		"  -d, --db, --database     A directory to store encrypted values in [database]\n"
		"\n"
		"Evaluation-specific options: (none)\n"
		"\n"
		"Files used:\n"
		"  <database>/<uid>/<position>/*.bin\n"
		"                            R  binary  the encrypted records\n"
		"  <public>/template.json    R  JSON    a description of the function being\n"
		"                                       evaluated\n"
		"  <public>/mife.pub         R  custom  public parameters for evaluating\n"
		);
	exit(code);
}

void parse_cmdline(int argc, char **argv, eval_inputs *const ins) {
	int i;

	/* set defaults */
	location public_location = { "public", true };
	/* since ins->database_location gets returned to the caller, we can't
	 * allocate it on our stack */
	if(ALLOC_FAILS(ins->database_location.path, strlen("database")+1)) {
		fprintf(stderr, "%s: out of memory when setting default database location\n", *argv);
		exit(-1);
	}
	strcpy(ins->database_location.path, "database");
	ins->database_location.stack_allocated = false;

	bool done = false;
	struct option long_opts[] =
		{ {"db"      , required_argument, NULL, 'd'}
		, {"database", required_argument, NULL, 'd'}
		, {"help"    ,       no_argument, NULL, 'h'}
		, {"public"  , required_argument, NULL, 'u'}
		};

	while(!done) {
		int c = getopt_long(argc, argv, "d:hu:", long_opts, NULL);
		switch(c) {
			case  -1: done = true; break;
			case   0: break; /* a long option with non-NULL flag; should never happen */
			case '?': usage(1); break; /* braking is good defensive driving */
			case 'd':
				location_free(ins->database_location);
				ins->database_location.path = optarg;
				ins->database_location.stack_allocated = true;
				break;
			case 'h': usage(0); break;
			case 'u':
				location_free(public_location);
				public_location = (location) { optarg, true };
				break;
			default:
				fprintf(stderr, "The impossible happened! getopt returned %d (%c)\n", c, c);
				exit(-1);
				break;
		}
	}

	/* read the mapping */
	if(optind != argc-1) {
		fprintf(stderr, "%s: specify exactly one mapping (found %d)\n", *argv, argc-optind);
		usage(2);
	}
	if(!jsmn_parse_ciphertext_mapping_string(argv[optind], &ins->mapping)) {
		fprintf(stderr, "%s: could not parse mapping as JSON object with string values\n", *argv);
		usage(3);
	}

	/* read the template */
	location template_location = location_append(public_location, "template.json");
	template *template = NULL;
	template_stats *stats = NULL;
	if(template_location.path == NULL || ALLOC_FAILS(template, 1) || ALLOC_FAILS(stats, 1)) {
		fprintf(stderr, "%s: out of memory while loading template\n", *argv);
		exit(-1);
	}
	if(!jsmn_parse_template_location(template_location, template)) {
		fprintf(stderr, "%s: could not parse '%s' as a\nJSON representation of a matrix branching program template over the field F_2\n", *argv, template_location.path);
		usage(4);
	}
	if(!template_to_mife_pp(ins->pp, template, stats)) {
		fprintf(stderr, "%s: internal error while computing statistics for template\n", *argv);
		exit(-1);
	}
	location_free(template_location);

	/* read the public parameters */
	location pp_location = location_append(public_location, "mife.pub");
	if(pp_location.path == NULL) {
		fprintf(stderr, "%s: out of memory while loading public parameters\n", *argv);
		exit(-1);
	}
	/* TODO: some error-checking would be nice here */
	fread_mife_pp(ins->pp, pp_location.path);
	location_free(pp_location);

	/* sanity check: are all and only the necessary positions specified in the
	 * mapping? */
	if(stats->positions_len != ins->mapping.positions_len) {
		fprintf(stderr, "arity mismatch:\n"
		                "\tfunction template from public parameters expects %u arguments,\n"
		                "\tmapping provided on the command line specifies %u arguments\n"
		              , stats->positions_len, ins->mapping.positions_len);
		usage(5);
	}
	bool position_missing = false;
	for(i = 0; i < stats->positions_len; i++) {
		if(NULL == uid_from_position(ins->mapping, stats->positions[i])) {
			fprintf(stderr, "no mapping for position %s\n", stats->positions[i]);
			position_missing = true;
		}
	}
	if(position_missing) usage(6);
}

bool load_matrix(const eval_inputs ins, unsigned int global_index, gghlite_enc_mat_t out_m) {
	const template_stats *const stats    = ins.pp->mbp_params;
	const template       *const template = stats->template;
	const int local_index      = stats->local_index[global_index];
	const char *const position = stats->positions[stats->position_index[global_index]];
	const char *const uid      = uid_from_position(ins.mapping, position);
	const char *const db_path  = ins.database_location.path;
	const int path_len  = snprintf(NULL, 0, "%s/%s/%s/%d.bin", db_path, uid, position, local_index);
	const int path_size = path_len+1;
	bool result = false;
	char *path;

	if(ALLOC_FAILS(path, path_size)) {
		fprintf(stderr, "out of memory trying to construct path:\n\t%s/%s/%s/%d.bin\n", db_path, uid, position, local_index);
		goto done;
	}

	int tmp = snprintf(path, path_size, "%s/%s/%s/%d.bin", db_path, uid, position, local_index);
	if(tmp != path_len) {
		fprintf(stderr, "The impossible happened: snprintf produced strings of two different lengths on two calls with all the same arguments.\n\t(%d first time, %d second)\n", path_len, tmp);
		goto free_path;
	}

	FILE *file = fopen(path, "rb");
	if(NULL == file) {
		fprintf(stderr, "Could not open ciphertext chunk at location\n\t%s\n", path);
		goto free_path;
	}

	/* TODO: would be nice to have some error checking here */
	fread_gghlite_enc_mat(ins.pp, out_m, file);
	result = true;

fclose:
	fclose(file);
free_path:
	free(path);
done:
	return result;
}

f2_matrix evaluate(const eval_inputs ins) {
	f2_matrix result = { .num_rows = 0, .num_cols = 0, .elems = NULL };
	const template *const template = template_from_eval_inputs(ins);
	gghlite_enc_mat_t product, multiplicand;
	int i;

	/* if there are no steps to evaluate, I guess we're done */
	if(template->steps_len < 1) goto done;
	if(!load_matrix(ins, 0, product)) goto done;
	for(i = 1; i < template->steps_len; i++) {
		if(!load_matrix(ins, i, multiplicand)) goto clear_product;
		gghlite_enc_mat_mul(*ins.pp->params_ref, product, product, multiplicand);
		gghlite_enc_mat_clear(multiplicand);
	}

	result = mife_zt_all(ins.pp, product);
clear_product:
	gghlite_enc_mat_clear(product);
done:
	return result;
}

unsigned int minui(const unsigned int l, const unsigned int r) {
	return l < r ? l : r;
}

void print_outputs(const template t, const f2_matrix m) {
	const unsigned int num_rows = minui(m.num_rows, t.outputs.num_rows);
	const unsigned int num_cols = minui(m.num_cols, t.outputs.num_cols);
	int i, j;
	for(i = 0; i < num_rows; i++)
		for(j = 0; j < num_cols; j++)
			if(m.elems[i][j])
				printf("%s\n", t.outputs.elems[i][j]);
}

void cleanup(eval_inputs ins, f2_matrix m) {
	template_stats *stats = ins.pp->mbp_params;
	template *templ = (template *)stats->template;
	template_stats_free(*stats); free(stats);
	template_free(*templ); free(templ);
	mife_clear_pp_read(ins.pp);
	location_free(ins.database_location);
	ciphertext_mapping_free(ins.mapping);
	f2_matrix_free(m);
}
