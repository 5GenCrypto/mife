#include <getopt.h>

#include <mife/mife.h>
#include <mmap/mmap_gghlite.h>
#include <mmap/mmap_clt.h>

#include "cmdline.h"
#include "mbp_types.h"
#include "mbp_glue.h"
#include "mife.h"
#include "parse.h"
#include "util.h"

typedef struct {
	mife_pp_t pp;
	location database_location;
	ciphertext_mapping mapping;
} eval_inputs;

static const mbp_template_stats *mbp_template_stats_from_eval_inputs(const eval_inputs ins) { return ins.pp->mbp_params; }
static const mbp_template       *mbp_template_from_eval_inputs      (const eval_inputs ins) { return mbp_template_stats_from_eval_inputs(ins)->template; }

void mife_eval_parse_cmdline(int argc, char **argv, eval_inputs *const ins, bool *use_clt);
f2_matrix mife_eval_evaluate(const_mmap_vtable mmap, const eval_inputs ins);
void mife_eval_print_outputs(const mbp_template t, const f2_matrix m);
void mife_eval_cleanup(const_mmap_vtable mmap, eval_inputs ins, f2_matrix m);

int main(int argc, char **argv) {
	eval_inputs ins;
	f2_matrix m;
	bool success;
    bool use_clt = false;

	mife_eval_parse_cmdline(argc, argv, &ins, &use_clt);

    const mmap_vtable *mmap;
    if (use_clt) {
        mmap = &clt_vtable;
        g_parallel = 1;
    } else {
        mmap = &gghlite_vtable;
    }

	m = mife_eval_evaluate(mmap, ins);
	success = NULL != m.elems;
	if(success) mife_eval_print_outputs(*mbp_template_from_eval_inputs(ins), m);
	mife_eval_cleanup(mmap, ins, m);

	return success ? 0 : -1;
}

static void mife_eval_usage(const int code) {
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
        "  -C, --clt13              Use CLT13 as the underlying multilinear map\n"
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

void mife_eval_parse_cmdline(int argc, char **argv, eval_inputs *const ins, bool *use_clt) {
	unsigned int i;

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
        , {"clt"     ,       no_argument, NULL, 'C'}
		};

	while(!done) {
		int c = getopt_long(argc, argv, "d:hu:C", long_opts, NULL);
		switch(c) {
			case  -1: done = true; break;
			case   0: break; /* a long option with non-NULL flag; should never happen */
			case '?': mife_eval_usage(1); break; /* braking is good defensive driving */
			case 'd':
				location_free(ins->database_location);
				ins->database_location.path = optarg;
				ins->database_location.stack_allocated = true;
				break;
			case 'h': mife_eval_usage(0); break;
			case 'u':
				location_free(public_location);
				public_location = (location) { optarg, true };
				break;
            case 'C':
                *use_clt = true;
                break;
			default:
				fprintf(stderr, "The impossible happened! getopt returned %d (%c)\n", c, c);
				exit(-1);
				break;
		}
	}

    const mmap_vtable *mmap;
    if (*use_clt) {
        mmap = &clt_vtable;
        g_parallel = 1;
    } else {
        mmap = &gghlite_vtable;
    }

	/* read the mapping */
	if(optind != argc-1) {
		fprintf(stderr, "%s: specify exactly one mapping (found %d)\n", *argv, argc-optind);
		mife_eval_usage(2);
	}
	if(!jsmn_parse_ciphertext_mapping_string(argv[optind], &ins->mapping)) {
		fprintf(stderr, "%s: could not parse mapping as JSON object with string values\n", *argv);
		mife_eval_usage(3);
	}

	/* read the template */
	location template_location = location_append(public_location, "template.json");
	mbp_template *template = NULL;
	mbp_template_stats *stats = NULL;
	if(template_location.path == NULL || ALLOC_FAILS(template, 1) || ALLOC_FAILS(stats, 1)) {
		fprintf(stderr, "%s: out of memory while loading template\n", *argv);
		exit(-1);
	}
	if(!jsmn_parse_mbp_template_location(template_location, template)) {
		fprintf(stderr, "%s: could not parse '%s' as a\nJSON representation of a matrix branching program template over the field F_2\n", *argv, template_location.path);
		mife_eval_usage(4);
	}
	if(!mbp_template_to_mife_pp(ins->pp, template, stats)) {
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
	fread_mife_pp(mmap, ins->pp, pp_location.path);
	location_free(pp_location);

	/* sanity check: are all and only the necessary positions specified in the
	 * mapping? */
	if(stats->positions_len != ins->mapping.positions_len) {
		fprintf(stderr, "arity mismatch:\n"
		                "\tfunction template from public parameters expects %u arguments,\n"
		                "\tmapping provided on the command line specifies %u arguments\n"
		              , stats->positions_len, ins->mapping.positions_len);
		mife_eval_usage(5);
	}
	bool position_missing = false;
	for(i = 0; i < stats->positions_len; i++) {
		if(NULL == uid_from_position(ins->mapping, stats->positions[i])) {
			fprintf(stderr, "no mapping for position %s\n", stats->positions[i]);
			position_missing = true;
		}
	}
	if(position_missing) mife_eval_usage(6);
}

static bool mife_eval_load_matrix(const_mmap_vtable mmap, const eval_inputs ins, unsigned int global_index, mmap_enc_mat_t out_m) {
	const mbp_template_stats *const stats    = ins.pp->mbp_params;
	/* const mbp_template       *const template = stats->template; */
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
	fread_mmap_enc_mat(mmap, out_m, file);
	result = true;

	fclose(file);
free_path:
	free(path);
done:
	return result;
}

f2_matrix mife_eval_evaluate(const_mmap_vtable mmap, const eval_inputs ins) {
	f2_matrix result = { .num_rows = 0, .num_cols = 0, .elems = NULL };
	const mbp_template *const template = mbp_template_from_eval_inputs(ins);
	mmap_enc_mat_t product, multiplicand;
	unsigned int i;

	/* if there are no steps to evaluate, I guess we're done */
	if(template->steps_len < 1) goto done;
	if(!mife_eval_load_matrix(mmap, ins, 0, product)) goto done;
	for(i = 1; i < template->steps_len; i++) {
		if(!mife_eval_load_matrix(mmap, ins, i, multiplicand)) goto clear_product;
        mmap_enc_mat_mul(mmap, ins.pp->params_ref, product, product, multiplicand);
		mmap_enc_mat_clear(mmap, multiplicand);
	}

	result = mife_zt_all(mmap, ins.pp, product);
clear_product:
	mmap_enc_mat_clear(mmap, product);
done:
	return result;
}

static unsigned int minui(const unsigned int l, const unsigned int r) {
	return l < r ? l : r;
}

void mife_eval_print_outputs(const mbp_template t, const f2_matrix m) {
	const unsigned int num_rows = minui(m.num_rows, t.outputs.num_rows);
	const unsigned int num_cols = minui(m.num_cols, t.outputs.num_cols);
	for(unsigned int i = 0; i < num_rows; i++) {
		for(unsigned int j = 0; j < num_cols; j++) {
			if(m.elems[i][j]) {
				printf("%s\n", t.outputs.elems[i][j]);
			}
		}
	}
}

void mife_eval_cleanup(const_mmap_vtable mmap, eval_inputs ins, f2_matrix m) {
	mbp_template_stats *stats = ins.pp->mbp_params;
	mbp_template *templ = (mbp_template *)stats->template;
	mbp_template_stats_free(*stats); free(stats);
	mbp_template_free(*templ); free(templ);
	mife_clear_pp_read(mmap, ins.pp);
	location_free(ins.database_location);
	ciphertext_mapping_free(ins.mapping);
	f2_matrix_free(m);
}
