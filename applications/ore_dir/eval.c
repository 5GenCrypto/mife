#include <getopt.h>

#include "cmdline.h"
#include "matrix.h"
#include "mife.h"
#include "mife_glue.h"
#include "parse.h"
#include "util.h"

typedef struct {
	mife_pp_t pp;
	location database_location;
	ciphertext_mapping mapping;
} eval_inputs;

void parse_cmdline(int argc, char **argv, eval_inputs *const ins);
void cleanup(eval_inputs *const ins);

int main(int argc, char **argv) {
	eval_inputs ins;

	parse_cmdline(argc, argv, &ins);
	/* TODO: evaluate */
	cleanup(&ins);

	return 0;
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

void cleanup(eval_inputs *const ins) {
	template_stats *stats = ins->pp->mbp_params;
	template *templ = (template *)stats->template;
	template_stats_free(*stats); free(stats);
	template_free(*templ); free(templ);
	mife_clear_pp_read(ins->pp);
	location_free(ins->database_location);
	ciphertext_mapping_free(ins->mapping);
}
