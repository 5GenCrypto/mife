#include <getopt.h>
#include <gghlite/gghlite.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "parse.h"

typedef struct {
	int sec_param, log_db_size;
	template template;
} keygen_inputs;

typedef struct {
	location public, private;
} keygen_locations;

void parse_cmdline(int argc, char **argv, keygen_inputs *const ins, keygen_locations *const outs);
void cleanup(const keygen_inputs *const ins, const keygen_locations *const outs);

int main(int argc, char **argv) {
	keygen_inputs ins;
	keygen_locations outs;
	parse_cmdline(argc, argv, &ins, &outs);

	fmpz_t order;
	gghlite_params_t public;
	fmpz_mat_t kilian[] = {};
	fmpz_mat_t kilian_inverse[] = {};

	/* TODO: do something interesting with this data */
	printf("security parameter %d\ndb size %d\npublic output to %s\nprivate output to %s\n"
	      , ins.sec_param, ins.log_db_size, outs.public.path, outs.private.path
	      );

	/* TODO: write to the output locations */

	cleanup(&ins, &outs);
	return 0;
}

void location_init(keygen_inputs *const ins, location *const loc, const char *const prefix, const int code) {
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

void usage(const int code) {
	printf(
		"The key generation phase initializes the parameters that are shared across\n"
		"an entire database. Brackets indicate default values for each argument.\n"
		"\n"
		"Common options:\n"
		"  -h, --help         Display this usage information\n"
		"  -r, --private      An output directory for private parameters [private-$secparam]\n"
		"  -u, --public       An input/output directory for public parameters [public-$secparam]\n"
		"\n"
		"Keygen-specific options:\n"
		"  -s, --secparam     Security parameter [80]\n"
		"  -n, --dbsize       Allow up to 2^n records [80]\n"
		);
	exit(code);
}

void parse_cmdline(int argc, char **argv, keygen_inputs *const ins, keygen_locations *const outs) {
	bool done = false;

	/* set defaults; the NULLs in outs will be overwritten */
	ins->sec_param = 80;
	ins->log_db_size = 80;
	*outs = (keygen_locations) { { NULL, false}, { NULL, false } };

	struct option long_opts[] =
		{ {"help"     ,       no_argument, NULL, 'h'}
		, {"dbsize"   , required_argument, NULL, 'n'}
		, {"private"  , required_argument, NULL, 'r'}
		, {"secparam" , required_argument, NULL, 's'}
		, {"public"   , required_argument, NULL, 'u'}
		, {}
		};

	while(!done) {
		int c = getopt_long(argc, argv, "hn:r:s:u:", long_opts, NULL);
		switch(c) {
			case  -1: done = true; break;
			case   0: break;
			case '?': usage(1); break; /* braking is good defensive driving */
			case 'h': usage(0); break;
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

	/* read template */
	location template_location = location_append(outs->public, "template.json");
	if(NULL == template_location.path) {
		fprintf(stderr, "%s: out of memory when trying to create path to template\n", *argv);
		exit(-1);
	}
	if(!jsmn_parse_template_location(template_location, &ins->template)) {
		fprintf(stderr, "%s: could not parse '%s' as a\nJSON representation of a matrix branching program template over the field F_2\n", *argv, template_location.path);
		usage(7);
	}
	location_free(template_location);
}

void cleanup(const keygen_inputs *const ins, const keygen_locations *const outs) {
	location_free(outs-> public);
	location_free(outs->private);
}
