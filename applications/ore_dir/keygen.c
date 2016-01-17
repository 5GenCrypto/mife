#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "cmdline.h"
#include "mife_glue.h"
#include "parse.h"
#include "util.h"

typedef struct {
	int sec_param, log_db_size;
	aes_randstate_t seed;
	template template;
} keygen_inputs;

typedef struct {
	location public, private;
} keygen_locations;

void parse_cmdline(int argc, char **argv, keygen_inputs *const ins, keygen_locations *const outs);
void print_outputs(keygen_locations outs, mife_pp_t pp, mife_sk_t sk);
void cleanup(keygen_inputs *const ins, keygen_locations *const outs);

const gghlite_flag_t ggh_flags = GGHLITE_FLAGS_QUIET | GGHLITE_FLAGS_GOOD_G_INV;

int main(int argc, char **argv) {
	keygen_inputs ins;
	keygen_locations outs;
	mife_pp_t pp;
	mife_sk_t sk;
	template_stats stats;

	parse_cmdline(argc, argv, &ins, &outs);
	if(!template_to_mife_pp(pp, &ins.template, &stats)) return -1;
	mife_setup(pp, sk, ins.log_db_size, ins.sec_param, ggh_flags, ins.seed);
	print_outputs(outs, pp, sk);

	cleanup(&ins, &outs);
	template_stats_free(stats);
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
	if(ALLOC_FAILS(loc->path, prefix_len+sec_param_len+1)) {
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

bool read_seed(location loc, char *dest) {
	FILE *src = fopen(loc.path, "rb");
	if(NULL == (src = fopen(loc.path      , "rb")) &&
	   NULL == (src = fopen("/dev/urandom", "rb")))
		return false;
	dest[SEED_SIZE] = '\0';
	return fread(dest, sizeof(*dest), SEED_SIZE, src) == SEED_SIZE;
}

void usage(const int code) {
	/* separate the diagnostic information from the usage information a little bit */
	if(0 != code) printf("\n\n");
	printf(
		"The key generation phase initializes the parameters that are shared across\n"
		"an entire database. Brackets indicate default values for each argument.\n"
		"\n"
		"Common options:\n"
		"  -h, --help         Display this usage information\n"
		"  -r, --private      A directory for private parameters [private-$secparam]\n"
		"  -u, --public       A directory for public parameters [public-$secparam]\n"
		"\n"
		"Keygen-specific options:\n"
		"  -s, --secparam     Security parameter [80]\n"
		"  -n, --dbsize       Allow up to 2^n records [80]\n"
		"\n"
		"Files used:\n"
		"  <public>/template.json  R  JSON    a description of the function being\n"
		"                                     encrypted\n"
		"  <public>/mife.pub        W custom  public parameters for evaluating\n"
		"  <private>/mife.priv      W custom  private parameters for encrypting\n"
		"  <private>/seed.bin      R  binary  32-byte seed for PRNG\n"
		"  /dev/urandom            R  binary  used in case above file is missing\n"
		);
	exit(code);
}

void parse_cmdline(int argc, char **argv, keygen_inputs *const ins, keygen_locations *const outs) {
	bool done = false;
	char seed[SEED_SIZE+1];

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
			case   0: break; /* a long option with non-NULL flag; should never happen */
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

	/* read seed */
	location seed_location = location_append(outs->private, "seed.bin");
	if(NULL == seed_location.path) {
		fprintf(stderr, "%s: out of memory when trying to create path to seed\n", *argv);
		exit(-1);
	}
	if(!read_seed(seed_location, seed)) {
		fprintf(stderr, "%s: could not read %d bytes of seed\n", *argv, SEED_SIZE);
		usage(8);
	}
	location_free(seed_location);
	aes_randinit_seed(ins->seed, seed, "keygen");
}

/* for now, use the mife library's custom format; would be good to upgrade this
 * to something a bit more redundant and human-readable to improve error
 * detection, error reporting, versioning and just generally make this tool a
 * bit more robust
 */
void print_outputs(keygen_locations outs, mife_pp_t pp, mife_sk_t sk) {
	/* the public directory has to exist -- we read template.json out of it! --
	 * but the private one might not yet */
	int err = mkdir(outs.private.path, S_IRWXU);
	if(0 != err && EEXIST != errno) {
		fprintf(stderr, "could not create output directory %s\n", outs.private.path);
		perror("print_outputs");
		return;
	}

	location  public_location = location_append(outs.public , "mife.pub" );
	location private_location = location_append(outs.private, "mife.priv");
	if( public_location.path == NULL ||
	   private_location.path == NULL) {
		location_free(public_location);
		fprintf(stderr, "out of memory when generating output paths\n");
		return;
	}

	fwrite_mife_pp(pp,  public_location.path);
	fwrite_mife_sk(sk, private_location.path);

	location_free( public_location);
	location_free(private_location);
}

void cleanup(keygen_inputs *const ins, keygen_locations *const outs) {
	location_free(outs-> public);
	location_free(outs->private);
	aes_randclear(ins->seed);
}
