#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "cmdline.h"
#include "matrix.h"
#include "mife.h"
#include "mife_glue.h"
#include "parse.h"
#include "util.h"

#define UID_TRIES 100
#define INT_STR_LEN (3 * sizeof(int))

typedef struct {
	plaintext pt;
	fmpz_t uid;
	mife_sk_t sk;
	mife_pp_t pp;
	aes_randstate_t seed;
} encrypt_inputs;

void parse_cmdline(int argc, char **argv, encrypt_inputs *const ins, location *const database_location);
bool print_output(mife_pp_t pp, int global_index, gghlite_enc_mat_t ct, location database_location, fmpz_t uid);
void cleanup(encrypt_inputs *const ins, location *const database_location);

int main(int argc, char **argv) {
	encrypt_inputs ins;
	location database_location;
	mife_mat_clr_t clr;
	int i, ***partitions;
	bool success = true;

	parse_cmdline(argc, argv, &ins, &database_location);
	mife_encrypt_setup(ins.pp, ins.uid, &ins.pt, clr, &partitions);
	for(i = 0; i < ((template_stats *)ins.pp->mbp_params)->template->steps_len; i++) {
		gghlite_enc_mat_t ct;
		mife_encrypt_single(ins.pp, ins.sk, ins.seed, i, clr, partitions, ct);
		success &= print_output(ins.pp, i, ct, database_location, ins.uid);
		gghlite_enc_mat_clear(ct);
	}
	mife_encrypt_cleanup(ins.pp, clr, partitions);

	cleanup(&ins, &database_location);
	return success ? 0 : -1;
}

void usage(const int code) {
	/* separate the diagnostic information from the usage information a little bit */
	if(0 != code) printf("\n\n");
	printf(
		"USAGE: encrypt [OPTIONS] PLAINTEXT\n"
		"The encryption operation hides the information in a single plaintext. The\n"
		"plaintext is should be represented as a JSON array containing strings naming\n"
		"symbols from the template available in the public parameters directory.\n"
		"Brackets indicate default values for each argument.\n"
		"\n"
		"Common options:\n"
		"  -h, --help               Display this usage information\n"
		"  -r, --private            A directory for private parameters [private]\n"
		"  -u, --public             A directory for public parameters [public]\n"
		"  -d, --db, --database     A directory to store encrypted values in [database]\n"
		"\n"
		"Encryption-specific options:\n"
		"  -i, --uid                A number that uniquely identifies this record and\n"
		"                           is at most the database size specified during\n"
		"                           the keygen phase [chosen at random]\n"
		"\n"
		"Files used:\n"
		"  <database>/<uid>/*/*.bin   W binary  the encrypted record\n"
		"  <public>/template.json    R  JSON    a description of the function being\n"
		"                                       encrypted\n"
		"  <public>/mife.pub         R  custom  public parameters for evaluating\n"
		"  <private>/mife.priv       R  custom  private parameters for encrypting\n"
		"  <private>/seed.bin        R  binary  32-byte seed for PRNG\n"
		"  /dev/urandom              R  binary  used in case above file is missing\n"
		);
	exit(code);
}

/* reads (num_bits/8 + 1)*8 bytes into n, mod 2^num_bits */
bool fmpz_read_bits(fmpz_t n, FILE *file, int num_bits) {
	fmpz_zero(n);
	while(num_bits >= 8) {
		int c = fgetc(file);
		if(EOF == c) return false;
		fmpz_mul_ui(n, n, 256);
		fmpz_add_ui(n, n, c);
		num_bits -= 8;
	}
	if(num_bits >= 0) {
		int c = fgetc(file);
		if(EOF == c) return false;
		c = c & (1 << num_bits);
		fmpz_mul_ui(n, n, 1 << num_bits);
		fmpz_add_ui(n, n, c);
	}
	return true;
}

void parse_cmdline(int argc, char **argv, encrypt_inputs *const ins, location *const database_location) {
	int i, j;
	bool done = false;
	bool have_uid = false;

	/* set defaults */
	location public_location = {  "public", true },
	        private_location = { "private", true };
	/* since database_location gets returned to the caller, we can't allocate
	 * it on our stack */
	if(ALLOC_FAILS(database_location->path, strlen("database")+1)) {
		fprintf(stderr, "%s: out of memory when setting default database location\n", *argv);
		exit(-1);
	}
	strcpy(database_location->path, "database");
	database_location->stack_allocated = false;

	struct option long_opts[] =
		{ {"db"      , required_argument, NULL, 'd'}
		, {"database", required_argument, NULL, 'd'}
		, {"help"    ,       no_argument, NULL, 'h'}
		, {"uid"     , required_argument, NULL, 'i'}
		, {"private" , required_argument, NULL, 'r'}
		, {"public"  , required_argument, NULL, 'u'}
		};

	while(!done) {
		int c = getopt_long(argc, argv, "d:hi:r:u:", long_opts, NULL);
		switch(c) {
			case  -1: done = true; break;
			case   0: break; /* a long option with non-NULL flag; should never happen */
			case '?': usage(1); break; /* braking is good defensive driving */
			case 'd':
				location_free(*database_location);
				database_location->path = optarg;
				database_location->stack_allocated = true;
				break;
			case 'h': usage(0); break;
			case 'i':
				if(0 != fmpz_set_str(ins->uid, optarg, 10)) {
					fprintf(stderr, "%s: could not read %s as a base-10 number\n", *argv, optarg);
					usage(2);
				}
				have_uid = true;
				break;
			case 'r':
				location_free(private_location);
				private_location = (location) { optarg, true };
				break;
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

	/* read the plaintext */
	if(optind != argc-1) {
		fprintf(stderr, "%s: specify exactly one plaintext (found %d)\n", *argv, argc-optind);
		usage(2);
	}
	if(!jsmn_parse_plaintext_string(argv[optind], &ins->pt)) {
		fprintf(stderr, "%s: could not parse plaintext as JSON array of strings\n", *argv);
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

	/* read the secret key */
	location sk_location = location_append(private_location, "mife.priv");
	if(sk_location.path == NULL) {
		fprintf(stderr, "%s: out of memory while loading private key\n", *argv);
		exit(-1);
	}
	/* TODO: error-checking */
	fread_mife_sk(ins->sk, sk_location.path);
	location_free(sk_location);

	/* initialize uid if !have_uid */
	/* it is not important that the uid be cryptographically random, so just
	 * use /dev/urandom */
	if(!have_uid) {
		char *record_path;
		const size_t offset = strlen(database_location->path) + 1;
		/* ought to be able to fit 3 bits per decimal character */
		const size_t len = offset + ins->pp->L / 3 + 1;
		FILE *urandom = fopen("/dev/urandom", "rb");
		struct stat file_stat;
		if(NULL == urandom) {
			fprintf(stderr, "no uid specified and could not open /dev/urandom to choose one\n");
			exit(-1);
		}
		if(ALLOC_FAILS(record_path, len)) {
			fclose(urandom);
			fprintf(stderr, "out of memory while finding unused uid\n");
			exit(-1);
		}

		memcpy(record_path, database_location->path, offset-1);
		record_path[offset-1] = '/';

		for(i = 0; i < UID_TRIES; i++) {
			fmpz_init(ins->uid);
			fmpz_read_bits(ins->uid, urandom, ins->pp->L - 1);
			fmpz_get_str(record_path + offset, 10, ins->uid);
			/* probably not totally foolproof, but is a decent quick check that
			 * a file doesn't exist
			 */
			if(stat(record_path, &file_stat)) break;
			fmpz_clear(ins->uid);
		}

		fclose(urandom);

		if(UID_TRIES == i) {
			fprintf(stderr, "tried %d random uids, but they were all in use\n", UID_TRIES);
			exit(-1);
		}
	}

	/* initialize the random seed */
	const char function_name[] = "encrypt";
	size_t context_len = fmpz_sizeinbase(ins->uid, 62) + sizeof(function_name);
	char *context;
	if(ALLOC_FAILS(context, context_len)) {
		fprintf(stderr, "%s: out of memory when generating context for RNG seed\n", *argv);
		exit(-1);
	}
	memcpy(context, function_name, sizeof(function_name)-1);
	fmpz_get_str(context + sizeof(function_name)-1, 62, ins->uid);
	check_parse_result(load_seed(private_location, context, ins->seed), usage, 5);

	/* from here on out, it's all sanity-checks */

	/* check that the uid is in range */
	if(fmpz_sizeinbase(ins->uid, 2) >= (size_t)ins->pp->L) {
		fprintf(stderr, "uid uses %zu bits, but the current key supports only up to %d bits\n", fmpz_sizeinbase(ins->uid, 2), ins->pp->L);
		usage(6);
	}

	/* check that the template and plaintext have the same length */
	if(template->steps_len != ins->pt.symbols_len) {
		fprintf(stderr, "the number of symbols in the plaintext (%d)\ndoes not match the number of steps in the template (%d)\n",
		        ins->pt.symbols_len, template->steps_len);
		usage(7);
	}

	/* check that the template and plaintext match up appropriately */
	bool match_everywhere = true;
	for(i = 0; i < template->steps_len; i++) {
		bool match_here = false;
		for(j = 0; j < template->steps[i].symbols_len; j++)
			match_here |= !strcmp(template->steps[i].symbols[j], ins->pt.symbols[i]);
		if(!match_here) {
			fprintf(stderr, "the plaintext symbol %s at index %d is unknown\n", ins->pt.symbols[i], i);
			fprintf(stderr, "\t(known symbols: ");
			for(j = 0; j < template->steps[i].symbols_len-1; j++)
				fprintf(stderr, "%s, ", template->steps[i].symbols[j]);
			fprintf(stderr, "%s)\n", template->steps[i].symbols[j]);
		}
		match_everywhere &= match_here;
	}
	if(!match_everywhere) usage(8);

	/* TODO: check that `ins->pp` and `stats` match up */
	/* TODO: check that the secret key is appropriately dimensioned */

	/* these do nothing for now, and are only here as a defensive measure
	 * against future refactorings */
	location_free(public_location);
	location_free(private_location);
}

FILE *fopen_bin(const location database_location, const char *const uid_str, const char *const position, const int local_index) {
	char *dir, *path;
	const size_t  dir_len = strlen(database_location.path) + 1 + strlen(uid_str) + 1 + strlen(position), dir_size = dir_len + 1;
	const size_t path_len = dir_len + 1 + INT_STR_LEN + 4, path_size = path_len + 1;
	FILE *result = NULL;

	if(ALLOC_FAILS( dir,  dir_size)) goto done;
	if(ALLOC_FAILS(path, path_size)) goto free_dir;
	assert(snprintf( dir,  dir_size, "%s/%s/%s" , database_location.path, uid_str, position) < dir_size);
	assert(snprintf(path, path_size, "%s/%d.bin", dir, local_index) < path_size);
	if(!create_directory_if_missing(dir)) goto free_path;
	result = fopen(path, "wb");

free_path:
	free(path);
free_dir:
	free(dir);
done:
	return result;
}

bool print_output(mife_pp_t pp, int global_index, gghlite_enc_mat_t ct, location database_location, fmpz_t uid) {
	bool success = false;
	const template_stats *const stats    = pp->mbp_params;
	const template       *const template = stats->template;
	      char           *const uid_str  = fmpz_get_str(NULL, 10, uid);

	if(NULL == uid_str) {
		fprintf(stderr, "could not produce string representing uid to build output directory name\n");
		goto done;
	}

	success = true;
	const int local_index = stats->local_index[global_index];
	const char *const position = template->steps[global_index].position;
	FILE *dest = fopen_bin(database_location, uid_str, position, local_index);
	if(NULL == dest) {
		fprintf(stderr, "could not write %s/%d.bin\n", position, local_index);
		success = false;
	}
	fwrite_gghlite_enc_mat(/* TODO: why is pp needed? */ pp, ct, dest);
	fclose(dest);

free_uid_str:
	free(uid_str);
done:
	return success;
}

void cleanup(encrypt_inputs *const ins, location *const database_location) {
	plaintext_free(ins->pt);
	fmpz_clear(ins->uid);
	mife_clear_sk(ins->sk);
	template_stats *stats = ins->pp->mbp_params;
	template *templ = (template *)stats->template;
	template_stats_free(*stats); free(stats);
	template_free(*templ); free(templ);
	mife_clear_pp_read(ins->pp);
	aes_randclear(ins->seed);
	location_free(*database_location);
}
