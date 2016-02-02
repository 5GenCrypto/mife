#include "cmdline.h"

const mife_flag_t mife_flags = MIFE_DEFAULT;

bool template_to_mife_pp(mife_pp_t pp, const template *const template, template_stats *const stats) {
	if(!template_to_template_stats(template, stats)) return false;
	mife_init_params(pp, mife_flags);
	mife_mbp_set(stats, pp, stats->positions_len,
		template_stats_to_params,
		template_stats_to_dimensions,
		template_stats_to_position,
		template_stats_to_cleartext,
		template_stats_to_result);
	return true;
}

parse_result load_seed(location private_location, char *context, aes_randstate_t seed) {
	FILE *src;
	char dest[SEED_SIZE];
	parse_result result;
	location seed_location = location_append(private_location, "seed.bin");

	if(NULL == seed_location.path) {
		fprintf(stderr, "out of memory when trying to create path to seed\n");
		result = PARSE_OUT_OF_MEMORY;
		goto fail_none;
	}
	if(NULL == (src = fopen(seed_location.path, "rb")) &&
	   NULL == (src = fopen("/dev/urandom"    , "rb"))) {
		fprintf(stderr, "could not open seed file for reading; attempted to read:\n");
		fprintf(stderr, "\t%s\n\t/dev/urandom\n", seed_location.path);
		result = PARSE_IO_ERROR;
		goto fail_free_location;
	}

	if(fread(dest, sizeof(*dest), SEED_SIZE, src) != SEED_SIZE) {
		fprintf(stderr, "could not read %d bytes of seed\n", SEED_SIZE);
		result = PARSE_INVALID;
		goto fail_close_src;
	}

	/* TODO: error checking */
	aes_randinit_seedn(seed, dest, SEED_SIZE, context, strlen(context));
	result = PARSE_SUCCESS;

fail_close_src:
	fclose(src);
fail_free_location:
	location_free(seed_location);
fail_none:
	return result;
}
