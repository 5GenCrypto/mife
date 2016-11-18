#include "cmdline.h"
#include <string.h>

const mife_flag_t mife_flags = MIFE_DEFAULT;

bool mbp_template_to_mife_pp(mife_pp_t pp, const mbp_template *const template, mbp_template_stats *const stats) {
	if(!mbp_template_to_mbp_template_stats(template, stats)) return false;
	mife_init_params(pp, mife_flags);
	mife_mbp_set(stats, pp, stats->positions_len,
		mbp_template_stats_to_params,
		mbp_template_stats_to_dimensions,
		mbp_template_stats_to_position,
		mbp_template_stats_to_cleartext,
		mbp_template_stats_to_result);
	return true;
}

parse_result load_seed(location private_location, char *context, aes_randstate_t seed) {
	FILE *src;
	char dest[AES_SEED_BYTE_SIZE];
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

	if(fread(dest, sizeof(*dest), AES_SEED_BYTE_SIZE, src) != AES_SEED_BYTE_SIZE) {
		fprintf(stderr, "could not read %d bytes of seed\n", AES_SEED_BYTE_SIZE);
		result = PARSE_INVALID;
		goto fail_close_src;
	}

	/* TODO: error checking */
	aes_randinit_seedn(seed, dest, AES_SEED_BYTE_SIZE, context, strlen(context));
	result = PARSE_SUCCESS;

fail_close_src:
	fclose(src);
fail_free_location:
	location_free(seed_location);
fail_none:
	return result;
}
