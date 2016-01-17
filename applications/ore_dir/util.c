#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

void location_free(location loc) { if(!loc.stack_allocated && NULL != loc.path) free(loc.path); }

location location_append(const location loc, const char *const path) {
	location result = { NULL, false };
	const unsigned int loc_len = strlen(loc.path), path_len = strlen(path);
	const unsigned int result_len = loc_len+1+path_len, result_size = result_len+1;
	if(ALLOC_FAILS(result.path, result_size)) return result;
	if(snprintf(result.path, result_size, "%s/%s", loc.path, path) != result_len) {
		free(result.path);
		result.path = NULL;
	}
	return result;
}

void check_parse_result(parse_result result, void usage(int), int problem) {
	switch(result) {
		case PARSE_SUCCESS: break;
		case PARSE_INVALID: usage(problem); break;
		case PARSE_OUT_OF_MEMORY: exit(-1); break;
		case PARSE_IO_ERROR: usage(problem); break;
	}
}
