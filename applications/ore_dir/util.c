#include <stdlib.h>
#include "util.h"

void check_parse_result(parse_result result, void usage(int), int problem) {
	switch(result) {
		case PARSE_SUCCESS: break;
		case PARSE_INVALID: usage(problem); break;
		case PARSE_OUT_OF_MEMORY: exit(-1); break;
		case PARSE_IO_ERROR: usage(problem); break;
	}
}
