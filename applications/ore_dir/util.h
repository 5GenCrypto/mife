#ifndef _ORE_UTILS_H
#define _ORE_UTILS_H

#include <stdbool.h>

#ifdef ALLOC_FAILS
#error "trying to redefine ALLOC_FAILS"
#else  /* ifdef ALLOC_FAILS */
#define ALLOC_FAILS(path, len) (NULL == ((path) = malloc((len) * sizeof(*(path)))))
#endif /* ifdef ALLOC_FAILS */

#ifdef SEED_SIZE
#error "trying to redefine SEED_SIZE"
#else /* ifdef SEED_SIZE */
#define SEED_SIZE 32
#endif /* ifdef SEED_SIZE */

typedef struct {
	char *path;
	bool stack_allocated;
} location;

void location_free(location loc);
location location_append(const location loc, const char *const path);

typedef enum {
	PARSE_SUCCESS,
	PARSE_INVALID,
	PARSE_OUT_OF_MEMORY,
	PARSE_IO_ERROR
} parse_result;

void check_parse_result(parse_result result, void usage(int), int problem);
bool create_directory_if_missing(char *dir);

#endif /* ifndef _ORE_UTILS_H */
