#ifndef _PARSE_H
#define _PARSE_H

#include <jsmn.h>
#include <stdbool.h>

#include "matrix.h"

typedef struct {
	char *path;
	bool stack_allocated;
} location;

void location_free(location loc);
location location_append(const location loc, const char *const path);

/* bool return type: whether the parse succeeded
 * json_string: the original string that was parsed
 * json_tokens: parser state (position in the token list produced by jsmn)
 * final argument: pointer to a C data structure where the result of the parse will be stored
 *
 * The parser for rows is a bit special: it takes an extra int and returns an
 * int. The input int tells how many columns to expect (if >=0); the returned
 * int tells how many columns were in the parsed row on a success or -1 on
 * failure.
 */
bool jsmn_parse_f2_elem  (const char *const json_string, const jsmntok_t **const json_tokens, bool      *const elem  );
int  jsmn_parse_f2_row   (const char *const json_string, const jsmntok_t **const json_tokens, bool     **const row   , int expected_num_cols);
bool jsmn_parse_f2_matrix(const char *const json_string, const jsmntok_t **const json_tokens, f2_matrix *const matrix);
bool jsmn_parse_f2_mbp   (const char *const json_string, const jsmntok_t **const json_tokens, f2_mbp    *const mbp   );
bool jsmn_parse_string   (const char *const json_string, const jsmntok_t **const json_tokens, char     **const string);
bool jsmn_parse_step     (const char *const json_string, const jsmntok_t **const json_tokens, step      *const step  );
bool jsmn_parse_template (const char *const json_string, const jsmntok_t **const json_tokens, template  *const templ );

/* top-level wrapper */
bool jsmn_parse_template_location(const location loc, template *const template);
#endif
