#ifndef _FLINT_RAW_IO_H_
#define _FLINT_RAW_IO_H_

#include "mife_defs.h"

/* raw reading and writing of fmpz types */
#define gghlite_enc_fprint_raw fmpz_mod_poly_fprint_raw
#define gghlite_enc_fread_raw fmpz_mod_poly_fread_raw
int fmpz_mod_poly_fprint_raw(FILE * file, const fmpz_mod_poly_t poly);
int gghlite_enc_fread_raw(FILE * f, gghlite_enc_t poly);
int fmpz_poly_fprint_raw(FILE * file, const fmpz_poly_t poly);
int fmpz_poly_fread_raw(FILE * file, fmpz_poly_t poly);
int fmpz_mat_fprint_raw(FILE * file, const fmpz_mat_t mat);
int fmpz_mat_fread_raw(FILE* file, fmpz_mat_t mat);

#endif /* _FLINT_RAW_IO_H_ */
