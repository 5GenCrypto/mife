#ifndef _FLINT_RAW_IO_H_
#define _FLINT_RAW_IO_H_

#include "mife_defs.h"

/* raw reading and writing of fmpz types */
int fmpz_mat_fprint_raw(FILE * file, const fmpz_mat_t mat);
int fmpz_mat_fread_raw(FILE* file, fmpz_mat_t mat);

#endif /* _FLINT_RAW_IO_H_ */
