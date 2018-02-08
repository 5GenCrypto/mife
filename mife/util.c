#include <assert.h>
#include <errno.h>
#include <libgen.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "util.h"

void location_free(location loc) { if(!loc.stack_allocated && NULL != loc.path) free(loc.path); }

uint64_t T = 0.0;
int PRINT_TIMERS = 0;

location location_append(const location loc, const char *const path) {
	location result = { NULL, false };
	const unsigned int loc_len = strlen(loc.path), path_len = strlen(path);
	const int result_len = loc_len+1+path_len, result_size = result_len+1;
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

bool create_directory_if_missing(char *dir) {
	if(!mkdir(dir, S_IRWXU)) return true;
	switch(errno) {
		case EEXIST: return true;
		case ENOENT: break; /* try to create the parent and go again */
		default: return false;
	}

	/* construct name of parent */
	char *dir_copy, *dir_parent;
	if(ALLOC_FAILS(dir_copy, strlen(dir)+1)) return false;
	strcpy(dir_copy, dir);
	dir_parent = dirname(dir_copy);

	/* if we just tried to create the root directory or the current directory,
	 * and THAT resulted in ENOENT, we are super-duper hosed */
	if(!strcmp(dir_parent, ".")) {
		free(dir_copy);
		return false;
	}

	bool success = create_directory_if_missing(dir_parent) && (0 == mkdir(dir, S_IRWXU));
	free(dir_copy);
	return success;
}


/**
 * Code to find the inverse of a matrix, adapted from:
 * https://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html
 */

// uses gaussian elimination to obtain the determinant of a matrix
static void fmpz_mat_det_modp(fmpz_t det, fmpz_mat_t a, int n, fmpz_t p) {
  assert(n >= 1);

  if(n == 1) {
    fmpz_set(det, fmpz_mat_entry(a, 0, 0));
    return;
  }

  if (n == 2) {
    fmpz_t tmp1;
    fmpz_init(tmp1);
    fmpz_mul(tmp1, fmpz_mat_entry(a,0,0), fmpz_mat_entry(a,1,1));
    fmpz_mod(tmp1, tmp1, p);
    fmpz_t tmp2;
    fmpz_init(tmp2);
    fmpz_mul(tmp2, fmpz_mat_entry(a,1,0), fmpz_mat_entry(a,0,1));
    fmpz_mod(tmp2, tmp2, p);
    fmpz_sub(det, tmp1, tmp2);
    fmpz_mod(det, det, p);
    fmpz_clear(tmp1);
    fmpz_clear(tmp2);
    return;
  }

  fmpz_mat_t m;
  fmpz_mat_init_set(m, a);

  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_t multfactor;
  fmpz_init(multfactor);

  int num_swaps = 0;

  for(int j = 0; j < n; j++) {
    for(int i = j+1; i < n; i++) {

      if(fmpz_is_zero(fmpz_mat_entry(m, j, j))) {
        // find first row that isn't a zero, and swap
        int was_swapped = 0;
        int h;
        for(h = j+1; h < n; h++) {
          if(fmpz_is_zero(fmpz_mat_entry(m, h, j))) {
            continue;
          }

          // swap row h with row j
          for(int k = 0; k < n; k++) {
            fmpz_set(tmp, fmpz_mat_entry(m, h, k));
            fmpz_set(fmpz_mat_entry(m, h, k), fmpz_mat_entry(m, j, k));
            fmpz_set(fmpz_mat_entry(m, j, k), tmp);
          }
          was_swapped = 1;
          break;
        }

        if(!was_swapped) {
          // matrix is not invertible!
          fmpz_set_ui(det, 0);
          fmpz_clear(multfactor);
          fmpz_clear(tmp);
          fmpz_mat_clear(m);
          return;
        }

        num_swaps++;
      }

      fmpz_invmod(multfactor, fmpz_mat_entry(m, j, j), p);
      fmpz_mul(multfactor, multfactor, fmpz_mat_entry(m, i, j));
      fmpz_mod(multfactor, multfactor, p);
      for(int k = j; k < n; k++) {
        fmpz_mul(tmp, fmpz_mat_entry(m, j, k), multfactor);
        fmpz_sub(fmpz_mat_entry(m, i, k), fmpz_mat_entry(m, i, k), tmp);
        fmpz_mod(fmpz_mat_entry(m, i, k), fmpz_mat_entry(m, i, k), p);
      }
    }
  }

  fmpz_clear(multfactor);
  fmpz_clear(tmp);

  fmpz_set_ui(det, 1);

  for(int j = 0; j < n; j++) {
    fmpz_mul(det, det, fmpz_mat_entry(m, j, j));
  }
  if(num_swaps % 2 == 1) {
    fmpz_neg(det, det);
  }
  fmpz_mod(det, det, p);
  fmpz_mat_clear(m);
}

/*
   Find the cofactor matrix of a square matrix
*/
static void fmpz_mat_cofactor_modp(fmpz_mat_t b, fmpz_mat_t a, int n, fmpz_t p) {
  int i,j,ii,jj,i1,j1;

  fmpz_t det;
  fmpz_init(det);

  fmpz_mat_t c;
  fmpz_mat_init(c, n-1, n-1);

  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      /* Form the adjoint a_ij */
      i1 = 0;
      for (ii=0;ii<n;ii++) {
        if (ii == i)
          continue;
        j1 = 0;
        for (jj=0;jj<n;jj++) {
          if (jj == j)
            continue;
          fmpz_set(fmpz_mat_entry(c, i1, j1), fmpz_mat_entry(a, ii, jj));
          j1++;
        }
        i1++;
      }

      /* Calculate the determinant */
      /* TODO: check if the determinant is zero */
      fmpz_mat_det_modp(det, c, n-1, p);

      /* Fill in the elements of the cofactor */
      if((i+j) % 2 == 1) {
        fmpz_negmod(det, det, p);
      }
      fmpz_mod(det, det, p);
      fmpz_set(fmpz_mat_entry(b, i, j), det);
    }
  }

  fmpz_clear(det);
  fmpz_mat_clear(c);
}

int fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p) {
  if (dim == 1) {
      fmpz_set(fmpz_mat_entry(inv,0,0), fmpz_mat_entry(a, 0, 0));
  } else {
      fmpz_t det;
      fmpz_init(det);
      /* TODO: check if the determinant is zero */
      fmpz_mat_det_modp(det, a, dim, p);
      if(0 == fmpz_sgn(det)) {
          fmpz_clear(det);
          return 1;
      }
      fmpz_mat_t cofactor;
      fmpz_mat_init(cofactor, dim, dim);
      fmpz_mat_cofactor_modp(cofactor, a, dim, p);

      fmpz_t invmod;
      fmpz_init(invmod);
      fmpz_invmod(invmod, det, p);
      fmpz_t tmp;
      fmpz_init(tmp);
      for(int i = 0; i < dim; i++) {
          for(int j = 0; j < dim; j++) {
              fmpz_mod(tmp, fmpz_mat_entry(cofactor, j, i), p);
              fmpz_mul(tmp, tmp, invmod);
              fmpz_mod(tmp, tmp, p);
              fmpz_set(fmpz_mat_entry(inv,i,j), tmp);
          }
      }
      fmpz_clear(tmp);
      fmpz_clear(invmod);
      fmpz_mat_clear(cofactor);
      fmpz_clear(det);
  }
  return 0;
}
