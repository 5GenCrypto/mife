#ifndef _MIFE_DEFS_H_
#define _MIFE_DEFS_H_

#include <aesrand/aesrand.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gghlite/misc.h>
#include "mbp_types.h"

#define debug_printf printf
#define CHECK(x, y) if((x) < (y)) { debug_printf( \
      "ERROR: fscanf() error encountered when trying to read from file\n" \
    ); }

typedef struct mmap_pp mmap_pp;

/* If we call init or fread, we will call clear. In particular, we will not
 * call clear on the mmap_pp we retrieve from an mmap_sk. */
typedef struct {
  void (*const clear)(mmap_pp *const);
  void (*const fread)(mmap_pp *const, FILE *const);
  void (*const fwrite)(const mmap_pp *const, FILE *const);
  const size_t size;
} mmap_pp_vtable;

typedef struct mmap_sk mmap_sk;

typedef struct {
  /* lambda: security parameter
   * kappa: how many multiplications we intend to do
   * gamma: the size of the universe that we will zero-test things at
   */
  void (*const init)(mmap_sk *const, size_t, size_t, size_t, aes_randstate_t);
  void (*const clear)(mmap_sk *const);
  void (*const fread)(mmap_sk *const, FILE *const);
  void (*const fwrite)(const mmap_sk *const, FILE *const);
  const size_t size;

  const mmap_pp *const (*const pp)(const mmap_sk *const);
  void (*const plaintext_field)(const mmap_sk *const, fmpz_t);
} mmap_sk_vtable;

typedef struct mmap_enc mmap_enc;

typedef struct {
  void (*const init)(mmap_enc *const, const mmap_pp *const);
  void (*const clear)(mmap_enc *const);
  void (*const fread)(mmap_enc *const, FILE *const);
  void (*const fwrite)(const mmap_enc *const, FILE *const);
  const size_t size;

  void (*const set)(mmap_enc *const, const mmap_enc *const);
  void (*const add)(mmap_enc *const, const mmap_pp *const, const mmap_enc *const, const mmap_enc *const);
  void (*const mul)(mmap_enc *const, const mmap_pp *const, const mmap_enc *const, const mmap_enc *const);
  bool (*const is_zero)(const mmap_enc *const, const mmap_pp *const);
  /* TODO: should this `int *` be `bool *`? */
  void (*const encode)(mmap_enc *const, const mmap_sk *const, const fmpz_t, int *, aes_randstate_t);
} mmap_enc_vtable;

typedef struct {
  const mmap_pp_vtable  *const pp;
  const mmap_sk_vtable  *const sk;
  const mmap_enc_vtable *const enc;
} mmap_vtable;

typedef const mmap_vtable *const const_mmap_vtable;

int NUM_ENCODINGS_GENERATED;
int PRINT_ENCODING_PROGRESS;
int NUM_ENCODINGS_TOTAL;

typedef enum {
  //!< default behaviour
  MIFE_DEFAULT = 0x00,
  
  //!< do not multiply kilian randomizers into the encodings
  MIFE_NO_KILIAN    = 0x01, 
  
  //!< do not multiply the scalar randomizers into the encodings
  MIFE_NO_RANDOMIZERS    = 0x02,

  //!< pick a simple partitioning (x[0] is encoded at the universe, all others 
  //are encoded at the empty set.)
  MIFE_SIMPLE_PARTITIONS  = 0x04,
} mife_flag_t;

struct _mmap_enc_mat_struct {
  int nrows; // number of rows in the matrix
  int ncols; // number of columns in the matrix
  mmap_enc ***m;
};

typedef struct _mmap_enc_mat_struct mmap_enc_mat_t[1];

struct _mife_mat_clr_struct {
  fmpz_mat_t **clr;
};

typedef struct _mife_mat_clr_struct mife_mat_clr_t[1];

struct _mife_ciphertext_struct {
  mmap_enc_mat_t **enc;
};

typedef struct _mife_ciphertext_struct mife_ciphertext_t[1];




struct _mife_pp_struct {
  int num_inputs; // the arity of the MBP (for comparisons, this is 2).
  int *n; // of length num_inputs
  int *gammas; // gamma for each input
  int L; // log # of plaintexts we can support
  int gamma; // should be sum of gammas[i] 
  int kappa; // the degree of multilinearity
  int numR; // number of kilian matrices. should be kappa-1
  mife_flag_t flags;
  fmpz_t p; // the prime, the order of the field
  mmap_pp *params_ref; // the underlying multilinear map's public parameters

  // MBP function pointers
  void *mbp_params; // additional parameters one can pass into the MBP setup
  int (*paramfn)  (struct _mife_pp_struct *, int); // for determining number of matrices per input
  void (*kilianfn)(struct _mife_pp_struct *, int *); // set kilian dimensions
  void (*orderfn) (struct _mife_pp_struct *, int, int *, int *); // function pointer for MBP ordering
  void (*setfn)   (struct _mife_pp_struct *, mife_mat_clr_t, void *);
  int (*parsefn)  (struct _mife_pp_struct *, f2_matrix); // function pointer for parsing output
};

typedef struct _mife_pp_struct mife_pp_t[1];

struct _mife_sk_struct {
  int numR;
  mmap_sk *self;
  fmpz_mat_t *R;
  fmpz_mat_t *R_inv;
};

typedef struct _mife_sk_struct mife_sk_t[1];




#endif /* _MIFE_IO_H_ */
