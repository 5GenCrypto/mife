#ifndef _ORE_H_
#define _ORE_H_

#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <gghlite/gghlite-defs.h>

#define MAXN 30 // the maximum length bitstring

int NUM_ENCODINGS_GENERATED;

/* these are used for get_matrix_bit() */
#define X_TYPE 0
#define Y_TYPE 1
#define NONZERO_VAL 1

typedef enum {
  //!< default behaviour
  ORE_DEFAULT    = 0x08,
  
  //!< do not multiply kilian randomizers into the encodings
  ORE_NO_KILIAN    = 0x01, 
  
  //!< do not multiply the scalar randomizers into the encodings
  ORE_NO_RANDOMIZERS    = 0x02,

  //!< pick a simple partitioning (x[0] is encoded at the universe, all others 
  //are encoded at the empty set.)
  ORE_SIMPLE_PARTITIONS  = 0x04,
  
  //!< the MBP produced is a series of 2n (d+3) x (d+3) matrices
  ORE_MBP_NORMAL = 0x08,

  //!< the MBP produced is degree-compressed by reading the input as: x0 (y0 y1) 
  //(x1 x2) (y2 y3) ...
  ORE_MBP_DC = 0x10,

  //!< the MBP produced is matrix-compressed, with dimensions 3 x d. This is not 
  //compatible with ORE_MBP_DC, though.
  ORE_MBP_MC = 0x20,
} ore_flag_t;

struct _gghlite_enc_mat_struct {
  gghlite_enc_t **m;
  int nrows; // number of rows in the matrix
  int ncols; // number of columns in the matrix
};

typedef struct _gghlite_enc_mat_struct gghlite_enc_mat_t[1];

struct _ore_mat_clr_struct {
  int dary_repr[MAXN]; // TODO make one for x and y
  fmpz_mat_t x_clr[MAXN];
  fmpz_mat_t y_clr[MAXN];
};

typedef struct _ore_mat_clr_struct ore_mat_clr_t[1];

struct _ore_ciphertext_struct {
  gghlite_enc_mat_t x_enc[MAXN];
  gghlite_enc_mat_t y_enc[MAXN];
};

typedef struct _ore_ciphertext_struct ore_ciphertext_t[1];

struct _ore_pp_struct {
  int d; // the base
  int nx; // number of x components
  int ny; // number of y components
  int L; // log # of plaintexts we can support
  int gammax; // number of indices needed for the x components
  int gammay; // number of indices needed for the y components
  int kappa; // the kappa for gghlite (degree of multilinearity)
  fmpz_t p; // the prime, the order of the field
  ore_flag_t flags;
  gghlite_params_t *params_ref; // gghlite's public parameters, by reference
};

struct _ore_sk_struct {
  // the following parameters actually need to be kept secret
  gghlite_sk_t self;
  fmpz_mat_t R[MAXN];
  fmpz_mat_t R_inv[MAXN];
  flint_rand_t randstate;
};

typedef struct _ore_pp_struct ore_pp_t[1];
typedef struct _ore_sk_struct ore_sk_t[1];


/* ORE interface */
void ore_setup(ore_pp_t pp, ore_sk_t sk, int bitstring_length, int base,
    int L, cmdline_params_t cmdline_params);
void ore_encrypt(ore_ciphertext_t ct, int message, ore_pp_t pp, ore_sk_t sk);
void compare(ore_pp_t pp, ore_ciphertext_t ct1, ore_ciphertext_t ct2);

/* functions dealing with fmpz types and matrix multiplications mod fmpz_t */
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);
int test_matrix_inv(int n, flint_rand_t randstate, fmpz_t modp);
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);
void fmpz_mat_scalar_mul_modp(fmpz_mat_t m, fmpz_t scalar, fmpz_t modp);
void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n,
    fmpz_t p);
void gghlite_enc_mat_init(gghlite_params_t params, gghlite_enc_mat_t m,
    int nrows, int ncols);
void gghlite_enc_mat_clear(gghlite_enc_mat_t m);

/* functions dealing with ORE challenge generation */
void apply_scalar_randomizers(ore_mat_clr_t met, ore_pp_t pp, ore_sk_t sk);
void message_to_dary(int dary[MAXN], int bitstring_len,
    int64_t message, int64_t d);
int get_matrix_bit_normal_mbp(int input, int i, int j, int type);
void set_matrices(ore_mat_clr_t met, int64_t message, ore_pp_t pp,
    ore_sk_t sk);
void gen_partitioning(int partitioning[GAMMA], int i, int L, int nu);
void mat_encode(ore_sk_t sk, gghlite_enc_mat_t enc, fmpz_mat_t m,
    int group[GAMMA]);
void gghlite_enc_mat_zeros_print(ore_pp_t pp, gghlite_enc_mat_t m);
void set_encodings(ore_ciphertext_t ct, ore_mat_clr_t met, int index,
    ore_pp_t pp, ore_sk_t sk);
void gghlite_enc_mat_mul(gghlite_params_t params, gghlite_enc_mat_t r,
    gghlite_enc_mat_t m1, gghlite_enc_mat_t m2);

/* functions for computing matrix inverse mod fmpz_t */
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);
void fmpz_mat_det_modp(fmpz_t det, fmpz_mat_t a, int n, fmpz_t p);
void fmpz_mat_cofactor_modp(fmpz_mat_t b, fmpz_mat_t a, int n, fmpz_t p);
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);

/* test functions */
int int_arrays_equal(int arr1[MAXN], int arr2[MAXN], int length);
void test_dary_conversion();
int test_matrix_inv(int n, flint_rand_t randstate, fmpz_t modp);

#endif /* _ORE_H_ */
