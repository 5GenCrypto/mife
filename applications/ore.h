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
  ORE_DEFAULT    = 0x00, //!< default behaviour
  ORE_NO_KILIAN    = 0x01, //!< do not multiply kilian randomizers into the encodings
  ORE_NO_RANDOMIZERS    = 0x02, //!< do not multiply the scalar randomizers into the encodings
  ORE_SIMPLE_PARTITIONS  = 0x04, //!< pick a simple partitioning (x[0] is encoded ta the universe, all others are encoded at the empty set.)
} ore_flag_t;

struct _gghlite_enc_mat_struct {
	gghlite_enc_t **m;
	int nrows; // number of rows in the matrix
	int ncols; // number of columns in the matrix
};

typedef struct _gghlite_enc_mat_struct gghlite_enc_mat_t[1];

struct _matrix_encodings_struct {
	int dary_repr[MAXN];
	int dim; // matrix dimension (FIXME: allow for rectangular matrices)
	int n; // the bitstring length
  gghlite_enc_mat_t x_enc[MAXN];
  gghlite_enc_mat_t y_enc[MAXN];
	fmpz_mat_t x_clr[MAXN];
	fmpz_mat_t y_clr[MAXN];
};

typedef struct _matrix_encodings_struct matrix_encodings_t[1];

struct _ore_sk_struct {
	gghlite_sk_t self;
	fmpz_mat_t R[MAXN];
	fmpz_mat_t R_inv[MAXN];
};

typedef struct _ore_sk_struct ore_sk_t[1];




/* functions dealing with fmpz types and matrix multiplications mod fmpz_t */
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);
int test_matrix_inv(int n, flint_rand_t randstate, fmpz_t modp);
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);
void fmpz_mat_scalar_mul_modp(fmpz_mat_t r, fmpz_mat_t m, fmpz_t scalar,
		fmpz_t modp, int dim);
void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n,
		fmpz_t p);
void gghlite_enc_mat_init(gghlite_params_t params, gghlite_enc_mat_t m,
		int nrows, int ncols);
void gghlite_enc_mat_clear(gghlite_enc_mat_t m);

/* functions dealing with ORE challenge generation */
void ore_sk_init(ore_sk_t sk, int n, int dim, flint_rand_t randstate,
		fmpz_t p);
void message_to_dary(int dary[MAXN], int bitstring_len,
		int64_t message, int64_t d);
int get_matrix_bit(int input, int i, int j, int type);
void set_matrices(matrix_encodings_t met, int64_t message, int d,
		int n, fmpz_t modp, flint_rand_t randstate, ore_flag_t flags);
void gen_partitioning(int partitioning[GAMMA], int i, int L, int d,
		int *len_p);
void mat_encode(gghlite_sk_t self, gghlite_enc_mat_t enc, fmpz_mat_t m,
		int group[GAMMA], flint_rand_t randstate);
void gghlite_enc_mat_zeros_print(gghlite_sk_t self, gghlite_enc_mat_t m);
void set_encodings(ore_sk_t sk, matrix_encodings_t met, int i, int L,
		flint_rand_t randstate, ore_flag_t flags);
void gghlite_enc_mat_mul(gghlite_params_t params, gghlite_enc_mat_t r,
		gghlite_enc_mat_t m1, gghlite_enc_mat_t m2);
void compare(gghlite_sk_t self, gghlite_enc_mat_t x[MAXN],
		gghlite_enc_mat_t y[MAXN], int n);

/* functions for computing matrix inverse mod fmpz_t */
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);
void fmpz_mat_det_modp(fmpz_t det, fmpz_mat_t a, int n, fmpz_t p);
void fmpz_mat_cofactor_modp(fmpz_mat_t b, fmpz_mat_t a, int n, fmpz_t p);
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);

/* test functions */
int int_arrays_equal(int arr1[MAXN], int arr2[MAXN], int length);
void test_gen_partitioning();
void test_dary_conversion();
int test_matrix_inv(int n, flint_rand_t randstate, fmpz_t modp);

// Function defines for tests
void test_gen_partitioning();
void test_dary_conversion();


#endif /* _ORE_H_ */
