#include <gghlite/gghlite.h>
#include "common.h"

#define MAXN 30 // the maximum length bitstring
#define MAXW 20 // maximum matrix dimension

void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);
int test_matrix_inv(int n, flint_rand_t randstate, fmpz_t modp);
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);


// Function defines for tests
static void test_gen_partitioning();
static void test_dary_conversion();

struct _enc_mat_struct {
	gghlite_enc_t m[MAXW][MAXW];
	int num_rows;
	int num_cols;
};

typedef struct _enc_mat_struct ent_mat_t[1];

struct _matrix_encodings_struct {
	int dary_repr[MAXN];
	int dim; // matrix dimension, dim <= MAXW
	int n; // the bitstring length
  gghlite_enc_t x_enc[MAXN][MAXW][MAXW];
  gghlite_enc_t y_enc[MAXN][MAXW][MAXW];
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

static void ore_sk_init(ore_sk_t sk, int n, int dim, flint_rand_t randstate, fmpz_t p) {
	// need 2 * n - 1 kilian randomizer matrices R
	int num_r = 2 * n - 1;
	for (int k = 0; k < num_r; k++) {
		fmpz_mat_init(sk->R[k], dim, dim);
		for (int i = 0; i < dim; i++) {
			for(int j = 0; j < dim; j++) {
				fmpz_randm(fmpz_mat_entry(sk->R[k], i, j), randstate, p);
			}
		}
		
		fmpz_mat_init(sk->R_inv[k], dim, dim);
		fmpz_modp_matrix_inverse(sk->R_inv[k], sk->R[k], dim, p);
	}
}

// message >= 0, d >= 2
static void message_to_dary(int dary[MAXN], int bitstring_len, int64_t message, int64_t d) {
	assert(message >= 0);
	assert(d >= 2);

	int i;
	for (i = bitstring_len - 1; i >= 0; i--) {
		dary[i] = message % d;
		message /= d;
	}
}

#define X_TYPE 0
#define Y_TYPE 1
#define NONZERO_VAL 1

/* input = the digit being read, (i,j) = coordinates of matrix, type = X or Y
 *
 * The DFA is defined as follows:
 * - state 0 is the equals state
 * - state 1 is the less than state
 * - state 2 is the greater than state
 *
 */
static int get_matrix_bit(int input, int i, int j, int type) {
	if (type == X_TYPE) {
		if ((i == 1 && j == 1) || (i == 2 && j == 2))
			return NONZERO_VAL;
		if (i == 0 && j == input+3)
			return NONZERO_VAL;
		return 0;
	} else {
		if ((i == 1 && j == 1) || (i == 2 && j == 2))
			return NONZERO_VAL;
		int input_state = i-3; // we just read this digit from x
		if (input_state < 0)
			return 0;
	
		// 'input_state' is the bit from x
		// 'input' is the bit from y

		if (j == 0 && input_state == input)
			return NONZERO_VAL;
		if (j == 1 && input_state < input)
			return NONZERO_VAL;
		if (j == 2 && input_state > input)
			return NONZERO_VAL;
		return 0;
	}
}

void fmpz_mat_scalar_mul_modp(fmpz_mat_t r, fmpz_mat_t m, fmpz_t scalar, fmpz_t modp, int dim) {
	for (int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			fmpz_mul(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), scalar);
			fmpz_mod(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), modp);
		}
	}	
}

/* sets the cleartext matrices x_clr and y_clr */
static void set_matrices(matrix_encodings_t met, int64_t message, int d, int n, flint_rand_t randstate, fmpz_t modp) {
	message_to_dary(met->dary_repr, n, message, d);
	met->n = n;
	met->dim = d+3;

	fmpz_t x_rand;
	fmpz_t y_rand;
	fmpz_init(x_rand);
	fmpz_init(y_rand);

	for (int k = 0; k < met->n; k++) {
		fmpz_mat_init(met->x_clr[k], met->dim, met->dim);
		fmpz_mat_init(met->y_clr[k], met->dim, met->dim);

		fmpz_randm(x_rand, randstate, modp);
		fmpz_randm(y_rand, randstate, modp);
		
		for (int i = 0; i < met->dim; i++) {
			for(int j = 0; j < met->dim; j++) {
				int x_digit = get_matrix_bit(met->dary_repr[k], i, j, X_TYPE);
				int y_digit = get_matrix_bit(met->dary_repr[k], i, j, Y_TYPE);
				fmpz_set_ui(fmpz_mat_entry(met->x_clr[k], i, j), x_digit);
				fmpz_set_ui(fmpz_mat_entry(met->y_clr[k], i, j), y_digit);

				/* apply scalar randomizers */
				fmpz_mat_scalar_mul_modp(met->x_clr[k], met->x_clr[k], x_rand, modp, met->dim);
				fmpz_mat_scalar_mul_modp(met->y_clr[k], met->y_clr[k], y_rand, modp, met->dim);
			}
		}
	}

	fmpz_clear(x_rand);
	fmpz_clear(y_rand);
}

/**
 * Generates the i^th member of the exclusive partition family for the index 
 * sets.
 *
 * @param partitioning The description of the partitioning, each entry is in 
 * [0,d-1] and it is of length (1 + (d-1)(L+1)).
 * @param i The i^th member of the partition family
 * @param L the log of the size of the partition family. So, i must be in the 
 * range [0,2^L-1]
 * @param d The number of total multiplications. partitioning[] will describe a 
 * d-partition of the universe set.
 * @param len_p A pointer to an integer to set the length of the partitioning 
 * array to. It should be (1 + (d-1)(L+1)).
 */ 
static void gen_partitioning(int partitioning[GAMMA], int i, int L, int d, int *len_p) {
	int j = 0;

	int bitstring[MAXN];
	message_to_dary(bitstring, L, i, 2);

	for(; j < d; j++) {
		partitioning[j] = j;
	}

	for(int k = 0; k < L; k++) {
		for(int j1 = 1; j1 < d; j1++) {
			partitioning[j] = (bitstring[k] == 1) ? j1 : 0;
			j++;
		}
	}
	*len_p = j;
}

void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n, fmpz_t p) {
	fmpz_mat_mul(a, b, c);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			fmpz_mod(fmpz_mat_entry(a, i, j), fmpz_mat_entry(a, i, j), p);
		}
	}
}

void mat_encode(gghlite_sk_t self, gghlite_enc_t enc[MAXW][MAXW], fmpz_mat_t m, int dim, int group[GAMMA], flint_rand_t randstate) {
	gghlite_clr_t e;
	gghlite_clr_init(e);
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			fmpz_poly_set_coeff_fmpz(e, 0, fmpz_mat_entry(m, i, j));
			gghlite_enc_init(enc[i][j], self->params);
			gghlite_enc_set_gghlite_clr(enc[i][j], self, e, 1, group, 1, randstate);
		}
	}
	gghlite_clr_clear(e);
}

void print_enc_matrix_zeros(gghlite_sk_t self, gghlite_enc_t m[MAXW][MAXW], int dim) {
	for(int i = 0; i < dim; i++) {
		printf("[");
		for(int j = 0; j < dim; j++) {
			printf(gghlite_enc_is_zero(self->params, m[i][j]) ? "0 " : "x " );
		}
		printf("]\n");
	}
}

void set_encodings(ore_sk_t sk, matrix_encodings_t met, int i, int L, flint_rand_t randstate) {
	int n = met->n;
	int ptnx[GAMMA];
	int lenx;
	gen_partitioning(ptnx, i, L, n, &lenx);
	int ptny[GAMMA];
	int leny;
	gen_partitioning(ptny, i, L, n, &leny);

	int partition_offset = 1 + (L+1) * (n-1);

	int group_x[MAXN][GAMMA];
	int group_y[MAXN][GAMMA];
	for(int j = 0; j < n; j++) {
		memset(group_x[j], 0, GAMMA * sizeof(int));
		memset(group_y[j], 0, GAMMA * sizeof(int));
		for(int k = 0; k < lenx; k++) {
			if(ptnx[k] == j) {
				group_x[j][k] = 1;
			}
		}
		for(int k = 0; k < leny; k++) {
			if(ptny[k] == j) {
				group_y[j][k+partition_offset] = 1;
			}
		}
	}

	// construct trivial partitioning (for testing correctness)
	/*
	for(int j = 0; j < n; j++) {
		memset(group_x[j], 0, GAMMA * sizeof(int));
		memset(group_y[j], 0, GAMMA * sizeof(int));
	}
	for(int k = 0; k < 2 * partition_offset; k++) {
		group_x[0][k] = 1;
	}
	*/


	// apply kilian and encode
	
	fmpz_mat_t tmp;
	fmpz_mat_init(tmp, met->dim, met->dim);

	fmpz_mat_mul(tmp, met->x_clr[0], sk->R[0]);
	mat_encode(sk->self, met->x_enc[0], tmp, met->dim, group_x[0], randstate);
	//mat_encode(sk->self, met->x_enc[0], met->x_clr[0], met->dim, group_x[0], 
	//randstate);


	for(int j = 1; j < n; j++) {
		fmpz_mat_mul(tmp, sk->R_inv[2 * j - 1], met->x_clr[j]);
		fmpz_mat_mul(tmp, tmp, sk->R[2 * j]);
		mat_encode(sk->self, met->x_enc[j], tmp, met->dim, group_x[j], randstate);
		//mat_encode(sk->self, met->x_enc[j], met->x_clr[j], met->dim, group_x[j], 
		//randstate);
	}

	for(int j = 0; j < n-1; j++) {
		fmpz_mat_mul(tmp, sk->R_inv[2 * j], met->y_clr[j]);
		fmpz_mat_mul(tmp, tmp, sk->R[2 * j + 1]);
		mat_encode(sk->self, met->y_enc[j], tmp, met->dim, group_y[j], randstate);
		//mat_encode(sk->self, met->y_enc[j], met->y_clr[j], met->dim, group_y[j], 
		//randstate);
	}

	fmpz_mat_mul(tmp, sk->R_inv[2 * (n-1)], met->y_clr[n-1]);
	mat_encode(sk->self, met->y_enc[n-1], tmp, met->dim, group_y[n-1], 
	randstate);
	//mat_encode(sk->self, met->y_enc[n-1], met->y_clr[n-1], met->dim, 
	//group_y[n-1], randstate);


	fmpz_mat_clear(tmp);
}

void enc_mul_init(gghlite_params_t params, gghlite_enc_t m[MAXW][MAXW], int dim) {
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			gghlite_enc_init(m[i][j], params);
		}
	}
}

void enc_mul_clear(gghlite_enc_t m[MAXW][MAXW], int dim) {
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			gghlite_enc_clear(m[i][j]);
		}
	}
}



void gghlite_enc_mmult(gghlite_params_t params, gghlite_enc_t r[MAXW][MAXW], gghlite_enc_t m1[MAXW][MAXW], gghlite_enc_t m2[MAXW][MAXW], int dim) {
	gghlite_enc_t tmp;
	gghlite_enc_init(tmp, params);

	gghlite_enc_t tmp_mat[MAXW][MAXW];
	enc_mul_init(params, tmp_mat, dim);
	
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			for(int k = 0; k < dim; k++) {
				gghlite_enc_mul(tmp, params, m1[i][k], m2[k][j]);
				gghlite_enc_add(tmp_mat[i][j], params, tmp_mat[i][j], tmp);
			}
		}
	}

	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			gghlite_enc_set(r[i][j], tmp_mat[i][j]);
		}
	}

	enc_mul_clear(tmp_mat, dim);
	gghlite_enc_clear(tmp);
}


void compare(gghlite_sk_t self, gghlite_enc_t x[MAXN][MAXW][MAXW], gghlite_enc_t y[MAXN][MAXW][MAXW], int n, int dim) {
	gghlite_enc_t tmp[MAXN][MAXW][MAXW];
	for(int i = 0; i < n; i++) {
		enc_mul_init(self->params, tmp[i], dim);
	}
	
	for(int i = 0; i < n; i++) {
		gghlite_enc_mmult(self->params, tmp[i], x[i], y[i], dim);
	}

	for(int i = 1; i < n; i++) {
		gghlite_enc_mmult(self->params, tmp[0], tmp[0], tmp[i], dim);
	}
	
	int equals = 1 - gghlite_enc_is_zero(self->params, tmp[0][0][0]);
	int lessthan = 1 - gghlite_enc_is_zero(self->params, tmp[0][0][1]);
	int greaterthan = 1 - gghlite_enc_is_zero(self->params, tmp[0][0][2]);

	printf("equals: %d, lessthan: %d, greaterthan: %d\n", equals, lessthan, greaterthan);

	print_enc_matrix_zeros(self, tmp[0], 5);
	assert(equals + lessthan + greaterthan == 1);

	if(equals) {
		printf("EQUALS\n");
	}

	if(lessthan) {
		printf("LESS THAN\n");
	}

	if(greaterthan) {
		printf("GREATER THAN\n");
	}

}


int main(int argc, char *argv[]) {
	cmdline_params_t cmdline_params;

  const char *name =  "Order Revealing Encryption";
  parse_cmdline(cmdline_params, argc, argv, name, NULL);

  print_header(name, cmdline_params);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, cmdline_params->seed, 1);

  uint64_t t = ggh_walltime(0);
  uint64_t t_total = ggh_walltime(0);

  uint64_t t_gen = 0;

	int n = 4;
	int d = 2;
	int dim = d+3;
	int L = 3; // 2^L = # of total messages we can encrypt

	int kappa = 2 * n;
	int gamma = 2 * (1 + (n-1) * (L+1));

	ore_sk_t sk;
	
  gghlite_jigsaw_init_gamma(sk->self,
                      cmdline_params->lambda,
                      kappa,
											gamma,
                      cmdline_params->flags,
                      randstate);

  printf("\n");
  gghlite_params_print(sk->self->params);
  printf("\n---\n\n");

  t_gen = ggh_walltime(t);
  printf("1. GGH InstGen wall time:                 %8.2f s\n", ggh_seconds(t_gen));

  t = ggh_walltime(0);
  fmpz_t p; fmpz_init(p);
  fmpz_poly_oz_ideal_norm(p, sk->self->g, sk->self->params->n, 0);

	ore_sk_init(sk, n, dim, randstate, p);


	//test_gen_partitioning(); // TODO: add a test for the partition selection
	test_matrix_inv(6, randstate, p);
	test_dary_conversion();


	matrix_encodings_t met1;
	set_matrices(met1, 9, d, n, randstate, p);
	set_encodings(sk, met1, 0, 3, randstate);

	matrix_encodings_t met2;
	set_matrices(met2, 10, d, n, randstate, p);
	set_encodings(sk, met2, 4, 3, randstate);

	/*
	gghlite_enc_t r[MAXW][MAXW];
	enc_mul_init(sk->self->params, r, 5);
	gghlite_enc_mmult(sk->self->params, r, met1->x_enc[0], met2->y_enc[0], 5);
	print_enc_matrix_zeros(sk->self, r, 5);
	*/

	compare(sk->self, met1->x_enc, met2->y_enc, 4, 5);


}

void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p) {
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			fmpz_mod(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), p);
		}
	}	
}

/**
 * Test code
 */

static int int_arrays_equal(int arr1[MAXN], int arr2[MAXN], int length) {
	for (int i = 0; i < length; i++) {
		if (arr1[i] != arr2[i])
			return 1;
	}
	return 0;
}

static void test_gen_partitioning() {
	int ptn[GAMMA];
	int L = 4;
	int i = 5;
	int d = 3;
	int len;

	for(int i = 0; i < 16; i++) {
		gen_partitioning(ptn, i, L, d, &len);

		printf("Partition %d: ", i);
		for(int j = 0; j < len; j++) {
			printf("%d ", ptn[j]);
		}
		printf("\n");
	}
}

static void test_dary_conversion() {
	printf("Testing d-ary conversion function...                          ");
	int dary1[MAXN];
	int dary2[MAXN];
	int dary3[MAXN];
	int correct1[] = {1,0,1,0};
	int correct2[] = {0,0,0,0,5,4,1,4};
	int correct3[] = {0,0,0,2};


	message_to_dary(dary1, 4, 10, 2);
	message_to_dary(dary2, 8, 1234, 6);
	message_to_dary(dary3, 4, 2, 11);


	int status = 0;
	status += int_arrays_equal(dary1, correct1, 4);
	status += int_arrays_equal(dary2, correct2, 8);
	status += int_arrays_equal(dary3, correct3, 4);



	if (status == 0)
		printf("SUCCESS\n");
	else
		printf("FAIL\n");	

}

int test_matrix_inv(int n, flint_rand_t randstate, fmpz_t modp) {
	printf("\nTesting matrix_inv function...                                ");
	fmpz_mat_t a;
	fmpz_mat_init(a, n, n);

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			fmpz_randm(fmpz_mat_entry(a, i, j), randstate, modp);
		}
	}

	fmpz_mat_t inv;
	fmpz_mat_init(inv, n, n);

	fmpz_mat_t prod;
	fmpz_mat_init(prod, n, n);

	fmpz_modp_matrix_inverse(inv, a, n, modp);
	fmpz_mat_mul_modp(prod, a, inv, n, modp);

	fmpz_mat_t identity;
	fmpz_mat_init(identity, n, n);
	fmpz_mat_one(identity);

	int status = fmpz_mat_equal(prod, identity);
	if (status != 0)
		printf("SUCCESS\n");
	else
		printf("FAIL\n");	
}

/**
 * Code to find the inverse of a matrix, adapted from:
 * https://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html
 */

/*
   Recursive definition of determinate using expansion by minors.
*/
void fmpz_mat_det_modp(fmpz_t det, fmpz_mat_t a, int n, fmpz_t p) {
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
	fmpz_mat_init(m, n-1, n-1);

  fmpz_set_ui(det, 0);
	for(int j1=0;j1<n;j1++) {
    for (int i=1;i<n;i++) {
			int j2 = 0;
      for (int j=0;j<n;j++) {
				if (j == j1)
					continue;
				fmpz_set(fmpz_mat_entry(m,i-1,j2), fmpz_mat_entry(a,i,j));
				j2++;
      }
    }
		fmpz_t det2;
		fmpz_init(det2);
		fmpz_mat_det_modp(det2, m, n-1, p);
		fmpz_mul(det2, det2, fmpz_mat_entry(a,0,j1));
		fmpz_mod(det2, det2, p);
		if(j1 % 2 == 1) {
			fmpz_negmod(det2, det2, p);
		}
		fmpz_add(det, det, det2);
		fmpz_clear(det2);
	}

	fmpz_mat_clear(m);
}

/*
   Find the cofactor matrix of a square matrix
*/
void fmpz_mat_cofactor_modp(fmpz_mat_t b, fmpz_mat_t a, int n, fmpz_t p) {
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

void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p) {
	fmpz_t det;
	fmpz_init(det);
	fmpz_mat_det_modp(det, a, dim, p);
	fmpz_mat_t cofactor;
	fmpz_mat_init(cofactor, dim, dim);
	fmpz_mat_cofactor_modp(cofactor, a, dim, p);

	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			fmpz_t invmod;
			fmpz_init(invmod);
			fmpz_invmod(invmod, det, p);
			fmpz_t tmp;
			fmpz_init(tmp);
			fmpz_mod(tmp, fmpz_mat_entry(cofactor, j, i), p);
			fmpz_mul(tmp, tmp, invmod);
			fmpz_mod(tmp, tmp, p);
			fmpz_set(fmpz_mat_entry(inv,i,j), tmp);
			fmpz_clear(invmod);
			fmpz_clear(tmp);
		}
	}

	fmpz_clear(det);
	fmpz_mat_clear(cofactor);
}
