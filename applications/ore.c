#include <gghlite/gghlite.h>
#include "common.h"

#define MAXN 30 // the maximum length bitstring
#define MAXW 20 // maximum matrix dimension

void print_int_matrix(int **inv, int n);
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int n, fmpz_t p);
int test_matrix_inv(int n, flint_rand_t randstate, fmpz_t modp);


// Function defines for tests
static void test_dary_conversion();
static void fmpz_print_matrix(fmpz_t matrix[MAXW][MAXW], int d);

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
	fmpz_mat_t kilian[MAXN];
};

typedef struct _ore_sk_struct ore_sk_t[1];

// initializes a matrix of fmpz_t of square dimension d
static void fmpz_init_matrix(fmpz_t m[MAXW][MAXW], int d) {
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			fmpz_init(m[i][j]);
		}
	}	
}

// clears a matrix of fmpz_t of square dimension d
static void fmpz_clear_matrix(fmpz_t m[MAXW][MAXW], int d) {
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			fmpz_clear(m[i][j]);
		}
	}	
}

// multiply two square dimension d matrices of fmpz_t types
static void fmpz_mmult(fmpz_t ret[MAXW][MAXW], fmpz_t m1[MAXW][MAXW], fmpz_t m2[MAXW][MAXW], int d) {
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			fmpz_zero(ret[i][j]);
			for (int k = 0; k < d; k++) {
				fmpz_t mul;
				fmpz_init(mul);
				fmpz_mul(mul, m1[i][k], m2[k][j]);
				fmpz_add(ret[i][j], ret[i][j], mul);
				fmpz_clear(mul);
			}
		}
	}
}


// multiply two square dimension d matrices of ints 
static void int_mmult(int **ret, int **m1, int **m2, int d) {
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			for (int k = 0; k < d; k++) {
				ret[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}

// multiply two square dimension d matrices of ints 
static void modp_mmult(int **ret, int **m1, int **m2, int d, int p) {
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			for (int k = 0; k < d; k++) {
				ret[i][j] += (m1[i][k] * m2[k][j]) % p;
			}
			ret[i][j] = ret[i][j] % p;
		}
	}
}



static void ore_sk_init(ore_sk_t sk, int n, int dim, flint_rand_t randstate, fmpz_t p) {
	for (int k = 0; k < n; k++) {
		fmpz_mat_init(sk->kilian[k], dim, dim);
		for (int i = 0; i < dim; i++) {
			for(int j = 0; j < dim; j++) {
				fmpz_randm(fmpz_mat_entry(sk->kilian[k], i, j), randstate, p);
			}
		}
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

/* sets the cleartext matrices x_clr and y_clr */
static void set_matrices(matrix_encodings_t met, int64_t message, int d, int n) {
	message_to_dary(met->dary_repr, n, message, d);
	met->n = n;
	met->dim = d+3;

	for (int k = 0; k < met->n; k++) {
		fmpz_mat_init(met->x_clr[k], met->dim, met->dim);
		fmpz_mat_init(met->y_clr[k], met->dim, met->dim);
		for (int i = 0; i < met->dim; i++) {
			for(int j = 0; j < met->dim; j++) {
				int x_digit = get_matrix_bit(met->dary_repr[k], i, j, X_TYPE);
				int y_digit = get_matrix_bit(met->dary_repr[k], i, j, Y_TYPE);
				fmpz_set_ui(fmpz_mat_entry(met->x_clr[k], i, j), x_digit);
				fmpz_set_ui(fmpz_mat_entry(met->y_clr[k], i, j), y_digit);
				/* TODO multiply by randomizer */
			}
		}
	}
}

void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n, fmpz_t p) {
	fmpz_mat_mul(a, b, c);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			fmpz_mod(fmpz_mat_entry(a, i, j), fmpz_mat_entry(a, i, j), p);
		}
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

  gghlite_sk_t self;


  gghlite_jigsaw_init_gamma(self,
                      cmdline_params->lambda,
                      cmdline_params->kappa,
											cmdline_params->gamma,
                      cmdline_params->flags,
                      randstate);

  printf("\n");
  gghlite_params_print(self->params);
  printf("\n---\n\n");

  t_gen = ggh_walltime(t);
  printf("1. GGH InstGen wall time:                 %8.2f s\n", ggh_seconds(t_gen));

  t = ggh_walltime(0);
  fmpz_t p; fmpz_init(p);
  fmpz_poly_oz_ideal_norm(p, self->g, self->params->n, 0);

	test_matrix_inv(6, randstate, p);
	test_dary_conversion();

	int n = 4;
	int d = 2;
	int dim = d+3;

	ore_sk_t sk;
	ore_sk_init(sk, n, dim, randstate, p);


/*
	matrix_encodings_t met;
	set_matrices(met, 10, d, n);

	printf("THE X0 MATRIX:\n");
	fmpz_mat_print_pretty(met->x_clr[0]);
	
	printf("THE Y0 MATRIX:\n");
	fmpz_mat_print_pretty(met->y_clr[0]);

	fmpz_mat_t C;
	fmpz_mat_init(C, dim, dim);
	fmpz_mat_mul(C, met->x_clr[0], sk->kilian[1]);
	fmpz_mat_print_pretty(C);
*/

/*
	fmpz_t result[MAXW][MAXW];
	fmpz_init_matrix(result, met->dim);
	fmpz_mmult(result, met->x_clr[0], met->y_clr[1], met->dim);

	printf("THE PRODUCT MATRIX:\n");
	fmpz_print_matrix(result, met->dim);
	fmpz_clear_matrix(result, met->dim);
*/

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

static void fmpz_print_matrix(fmpz_t matrix[MAXW][MAXW], int d) {
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			fmpz_print(matrix[i][j]);
			printf(" ");
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


static void test_clr_matrix_gen() {
	printf("\nTesting clr_matrix_gen function...                          ");

}

void print_int_matrix(int **inv, int n) {
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			printf("%d ", inv[i][j]);
		}
		printf("\n");
	}
}

// From Wikipedia: a^-1 (mod b)
int mod_inv(int a, int b)
{
	int b0 = b, t, q;
	int x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

// Computes a - b (mod p)
int mod_sub(int a, int b, int p)
{
	int tmp = a-b;
		return p+tmp;
	return tmp;
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

void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int n, fmpz_t p) {
	fmpz_t det;
	fmpz_init(det);
	fmpz_mat_det_modp(det, a, n, p);
	fmpz_mat_t cofactor;
	fmpz_mat_init(cofactor, n, n);
	fmpz_mat_cofactor_modp(cofactor, a, n, p);

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
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
