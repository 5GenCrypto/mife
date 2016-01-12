/**
 * do the following to generate a seed:
 * $ cat /dev/urandom | head -c 1024 | sha256sum | head -c 32 && echo
 *
 * And then use -h ${seed}
 *
 */

#include "ore.h"

#define CHECK(x) if(x < 0) { assert(0); }

int main(int argc, char *argv[]) {
  //printf("Last change from 2AM on 1/7/2016\n");
  run_tests();
  //ore_pp_gen(argc, argv);
  //ore_challenge_gen(argc, argv);
  //generate_plaintexts(argc, argv);
  //print_random_matrices_with_inverse(atoi(argv[1]), argv[2], argv[3], 
  //argv[4]);
}

int test_ciphertexts(char *pp_file, char *ct1_file, char *ct2_file) {
  mife_pp_t pp;
  mife_ciphertext_t *cts = malloc(2 * sizeof(mife_ciphertext_t));
  fread_mife_pp(pp, pp_file);
  fread_mife_ciphertext(pp, cts[0], ct1_file);
  fread_mife_ciphertext(pp, cts[1], ct2_file);
  int ret = mife_evaluate(pp, cts);
  mife_ciphertext_clear(pp, cts[0]);
  mife_ciphertext_clear(pp, cts[1]);
  free(cts);
  mife_clear_pp_read(pp);
  return ret;
}

void test_all_ciphertexts(char *pp_file, int n) {
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      char *str1 = malloc(20);
      char *str2 = malloc(20);
      sprintf(str1, "ct%d.out", i+1);
      sprintf(str2, "ct%d.out", j+1);
      printf("Checking ciphertext %d with ciphertext %d...           ", i, j);
      int result = test_ciphertexts(pp_file, str1, str2);
      if((i == j && result == 0) ||
          (i < j && result == 1) ||
          (i > j && result == 2)) {
        printf("SUCCESS\n");
      } else {
        printf("FAIL\n");
      }
    }
  }
}

void generate_plaintexts(int argc, char *argv[]) {
  cmdline_params_t cmdline_params;

  const char *name =  "Order Revealing Encryption";
  parse_cmdline(cmdline_params, argc, argv, name, NULL);

  int lambda = 80;
  int num_messages = 20;

  /* message space size will be base_d ^ base_n */
  int base_d = 10;
  int base_n = 10;

  mife_pp_t pp;

  fmpz_t message_space_size;
  fmpz_init_exp(message_space_size, base_d, base_n);
  ore_params_t *params = malloc(sizeof(ore_params_t));
  ore_set_best_params(pp, lambda, message_space_size, params);

  printf("Using SHA256 seed: %s\n", cmdline_params->shaseed);  

  /* Generate the messages first with a copied randstate */
  aes_randstate_t gen_randstate;
  aes_randinit_seed(gen_randstate, cmdline_params->shaseed);
  fmpz_t *messages = malloc(num_messages * sizeof(fmpz_t));
  printf("The plaintexts:\n");
  for(int i = 0; i < num_messages; i++) {
    fmpz_init(messages[i]);
    fmpz_randm_aes(messages[i], gen_randstate, message_space_size);
    printf("%lu\n", fmpz_get_ui(messages[i]));
  }
  fmpz_clear(message_space_size);
  aes_randclear(gen_randstate);
  free(params);
  printf("\n");
}

void ore_pp_gen(int argc, char *argv[]) {
  cmdline_params_t cmdline_params;

  const char *name =  "Order Revealing Encryption";
  parse_cmdline(cmdline_params, argc, argv, name, NULL);

  int L = 80; // 2^L = # of total messages we can encrypt
  int lambda = cmdline_params->lambda; // security parameter

  /* message space size will be base_d ^ base_n */
  int base_d = 10;
  int base_n = 10;

  mife_pp_t pp;
  mife_sk_t sk;

  fmpz_t message_space_size;
  fmpz_init_exp(message_space_size, base_d, base_n);
  ore_params_t *params = malloc(sizeof(ore_params_t));
  ore_set_best_params(pp, lambda, message_space_size, params);

  printf("Estimated ciphertext size for lambda = 80: ");
  const char *units[3] = {"KB","MB","GB"};
  double sd = CT_SIZE / 8.0;
  int i;
  for(i=0; i<3; i++) {
    if (sd < 1024.0)
      break;
    sd = sd/1024;
  }
  printf("%6.2f %s\n\n", sd, units[i-1]);

  mife_mbp_set(params, pp, 2, &ore_mbp_param, &ore_mbp_kilian,
      &ore_mbp_ordering, &ore_mbp_set_matrices, &ore_mbp_parse);

  gghlite_flag_t ggh_flags = GGHLITE_FLAGS_DEFAULT | GGHLITE_FLAGS_GOOD_G_INV;
  reset_T();
  mife_setup(pp, sk, L, lambda, ggh_flags, cmdline_params->shaseed);
  printf("Completed MIFE setup. (Time elapsed: %8.2f s)\n", get_T());

  fwrite_mife_pp(pp, "pp.out");
  mife_clear_pp(pp);
  mife_clear_sk(sk);
  mpfr_free_cache();
  flint_cleanup();
  free(params);
  printf("\nSUCCESS\n");

}

void ore_challenge_gen(int argc, char *argv[]) {
  cmdline_params_t cmdline_params;

  const char *name =  "Order Revealing Encryption";
  parse_cmdline(cmdline_params, argc, argv, name, NULL);

  int L = 80; // 2^L = # of total messages we can encrypt
  int lambda = cmdline_params->lambda; // security parameter
  int num_messages = 20;

  /* message space size will be base_d ^ base_n */
  int base_d = 10;
  int base_n = 10;

  mife_pp_t pp;
  mife_sk_t sk;

  fmpz_t message_space_size;
  fmpz_init_exp(message_space_size, base_d, base_n);
  ore_params_t *params = malloc(sizeof(ore_params_t));
  ore_set_best_params(pp, lambda, message_space_size, params);

  /* Read in the messages */
  fmpz_t *messages = malloc(num_messages * sizeof(fmpz_t));
  FILE *fp = fopen("plaintexts.secret", "r");
  if(!fp) {
    printf("ERROR: Could not find file plaintexts.secret\n");
    exit(1);
  }
  for(int i = 0; i < num_messages; i++) {
    unsigned long m;
    CHECK(fscanf(fp, "%lu\n", &m));
    fmpz_init_set_ui(messages[i], m);
  }
  fclose(fp);
  fmpz_clear(message_space_size);


  printf("Estimated ciphertext size for lambda = 80: ");
  const char *units[3] = {"KB","MB","GB"};
  double sd = CT_SIZE / 8.0;
  int i;
  for(i=0; i<3; i++) {
    if (sd < 1024.0)
      break;
    sd = sd/1024;
  }
  printf("%6.2f %s\n\n", sd, units[i-1]);

  mife_mbp_set(params, pp, 2, &ore_mbp_param, &ore_mbp_kilian,
      &ore_mbp_ordering, &ore_mbp_set_matrices, &ore_mbp_parse);

  gghlite_flag_t ggh_flags = GGHLITE_FLAGS_DEFAULT | GGHLITE_FLAGS_GOOD_G_INV;
  reset_T();
  mife_setup(pp, sk, L, lambda, ggh_flags, cmdline_params->shaseed);
  printf("Completed MIFE setup. (Time elapsed: %8.2f s)\n", get_T());

  fwrite_mife_pp(pp, "pp.out");

  int ci = cmdline_params->challenge_index;
 
  PRINT_ENCODING_PROGRESS = 1;
  mife_ciphertext_t *ciphertexts =
    malloc(num_messages * sizeof(mife_ciphertext_t));
  for(int i = 0; i < num_messages; i++) { 
    if(ci == 0 || i == ci-1) {
      printf("Encrypting message %d\n", i+1);
      NUM_ENCODINGS_GENERATED = 0;
      reset_T();
      mife_encrypt(ciphertexts[i], messages[i], pp, sk);
      char *str = malloc(10 * sizeof(char));
      sprintf(str, "ct%d.out", i+1);
      fwrite_mife_ciphertext(pp, ciphertexts[i], str);
      free(str);
      mife_ciphertext_clear(pp, ciphertexts[i]);

    }
  }
  PRINT_ENCODING_PROGRESS = 0;

  free(ciphertexts);
  free(messages);

  mife_clear_pp(pp);
  mife_clear_sk(sk);
  mpfr_free_cache();
  flint_cleanup();
  free(params);
  printf("\nSUCCESS\n");
}

long dc_enc_size(int n) {
  int kappa = n+1;
  if(kappa < 2 || kappa > MAX_KAPPA_BENCH) {
    return -1;
  }
  return KAPPA_BENCH[kappa];
}

long mc_enc_size(int n) {
  int kappa = 2*n;
  if(kappa < 2 || kappa > MAX_KAPPA_BENCH) {
    return -1;
  }
  return KAPPA_BENCH[kappa];
}

int dc_num_enc(int d, int n) {
  return d*d*(n-1) + (d+1)*(4*n-2);
}

int mc_num_enc(int d, int n) {
  return 6*(n-1)*(d+2) + 4*d;
}

/**
 * The message space size is d^n.
 */
void ore_set_best_params(mife_pp_t pp, int lambda, fmpz_t message_space_size,
    ore_params_t *params) {
  int max_base = 5;

  mpfr_t total;
  mpfr_init_set_ui(total, fmpz_get_ui(message_space_size), MPFR_RNDN);
  
  mpfr_t tmp1;
  mpfr_t tmp2;
  mpfr_init(tmp1);
  mpfr_init(tmp2);
  mpfr_t bt;
  mpfr_init(bt);

  long dc_vals[max_base+1];
  long mc_vals[max_base+1];
  long nmap[max_base+1];

  for(int base = 2; base <= max_base; base++) {
    // compute minimum n
    mpfr_set_ui(bt, base, MPFR_RNDN);
    mpfr_log(tmp1, total, MPFR_RNDN);
    mpfr_log(tmp2, bt, MPFR_RNDN);
    mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);
    mpfr_ceil(tmp1, tmp1);
    unsigned long n = mpfr_get_ui(tmp1, MPFR_RNDN);
    nmap[base] = n;
    
    // num encodings #1, #2
    int dc_enc = dc_num_enc(base,n);
    int mc_enc = mc_num_enc(base,n);
   
    // use n to compute enc size
    dc_vals[base] = dc_enc * dc_enc_size(n);
    mc_vals[base] = mc_enc * mc_enc_size(n);
  }

  mpfr_clear(bt);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  mpfr_clear(total);

  long min_enc = dc_vals[5]; 
  int min_type = 0;
  int min_base;
  int i = 2;
  for(; i <= max_base; i++) {
    if(dc_vals[i] > 0) {
      min_type = 0;
      min_enc = dc_vals[i];
      break;
    }
    if(mc_vals[i] > 0) {
      min_type = 1;
      min_enc = mc_vals[i];
      break;
    }
  }
  min_base = i;

  for(; i <= max_base; i++) {
    if(dc_vals[i] > 0 && dc_vals[i] < min_enc) {
      min_type = 0;
      min_base = i;
      min_enc = dc_vals[i];
    }
    if(mc_vals[i] > 0 && mc_vals[i] < min_enc) {
      min_type = 1;
      min_base = i;
      min_enc = mc_vals[i];
    }
  }

  if(min_type == 0) {
    ORE_GLOBAL_FLAGS = ORE_MBP_DC;
  } else {
    ORE_GLOBAL_FLAGS = ORE_MBP_MC;
  }

  params->d = min_base;  
  mife_init_params(pp, nmap[min_base], MIFE_DEFAULT);

  if(ORE_GLOBAL_FLAGS & ORE_MBP_DC) {
    set_NUM_ENC(dc_num_enc(params->d, pp->bitstr_len));
  } else if(ORE_GLOBAL_FLAGS & ORE_MBP_MC) {
    set_NUM_ENC(mc_num_enc(params->d, pp->bitstr_len));
  }
  CT_SIZE = min_enc;
}

/**
 * Creates the first x-matrix in the matrix-compressed version of ORE.
 *
 * This matrix has a single row (so it's actually a vector), and d columns, 
 * which each represent the bit that is being read.
 *
 * @param m The matrix
 * @param input1 A number in [0,d-1]
 * @param input2 A number in [0,d-1]
 * @param d The base
 */
void ore_mc_clrmat_init_XFIRST(fmpz_mat_t m, int input, int d) {
  fmpz_mat_init(m, 1, d);
  fmpz_mat_zero(m);
  fmpz_set_ui(fmpz_mat_entry(m, 0, input), NONZERO_VAL);
}

/**
 * Creates the first y-matrix in the matrix-compressed version of ORE.
 *
 * This matrix is d x 3.
 *
 * @param m The matrix
 * @param input1 A number in [0,d-1]
 * @param input2 A number in [0,d-1]
 * @param d The base
 */
void ore_mc_clrmat_init_YFIRST(fmpz_mat_t m, int input, int d) {
  fmpz_mat_init(m, d, 3);
  fmpz_mat_zero(m);
 
  for(int i = 0; i < d; i++) {
    int j;
    if(i == input) {
      j = 0;
    } else if(i < input) {
      j = 1;
    } else { // i > input
      j = 2;
    }
    fmpz_set_ui(fmpz_mat_entry(m, i, j), NONZERO_VAL);
  }
}

/**
 * Creates the x matrices (not the first) in the matrix-compressed version of 
 * ORE.
 *
 * This matrix is 3 x (d+2).
 *
 * @param m The matrix
 * @param input1 A number in [0,d-1]
 * @param input2 A number in [0,d-1]
 * @param d The base
 */
void ore_mc_clrmat_init_XREST(fmpz_mat_t m, int input, int d) {
  fmpz_mat_init(m, 3, d+2);
  fmpz_mat_zero(m);
 
  for(int i = 0; i < 3; i++) {
    int j;
    if(i == 0) {
      j = input;
    } else if(i == 1) {
      j = d;
    } else { // i == 2
      j = d+1;
    }
    fmpz_set_ui(fmpz_mat_entry(m, i, j), NONZERO_VAL);
  }
}

/**
 * Creates the y matrices (not the first) in the matrix-compressed version of 
 * ORE.
 *
 * This matrix is (d+2) x 3.
 *
 * @param m The matrix
 * @param input1 A number in [0,d-1]
 * @param input2 A number in [0,d-1]
 * @param d The base
 */
void ore_mc_clrmat_init_YREST(fmpz_mat_t m, int input, int d) {
  fmpz_mat_init(m, d+2, 3);
  fmpz_mat_zero(m);
 
  for(int i = 0; i < d+2; i++) {
    int j;
    if(i == d) {
      j = 1;
    } else if(i == d+1) {
      j = 2;
    } else if(i == input) {
      j = 0;
    } else if(i < input) {
      j = 1;
    } else { // i > input
      j = 2;
    }
    fmpz_set_ui(fmpz_mat_entry(m, i, j), NONZERO_VAL);
  }
}



/**
 * Creates the first x-matrix in the degree-compressed version of ORE.
 *
 * This matrix has a single row (so it's actually a vector), and d columns, 
 * which each represent the bit that is being read.
 *
 * @param m The matrix
 * @param input1 A number in [0,d-1]
 * @param input2 A number in [0,d-1]
 * @param d The base
 */
void ore_dc_clrmat_init_FIRST(fmpz_mat_t m, int input, int d) {
  fmpz_mat_init(m, 1, d);
  fmpz_mat_zero(m);
  fmpz_set_ui(fmpz_mat_entry(m, 0, input), NONZERO_VAL);
}

/**
 * Creates the first y-matrix in the degree-compressed version of ORE.
 *
 * We use columns [0,d-1] to represent the bit being read, column d to represent 
 * '<', and d+1 to represent '>'.
 *
 * @param m The matrix
 * @param input1 A number in [0,d-1]
 * @param input2 A number in [0,d-1]
 * @param d The base
 */
void ore_dc_clrmat_init_SECOND(fmpz_mat_t m, int input1, int input2, int d) {
  fmpz_mat_init(m, d, d+2);
  fmpz_mat_zero(m);
  
  for(int i = 0; i < d; i++) {
    int j;
    if(i < input1) {
      j = d;
    } else if(i > input1) {
      j = d+1;
    } else { // i == input1
      j = input2;
    }
    fmpz_set_ui(fmpz_mat_entry(m, i, j), NONZERO_VAL);
  }
}

/**
 * Creates the "middle" matrices in the degree-compressed version of ORE.
 *
 * We use columns [0,d-1] to represent the bit being read, column d to represent 
 * '<', and d+1 to represent '>'. This is the same as the second matrix, except 
 * we have two more rows, and the last two rows just follow the identity matrix 
 * (since if we were in the '<' state, we stay in that state, and the same is 
 * true for the '>' state).
 *
 * @param m The matrix
 * @param input1 A number in [0,d-1]
 * @param input2 A number in [0,d-1]
 * @param d The base
 * @param type Either X_TYPE or Y_TYPE
 */
void ore_dc_clrmat_init_MIDDLE(fmpz_mat_t m, int input1, int input2, int d,
    int type) {
  fmpz_mat_init(m, d+2, d+2);
  fmpz_mat_zero(m);
  
  for(int i = 0; i < d; i++) {
    int j = -1;
    if(i < input1) {
      j = (type == X_TYPE) ? d+1 : d;
    } else if(i > input1) {
      j = (type == X_TYPE) ? d : d+1;
    } else { // i == input1
      j = input2;
    }
    fmpz_set_ui(fmpz_mat_entry(m, i, j), NONZERO_VAL);
  }
  fmpz_set_ui(fmpz_mat_entry(m, d, d), NONZERO_VAL);
  fmpz_set_ui(fmpz_mat_entry(m, d+1, d+1), NONZERO_VAL);
}

/**
 * Creates the last matrix in the degree-compressed version of ORE.
 *
 * We use column 0 to represent '=', column 1 for '<', and column 2 for '>'.
 *
 * @param m The matrix
 * @param input A number in [0,d-1]
 * @param d The base
 * @param type Either X_TYPE or Y_TYPE
 */
void ore_dc_clrmat_init_LAST(fmpz_mat_t m, int input, int d, int type) {
  fmpz_mat_init(m, d+2, 3);
  fmpz_mat_zero(m);
  
  for(int i = 0; i < d+2; i++) {
    int j;
    if(i == d) {
      j = 1;
    } else if(i == d+1) {
      j = 2;
    } else if(i == input) {
      j = 0;
    } else if(i < input) {
      j = (type == X_TYPE) ? 2 : 1;
    } else {
      assert(i > input);
      j = (type == X_TYPE) ? 1 : 2;
    }
    fmpz_set_ui(fmpz_mat_entry(m, i, j), NONZERO_VAL);
  }
}


/**
 * Creates the (special case) second matrix in the degree-compressed version of 
 * ORE.
 *
 * This only is used when the second matrix is also the last matrix, in which 
 * case the dimensions are d x 3.
 *
 * We use column 0 to represent '=', column 1 for '<', and column 2 for '>'.
 *
 * @param m The matrix
 * @param input A number in [0,d-1]
 * @param d The base
 */
void ore_dc_clrmat_init_SECONDANDLAST(fmpz_mat_t m, int input, int d) {
  fmpz_mat_init(m, d, 3);
  fmpz_mat_zero(m);
  
  for(int i = 0; i < d; i++) {
    int j = -1;
    if(i == input) {
      j = 0;
    } else if(i < input) {
      j = 1;
    } else if(i > input) {
      j = 2;
    }
    fmpz_set_ui(fmpz_mat_entry(m, i, j), NONZERO_VAL);
  }
}


/* input = the digit being read, (i,j) = coordinates of matrix, type = X or Y
 *
 * The DFA is defined as follows:
 * - state 0 is the equals state
 * - state 1 is the less than state
 * - state 2 is the greater than state
 *
 */
int ore_get_matrix_bit_normal_mbp(int input, int i, int j, int type) {
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

int ore_mbp_param(int bitstr_len, int index) {
  if(ORE_GLOBAL_FLAGS & ORE_MBP_NORMAL) {
    return bitstr_len;
  }

  if(ORE_GLOBAL_FLAGS & ORE_MBP_DC) {
    if(index == 0) {
      return bitstr_len / 2 + 1;
    }
    if(index == 1) {
      return (bitstr_len + 1) / 2;
    }
  }
  
  if(ORE_GLOBAL_FLAGS & ORE_MBP_MC) {
    return bitstr_len;
  }

  return -1; // error
}

void ore_mbp_kilian(mife_pp_t pp, int *dims) {
  // set the kilian dimensions
  int d = ((ore_params_t *) pp->mbp_params)->d;
  if(ORE_GLOBAL_FLAGS & ORE_MBP_NORMAL) {
    for(int i = 0; i < pp->numR; i++) {
      dims[i] = d+3;
    }
  } else if(ORE_GLOBAL_FLAGS & ORE_MBP_DC) {
    dims[0] = d;
    for(int i = 1; i < pp->numR; i++) {
      dims[i] = d+2;
    }
  } else if(ORE_GLOBAL_FLAGS & ORE_MBP_MC) {
    dims[0] = d;
    for(int i = 1; i < pp->numR; i++) {
      if(i % 2 == 1) {
        dims[i] = 3;
      } else {
        dims[i] = d+2;
      }
    }
  } else {
    assert(0);
  }
}


/* sets the cleartext matrices */
void ore_mbp_set_matrices(mife_mat_clr_t met, fmpz_t message, mife_pp_t pp,
    mife_sk_t sk) {
  int d = ((ore_params_t *) pp->mbp_params)->d;
  ulong *dary_repr = malloc(pp->bitstr_len * sizeof(ulong));
  message_to_dary(dary_repr, pp->bitstr_len, message, d);

  assert(pp->num_inputs == 2);

  met->clr = malloc(pp->num_inputs * sizeof(fmpz_mat_t *));
  for(int i = 0; i < pp->num_inputs; i++) {
    met->clr[i] = malloc(pp->n[i] * sizeof(fmpz_mat_t));
  }

  if(ORE_GLOBAL_FLAGS & ORE_MBP_NORMAL) {
    assert(pp->n[0] == pp->n[1]);
    for(int k = 0; k < pp->n[0]; k++) {
      int dim = d+3;
      for(int i = 0; i < pp->num_inputs; i++) {
        fmpz_mat_init(met->clr[i][k], dim, dim);
      }

      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          int x_digit = ore_get_matrix_bit_normal_mbp(dary_repr[k],
              i, j, X_TYPE);
          int y_digit = ore_get_matrix_bit_normal_mbp(dary_repr[k],
              i, j, Y_TYPE);
          fmpz_set_ui(fmpz_mat_entry(met->clr[0][k], i, j), x_digit);
          fmpz_set_ui(fmpz_mat_entry(met->clr[1][k], i, j), y_digit);
        }
      }
    }
  } else if(ORE_GLOBAL_FLAGS & ORE_MBP_DC) {
    for(int k = 0, bc = 0; k < pp->n[0]; k++, bc++) {
      if(bc == 0) {
        ore_dc_clrmat_init_FIRST(met->clr[0][k], dary_repr[bc], d);
      } else if(bc == pp->bitstr_len - 1) {
        ore_dc_clrmat_init_LAST(met->clr[0][k], dary_repr[bc], d,
            X_TYPE);
      } else {
        ore_dc_clrmat_init_MIDDLE(met->clr[0][k], dary_repr[bc],
            dary_repr[bc+1], d, X_TYPE);
        bc++;
      }
    }
    
    for(int k = 0, bc = 0; k < pp->n[1]; k++, bc++) {
      if(k == 0 && pp->n[1] > 1) {
        ore_dc_clrmat_init_SECOND(met->clr[1][k], dary_repr[bc],
            dary_repr[bc+1], d);
        bc++;
      } else if(k == 0 && pp->n[1] == 1) {
        ore_dc_clrmat_init_SECONDANDLAST(met->clr[1][k], dary_repr[bc],
            d);
      } else if((pp->bitstr_len % 2 == 1) && (k == pp->n[1]-1)) {
        ore_dc_clrmat_init_LAST(met->clr[1][k], dary_repr[bc], d,
            Y_TYPE);
      } else {
        ore_dc_clrmat_init_MIDDLE(met->clr[1][k], dary_repr[bc],
            dary_repr[bc+1], d, Y_TYPE);
        bc++;
      }
    }
  } else if(ORE_GLOBAL_FLAGS & ORE_MBP_MC) {
    for(int k = 0; k < pp->n[0]; k++) {
      if(k == 0) {
        ore_mc_clrmat_init_XFIRST(met->clr[0][k], dary_repr[k], d);
      } else {
        ore_mc_clrmat_init_XREST(met->clr[0][k], dary_repr[k], d);
      }
    }
    for(int k = 0; k < pp->n[1]; k++) {
      if(k == 0) {
        ore_mc_clrmat_init_YFIRST(met->clr[1][k], dary_repr[k], d);
      } else {
        ore_mc_clrmat_init_YREST(met->clr[1][k], dary_repr[k], d);
      }
    }
  } else {
    assert(0);
  }
  free(dary_repr);    
}

// index ranges from [0, kappa - 1]
void ore_mbp_ordering(int index, int *ip, int *jp) {
  *ip = index % 2;
  *jp = index / 2;
}

int ore_mbp_parse(char **m) {
  if(m[0][0] + m[0][1] + m[0][2] != 1)
    return -1;

  for(int i = 0; i < 3; i++) {
    if(m[0][i] == 1)
      return i;
  }

  return -1;
}

int test_ore(int lambda, int mspace_size, int num_messages, int d,
    int bitstr_len, ore_flag_t flags, int verbose) {
  printf("Testing ORE...                          ");
  ORE_GLOBAL_FLAGS = flags;
  int status = 0;
  int L = 80;
  mife_pp_t pp;
  mife_sk_t sk;

  mife_init_params(pp, bitstr_len, flags);
  
  ore_params_t *params = malloc(sizeof(ore_params_t));
  params->d = d;
  mife_mbp_set(params, pp, 2, &ore_mbp_param, &ore_mbp_kilian,
      &ore_mbp_ordering, &ore_mbp_set_matrices, &ore_mbp_parse);

  gghlite_flag_t ggh_flags = GGHLITE_FLAGS_QUIET | GGHLITE_FLAGS_GOOD_G_INV;
  mife_setup(pp, sk, L, lambda, ggh_flags, DEFAULT_SHA_SEED);

  fwrite_mife_pp(pp, "pp_test.out");
  mife_clear_pp(pp);
  mife_pp_t pp2;
  fread_mife_pp(pp2, "pp_test.out");

  fwrite_mife_sk(sk, "sk_test.out");
  mife_clear_sk(sk);
  mife_sk_t sk2;
  fread_mife_sk(sk2, "sk_test.out");
  aes_randinit(sk2->randstate);

  mife_mbp_set(params, pp2, 2, &ore_mbp_param, &ore_mbp_kilian,
      &ore_mbp_ordering, &ore_mbp_set_matrices, &ore_mbp_parse);

  fmpz_t message_space_size;
  fmpz_init_set_ui(message_space_size, mspace_size);

  fmpz_t *messages = malloc(num_messages * sizeof(fmpz_t));
  for(int i = 0; i < num_messages; i++) {
    fmpz_init(messages[i]);
    fmpz_randm_aes(messages[i], sk2->randstate, message_space_size);
    if(verbose) {
      fmpz_print(messages[i]);
      printf("\n");
    }
  }

  mife_ciphertext_t *ciphertexts = malloc(num_messages * sizeof(mife_ciphertext_t));
  for(int i = 0; i < num_messages; i++) {
    mife_encrypt(ciphertexts[i], messages[i], pp2, sk2);
  }

  for(int i = 0; i < num_messages; i++) {
    for(int j = 0; j < num_messages; j++) {
      mife_ciphertext_t *cts = malloc(2 * sizeof(mife_ciphertext_t));
      memcpy(cts[0], ciphertexts[i], sizeof(mife_ciphertext_t));
      memcpy(cts[1], ciphertexts[j], sizeof(mife_ciphertext_t));

      int compare = mife_evaluate(pp2, cts);
      fmpz_t modi, modj, true_mspace_size, fmpzd;
      fmpz_init(modi);
      fmpz_init(modj);
      fmpz_init(true_mspace_size);
      fmpz_init_set_ui(fmpzd, d);
      fmpz_pow_ui(true_mspace_size, fmpzd, bitstr_len);
      fmpz_mod(modi, messages[i], true_mspace_size);
      fmpz_mod(modj, messages[j], true_mspace_size);
      int true_compare = fmpz_cmp(modi, modj);
      fmpz_clear(modi);
      fmpz_clear(modj);
      fmpz_clear(true_mspace_size);
      fmpz_clear(fmpzd);
      if(true_compare < 0) {
        true_compare = 1;
      } else if(true_compare > 0) {
        true_compare = 2;
      }
      if(verbose) {
        printf(
          "messages #1: %lu, message #2: %lu, compare: %d, true_compare: %d\n",
          fmpz_get_ui(messages[i]), fmpz_get_ui(messages[j]),
          compare, true_compare
        );
      }
      if(compare != true_compare) {
        status += 1;
      }
      free(cts);
    }
  }

  for(int i = 0; i < num_messages; i++) {
    mife_ciphertext_clear(pp2, ciphertexts[i]);
  }
  free(ciphertexts);
  free(messages);

  mife_clear_pp_read(pp2);
  mife_clear_sk(sk2);
  free(params);
  
  if(status == 0) {
    printf("SUCCESS\n");
  } else {
    printf("FAIL\n");
  }
  return status;
}

void test_rand() {
  aes_randstate_t gen_randstate;
  aes_randinit(gen_randstate);
  fmpz_t a;
  fmpz_t b;
  fmpz_init(a);
  fmpz_init_set_ui(b, 100008);
  fmpz_randm_aes(a, gen_randstate, b);
  //fmpz_print(a); printf("\n");
  fmpz_randm_aes(a, gen_randstate, b);
  //fmpz_print(a); printf("\n");
  fmpz_randm_aes(a, gen_randstate, b);
  //fmpz_print(a); printf("\n");
  aes_randclear(gen_randstate);
}

void run_tests() {
  fmpz_t p;
  fmpz_init_set_ui(p, 19);
  aes_randstate_t randstate;
  aes_randinit(randstate);
  test_matrix_inv(30, randstate, p);
  aes_randclear(randstate);
  
  test_rand();

  test_dary_conversion();
  //test_ore(5, 16, 5, 2, 2, ORE_MBP_NORMAL, 0);
  test_ore(5, 16, 5, 2, 3, ORE_MBP_DC, 0);
  // FIXME: test_ore(5, 16, 5, 2, 2, ORE_MBP_DC, 0); FAILS!
  /* 
  test_ore(5, 16, 5, 2, 4, ORE_MBP_MC, 0);
  test_ore(5, 1000, 10, 5, 5, ORE_MBP_NORMAL, 0);
  test_ore(5, 1000, 10, 5, 5, ORE_MBP_DC, 0);
  test_ore(5, 1000, 10, 5, 5, ORE_MBP_MC, 0);
  */

  mpfr_free_cache();
  flint_cleanup();

}

void print_random_matrices_with_inverse(int n, char *p_file,
    char *a_file, char *inv_file) {
  aes_randstate_t randstate;
  aes_randinit_seed(randstate, "62f1b2d2733a5e5fd63de8eff7fcf965");
  fmpz_t modp;
  fmpz_init(modp);
  FILE *fp = fopen(p_file, "r");
  fmpz_fread(fp, modp);
  fclose(fp);

  fmpz_mat_t a;
  fmpz_mat_init(a, n, n);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      fmpz_randm_aes(fmpz_mat_entry(a, i, j), randstate, modp);
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
  if (status == 0)
    printf("ERROR: Failed to produce the proper matrices.\n");

  fp = fopen(a_file, "w");
  fmpz_mat_fprint(fp, a);
  fclose(fp);

  fp = fopen(inv_file, "w");
  fmpz_mat_transpose(inv, inv);
  fmpz_mat_fprint(fp, inv);
  fclose(fp);

  fmpz_mat_clear(a);
  fmpz_mat_clear(inv);
  fmpz_mat_clear(prod);
  fmpz_mat_clear(identity);

  aes_randclear(randstate);

}

