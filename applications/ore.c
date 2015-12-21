#include "ore.h"

/**
 * TODO:
 * - add tests
 *
 *
 */

int main(int argc, char *argv[]) {
  cmdline_params_t cmdline_params;

  const char *name =  "Order Revealing Encryption";
  parse_cmdline(cmdline_params, argc, argv, name, NULL);
  flint_rand_t randstate;

/*
  get_best_params(80, 10, 11);
  exit(0);

  FILE *fp = fopen("enc_sizes.out", "w"); 
  gghlite_params_test_kappa_enc_size(80, 40, fp);
  fclose(fp);
  exit(0);
*/

  int bitstr_len = 4;
  int base = 2;
  int L = 80; // 2^L = # of total messages we can encrypt

  ore_pp_t pp;
  ore_sk_t sk;

  T = ggh_walltime(0);
  ore_setup(pp, sk, bitstr_len, base, L, cmdline_params);

  //test_matrix_inv(6, randstate, p);
  //test_dary_conversion();

  int num1 = 5;
  int num2 = 9;

  NUM_ENCODINGS_GENERATED = 0;

  ore_ciphertext_t ct1;
  T = ggh_walltime(0);
  ore_encrypt(ct1, num1, pp, sk);
  printf("3. Time it takes to create a single ciphertext: %8.2f s\n",
    ggh_seconds(ggh_walltime(T)));
  int num_encodings_generated = NUM_ENCODINGS_GENERATED; 

  T = ggh_walltime(0);
  fwrite_ore_ciphertext(pp, ct1, "ct1.out");
  printf("Time it takes to write the ciphertext to a file: %8.2f s\n",
    ggh_seconds(ggh_walltime(T)));
  ore_ciphertext_clear(pp, ct1);

  ore_ciphertext_t ct1_read;
  T = ggh_walltime(0);
  fread_ore_ciphertext(pp, ct1_read, "ct1.out");
  printf("Time it takes to read a ciphertext from a file: %8.2f s\n",
    ggh_seconds(ggh_walltime(T)));

  ore_ciphertext_t ct2, ct2_read;
  ore_encrypt(ct2, num2, pp, sk);
  fwrite_ore_ciphertext(pp, ct2, "ct2.out");
  ore_ciphertext_clear(pp, ct2);
  fread_ore_ciphertext(pp, ct2_read, "ct2.out");

  fwrite_ore_pp(pp, "pp.out");
  ore_clear_pp(pp);
  ore_pp_t pp_read;
  fread_ore_pp(pp_read, "pp.out");

  printf("Comparing %d with %d: ", num1, num2);
  
  T = ggh_walltime(0);
  compare(pp_read, ct1_read, ct2_read);

  printf("4. Time it takes to run comparison: %8.2f s\n",
      ggh_seconds(ggh_walltime(T)));

  printf("Number of encodings generated per ciphertext: %d\n",
      num_encodings_generated);

  ore_ciphertext_clear(pp, ct1_read);
  ore_ciphertext_clear(pp, ct2_read);
  ore_clear_pp_read(pp_read);
  ore_clear_sk(sk);
  mpfr_free_cache();
  flint_cleanup();
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
void get_best_params(int lambda, int message_d, int message_n) {
  mpfr_t dt; 
  mpfr_t nt;
  mpfr_t total;
  mpfr_init(dt);
  mpfr_init(nt);
  mpfr_init(total);

  mpfr_t tmp1;
  mpfr_t tmp2;
  mpfr_init(tmp1);
  mpfr_init(tmp2);
  mpfr_t bt;
  mpfr_init(bt);

  mpfr_set_ui(dt, message_d, MPFR_RNDN);
  mpfr_set_ui(nt, message_n, MPFR_RNDN);
  mpfr_pow(total, dt, nt, MPFR_RNDN);
  
  int max_base = 60;

  long dc_vals[max_base];
  long mc_vals[max_base];

  for(int base = 2; base < max_base; base++) {
    // compute minimum n
    mpfr_set_ui(bt, base, MPFR_RNDN);
    mpfr_log(tmp1, total, MPFR_RNDN);
    mpfr_log(tmp2, bt, MPFR_RNDN);
    mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);
    mpfr_ceil(tmp1, tmp1);
    unsigned long n = mpfr_get_ui(tmp1, MPFR_RNDN);
    
    // num encodings #1, #2
    int dc_enc = dc_num_enc(base,n);
    int mc_enc = mc_num_enc(base,n);
   
    // use n to compute enc size
    if(base == 4) {
      printf("Base 4 dc num enc: %d\n", dc_num_enc(4,n));
    }
    dc_vals[base] = dc_enc * dc_enc_size(n);
    mc_vals[base] = mc_enc * mc_enc_size(n);
  }

  mpfr_clear(dt);
  mpfr_clear(bt);
  mpfr_clear(nt);
  mpfr_clear(total);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);


  for(int i = 2; i < max_base; i++) {
    printf("%d: %ld %ld\n", i, dc_vals[i] / 8 / 1024 / 1024, mc_vals[i] / 8 / 1024 / 1024);
  }

  unsigned long min_enc = dc_vals[2]; 
  int min_type = 0;
  int min_base = 2;
  for(int i = 2; i < max_base; i++) {
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
    printf("Using ORE_MBP_DC\n");
  } else {
    printf("Using ORE_MBP_MC\n");
  }
  printf("min_enc: %lu, min_base: %d,\n", min_enc, min_base);
}

void gghlite_params_clear_read(gghlite_params_t self) {
  for(int i = 0; i < self->gamma; i++) {
    for(int j = 0; j < KAPPA; j++) {
      free(self->x[i][j]);
    }
    free(self->x[i]);
  }
  free(self->x);
  free(self->y);

  fmpz_mod_poly_clear(self->pzt);
  mpfr_clear(self->xi);
  mpfr_clear(self->sigma_s);
  mpfr_clear(self->ell_b);
  mpfr_clear(self->sigma_p);
  mpfr_clear(self->ell_g);
  mpfr_clear(self->sigma);
  fmpz_mod_poly_oz_ntt_precomp_clear(self->ntt);
  fmpz_clear(self->q);
  free(self);
}

/**
 * 
 * Members of pp that are not currently transferred:
 * pp->flags
 * params->rerand_mask
 * params->flags
 * params->x
 * params->y
 * params->D_sigma_p
 * params->D_sigma_s
 */
void fwrite_ore_pp(ore_pp_t pp, char *filepath) {
  int mpfr_base = 10;
  FILE *fp = fopen(filepath, "w");
  fprintf(fp, "%d %d %d %d %d %d %d %d %d\n",
    pp->bitstr_len,
    pp->d,
    pp->nx,
    pp->ny,
    pp->L,
    pp->gammax,
    pp->gammay,
    pp->kappa,
    pp->numR
  );
  fmpz_fprint(fp, pp->p);
  fprintf(fp, "\n");
  fprintf(fp, "%zd %zd %zd %ld %ld\n",
    (*pp->params_ref)->lambda,
    (*pp->params_ref)->gamma,
    (*pp->params_ref)->kappa,
    (*pp->params_ref)->n,
    (*pp->params_ref)->ell
  );
  fmpz_fprint(fp, (*pp->params_ref)->q);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, (*pp->params_ref)->sigma, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, (*pp->params_ref)->sigma_p, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, (*pp->params_ref)->sigma_s, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, (*pp->params_ref)->ell_b, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, (*pp->params_ref)->ell_g, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, (*pp->params_ref)->xi, MPFR_RNDN);
  fprintf(fp, "\n");
  gghlite_enc_fprint(fp, (*pp->params_ref)->pzt);
  fprintf(fp, "\n");
  fprintf(fp, "%zd\n", (*pp->params_ref)->ntt->n);
  fmpz_mod_poly_fprint(fp, (*pp->params_ref)->ntt->w);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint(fp, (*pp->params_ref)->ntt->w_inv);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint(fp, (*pp->params_ref)->ntt->phi);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint(fp, (*pp->params_ref)->ntt->phi_inv);
  fclose(fp);
}

void fread_ore_pp(ore_pp_t pp, char *filepath) {
  int mpfr_base = 10;
  FILE *fp = fopen(filepath, "r");
  fscanf(fp, "%d %d %d %d %d %d %d %d %d\n",
    &pp->bitstr_len,
    &pp->d,
    &pp->nx,
    &pp->ny,
    &pp->L,
    &pp->gammax,
    &pp->gammay,
    &pp->kappa,
    &pp->numR
  ) > 0;

  fmpz_init(pp->p);
  fmpz_fread(fp, pp->p);
  fscanf(fp, "\n") == 1;

  pp->params_ref = malloc(sizeof(gghlite_params_t));
  size_t lambda, kappa, gamma, n, ell;
  fscanf(fp, "%zd %zd %zd %ld %ld\n",
    &lambda,
    &gamma,
    &kappa,
    &n,
    &ell
  ) > 0;
  
  gghlite_params_initzero(*pp->params_ref, lambda, kappa, gamma);
  (*pp->params_ref)->n = n;
  (*pp->params_ref)->ell = ell;

  fmpz_fread(fp, (*pp->params_ref)->q);
  fscanf(fp, "\n") == 1;
  mpfr_inp_str((*pp->params_ref)->sigma, fp, mpfr_base, MPFR_RNDN);
  fscanf(fp, "\n") == 1;
  mpfr_inp_str((*pp->params_ref)->sigma_p, fp, mpfr_base, MPFR_RNDN);
  fscanf(fp, "\n") == 1;
  mpfr_inp_str((*pp->params_ref)->sigma_s, fp, mpfr_base, MPFR_RNDN);
  fscanf(fp, "\n") == 1;
  mpfr_inp_str((*pp->params_ref)->ell_b, fp, mpfr_base, MPFR_RNDN);
  fscanf(fp, "\n") == 1;
  mpfr_inp_str((*pp->params_ref)->ell_g, fp, mpfr_base, MPFR_RNDN);
  fscanf(fp, "\n") == 1;
  mpfr_inp_str((*pp->params_ref)->xi, fp, mpfr_base, MPFR_RNDN);
  fscanf(fp, "\n") == 1;
  
  gghlite_enc_init((*pp->params_ref)->pzt, *pp->params_ref);
  gghlite_enc_init((*pp->params_ref)->ntt->w, *pp->params_ref);
  gghlite_enc_init((*pp->params_ref)->ntt->w_inv, *pp->params_ref);
  gghlite_enc_init((*pp->params_ref)->ntt->phi, *pp->params_ref);
  gghlite_enc_init((*pp->params_ref)->ntt->phi_inv, *pp->params_ref);

  gghlite_enc_fread(fp, (*pp->params_ref)->pzt);
  fscanf(fp, "\n") == 1;
  fscanf(fp, "%zd\n", &(*pp->params_ref)->ntt->n) == 1;
  gghlite_enc_fread(fp, (*pp->params_ref)->ntt->w);
  fscanf(fp, "\n") == 1;
  gghlite_enc_fread(fp, (*pp->params_ref)->ntt->w_inv);
  fscanf(fp, "\n") == 1;
  gghlite_enc_fread(fp, (*pp->params_ref)->ntt->phi);
  fscanf(fp, "\n") == 1;
  gghlite_enc_fread(fp, (*pp->params_ref)->ntt->phi_inv);
  fclose(fp);
}

void fwrite_ore_ciphertext(ore_pp_t pp, ore_ciphertext_t ct, char *filepath) {
  FILE *fp = fopen(filepath, "w");
  for(int i = 0; i < pp->nx; i++) {
    fwrite_gghlite_enc_mat(pp, ct->x_enc[i], fp);
  }
  for(int i = 0; i < pp->ny; i++) {
    fwrite_gghlite_enc_mat(pp, ct->y_enc[i], fp);
  }
  fclose(fp);
}

void fwrite_gghlite_enc_mat(ore_pp_t pp, gghlite_enc_mat_t m, FILE *fp) {
  fprintf(fp, " %d ", m->nrows);
  fprintf(fp, " %d ", m->ncols);
  for(int i = 0; i < m->nrows; i++) {
    for(int j = 0; j < m->ncols; j++) {
      gghlite_enc_fprint(fp, m->m[i][j]);
      fprintf(fp, "\n");
    }
  }
}

void fread_ore_ciphertext(ore_pp_t pp, ore_ciphertext_t ct, char *filepath) {
  FILE *fp = fopen(filepath, "r");
  ct->x_enc = malloc(pp->nx * sizeof(gghlite_enc_mat_t));
  ct->y_enc = malloc(pp->ny * sizeof(gghlite_enc_mat_t));
  assert(ct->x_enc && ct->y_enc);
  for(int i = 0; i < pp->nx; i++) {
    fread_gghlite_enc_mat(pp, ct->x_enc[i], fp);
  }
  for(int i = 0; i < pp->ny; i++) {
    fread_gghlite_enc_mat(pp, ct->y_enc[i], fp);
  }
  fclose(fp);
}

void fread_gghlite_enc_mat(ore_pp_t pp, gghlite_enc_mat_t m, FILE *fp) {
  int check1 = fscanf(fp, " %d ", &m->nrows);
  int check2 = fscanf(fp, " %d ", &m->ncols);
  assert(check1 == 1 && check2 == 1);
  m->m = malloc(m->nrows * sizeof(gghlite_enc_t *));
  assert(m->m);
  for(int i = 0; i < m->nrows; i++) {
    m->m[i] = malloc(m->ncols * sizeof(gghlite_enc_t));
    assert(m->m[i]);
    for(int j = 0; j < m->ncols; j++) {
      gghlite_enc_init(m->m[i][j], *pp->params_ref);
      int check3 = gghlite_enc_fread(fp, m->m[i][j]);
      assert(check3 == 1);
      fscanf(fp, "\n") == 1;
    }
  }
}

void ore_ciphertext_clear(ore_pp_t pp, ore_ciphertext_t ct) {
  for(int i = 0; i < pp->nx; i++) {
    gghlite_enc_mat_clear(ct->x_enc[i]);
  }
  for(int i = 0; i < pp->ny; i++) {
    gghlite_enc_mat_clear(ct->y_enc[i]);
  }
  free(ct->x_enc);
  free(ct->y_enc);
}

void ore_clear_pp_read(ore_pp_t pp) {
  fmpz_clear(pp->p);
  gghlite_params_clear_read(*pp->params_ref);
}

void ore_clear_pp(ore_pp_t pp) {
  fmpz_clear(pp->p);
}

void ore_clear_sk(ore_sk_t sk) {
  gghlite_sk_clear(sk->self, 1);
  for(int i = 0; i < sk->numR; i++) {
    fmpz_mat_clear(sk->R[i]);
    fmpz_mat_clear(sk->R_inv[i]);
  }
  free(sk->R);
  free(sk->R_inv);
  flint_randclear(sk->randstate);
}

void ore_mat_clr_clear(ore_pp_t pp, ore_mat_clr_t met) {
  free(met->dary_repr);
  for(int i = 0; i < pp->nx; i++) {
    fmpz_mat_clear(met->x_clr[i]);
  }
  free(met->x_clr);
  for(int i = 0; i < pp->ny; i++) {
    fmpz_mat_clear(met->y_clr[i]);
  }
  free(met->y_clr);
}

void ore_encrypt(ore_ciphertext_t ct, int message, ore_pp_t pp, ore_sk_t sk) {
  int index = 0; // FIXME choose it randomly using sk->randstate
  ore_mat_clr_t met;
  set_matrices(met, message, pp, sk);
  if(! (pp->flags & ORE_NO_RANDOMIZERS)) {
    apply_scalar_randomizers(met, pp, sk);     
  }
  set_encodings(ct, met, index, pp, sk);
  ore_mat_clr_clear(pp, met);
}

void ore_setup(ore_pp_t pp, ore_sk_t sk, int bitstr_len, int base,
    int L, cmdline_params_t cmdline_params) {
  flint_randinit_seed(sk->randstate, cmdline_params->seed, 1);
  pp->flags = ORE_DEFAULT;
  pp->d = base;
  pp->L = L;
  pp->bitstr_len = bitstr_len;

  if(pp->flags & ORE_MBP_NORMAL) {
    pp->kappa = 2 * pp->bitstr_len;
    pp->nx = pp->bitstr_len;
    pp->ny = pp->bitstr_len;
  } else if(pp->flags & ORE_MBP_DC) {
    pp->kappa = 1 + pp->bitstr_len;
    pp->nx = pp->bitstr_len / 2 + 1;
    pp->ny = (pp->bitstr_len+1) / 2;
  } else if(pp->flags & ORE_MBP_MC) {
    pp->kappa = 2 * pp->bitstr_len;
    pp->nx = pp->bitstr_len;
    pp->ny = pp->bitstr_len;
  }

  pp->gammax = 1 + (pp->nx-1) * (pp->L+1);
  pp->gammay = 1 + (pp->ny-1) * (pp->L+1);
  pp->gamma = pp->gammax + pp->gammay;

  T = ggh_walltime(0);
  gghlite_jigsaw_init_gamma(sk->self,
                      cmdline_params->lambda,
                      pp->kappa,
                      pp->gamma,
                      cmdline_params->flags,
                      sk->randstate);
  printf("\n");
  gghlite_params_print(sk->self->params);
  printf("\n---\n\n");

  pp->params_ref = &(sk->self->params);

  printf("Supporting at most 2^%d plaintexts, each in base %d,\n", pp->L,
      pp->d);
  printf("of length %d, with gamma = %d\n\n", pp->bitstr_len, pp->gamma);

  fmpz_init(pp->p);
  fmpz_poly_oz_ideal_norm(pp->p, sk->self->g, sk->self->params->n, 0);

  printf("1. GGH InstGen wall time:                 %8.2f s\n",
      ggh_seconds(ggh_walltime(T)));

  T = ggh_walltime(0);

  // set the kilian randomizers in sk
  pp->numR = pp->kappa - 1;
  sk->numR = pp->numR;
  int *dims = malloc(pp->numR * sizeof(int));
  assert(dims);

  // set the kilian dimensions
  if(pp->flags & ORE_MBP_NORMAL) {
    for(int i = 0; i < pp->numR; i++) {
      dims[i] = pp->d+3;
    }
  } else if(pp->flags & ORE_MBP_DC) {
    dims[0] = pp->d;
    for(int i = 1; i < pp->numR; i++) {
      dims[i] = pp->d+2;
    }
  } else if(pp->flags & ORE_MBP_MC) {
    dims[0] = pp->d;
    for(int i = 1; i < pp->numR; i++) {
      if(i % 2 == 1) {
        dims[i] = 3;
      } else {
        dims[i] = pp->d+2;
      }
    }
  } else {
    assert(0);
  }

  sk->R = malloc(sk->numR * sizeof(fmpz_mat_t));
  sk->R_inv = malloc(sk->numR * sizeof(fmpz_mat_t));

  for (int k = 0; k < pp->numR; k++) {
    fmpz_mat_init(sk->R[k], dims[k], dims[k]);
    for (int i = 0; i < dims[k]; i++) {
      for(int j = 0; j < dims[k]; j++) {
        fmpz_randm(fmpz_mat_entry(sk->R[k], i, j), sk->randstate, pp->p);
      }
    }
	
    fmpz_mat_init(sk->R_inv[k], dims[k], dims[k]);
    fmpz_modp_matrix_inverse(sk->R_inv[k], sk->R[k], dims[k], pp->p);
  }

  free(dims);

  printf("2. Time it takes to generate the ORE secret key:    %8.2f s\n",
      ggh_seconds(ggh_walltime(T)));
}

// message >= 0, d >= 2
void message_to_dary(int *dary, int bitstring_len, int64_t message, int64_t d) {
  assert(message >= 0);
  assert(d >= 2);

  int i;
  for (i = bitstring_len - 1; i >= 0; i--) {
    dary[i] = message % d;
    message /= d;
  }
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
 
  for(int i = 0; i < 3; i++) {
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
 * Creates the first matrix in the degree-compressed version of ORE.
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
 * Creates the second matrix in the degree-compressed version of ORE.
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
 */
void ore_dc_clrmat_init_MIDDLE(fmpz_mat_t m, int input1, int input2, int d) {
  fmpz_mat_init(m, d+2, d+2);
  fmpz_mat_zero(m);
  
  for(int i = 0; i < d; i++) {
    int j = -1;
    if(i < input1) {
      j = d;
    } else if(i > input1) {
      j = d+1;
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
 */
void ore_dc_clrmat_init_LAST(fmpz_mat_t m, int input, int d) {
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
    } else {
      assert(i > input);
      j = 2;
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
    int j;
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
int get_matrix_bit_normal_mbp(int input, int i, int j, int type) {
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

void fmpz_mat_scalar_mul_modp(fmpz_mat_t m, fmpz_t scalar, fmpz_t modp) {
  for (int i = 0; i < m->r; i++) {
    for(int j = 0; j < m->c; j++) {
      fmpz_mul(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), scalar);
      fmpz_mod(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), modp);
    }
  }	}

/* sets the cleartext matrices x_clr and y_clr */
void set_matrices(ore_mat_clr_t met, int64_t message, ore_pp_t pp,
    ore_sk_t sk) {
  met->dary_repr = malloc(pp->bitstr_len * sizeof(int));
  message_to_dary(met->dary_repr, pp->bitstr_len, message, pp->d);

  met->x_clr = malloc(pp->nx * sizeof(fmpz_mat_t));
  met->y_clr = malloc(pp->ny * sizeof(fmpz_mat_t));

  if(pp->flags & ORE_MBP_NORMAL) {
    assert(pp->nx == pp->ny);
    for (int k = 0; k < pp->nx; k++) {
      int dim = pp->d+3;
      fmpz_mat_init(met->x_clr[k], dim, dim);
      fmpz_mat_init(met->y_clr[k], dim, dim);

      for (int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          int x_digit = get_matrix_bit_normal_mbp(met->dary_repr[k],
              i, j, X_TYPE);
          int y_digit = get_matrix_bit_normal_mbp(met->dary_repr[k],
              i, j, Y_TYPE);
          fmpz_set_ui(fmpz_mat_entry(met->x_clr[k], i, j), x_digit);
          fmpz_set_ui(fmpz_mat_entry(met->y_clr[k], i, j), y_digit);
        }
      }
    }
  } else if(pp->flags & ORE_MBP_DC) {

    for(int k = 0, bc = 0; k < pp->nx; k++, bc++) {
      if(bc == 0) {
        ore_dc_clrmat_init_FIRST(met->x_clr[k], met->dary_repr[bc], pp->d);
      } else if(bc == pp->bitstr_len - 1) {
        ore_dc_clrmat_init_LAST(met->x_clr[k], met->dary_repr[bc], pp->d);
      } else {
        ore_dc_clrmat_init_MIDDLE(met->x_clr[k], met->dary_repr[bc],
            met->dary_repr[bc+1], pp->d);
        bc++;
      }
    }
    
    for(int k = 0, bc = 0; k < pp->ny; k++, bc++) {
      if(k == 0 && pp->ny > 1) {
        ore_dc_clrmat_init_SECOND(met->y_clr[k], met->dary_repr[bc],
            met->dary_repr[bc+1], pp->d);
        bc++;
      } else if(k == 0 && pp->ny == 1) {
        ore_dc_clrmat_init_SECONDANDLAST(met->y_clr[k], met->dary_repr[bc],
            pp->d);
      } else if((pp->bitstr_len % 2 == 1) && (k == pp->ny-1)) {
        ore_dc_clrmat_init_LAST(met->x_clr[k], met->dary_repr[bc], pp->d);
      } else {
        ore_dc_clrmat_init_MIDDLE(met->y_clr[k], met->dary_repr[bc],
            met->dary_repr[bc+1], pp->d);
        bc++;
      }
    }
  } else if(pp->flags & ORE_MBP_MC) {
    for(int k = 0; k < pp->nx; k++) {
      if(k == 0) {
        ore_mc_clrmat_init_XFIRST(met->x_clr[k], met->dary_repr[k], pp->d);
      } else {
        ore_mc_clrmat_init_XREST(met->x_clr[k], met->dary_repr[k], pp->d);
      }
    }
    for(int k = 0; k < pp->ny; k++) {
      if(k == 0) {
        ore_mc_clrmat_init_YFIRST(met->y_clr[k], met->dary_repr[k], pp->d);
      } else {
        ore_mc_clrmat_init_YREST(met->y_clr[k], met->dary_repr[k], pp->d);
      }
    }
  } else {
    assert(0);
  }
    
}

void apply_scalar_randomizers(ore_mat_clr_t met, ore_pp_t pp, ore_sk_t sk) {
  for(int k = 0; k < pp->nx; k++) {
    fmpz_t x_rand;
    fmpz_init(x_rand);
    fmpz_randm(x_rand, sk->randstate, pp->p);
    fmpz_mat_scalar_mul_modp(met->x_clr[k], x_rand, pp->p);
    fmpz_clear(x_rand);
  }

  for(int k = 0; k < pp->ny; k++) {
    fmpz_t y_rand;
    fmpz_init(y_rand);
    fmpz_randm(y_rand, sk->randstate, pp->p);
    fmpz_mat_scalar_mul_modp(met->y_clr[k], y_rand, pp->p);
    fmpz_clear(y_rand);
  }
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
 * @param nu The number of total elements to be multiplied. partitioning[] 
 * will describe a nu-partition of the universe set.
 */ 
void gen_partitioning(int *partitioning, int i, int L, int nu) {
  int j = 0;

  int bitstring[L];
  message_to_dary(bitstring, L, i, 2);

  for(; j < nu; j++) {
    partitioning[j] = j;
  }

  for(int k = 0; k < L; k++) {
    for(int j1 = 1; j1 < nu; j1++) {
      partitioning[j] = (bitstring[k] == 1) ? j1 : 0;
      j++;
    }
  }
}

void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n,
    fmpz_t p) {
  fmpz_mat_mul(a, b, c);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      fmpz_mod(fmpz_mat_entry(a, i, j), fmpz_mat_entry(a, i, j), p);
    }
  }
}

void gghlite_enc_mat_init(gghlite_params_t params, gghlite_enc_mat_t m,
    int nrows, int ncols) {
  m->nrows = nrows;
  m->ncols = ncols;
  m->m = malloc(nrows * sizeof(gghlite_enc_t *));
  assert(m->m);
  for(int i = 0; i < m->nrows; i++) {
    m->m[i] = malloc(m->ncols * sizeof(gghlite_enc_t));
    assert(m->m[i]);
    for(int j = 0; j < m->ncols; j++) {
      gghlite_enc_init(m->m[i][j], params);
    }
  }
}

void gghlite_enc_mat_clear(gghlite_enc_mat_t m) {
  for(int i = 0; i < m->nrows; i++) {
    for(int j = 0; j < m->ncols; j++) {
      gghlite_enc_clear(m->m[i][j]);
    }
    free(m->m[i]);
  }
  free(m->m);
}


void mat_encode(ore_sk_t sk, gghlite_enc_mat_t enc, fmpz_mat_t m, int *group) {
  gghlite_clr_t e;
  gghlite_clr_init(e);
  for(int i = 0; i < enc->nrows; i++) {
    for(int j = 0; j < enc->ncols; j++) {
      fmpz_poly_set_coeff_fmpz(e, 0, fmpz_mat_entry(m, i, j));
      gghlite_enc_set_gghlite_clr(enc->m[i][j], sk->self, e, 1, group, 1,
          sk->randstate);
      NUM_ENCODINGS_GENERATED++;
    }
  }
  gghlite_clr_clear(e);
}

void gghlite_enc_mat_zeros_print(ore_pp_t pp, gghlite_enc_mat_t m) {
  for(int i = 0; i < m->nrows; i++) {
    printf("[");
    for(int j = 0; j < m->ncols; j++) {
      printf(gghlite_enc_is_zero(*pp->params_ref, m->m[i][j]) ? "0 " : "x " );
    }
    printf("]\n");
  }
}

void set_encodings(ore_ciphertext_t ct, ore_mat_clr_t met, int index,
    ore_pp_t pp, ore_sk_t sk) {

  int *ptnx = malloc(pp->gammax * sizeof(int));
  int *ptny = malloc(pp->gammay * sizeof(int));
  gen_partitioning(ptnx, index, pp->L, pp->nx);
  gen_partitioning(ptny, index, pp->L, pp->ny);


  /* construct the partitions in the group array form */
  int **group_x = malloc(pp->nx * sizeof(int *));
  int **group_y = malloc(pp->ny * sizeof(int *));
  for(int j = 0; j < pp->nx; j++) {
    group_x[j] = malloc(pp->gamma * sizeof(int));
    memset(group_x[j], 0, pp->gamma * sizeof(int));
    for(int k = 0; k < pp->gammax; k++) {
      if(ptnx[k] == j) {
        group_x[j][k] = 1;
      }
    }
  }

  for(int j = 0; j < pp->ny; j++) {
    group_y[j] = malloc(pp->gamma * sizeof(int));
    memset(group_y[j], 0, pp->gamma * sizeof(int));
    for(int k = 0; k < pp->gammay; k++) {
      if(ptny[k] == j) {
        group_y[j][k + pp->gammax] = 1;
      }
    }
  }

  free(ptnx);
  free(ptny);

  if(pp->flags & ORE_SIMPLE_PARTITIONS) {
    // override group arrays with trivial partitioning
    for(int j = 0; j < pp->nx; j++) {
      memset(group_x[j], 0, pp->gamma * sizeof(int));
    }
    for(int j = 0; j < pp->ny; j++) {
      memset(group_y[j], 0, pp->gamma * sizeof(int));
    }
    for(int k = 0; k < pp->gamma; k++) {
      group_x[0][k] = 1;
    }
  }

  if(! (pp->flags & ORE_NO_KILIAN)) {
    // apply kilian to the cleartext matrices (overwriting them in the process)
    fmpz_mat_t tmp;
    fmpz_mat_init(tmp, met->x_clr[0]->r, sk->R[0]->c);
    fmpz_mat_mul(tmp, met->x_clr[0], sk->R[0]);
    fmpz_mat_set(met->x_clr[0], tmp);
    fmpz_mat_clear(tmp);
   
    for(int j = 1; j < pp->nx; j++) {
      fmpz_mat_init(tmp, sk->R_inv[2 * j - 1]->r, met->x_clr[j]->c);
      fmpz_mat_mul(tmp, sk->R_inv[2 * j - 1], met->x_clr[j]);
      if(2 * j < pp->numR) {
        fmpz_mat_mul(tmp, tmp, sk->R[2 * j]);
      }
      fmpz_mat_set(met->x_clr[j], tmp);
      fmpz_mat_clear(tmp);
    }

    for(int j = 0; j < pp->ny; j++) {
      fmpz_mat_init(tmp, sk->R_inv[2 * j]->r, met->y_clr[j]->c);
      fmpz_mat_mul(tmp, sk->R_inv[2 * j], met->y_clr[j]);
      if(2 * j + 1 < pp->numR) {
        fmpz_mat_mul(tmp, tmp, sk->R[2 * j + 1]);
      }
      fmpz_mat_set(met->y_clr[j], tmp);
      fmpz_mat_clear(tmp);
    }
  }

  // encode
  ct->x_enc = malloc(pp->nx * sizeof(gghlite_enc_mat_t));
  ct->y_enc = malloc(pp->ny * sizeof(gghlite_enc_mat_t));
  assert(ct->x_enc && ct->y_enc);
  for(int j = 0; j < pp->nx; j++) {
    gghlite_enc_mat_init(sk->self->params, ct->x_enc[j],
        met->x_clr[j]->r, met->x_clr[j]->c);
    mat_encode(sk, ct->x_enc[j], met->x_clr[j], group_x[j]);
  }
  for(int j = 0; j < pp->ny; j++) {
    gghlite_enc_mat_init(sk->self->params, ct->y_enc[j],
        met->y_clr[j]->r, met->y_clr[j]->c);
    mat_encode(sk, ct->y_enc[j], met->y_clr[j], group_y[j]);
  }
  
  // free group arrays
  for(int j = 0; j < pp->nx; j++) {
    free(group_x[j]);
  }
  free(group_x);
  for(int j = 0; j < pp->ny; j++) {
    free(group_y[j]);
  }
  free(group_y);
}

void gghlite_enc_mat_mul(gghlite_params_t params, gghlite_enc_mat_t r,
    gghlite_enc_mat_t m1, gghlite_enc_mat_t m2) {
  gghlite_enc_t tmp;
  gghlite_enc_init(tmp, params);

  gghlite_enc_mat_t tmp_mat;
  gghlite_enc_mat_init(params, tmp_mat, m1->nrows, m2->ncols);

  assert(m1->ncols == m2->nrows);

  for(int i = 0; i < m1->nrows; i++) {
    for(int j = 0; j < m2->ncols; j++) {
      for(int k = 0; k < m1->ncols; k++) {
        gghlite_enc_mul(tmp, params, m1->m[i][k], m2->m[k][j]);
        gghlite_enc_add(tmp_mat->m[i][j], params, tmp_mat->m[i][j], tmp);
      }
    }
  }

  gghlite_enc_mat_clear(r);
  gghlite_enc_mat_init(params, r, m1->nrows, m2->ncols);

  for(int i = 0; i < r->nrows; i++) {
    for(int j = 0; j < r->ncols; j++) {
      gghlite_enc_set(r->m[i][j], tmp_mat->m[i][j]);
    }
  }

  gghlite_enc_mat_clear(tmp_mat);
  gghlite_enc_clear(tmp);
}

void compare(ore_pp_t pp, ore_ciphertext_t ct1, ore_ciphertext_t ct2) {
  gghlite_enc_mat_t tmp;
  gghlite_enc_mat_init(*pp->params_ref, tmp,
      ct1->x_enc[0]->nrows, ct2->y_enc[0]->ncols);
  gghlite_enc_mat_mul(*pp->params_ref, tmp, ct1->x_enc[0], ct2->y_enc[0]);

  for(int i = 0, xc = 1, yc = 1; i < pp->nx + pp->ny - 2; i++) {
    if(i % 2 == 0) {
      // multiply tmp by x
      gghlite_enc_mat_mul(*pp->params_ref, tmp, tmp, ct1->x_enc[xc]);
      xc++;
    } else {
      // multiply tmp by y
      gghlite_enc_mat_mul(*pp->params_ref, tmp, tmp, ct2->y_enc[yc]);
      yc++;
    }
  } 
	
  int equals = 1 - gghlite_enc_is_zero(*pp->params_ref, tmp->m[0][0]);
  int lessthan = 1 - gghlite_enc_is_zero(*pp->params_ref, tmp->m[0][1]);
  int greaterthan = 1 - gghlite_enc_is_zero(*pp->params_ref, tmp->m[0][2]);
  gghlite_enc_mat_zeros_print(pp, tmp);
  gghlite_enc_mat_clear(tmp);

  assert((equals + lessthan + greaterthan) == 1);

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

int int_arrays_equal(int *arr1, int *arr2, int length) {
  for (int i = 0; i < length; i++) {
    if (arr1[i] != arr2[i])
      return 1;
  }
  return 0;
}

void test_dary_conversion() {
  printf("Testing d-ary conversion function...                          ");
  int dary1[4];
  int dary2[8];
  int dary3[4];
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

  fmpz_mat_clear(a);
  fmpz_mat_clear(inv);
  fmpz_mat_clear(prod);
  fmpz_mat_clear(identity);
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

/* copied from fmpz_mod_poly_fread, mostly (but fixed) */
int gghlite_enc_fread(FILE * f, fmpz_mod_poly_t poly)
{
    slong i, length;
    fmpz_t coeff;
    ulong res;

    fmpz_init(coeff);
    if (flint_fscanf(f, "%wd", &length) != 1) {
        fmpz_clear(coeff);
        return 0;
    }

    fmpz_fread(f,coeff);
    fmpz_mod_poly_clear(poly);
    fmpz_mod_poly_init(poly, coeff);
    fmpz_mod_poly_fit_length(poly, length);

    poly->length = length;
    flint_fscanf(f, " ");
  
    for (i = 0; i < length; i++)
    {
        flint_fscanf(f, " ");
        res = fmpz_fread(f, coeff);
        fmpz_mod_poly_set_coeff_fmpz(poly,i,coeff);

        if (!res)
        {
            poly->length = i;
            fmpz_clear(coeff);
            return 0;
        }
    }

    fmpz_clear(coeff);
    _fmpz_mod_poly_normalise(poly);

    return 1;
}

