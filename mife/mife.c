#include "mife.h"

#define CHECK(x) if(x < 0) { assert(0); }

/**
 * sets exp = base^n, where exp is an mpfr_t
 */
void fmpz_init_exp(fmpz_t exp, int base, int n) {
  fmpz_init(exp);
  fmpz_t tmp;
  fmpz_init_set_ui(tmp, base);
  fmpz_pow_ui(exp, tmp, n);
  fmpz_clear(tmp);
}

void mife_init_params(mife_pp_t pp, int bitstr_len, mife_flag_t flags) {
  pp->flags = flags;
  pp->bitstr_len = bitstr_len;
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

void fwrite_gghlite_params(FILE *fp, gghlite_params_t params) {
  int mpfr_base = 10;
  fprintf(fp, "%zd %zd %zd %ld %ld %lu %d\n",
    params->lambda,
    params->gamma,
    params->kappa,
    params->n,
    params->ell,
    params->rerand_mask,
    params->flags
  );
  fmpz_fprint(fp, params->q);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, params->sigma, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, params->sigma_p, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, params->sigma_s, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, params->ell_b, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, params->ell_g, MPFR_RNDN);
  fprintf(fp, "\n");
  mpfr_out_str(fp, mpfr_base, 0, params->xi, MPFR_RNDN);
  fprintf(fp, "\n");
  gghlite_enc_fprint(fp, params->pzt);
  fprintf(fp, "\n");
  fprintf(fp, "%zd\n", params->ntt->n);
  fmpz_mod_poly_fprint(fp, params->ntt->w);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint(fp, params->ntt->w_inv);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint(fp, params->ntt->phi);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint(fp, params->ntt->phi_inv);
}

/**
 * 
 * Members of pp that are not currently transferred:
 * params->x
 * params->y
 * params->D_sigma_p
 * params->D_sigma_s
 */
void fwrite_mife_pp(mife_pp_t pp, char *filepath) {
  FILE *fp = fopen(filepath, "w");
  fprintf(fp, "%d %d %d %d %d %d %d\n",
    pp->num_inputs,
    pp->bitstr_len,
    pp->L,
    pp->gamma,
    pp->kappa,
    pp->numR,
    pp->flags
  );
  for(int i = 0; i < pp->num_inputs; i++) {
    fprintf(fp, "%d ", pp->n[i]);
  }
  fprintf(fp, "\n");
  for(int i = 0; i < pp->num_inputs; i++) {
    fprintf(fp, "%d ", pp->gammas[i]);
  }
  fprintf(fp, "\n");
  fmpz_fprint(fp, pp->p);
  fprintf(fp, "\n");

  fwrite_gghlite_params(fp, *pp->params_ref);
  fclose(fp);
}

void fwrite_mife_sk(mife_sk_t sk, char *filepath) {
  FILE *fp = fopen(filepath, "w");
  fprintf(fp, "%d\n", sk->numR);
  for(int i = 0; i < sk->numR; i++) {
    fprintf(fp, "%ld %ld\n", sk->R[i]->r, sk->R[i]->c);
    fmpz_mat_fprint(fp, sk->R[i]);
    fprintf(fp, "\n");
    fprintf(fp, "%ld %ld\n", sk->R_inv[i]->r, sk->R_inv[i]->c);
    fmpz_mat_fprint(fp, sk->R_inv[i]);
    fprintf(fp, "\n");
  }
  fwrite_gghlite_params(fp, sk->self->params);
  fprintf(fp, "\n");
  fmpz_poly_fprint(fp, sk->self->g);
  fprintf(fp, "\n");
  fmpq_poly_fprint(fp, sk->self->g_inv);
  fprintf(fp, "\n");
  fmpz_poly_fprint(fp, sk->self->h);
  fprintf(fp, "\n");
  for(int i = 0; i < sk->self->params->gamma; i++) {
    fmpz_fprint(fp, fmpz_mod_poly_modulus(sk->self->z[i]));
    fprintf(fp, "\n");
    fmpz_mod_poly_fprint(fp, sk->self->z[i]);
    fprintf(fp, "\n");
    fmpz_fprint(fp, fmpz_mod_poly_modulus(sk->self->z_inv[i]));
    fprintf(fp, "\n");
    fmpz_mod_poly_fprint(fp, sk->self->z_inv[i]);
    fprintf(fp, "\n");
    fmpz_poly_fprint(fp, sk->self->a[i]);
    fprintf(fp, "\n");
    for(int j = 0; j < sk->self->params->kappa; j++) {
      fmpz_poly_fprint(fp, sk->self->b[i][j][0]);
      fprintf(fp, "\n");
      fmpz_poly_fprint(fp, sk->self->b[i][j][1]);
      fprintf(fp, "\n");
    }
  }

  fclose(fp);

}

void fread_mife_sk(mife_sk_t sk, char *filepath) {
  FILE *fp = fopen(filepath, "r");
  CHECK(fscanf(fp, "%d\n", &sk->numR));
  sk->R = malloc(sk->numR * sizeof(fmpz_mat_t));
  sk->R_inv = malloc(sk->numR * sizeof(fmpz_mat_t));

  for(int i = 0; i < sk->numR; i++) {
    unsigned long r1, c1, r2, c2;
    CHECK(fscanf(fp, "%ld %ld\n", &r1, &c1));
    fmpz_mat_init(sk->R[i], r1, c1);
    fmpz_mat_fread(fp, sk->R[i]);
    CHECK(fscanf(fp, "\n"));
    CHECK(fscanf(fp, "%ld %ld\n", &r2, &c2));
    fmpz_mat_init(sk->R_inv[i], r2, c2);
    fmpz_mat_fread(fp, sk->R_inv[i]);
    CHECK(fscanf(fp, "\n"));
  }
  fread_gghlite_params(fp, sk->self->params);
  CHECK(fscanf(fp, "\n"));
  fmpz_poly_init(sk->self->g);
  fmpz_poly_fread(fp, sk->self->g);
  CHECK(fscanf(fp, "\n"));
  fmpq_poly_init(sk->self->g_inv);
  fmpq_poly_fread(fp, sk->self->g_inv);
  CHECK(fscanf(fp, "\n"));
  fmpz_poly_init(sk->self->h);
  fmpz_poly_fread(fp, sk->self->h);
  CHECK(fscanf(fp, "\n"));
  
  sk->self->z = malloc(sk->self->params->gamma * sizeof(gghlite_enc_t));
  sk->self->z_inv = malloc(sk->self->params->gamma * sizeof(gghlite_enc_t));
  sk->self->a = malloc(sk->self->params->gamma * sizeof(gghlite_clr_t));
  sk->self->b = malloc(sk->self->params->gamma * sizeof(gghlite_clr_t **));
  for(int i = 0; i < sk->self->params->gamma; i++) {
    fmpz_t p1, p2;
    fmpz_init(p1);
    fmpz_init(p2);
    fmpz_fread(fp, p1);
    CHECK(fscanf(fp, "\n"));
    fmpz_mod_poly_init(sk->self->z[i], p1);
    fmpz_mod_poly_fread(fp, sk->self->z[i]);
    CHECK(fscanf(fp, "\n"));
    fmpz_fread(fp, p2);
    CHECK(fscanf(fp, "\n"));
    fmpz_mod_poly_init(sk->self->z_inv[i], p2);
    fmpz_mod_poly_fread(fp, sk->self->z_inv[i]);
    CHECK(fscanf(fp, "\n"));
    fmpz_poly_init(sk->self->a[i]);
    fmpz_poly_fread(fp, sk->self->a[i]);
    CHECK(fscanf(fp, "\n"));
    sk->self->b[i] = malloc(sk->self->params->kappa * sizeof(gghlite_clr_t *));
    for(int j = 0; j < sk->self->params->kappa; j++) {
      sk->self->b[i][j] = malloc(2 * sizeof(gghlite_clr_t));
      fmpz_poly_init(sk->self->b[i][j][0]);
      fmpz_poly_fread(fp, sk->self->b[i][j][0]);
      CHECK(fscanf(fp, "\n"));
      fmpz_poly_init(sk->self->b[i][j][1]);
      fmpz_poly_fread(fp, sk->self->b[i][j][1]);
      CHECK(fscanf(fp, "\n"));
    }
    fmpz_clear(p1);
    fmpz_clear(p2);
  }

  gghlite_sk_set_D_g(sk->self);
  fclose(fp);
}


void fread_gghlite_params(FILE *fp, gghlite_params_t params) {
  int mpfr_base = 10;
  size_t lambda, kappa, gamma, n, ell;
  uint64_t rerand_mask;
  int gghlite_flag_int;
  CHECK(fscanf(fp, "%zd %zd %zd %ld %ld %lu %d\n",
    &lambda,
    &gamma,
    &kappa,
    &n,
    &ell,
    &rerand_mask,
    &gghlite_flag_int
  ));
  
  gghlite_params_initzero(params, lambda, kappa, gamma);
  params->n = n;
  params->ell = ell;
  params->rerand_mask = rerand_mask;
  params->flags = gghlite_flag_int;

  fmpz_fread(fp, params->q);
  CHECK(fscanf(fp, "\n"));
  mpfr_inp_str(params->sigma, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"));
  mpfr_inp_str(params->sigma_p, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"));
  mpfr_inp_str(params->sigma_s, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"));
  mpfr_inp_str(params->ell_b, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"));
  mpfr_inp_str(params->ell_g, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"));
  mpfr_inp_str(params->xi, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"));
  
  gghlite_enc_init(params->pzt, params);
  gghlite_enc_init(params->ntt->w, params);
  gghlite_enc_init(params->ntt->w_inv, params);
  gghlite_enc_init(params->ntt->phi, params);
  gghlite_enc_init(params->ntt->phi_inv, params);

  gghlite_enc_fread(fp, params->pzt);
  CHECK(fscanf(fp, "\n"));
  CHECK(fscanf(fp, "%zd\n", &params->ntt->n));
  gghlite_enc_fread(fp, params->ntt->w);
  CHECK(fscanf(fp, "\n"));
  gghlite_enc_fread(fp, params->ntt->w_inv);
  CHECK(fscanf(fp, "\n"));
  gghlite_enc_fread(fp, params->ntt->phi);
  CHECK(fscanf(fp, "\n"));
  gghlite_enc_fread(fp, params->ntt->phi_inv);

  gghlite_params_set_D_sigmas(params);
}

void fread_mife_pp(mife_pp_t pp, char *filepath) {
  FILE *fp = fopen(filepath, "r");
  int flag_int;
  CHECK(fscanf(fp, "%d %d %d %d %d %d %d\n",
    &pp->num_inputs,
    &pp->bitstr_len,
    &pp->L,
    &pp->gamma,
    &pp->kappa,
    &pp->numR,
    &flag_int
  ));
  pp->n = malloc(pp->num_inputs * sizeof(int));
  for(int i = 0; i < pp->num_inputs; i++) {
    CHECK(fscanf(fp, "%d ", &pp->n[i]));
  }
  CHECK(fscanf(fp, "\n"));
  pp->gammas = malloc(pp->num_inputs * sizeof(int));
  for(int i = 0; i < pp->num_inputs; i++) {
    CHECK(fscanf(fp, "%d ", &pp->gammas[i]));
  }
  CHECK(fscanf(fp, "\n"));
  pp->flags = flag_int;
  fmpz_init(pp->p);
  fmpz_fread(fp, pp->p);
  CHECK(fscanf(fp, "\n"));

  pp->params_ref = malloc(sizeof(gghlite_params_t));
  fread_gghlite_params(fp, *pp->params_ref);
  fclose(fp);
}

void fwrite_mife_ciphertext(mife_pp_t pp, mife_ciphertext_t ct, char *filepath) {
  FILE *fp = fopen(filepath, "w");
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int j = 0; j < pp->n[i]; j++) {
      fwrite_gghlite_enc_mat(pp, ct->enc[i][j], fp);
    }
  }
  fclose(fp);
}

void fwrite_gghlite_enc_mat(mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp) {
  fprintf(fp, " %d ", m->nrows);
  fprintf(fp, " %d ", m->ncols);
  for(int i = 0; i < m->nrows; i++) {
    for(int j = 0; j < m->ncols; j++) {
      gghlite_enc_fprint(fp, m->m[i][j]);
      fprintf(fp, "\n");
    }
  }
}

void fread_mife_ciphertext(mife_pp_t pp, mife_ciphertext_t ct, char *filepath) {
  FILE *fp = fopen(filepath, "r");

  ct->enc = malloc(pp->num_inputs * sizeof(gghlite_enc_mat_t *));
  for(int i = 0; i < pp->num_inputs; i++) {
    ct->enc[i] = malloc(pp->n[i] * sizeof(gghlite_enc_mat_t));
    for(int j = 0; j < pp->n[i]; j++) {
      fread_gghlite_enc_mat(pp, ct->enc[i][j], fp);
    }
  }
  fclose(fp);
}

void fread_gghlite_enc_mat(mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp) {
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
      CHECK(fscanf(fp, "\n"));
    }
  }
}

void mife_ciphertext_clear(mife_pp_t pp, mife_ciphertext_t ct) {
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int j = 0; j < pp->n[i]; j++) {
      gghlite_enc_mat_clear(ct->enc[i][j]);
    }
    free(ct->enc[i]);
  }
  free(ct->enc);
}

void mife_clear_pp_read(mife_pp_t pp) {
  mife_clear_pp(pp);
  gghlite_params_clear_read(*pp->params_ref);
}

void mife_clear_pp(mife_pp_t pp) {
  fmpz_clear(pp->p);
  free(pp->n);
  free(pp->gammas);
}

void mife_clear_sk(mife_sk_t sk) {
  gghlite_sk_clear(sk->self, 1);
  for(int i = 0; i < sk->numR; i++) {
    fmpz_mat_clear(sk->R[i]);
    fmpz_mat_clear(sk->R_inv[i]);
  }
  free(sk->R);
  free(sk->R_inv);
  aes_randclear(sk->randstate);
}

void mife_mat_clr_clear(mife_pp_t pp, mife_mat_clr_t met) {
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int j = 0; j < pp->n[i]; j++) {
      fmpz_mat_clear(met->clr[i][j]);
    }
    free(met->clr[i]);
  }
  free(met->clr);
}

void print_mife_mat_clr(mife_pp_t pp, mife_mat_clr_t met) {
  for(int i = 0; i < pp->num_inputs; i++) {
    printf("clr[%d] matrices: \n", i);
    for(int j = 0; j < pp->n[i]; j++) {
      fmpz_mat_print_pretty(met->clr[i][j]);
      printf("\n\n");
    }
  }
  printf("\n");
}

void mife_encrypt(mife_ciphertext_t ct, fmpz_t message, mife_pp_t pp,
    mife_sk_t sk) {
  // compute a random index in the range [0,2^L]
  fmpz_t index, powL, two;
  fmpz_init(index);
  fmpz_init(powL);
  fmpz_init_set_ui(two, 2);
  fmpz_pow_ui(powL, two, pp->L); // computes powL = 2^L
  fmpz_set_ui(index, 0);
  fmpz_randm_aes(index, sk->randstate, powL);
  fmpz_clear(powL);
  fmpz_clear(two);
  
  mife_mat_clr_t met;
  pp->setfn(met, message, pp, sk);

  if(! (pp->flags & MIFE_NO_RANDOMIZERS)) {
    mife_apply_randomizers(met, pp, sk);     
  }
  mife_set_encodings(ct, met, index, pp, sk);
  fmpz_clear(index);
  mife_mat_clr_clear(pp, met);
}

void mife_mbp_set(
    void *mbp_params,
    mife_pp_t pp,
    int num_inputs,
    int (*paramfn)(int, int),
    void (*kilianfn)(struct _mife_pp_struct *, int *),
    void (*orderfn)(int, int *, int *),
    void (*setfn)(mife_mat_clr_t, fmpz_t, struct _mife_pp_struct *, mife_sk_t),
    int (*parsefn)(char **)
    ) {
  pp->mbp_params = mbp_params;
  pp->num_inputs = num_inputs;
  pp->paramfn = paramfn;
  pp->kilianfn = kilianfn;
  pp->orderfn = orderfn;
  pp->setfn = setfn;
  pp->parsefn = parsefn;
} 

void mife_setup(mife_pp_t pp, mife_sk_t sk, int L, int lambda,
    gghlite_flag_t ggh_flags, char *shaseed) {
  aes_randinit_seed(sk->randstate, shaseed);

  pp->n = malloc(pp->num_inputs * sizeof(int));
  for(int index = 0; index < pp->num_inputs; index++) { 
    pp->n[index] = pp->paramfn(pp->bitstr_len, index);
  }
  
  pp->kappa = pp->n[0] + pp->n[1];
  pp->L = L;

  pp->gamma = 0;
  pp->gammas = malloc(pp->num_inputs * sizeof(int));
  for(int i = 0 ; i < pp->num_inputs; i++) {
    pp->gammas[i] = 1 + (pp->n[i]-1) * (pp->L+1);
    pp->gamma += pp->gammas[i];
  }

  //printf("kappa: %d\n", pp->kappa);

  gghlite_jigsaw_init_gamma(sk->self,
                      lambda,
                      pp->kappa,
                      pp->gamma,
                      ggh_flags,
                      sk->randstate);

  pp->params_ref = &(sk->self->params);

  /*
  printf("Supporting at most 2^%d plaintexts \n", pp->L);
  printf("of length %d, with gamma = %d\n\n", pp->bitstr_len, pp->gamma);
  */
  
  fmpz_init(pp->p);
  fmpz_poly_oz_ideal_norm(pp->p, sk->self->g, sk->self->params->n, 0);

  // set the kilian randomizers in sk
  pp->numR = pp->kappa - 1;
  sk->numR = pp->numR;
  int *dims = malloc(pp->numR * sizeof(int));
  pp->kilianfn(pp, dims);

  sk->R = malloc(sk->numR * sizeof(fmpz_mat_t));
  sk->R_inv = malloc(sk->numR * sizeof(fmpz_mat_t));

  for (int k = 0; k < pp->numR; k++) {
    fmpz_mat_init(sk->R[k], dims[k], dims[k]);
    for (int i = 0; i < dims[k]; i++) {
      for(int j = 0; j < dims[k]; j++) {
        fmpz_randm_aes(fmpz_mat_entry(sk->R[k], i, j), sk->randstate, pp->p);
      }
    }
	
    fmpz_mat_init(sk->R_inv[k], dims[k], dims[k]);
    fmpz_modp_matrix_inverse(sk->R_inv[k], sk->R[k], dims[k], pp->p);
  }

  free(dims);
}

void fmpz_mat_scalar_mul_modp(fmpz_mat_t m, fmpz_t scalar, fmpz_t modp) {
  for (int i = 0; i < m->r; i++) {
    for(int j = 0; j < m->c; j++) {
      fmpz_mul(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), scalar);
      fmpz_mod(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), modp);
    }
  }
}

void mife_apply_randomizers(mife_mat_clr_t met, mife_pp_t pp, mife_sk_t sk) {
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int k = 0; k < pp->n[i]; k++) {
      fmpz_t rand;
      fmpz_init(rand);
      fmpz_randm_aes(rand, sk->randstate, pp->p);
      fmpz_mat_scalar_mul_modp(met->clr[i][k], rand, pp->p);
      fmpz_clear(rand);
    }
  }
}

// message >= 0, d >= 2
void message_to_dary(ulong *dary, int bitstring_len, fmpz_t message, int d) {
  assert(d >= 2);
  fmpz_t message2;
  fmpz_init_set(message2, message);
  fmpz_t modresult;
  fmpz_init(modresult);
  fmpz_t modd;
  fmpz_init_set_ui(modd, d);

  int i;
  for (i = bitstring_len - 1; i >= 0; i--) {
    fmpz_tdiv_qr(message2, modresult, message2, modd);
    dary[i] = fmpz_get_ui(modresult);
  }
  
  fmpz_clear(message2);
  fmpz_clear(modd);
  fmpz_clear(modresult);
}



/**
 * Generates the i^th member of the exclusive partition family for the index 
 * sets.
 *
 * @param partitioning The description of the partitioning, each entry is in 
 * [0,d-1] and it is of length (1 + (d-1)(L+1)).
 * @param index The index being the i^th member of the partition family
 * @param L the log of the size of the partition family. So, i must be in the 
 * range [0,2^L-1]
 * @param nu The number of total elements to be multiplied. partitioning[] 
 * will describe a nu-partition of the universe set.
 */ 
void mife_gen_partitioning(int *partitioning, fmpz_t index, int L, int nu) {
  int j = 0;

  ulong *bitstring = malloc(L * sizeof(ulong));
  memset(bitstring, 0, L * sizeof(ulong));
  message_to_dary(bitstring, L, index, 2);

  for(; j < nu; j++) {
    partitioning[j] = j;
  }

  for(int k = 0; k < L; k++) {
    for(int j1 = 1; j1 < nu; j1++) {
      partitioning[j] = (bitstring[k] == 1) ? j1 : 0;
      j++;
    }
  }
  free(bitstring);
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


void mife_mat_encode(mife_pp_t pp, mife_sk_t sk, gghlite_enc_mat_t enc,
    fmpz_mat_t m, int *group) {
  gghlite_clr_t e;
  gghlite_clr_init(e);
  for(int i = 0; i < enc->nrows; i++) {
    for(int j = 0; j < enc->ncols; j++) {
      fmpz_poly_set_coeff_fmpz(e, 0, fmpz_mat_entry(m, i, j));
      gghlite_enc_set_gghlite_clr(enc->m[i][j], sk->self, e, 1, group, 1,
          sk->randstate);
      NUM_ENCODINGS_GENERATED++;
      if(PRINT_ENCODING_PROGRESS) {
        printf("Generated encoding %d / %d (Time elapsed: %8.2f s)\n",
            NUM_ENCODINGS_GENERATED,
            get_NUM_ENC(),
            get_T());
      }
    }
  }
  gghlite_clr_clear(e);
}

void gghlite_enc_mat_zeros_print(mife_pp_t pp, gghlite_enc_mat_t m) {
  for(int i = 0; i < m->nrows; i++) {
    printf("[");
    for(int j = 0; j < m->ncols; j++) {
      printf(gghlite_enc_is_zero(*pp->params_ref, m->m[i][j]) ? "0 " : "x " );
    }
    printf("]\n");
  }
}

void mife_set_encodings(mife_ciphertext_t ct, mife_mat_clr_t met, fmpz_t index,
    mife_pp_t pp, mife_sk_t sk) {

  int **ptns = malloc(pp->num_inputs * sizeof(int *));
  for(int i = 0; i < pp->num_inputs; i++) {
    ptns[i] = malloc(pp->gammas[i] * sizeof(int));
    mife_gen_partitioning(ptns[i], index, pp->L, pp->n[i]);
  }

  /* construct the partitions in the group array form */
  int ***groups = malloc(pp->num_inputs * sizeof(int **));
  for(int i = 0; i < pp->num_inputs; i++) {
    int gamma_offset = 0;
    for(int h = 0; h < i; h++) {
      gamma_offset += pp->gammas[h];
    }

    groups[i] = malloc(pp->n[i] * sizeof(int *));
    for(int j = 0; j < pp->n[i]; j++) {
      groups[i][j] = malloc(pp->gamma * sizeof(int));
      memset(groups[i][j], 0, pp->gamma * sizeof(int));
      for(int k = 0; k < pp->gammas[i]; k++) {
        if(ptns[i][k] == j) {
          groups[i][j][k + i * pp->gammas[0]] = 1;
        }
      }
    }
  }

  for(int i = 0; i < pp->num_inputs; i++) {
    free(ptns[i]);
  }
  free(ptns);

  if(pp->flags & MIFE_SIMPLE_PARTITIONS) {
    // override group arrays with trivial partitioning
    for(int i = 0; i < pp->num_inputs; i++) {
      for(int j = 0; j < pp->n[i]; j++) {
        memset(groups[i][j], 0, pp->gamma * sizeof(int));
      }
    }
    
    for(int k = 0; k < pp->gamma; k++) {
      groups[0][0][k] = 1;
    }
  }

  if(! (pp->flags & MIFE_NO_KILIAN)) {
    // apply kilian to the cleartext matrices (overwriting them in the process)

    fmpz_mat_t tmp;

    for(int index = 0; index < pp->kappa; index++) {
      int i, j;
      pp->orderfn(index, &i, &j);

      // first one
      if(index == 0) {
        fmpz_mat_init(tmp, met->clr[i][j]->r, sk->R[0]->c);
        fmpz_mat_mul(tmp, met->clr[i][j], sk->R[0]);
        fmpz_mat_set(met->clr[i][j], tmp);
        fmpz_mat_clear(tmp);
        continue;
      }

      // last one
      if(index == pp->kappa - 1) {
        fmpz_mat_init(tmp, sk->R_inv[pp->numR-1]->r, met->clr[i][j]->c);
        fmpz_mat_mul(tmp, sk->R_inv[pp->numR-1], met->clr[i][j]);
        fmpz_mat_set(met->clr[i][j], tmp);
        fmpz_mat_clear(tmp);
        continue;
      }
      
      // all others
      fmpz_mat_init(tmp, sk->R_inv[index-1]->r, met->clr[i][j]->c);
      fmpz_mat_mul(tmp, sk->R_inv[index-1], met->clr[i][j]);
      fmpz_mat_mul(tmp, tmp, sk->R[index]);
      fmpz_mat_set(met->clr[i][j], tmp);
      fmpz_mat_clear(tmp);
    }
  }

  // encode
  ct->enc = malloc(pp->num_inputs * sizeof(gghlite_enc_mat_t *));
  for(int i = 0; i < pp->num_inputs; i++) {
    ct->enc[i] = malloc(pp->n[i] * sizeof(gghlite_enc_mat_t));
    for(int j = 0; j < pp->n[i]; j++) {
      gghlite_enc_mat_init(sk->self->params, ct->enc[i][j],
          met->clr[i][j]->r, met->clr[i][j]->c);
      mife_mat_encode(pp, sk, ct->enc[i][j], met->clr[i][j], groups[i][j]);
    }
  }
  
  // free group arrays
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int j = 0; j < pp->n[i]; j++) {
      free(groups[i][j]);
    }
    free(groups[i]);
  }
  free(groups);
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

int mife_evaluate(mife_pp_t pp, mife_ciphertext_t *cts) {
  gghlite_enc_mat_t tmp;

  for(int index = 1; index < pp->kappa; index++) {
    int i, j;
    pp->orderfn(index, &i, &j);
    
    if(index == 1) {
      // multiply the 0th index with the 1st index
      int i0, j0;
      pp->orderfn(0, &i0, &j0);
      gghlite_enc_mat_init(*pp->params_ref, tmp,
        cts[i0]->enc[i0][j0]->nrows, cts[i]->enc[i][j]->ncols);
      gghlite_enc_mat_mul(*pp->params_ref, tmp,
        cts[i0]->enc[i0][j0], cts[i]->enc[i][j]);
      continue;
    }

    gghlite_enc_mat_mul(*pp->params_ref, tmp, tmp, cts[i]->enc[i][j]);
  }

  char **result = malloc(tmp->nrows * sizeof(char *));
  for(int i = 0; i < tmp->nrows; i++) {
    result[i] = malloc(tmp->ncols * sizeof(char));
    for(int j = 0; j < tmp->ncols; j++) {
      result[i][j] = 1 - gghlite_enc_is_zero(*pp->params_ref, tmp->m[i][j]);
    }
  }
  gghlite_enc_mat_clear(tmp);

  int ret = pp->parsefn(result);

  for(int i = 0; i < tmp->nrows; i++) {
    free(result[i]);
  }
  free(result);

  return ret;
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

int int_arrays_equal(ulong *arr1, ulong *arr2, int length) {
  for (int i = 0; i < length; i++) {
    if (arr1[i] != arr2[i])
      return 1;
  }
  return 0;
}

void test_dary_conversion() {
  printf("Testing d-ary conversion function...                          ");
  ulong dary1[4];
  ulong dary2[8];
  ulong dary3[4];
  ulong correct1[] = {1,0,1,0};
  ulong correct2[] = {0,0,0,0,5,4,1,4};
  ulong correct3[] = {0,0,0,2};
  fmpz_t num1, num2, num3;
  fmpz_init_set_ui(num1, 10);
  fmpz_init_set_ui(num2, 1234);
  fmpz_init_set_ui(num3, 2);

  message_to_dary(dary1, 4, num1, 2);
  message_to_dary(dary2, 8, num2, 6);
  message_to_dary(dary3, 4, num3, 11);


  int status = 0;
  status += int_arrays_equal(dary1, correct1, 4);
  status += int_arrays_equal(dary2, correct2, 8);
  status += int_arrays_equal(dary3, correct3, 4);



  if (status == 0)
    printf("SUCCESS\n");
  else
    printf("FAIL\n");	

}

void print_matrix_sage(fmpz_mat_t a) {
  printf("\n[");
  for(int i = 0; i < a->r; i++) {
    printf("[");
    for(int j = 0; j < a->c; j++) {
      fmpz_print(fmpz_mat_entry(a, i, j));
      if(j != a->c-1) {
        printf(", ");
      }
    }
    printf("]");
    if(i != a->r-1) {
      printf(",");
    }
    printf("\n");
  }
  printf("]\n");
}


int test_matrix_inv(int n, aes_randstate_t randstate, fmpz_t modp) {
  printf("\nTesting matrix_inv function...                                ");
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
  if (status != 0)
    printf("SUCCESS\n");
  else
    printf("FAIL\n");

  fmpz_mat_clear(a);
  fmpz_mat_clear(inv);
  fmpz_mat_clear(prod);
  fmpz_mat_clear(identity);
  return status;
}

/**
 * Code to find the inverse of a matrix, adapted from:
 * https://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html
 */

// uses gaussian elimination to obtain the determinant of a matrix
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

void set_NUM_ENC(int val) {
  NUM_ENCODINGS_TOTAL = val;
}

int get_NUM_ENC() {
  return NUM_ENCODINGS_TOTAL;
}

void reset_T() {
  T = ggh_walltime(0);
}

float get_T() {
  return ggh_seconds(ggh_walltime(T));
}

