#include "mife_internals.h"

int g_parallel;

void mife_ciphertext_clear(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t ct) {
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int j = 0; j < pp->n[i]; j++) {
      mmap_enc_mat_clear(mmap, ct->enc[i][j]);
    }
    free(ct->enc[i]);
  }
  free(ct->enc);
}

void mife_clear_pp_read(const_mmap_vtable mmap, mife_pp_t pp) {
  mmap->pp->clear(pp->params_ref);
  free(pp->params_ref);
  mife_clear_pp(pp);
}

void mife_clear_pp(mife_pp_t pp) {
  fmpz_clear(pp->p);
  free(pp->n);
  free(pp->gammas);
}

void mife_clear_sk(const_mmap_vtable mmap, mife_sk_t sk) {
  mmap->sk->clear(sk->self);
  free(sk->self);
  for(int i = 0; i < sk->numR; i++) {
    fmpz_mat_clear(sk->R[i]);
    fmpz_mat_clear(sk->R_inv[i]);
  }
  free(sk->R);
  free(sk->R_inv);
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

/* static void print_mife_mat_clr(mife_pp_t pp, mife_mat_clr_t met) { */
/*   for(int i = 0; i < pp->num_inputs; i++) { */
/*     printf("clr[%d] matrices: \n", i); */
/*     for(int j = 0; j < pp->n[i]; j++) { */
/*       fmpz_mat_print_pretty(met->clr[i][j]); */
/*       printf("\n\n"); */
/*     } */
/*   } */
/*   printf("\n"); */
/* } */


void fmpz_mat_scalar_mul_modp(fmpz_mat_t m, fmpz_t scalar, fmpz_t modp) {
  for (int i = 0; i < m->r; i++) {
    for(int j = 0; j < m->c; j++) {
      fmpz_mul(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), scalar);
      fmpz_mod(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), modp);
    }
  }
}

void mife_apply_randomizer(mife_pp_t pp, aes_randstate_t randstate, fmpz_mat_t m) {
  fmpz_t rand;
  fmpz_init(rand);
  do {
    fmpz_randm_aes(rand, randstate, pp->p);
  } while (fmpz_is_zero(rand) || fmpz_is_one(rand));
  fmpz_mat_scalar_mul_modp(m, rand, pp->p);
  fmpz_clear(rand);
}

void mife_apply_randomizers(mife_mat_clr_t met, mife_pp_t pp, mife_sk_t sk,
    aes_randstate_t randstate) {
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int k = 0; k < pp->n[i]; k++) {
      mife_apply_randomizer(pp, randstate, met->clr[i][k]);
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

void mife_mat_encode(const_mmap_vtable mmap, mife_pp_t pp, mife_sk_t sk, mmap_enc_mat_t enc,
    fmpz_mat_t m, int *group, aes_randstate_t randstate) {
  if (g_parallel) {
#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for(int i = 0; i < enc->nrows; i++) {
      for(int j = 0; j < enc->ncols; j++) {
          fmpz_t *plaintext = malloc(sizeof(fmpz_t));
          memcpy(plaintext, fmpz_mat_entry(m, i, j), sizeof(fmpz_t));
          mmap->enc->encode(enc->m[i][j], sk->self, 1, plaintext, group, randstate);
          NUM_ENCODINGS_GENERATED++;
          timer_printf("\r    Generated encoding [%d / %d] (Time elapsed: %8.2f s)",
              NUM_ENCODINGS_GENERATED,
              get_NUM_ENC(),
              get_T());
          free(plaintext);
      }
    }
  } else {
    for(int i = 0; i < enc->nrows; i++) {
      for(int j = 0; j < enc->ncols; j++) {
          fmpz_t *plaintext = malloc(sizeof(fmpz_t));
          memcpy(plaintext, fmpz_mat_entry(m, i, j), sizeof(fmpz_t));
          mmap->enc->encode(enc->m[i][j], sk->self, 1, plaintext, group, randstate);
          NUM_ENCODINGS_GENERATED++;
          timer_printf("\r    Generated encoding [%d / %d] (Time elapsed: %8.2f s)",
              NUM_ENCODINGS_GENERATED,
              get_NUM_ENC(),
              get_T());
          free(plaintext);
      }
    }
  }
}

int ***mife_partitions(mife_pp_t pp, fmpz_t index) {
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
          groups[i][j][k + gamma_offset] = 1;
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

  return groups;
}

void mife_partitions_clear(mife_pp_t pp, int ***groups) {
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int j = 0; j < pp->n[i]; j++) {
      free(groups[i][j]);
    }
    free(groups[i]);
  }
  free(groups);
}

void mife_apply_kilian(mife_pp_t pp, mife_sk_t sk, fmpz_mat_t m, int global_index) {
  fmpz_mat_t tmp;

  // first one
  if(global_index == 0) {
    fmpz_mat_init(tmp, m->r, sk->R[0]->c);
    fmpz_mat_mul(tmp, m, sk->R[0]);
  }

  // last one
  else if(global_index == pp->kappa - 1) {
    fmpz_mat_init(tmp, sk->R_inv[pp->numR-1]->r, m->c);
    fmpz_mat_mul(tmp, sk->R_inv[pp->numR-1], m);
  }

  // all others
  else {
    fmpz_mat_init(tmp, sk->R_inv[global_index-1]->r, m->c);
    fmpz_mat_mul(tmp, sk->R_inv[global_index-1], m);
    fmpz_mat_mul(tmp, tmp, sk->R[global_index]);
  }

  fmpz_mat_set(m, tmp);
  fmpz_mat_clear(tmp);
}

void mife_set_encodings(const_mmap_vtable mmap, mife_ciphertext_t ct, mife_mat_clr_t met, fmpz_t index,
    mife_pp_t pp, mife_sk_t sk, aes_randstate_t randstate) {

  int ***groups = mife_partitions(pp, index);

  if(! (pp->flags & MIFE_NO_KILIAN)) {
    // apply kilian to the cleartext matrices (overwriting them in the process)
    for(int idx = 0; idx < pp->kappa; idx++) {
      int i, j;
      pp->orderfn(pp, idx, &i, &j);
      mife_apply_kilian(pp, sk, met->clr[i][j], idx);
    }
  }

  // encode
  ct->enc = malloc(pp->num_inputs * sizeof(mmap_enc_mat_t *));
  for(int i = 0; i < pp->num_inputs; i++) {
    ct->enc[i] = malloc(pp->n[i] * mmap->enc->size);
    for(int j = 0; j < pp->n[i]; j++) {
      mmap_enc_mat_init(mmap, pp->params_ref, ct->enc[i][j],
          met->clr[i][j]->r, met->clr[i][j]->c);
      mife_mat_encode(mmap, pp, sk, ct->enc[i][j], met->clr[i][j], groups[i][j],
          randstate);
    }
  }

  // free group arrays
  mife_partitions_clear(pp, groups);
}

void mmap_enc_mat_zeros_print(const_mmap_vtable mmap, mife_pp_t pp, mmap_enc_mat_t m) {
  for(int i = 0; i < m->nrows; i++) {
    printf("[");
    for(int j = 0; j < m->ncols; j++) {
      printf(mmap->enc->is_zero(m->m[i][j], pp->params_ref) ? "0 " : "x " );
    }
    printf("]\n");
  }
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

