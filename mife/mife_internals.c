#include "mife_internals.h"


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
  gghlite_params_clear_read(*pp->params_ref);
  free(pp->params_ref);
  mife_clear_pp(pp);
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
  fmpz_randm_aes(rand, randstate, pp->p);
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
    fmpz_mat_t m, int *group, aes_randstate_t randstate) {
  gghlite_clr_t e;
  gghlite_clr_init(e);
  for(int i = 0; i < enc->nrows; i++) {
    for(int j = 0; j < enc->ncols; j++) {
      fmpz_poly_set_coeff_fmpz(e, 0, fmpz_mat_entry(m, i, j));
      gghlite_enc_set_gghlite_clr(enc->m[i][j], sk->self, e, 1, group, 1,
          randstate);
      NUM_ENCODINGS_GENERATED++;
        timer_printf("\r    Generated encoding [%d / %d] (Time elapsed: %8.2f s)",
            NUM_ENCODINGS_GENERATED,
            get_NUM_ENC(),
            get_T());
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

void mife_set_encodings(mife_ciphertext_t ct, mife_mat_clr_t met, fmpz_t index,
    mife_pp_t pp, mife_sk_t sk, aes_randstate_t randstate) {

  int ***groups = mife_partitions(pp, index);

  if(! (pp->flags & MIFE_NO_KILIAN)) {
    // apply kilian to the cleartext matrices (overwriting them in the process)
    for(int index = 0; index < pp->kappa; index++) {
      int i, j;
      pp->orderfn(pp, index, &i, &j);
      mife_apply_kilian(pp, sk, met->clr[i][j], index);
    }
  }

  // encode
  ct->enc = malloc(pp->num_inputs * sizeof(gghlite_enc_mat_t *));
  for(int i = 0; i < pp->num_inputs; i++) {
    ct->enc[i] = malloc(pp->n[i] * sizeof(gghlite_enc_mat_t));
    for(int j = 0; j < pp->n[i]; j++) {
      gghlite_enc_mat_init(sk->self->params, ct->enc[i][j],
          met->clr[i][j]->r, met->clr[i][j]->c);
      mife_mat_encode(pp, sk, ct->enc[i][j], met->clr[i][j], groups[i][j],
          randstate);
    }
  }
  
  // free group arrays
  mife_partitions_clear(pp, groups);
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

void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p) {
  for(int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      fmpz_mod(fmpz_mat_entry(m, i, j), fmpz_mat_entry(m, i, j), p);
    }
  }
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

