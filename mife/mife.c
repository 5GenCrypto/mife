#include "mife.h"
#include "util.h"

int NUM_ENCODINGS_GENERATED;
int PRINT_ENCODING_PROGRESS;
int NUM_ENCODINGS_TOTAL;

f2_matrix mife_zt_all(const_mmap_vtable mmap, const mife_pp_t pp, mmap_enc_mat_t ct) {
  f2_matrix pt;
  if(!f2_matrix_zero(&pt, ct->nrows, ct->ncols))
    return pt;

  for(int i = 0; i < ct->nrows; i++) {
    for(int j = 0; j < ct->ncols; j++) {
      pt.elems[i][j] = !mmap->enc->is_zero(ct->m[i][j], pp->params_ref);
    }
  }

  return pt;
}

void mife_init_params(mife_pp_t pp, mife_flag_t flags) {
  pp->flags = flags;
}

void mife_mbp_set(
    void *mbp_params,
    mife_pp_t pp,
    int num_inputs,
    int (*paramfn)  (mife_pp_t, int),
    void (*kilianfn)(mife_pp_t, int *),
    void (*orderfn) (mife_pp_t, int, int *, int *),
    void (*setfn)   (mife_pp_t, mife_mat_clr_t, void *),
    int (*parsefn)  (mife_pp_t, f2_matrix)
    ) {
  pp->mbp_params = mbp_params;
  pp->num_inputs = num_inputs;
  pp->paramfn = paramfn;
  pp->kilianfn = kilianfn;
  pp->orderfn = orderfn;
  pp->setfn = setfn;
  pp->parsefn = parsefn;
} 

void fmpz_rand_mat_square_aes(fmpz_mat_t m, int dim, aes_randstate_t randstate, fmpz_t p) {
  for (int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      fmpz_randm_aes(fmpz_mat_entry(m, i, j), randstate, p);
    }
  }
}

void mife_setup(const_mmap_vtable mmap, mife_pp_t pp, mife_sk_t sk, int L, int lambda,
    aes_randstate_t randstate) {

  pp->n = malloc(pp->num_inputs * sizeof(int));
  pp->kappa = 0;
  for(int index = 0; index < pp->num_inputs; index++) {
    pp->n[index] = pp->paramfn(pp, index);
    pp->kappa += pp->n[index];
  }
  pp->L = L;

  pp->gamma = 0;
  pp->gammas = malloc(pp->num_inputs * sizeof(int));
  for(int i = 0 ; i < pp->num_inputs; i++) {
    pp->gammas[i] = 1 + (pp->n[i]-1) * (pp->L+1);
    pp->gamma += pp->gammas[i];
  }

  timer_printf("Starting MMAP secret key initialization: %d %d %d...\n",
      lambda, pp->kappa, pp->gamma);
  sk->self = malloc(mmap->sk->size);
  mmap->sk->init(sk->self, lambda, pp->kappa, pp->gamma, randstate, false);
  timer_printf("Finished MMAP secret key initialization\n");

  /* For const correctness, we should probably have two separate
   * _mife_pp_struct types, one for pp's read from disk and one for pp's
   * produced by picking a new keypair. Instead we keep only an informal
   * (non-compiler-checked) invariant that we do not call read or clear on
   * params_ref's produced from a new keypair generation.
   */
  pp->params_ref = (mmap_pp *)mmap->sk->pp(sk->self);

  timer_printf("Starting setting p...\n");
  start_timer();
  fmpz_init(pp->p);
  mmap->sk->plaintext_field(sk->self, pp->p);
  timer_printf("Finished setting p");
  print_timer();
  timer_printf("\n");

  // set the kilian randomizers in sk
  pp->numR = pp->kappa - 1;
  sk->numR = pp->numR;
  int *dims = malloc(pp->numR * sizeof(int));
  pp->kilianfn(pp, dims);

  sk->R = malloc(sk->numR * sizeof(fmpz_mat_t));
  sk->R_inv = malloc(sk->numR * sizeof(fmpz_mat_t));

  timer_printf("Starting setting Kilian matrices...\n");
  start_timer();

  /* do not parallelize calls to randomness generation! */  
  uint64_t t_init = ggh_walltime(0);
  for (int k = 0; k < pp->numR; k++) {
    fmpz_mat_init(sk->R[k], dims[k], dims[k]);
    fmpz_rand_mat_square_aes(sk->R[k], dims[k], randstate, pp->p);
    timer_printf("\r    Init Progress: [%d / %d] %8.2fs", pp->numR,
        k, ggh_seconds(ggh_walltime(t_init)));
  }
  timer_printf("\n");
  
  int progress_count_approx = 0;
  uint64_t t = ggh_walltime(0);
  int *non_invertible;
  if(ALLOC_FAILS(non_invertible, pp->numR)) assert(false);
#pragma omp parallel for
  for (int k = 0; k < pp->numR; k++) {
    fmpz_mat_init(sk->R_inv[k], dims[k], dims[k]);
    non_invertible[k] = fmpz_modp_matrix_inverse(sk->R_inv[k], sk->R[k], dims[k], pp->p);
    progress_count_approx++;
    timer_printf("\r    Inverse Computation Progress (Parallel): \
        [%lu / %lu] %8.2fs",
        progress_count_approx, pp->numR, ggh_seconds(ggh_walltime(t)));
  }
  for (int k = 0; k < pp->numR; k++) {
    while(non_invertible[k]) {
      timer_printf("Retrying matrix %d\n", k);
      fmpz_rand_mat_square_aes(sk->R[k], dims[k], randstate, pp->p);
      non_invertible[k] = fmpz_modp_matrix_inverse(sk->R_inv[k], sk->R[k], dims[k], pp->p);
    }
  }
  timer_printf("\n");
  timer_printf("Finished setting Kilian matrices");
  print_timer();
  timer_printf("\n");

  free(dims);
}

void mife_encrypt(const_mmap_vtable mmap, mife_ciphertext_t ct, void *message, mife_pp_t pp,
    mife_sk_t sk, aes_randstate_t randstate) {
  // compute a random index in the range [0,2^L]
  fmpz_t index, powL, two;
  fmpz_init(index);
  fmpz_init(powL);
  fmpz_init_set_ui(two, 2);
  fmpz_pow_ui(powL, two, pp->L); // computes powL = 2^L
  fmpz_set_ui(index, 0);
  fmpz_randm_aes(index, randstate, powL);
  fmpz_clear(powL);
  fmpz_clear(two);
  
  mife_mat_clr_t met;
  pp->setfn(pp, met, message);

  if(! (pp->flags & MIFE_NO_RANDOMIZERS)) {
    mife_apply_randomizers(met, pp, sk, randstate);
  }
  mife_set_encodings(mmap, ct, met, index, pp, sk, randstate);
  fmpz_clear(index);
  mife_mat_clr_clear(pp, met);
}

int mife_evaluate(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t *cts) {
  mmap_enc_mat_t tmp;

  for(int index = 1; index < pp->kappa; index++) {
    int i, j;
    pp->orderfn(pp, index, &i, &j);
    
    if(index == 1) {
      // multiply the 0th index with the 1st index
      int i0, j0;
      pp->orderfn(pp, 0, &i0, &j0);
      mmap_enc_mat_init(mmap, pp->params_ref, tmp,
        cts[i0]->enc[i0][j0]->nrows, cts[i]->enc[i][j]->ncols);
      mmap_enc_mat_mul(mmap, pp->params_ref, tmp,
        cts[i0]->enc[i0][j0], cts[i]->enc[i][j]);
      continue;
    }

    mmap_enc_mat_mul(mmap, pp->params_ref, tmp, tmp, cts[i]->enc[i][j]);
  }

  f2_matrix result = mife_zt_all(mmap, pp, tmp);
  mmap_enc_mat_clear(mmap, tmp);
  int ret = pp->parsefn(pp, result);
  f2_matrix_free(result);

  return ret;
}

void mife_encrypt_setup(mife_pp_t pp, fmpz_t uid, void *message,
    mife_mat_clr_t out_clr, int ****out_partitions) {
  pp->setfn(pp, out_clr, message);
  *out_partitions = mife_partitions(pp, uid);
}

void mife_encrypt_single(const_mmap_vtable mmap, mife_pp_t pp, mife_sk_t sk, aes_randstate_t randstate,
    int global_index, mife_mat_clr_t clr, int ***partitions,
    mmap_enc_mat_t dest) {
  int position_index, local_index;
  fmpz_mat_t src;

  pp->orderfn(pp, global_index, &position_index, &local_index);
  fmpz_mat_init_set(src, clr->clr[position_index][local_index]);

  if(!(pp->flags & MIFE_NO_RANDOMIZERS))
    mife_apply_randomizer(pp, randstate, src);

  if(!(pp->flags & MIFE_NO_KILIAN))
    mife_apply_kilian(pp, sk, src, global_index);

  mmap_enc_mat_init(mmap, pp->params_ref, dest, src->r, src->c);
  mife_mat_encode(mmap, pp, sk, dest, src, partitions[position_index][local_index], randstate);
  fmpz_mat_clear(src);
}

void mife_encrypt_clear(mife_pp_t pp, mife_mat_clr_t clr, int ***partitions) {
  mife_mat_clr_clear(pp, clr);
  mife_partitions_clear(pp, partitions);
}

