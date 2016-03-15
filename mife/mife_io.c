#include "mife_io.h"

/**
 * 
 * Members of pp that are not transferred:
 * params->x
 * params->y
 * params->D_sigma_p
 * params->D_sigma_s
 */
void fwrite_mife_pp(const_mmap_vtable mmap, mife_pp_t pp, char *filepath) {
  FILE *fp = fopen(filepath, "wb");
  fprintf(fp, "%d %d %d %d %d %d\n",
    pp->num_inputs,
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
  fmpz_out_raw(fp, pp->p);
  fprintf(fp, "\n");

  mmap->pp->fwrite(pp->params_ref, fp);
  fclose(fp);
}

void fread_mife_pp(const_mmap_vtable mmap, mife_pp_t pp, char *filepath) {
  FILE *fp = fopen(filepath, "rb");
  int flag_int;
  CHECK(fscanf(fp, "%d %d %d %d %d %d\n",
    &pp->num_inputs,
    &pp->L,
    &pp->gamma,
    &pp->kappa,
    &pp->numR,
    &flag_int
  ), 6);
  pp->n = malloc(pp->num_inputs * sizeof(int));
  for(int i = 0; i < pp->num_inputs; i++) {
    CHECK(fscanf(fp, "%d ", &pp->n[i]), 1);
  }
  CHECK(fscanf(fp, "\n"), 0);
  pp->gammas = malloc(pp->num_inputs * sizeof(int));
  for(int i = 0; i < pp->num_inputs; i++) {
    CHECK(fscanf(fp, "%d ", &pp->gammas[i]), 1);
  }
  CHECK(fscanf(fp, "\n"), 0);
  pp->flags = flag_int;
  fmpz_init(pp->p);
  fmpz_inp_raw(pp->p, fp);
  CHECK(fscanf(fp, "\n"), 0);

  pp->params_ref = malloc(mmap->pp->size);
  mmap->pp->fread(pp->params_ref, fp);
  fclose(fp);
}

void fwrite_mife_sk(const_mmap_vtable mmap, mife_sk_t sk, char *filepath) {
	uint64_t t = ggh_walltime(0);
  FILE *fp = fopen(filepath, "wb");
	timer_printf("Starting writing Kilian matrices...\n");
  fprintf(fp, "%d\n", sk->numR);
  for(int i = 0; i < sk->numR; i++) {
    fprintf(fp, "%ld %ld\n", sk->R[i]->r, sk->R[i]->c);
    fmpz_mat_fprint_raw(fp, sk->R[i]);
    fprintf(fp, "\n");
    fprintf(fp, "%ld %ld\n", sk->R_inv[i]->r, sk->R_inv[i]->c);
    fmpz_mat_fprint_raw(fp, sk->R_inv[i]);
    fprintf(fp, "\n");
	timer_printf("\r    Progress: [%lu / %lu] %8.2fs",
		i, sk->numR, ggh_seconds(ggh_walltime(t)));

  }
	timer_printf("\n");
	timer_printf("Finished writing Kilian matrices %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	mmap->sk->fwrite(sk->self, fp);

  fclose(fp);

}

void fread_mife_sk(const_mmap_vtable mmap, mife_sk_t sk, char *filepath) {
	uint64_t t = ggh_walltime(0);
  FILE *fp = fopen(filepath, "rb");
	timer_printf("Starting reading Kilian matrices...\n");
  CHECK(fscanf(fp, "%d\n", &sk->numR), 1);
  sk->R = malloc(sk->numR * sizeof(fmpz_mat_t));
  sk->R_inv = malloc(sk->numR * sizeof(fmpz_mat_t));

  for(int i = 0; i < sk->numR; i++) {
    unsigned long r1, c1, r2, c2;
    CHECK(fscanf(fp, "%ld %ld\n", &r1, &c1), 2);
    fmpz_mat_init(sk->R[i], r1, c1);
    fmpz_mat_fread_raw(fp, sk->R[i]);
    CHECK(fscanf(fp, "\n"), 0);
    CHECK(fscanf(fp, "%ld %ld\n", &r2, &c2), 2);
    fmpz_mat_init(sk->R_inv[i], r2, c2);
    fmpz_mat_fread_raw(fp, sk->R_inv[i]);
    CHECK(fscanf(fp, "\n"), 0);
	timer_printf("\r    Progress: [%lu / %lu] %8.2fs",
		i, sk->numR, ggh_seconds(ggh_walltime(t)));
  }
	timer_printf("\n");
	timer_printf("Finished reading Kilian matrices %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	sk->self = malloc(mmap->sk->size);
	mmap->sk->fread(sk->self, fp);

  fclose(fp);
}

void fwrite_mife_ciphertext(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t ct, char *filepath) {
  FILE *fp = fopen(filepath, "wb");
  for(int i = 0; i < pp->num_inputs; i++) {
    for(int j = 0; j < pp->n[i]; j++) {
      fwrite_gghlite_enc_mat(mmap, pp, ct->enc[i][j], fp);
    }
  }
  fclose(fp);
}

void fwrite_gghlite_enc_mat(const_mmap_vtable mmap, mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp) {
  fprintf(fp, " %d ", m->nrows);
  fprintf(fp, " %d ", m->ncols);
  for(int i = 0; i < m->nrows; i++) {
    for(int j = 0; j < m->ncols; j++) {
      mmap->enc->fwrite(m->m[i][j], fp);
      fprintf(fp, "\n");
    }
  }
}

void fread_mife_ciphertext(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t ct, char *filepath) {
  FILE *fp = fopen(filepath, "rb");

  ct->enc = malloc(pp->num_inputs * sizeof(gghlite_enc_mat_t *));
  for(int i = 0; i < pp->num_inputs; i++) {
    ct->enc[i] = malloc(pp->n[i] * sizeof(gghlite_enc_mat_t));
    for(int j = 0; j < pp->n[i]; j++) {
      fread_gghlite_enc_mat(mmap, pp, ct->enc[i][j], fp);
    }
  }
  fclose(fp);
}

void fread_gghlite_enc_mat(const_mmap_vtable mmap, const mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp) {
  int check1 = fscanf(fp, " %d ", &m->nrows);
  int check2 = fscanf(fp, " %d ", &m->ncols);
  assert(check1 == 1 && check2 == 1);
  m->m = malloc(m->nrows * sizeof(mmap_enc **));
  assert(m->m);
  for(int i = 0; i < m->nrows; i++) {
    m->m[i] = malloc(m->ncols * sizeof(mmap_enc *));
    assert(m->m[i]);
    for(int j = 0; j < m->ncols; j++) {
      m->m[i][j] = malloc(mmap->enc->size);
      assert(m->m[i][j]);
      mmap->enc->fread(m->m[i][j], fp);
      CHECK(fscanf(fp, "\n"), 0);
    }
  }
}

