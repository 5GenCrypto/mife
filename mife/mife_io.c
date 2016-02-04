#include "mife_io.h"

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
  ), 7);
  
  gghlite_params_initzero(params, lambda, kappa, gamma);
  params->n = n;
  params->ell = ell;
  params->rerand_mask = rerand_mask;
  params->flags = gghlite_flag_int;

  fmpz_fread(fp, params->q);
  CHECK(fscanf(fp, "\n"), 0);
  mpfr_inp_str(params->sigma, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"), 0);
  mpfr_inp_str(params->sigma_p, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"), 0);
  mpfr_inp_str(params->sigma_s, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"), 0);
  mpfr_inp_str(params->ell_b, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"), 0);
  mpfr_inp_str(params->ell_g, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"), 0);
  mpfr_inp_str(params->xi, fp, mpfr_base, MPFR_RNDN);
  CHECK(fscanf(fp, "\n"), 0);
  
  gghlite_enc_init(params->pzt, params);
  gghlite_enc_init(params->ntt->w, params);
  gghlite_enc_init(params->ntt->w_inv, params);
  gghlite_enc_init(params->ntt->phi, params);
  gghlite_enc_init(params->ntt->phi_inv, params);

  gghlite_enc_fread(fp, params->pzt);
  CHECK(fscanf(fp, "\n"), 0);
  CHECK(fscanf(fp, "%zd\n", &params->ntt->n), 1);
  gghlite_enc_fread(fp, params->ntt->w);
  CHECK(fscanf(fp, "\n"), 0);
  gghlite_enc_fread(fp, params->ntt->w_inv);
  CHECK(fscanf(fp, "\n"), 0);
  gghlite_enc_fread(fp, params->ntt->phi);
  CHECK(fscanf(fp, "\n"), 0);
  gghlite_enc_fread(fp, params->ntt->phi_inv);

  gghlite_params_set_D_sigmas(params);
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

void gghlite_params_clear_read(gghlite_params_t self) {
  for(int i = 0; i < self->gamma; i++) {
    for(int j = 0; j < self->kappa; j++) {
      free(self->x[i][j]);
    }
    free(self->x[i]);
  }
  free(self->x);
  free(self->y);

  dgsl_rot_mp_clear(self->D_sigma_s);
  dgsl_rot_mp_clear(self->D_sigma_p);

  fmpz_mod_poly_clear(self->pzt);
  mpfr_clear(self->xi);
  mpfr_clear(self->sigma_s);
  mpfr_clear(self->ell_b);
  mpfr_clear(self->sigma_p);
  mpfr_clear(self->ell_g);
  mpfr_clear(self->sigma);
  fmpz_mod_poly_oz_ntt_precomp_clear(self->ntt);
  fmpz_clear(self->q);
}

/**
 * 
 * Members of pp that are not transferred:
 * params->x
 * params->y
 * params->D_sigma_p
 * params->D_sigma_s
 */
void fwrite_mife_pp(mife_pp_t pp, char *filepath) {
  FILE *fp = fopen(filepath, "w");
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
  fmpz_fprint(fp, pp->p);
  fprintf(fp, "\n");

  fwrite_gghlite_params(fp, *pp->params_ref);
  fclose(fp);
}

void fread_mife_pp(mife_pp_t pp, char *filepath) {
  FILE *fp = fopen(filepath, "r");
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
  fmpz_fread(fp, pp->p);
  CHECK(fscanf(fp, "\n"), 0);

  pp->params_ref = malloc(sizeof(gghlite_params_t));
  fread_gghlite_params(fp, *pp->params_ref);
  fclose(fp);
}

void fwrite_mife_sk(mife_sk_t sk, char *filepath) {
	uint64_t t = ggh_walltime(0);
  FILE *fp = fopen(filepath, "w");
	timer_printf("Starting writing Kilian matrices...\n");
  fprintf(fp, "%d\n", sk->numR);
  for(int i = 0; i < sk->numR; i++) {
    fprintf(fp, "%ld %ld\n", sk->R[i]->r, sk->R[i]->c);
    fmpz_mat_fprint(fp, sk->R[i]);
    fprintf(fp, "\n");
    fprintf(fp, "%ld %ld\n", sk->R_inv[i]->r, sk->R_inv[i]->c);
    fmpz_mat_fprint(fp, sk->R_inv[i]);
    fprintf(fp, "\n");
	timer_printf("\r    Progress: [%lu / %lu] %8.2fs",
		i, sk->numR, ggh_seconds(ggh_walltime(t)));

  }
	timer_printf("\n");
	timer_printf("Finished writing Kilian matrices %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	t = ggh_walltime(0);
	timer_printf("Starting writing gghlite params...\n");
  fwrite_gghlite_params(fp, sk->self->params);
  fprintf(fp, "\n");
	timer_printf("Finished writing gghlite params %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	t = ggh_walltime(0);
	timer_printf("Starting writing g, g_inv, h...\n");
  fmpz_poly_fprint(fp, sk->self->g);
  fprintf(fp, "\n");
  fmpq_poly_fprint(fp, sk->self->g_inv);
  fprintf(fp, "\n");
  fmpz_poly_fprint(fp, sk->self->h);
  fprintf(fp, "\n");
	timer_printf("Finished writing g, g_inv, h %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	t = ggh_walltime(0);
	timer_printf("Starting writing z, z_inv, a, b...\n");
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
	timer_printf("\r    Progress: [%lu / %lu] %8.2fs",
		i, sk->self->params->gamma, ggh_seconds(ggh_walltime(t)));
  }
	timer_printf("\n");
	timer_printf("Finished writing z, z_inv, a, b %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

  fclose(fp);

}

void fread_mife_sk(mife_sk_t sk, char *filepath) {
	uint64_t t = ggh_walltime(0);
  FILE *fp = fopen(filepath, "r");
	timer_printf("Starting reading Kilian matrices...\n");
  CHECK(fscanf(fp, "%d\n", &sk->numR), 1);
  sk->R = malloc(sk->numR * sizeof(fmpz_mat_t));
  sk->R_inv = malloc(sk->numR * sizeof(fmpz_mat_t));

  for(int i = 0; i < sk->numR; i++) {
    unsigned long r1, c1, r2, c2;
    CHECK(fscanf(fp, "%ld %ld\n", &r1, &c1), 2);
    fmpz_mat_init(sk->R[i], r1, c1);
    fmpz_mat_fread(fp, sk->R[i]);
    CHECK(fscanf(fp, "\n"), 0);
    CHECK(fscanf(fp, "%ld %ld\n", &r2, &c2), 2);
    fmpz_mat_init(sk->R_inv[i], r2, c2);
    fmpz_mat_fread(fp, sk->R_inv[i]);
    CHECK(fscanf(fp, "\n"), 0);
	timer_printf("\r    Progress: [%lu / %lu] %8.2fs",
		i, sk->numR, ggh_seconds(ggh_walltime(t)));
  }
	timer_printf("\n");
	timer_printf("Finished reading Kilian matrices %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	t = ggh_walltime(0);
	timer_printf("Starting reading gghlite params...\n");
  fread_gghlite_params(fp, sk->self->params);
  CHECK(fscanf(fp, "\n"), 0);
	timer_printf("Finished reading gghlite params %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	t = ggh_walltime(0);
	timer_printf("Starting reading g, g_inv, h...\n");
  fmpz_poly_init(sk->self->g);
  fmpz_poly_fread(fp, sk->self->g);
  CHECK(fscanf(fp, "\n"), 0);
  fmpq_poly_init(sk->self->g_inv);
  fmpq_poly_fread(fp, sk->self->g_inv);
  CHECK(fscanf(fp, "\n"), 0);
  fmpz_poly_init(sk->self->h);
  fmpz_poly_fread(fp, sk->self->h);
  CHECK(fscanf(fp, "\n"), 0);
  	timer_printf("Finished reading g, g_inv, h %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	t = ggh_walltime(0);
	timer_printf("Starting reading z, z_inv, a, b...\n");
  sk->self->z = malloc(sk->self->params->gamma * sizeof(gghlite_enc_t));
  sk->self->z_inv = malloc(sk->self->params->gamma * sizeof(gghlite_enc_t));
  sk->self->a = malloc(sk->self->params->gamma * sizeof(gghlite_clr_t));
  sk->self->b = malloc(sk->self->params->gamma * sizeof(gghlite_clr_t **));
  for(int i = 0; i < sk->self->params->gamma; i++) {
    fmpz_t p1, p2;
    fmpz_init(p1);
    fmpz_init(p2);
    fmpz_fread(fp, p1);
    CHECK(fscanf(fp, "\n"), 0);
    fmpz_mod_poly_init(sk->self->z[i], p1);
    fmpz_mod_poly_fread(fp, sk->self->z[i]);
    CHECK(fscanf(fp, "\n"), 0);
    fmpz_fread(fp, p2);
    CHECK(fscanf(fp, "\n"), 0);
    fmpz_mod_poly_init(sk->self->z_inv[i], p2);
    fmpz_mod_poly_fread(fp, sk->self->z_inv[i]);
    CHECK(fscanf(fp, "\n"), 0);
    fmpz_poly_init(sk->self->a[i]);
    fmpz_poly_fread(fp, sk->self->a[i]);
    CHECK(fscanf(fp, "\n"), 0);
    sk->self->b[i] = malloc(sk->self->params->kappa * sizeof(gghlite_clr_t *));
    for(int j = 0; j < sk->self->params->kappa; j++) {
      sk->self->b[i][j] = malloc(2 * sizeof(gghlite_clr_t));
      fmpz_poly_init(sk->self->b[i][j][0]);
      fmpz_poly_fread(fp, sk->self->b[i][j][0]);
      CHECK(fscanf(fp, "\n"), 0);
      fmpz_poly_init(sk->self->b[i][j][1]);
      fmpz_poly_fread(fp, sk->self->b[i][j][1]);
      CHECK(fscanf(fp, "\n"), 0);
    }
    fmpz_clear(p1);
    fmpz_clear(p2);
	timer_printf("\r    Progress: [%lu / %lu] %8.2fs",
		i, sk->self->params->gamma, ggh_seconds(ggh_walltime(t)));
  }
	timer_printf("\n");
	timer_printf("Finished reading z, z_inv, a, b %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));

	t = ggh_walltime(0);
	timer_printf("Starting setting D_g...\n");
  gghlite_sk_set_D_g(sk->self);
	timer_printf("Finished setting D_g %8.2fs\n",
		ggh_seconds(ggh_walltime(t)));
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

void fread_gghlite_enc_mat(const mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp) {
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
      CHECK(fscanf(fp, "\n"), 0);
    }
  }
}

