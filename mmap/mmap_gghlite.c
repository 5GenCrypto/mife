#include "mmap_gghlite.h"

const mmap_pp_vtable gghlite_pp_vtable =
  { .clear = gghlite_params_clear_read_wrapper
  , .fread = fread_gghlite_params_wrapper
  , .fwrite = fwrite_gghlite_params_wrapper
  , .size = sizeof(mmap_pp)
  };

const mmap_sk_vtable gghlite_sk_vtable =
  { .init = gghlite_jigsaw_init_gamma_wrapper
  , .clear = gghlite_sk_clear_wrapper
  , .fread = fread_gghlite_sk_wrapper
  , .fwrite = fwrite_gghlite_sk_wrapper
  , .size = sizeof(mmap_sk)
  , .pp = gghlite_sk_to_pp
  , .plaintext_field = fmpz_poly_oz_ideal_norm_wrapper
  };

const mmap_enc_vtable gghlite_enc_vtable =
  { .init = gghlite_enc_init_wrapper
  , .clear = gghlite_enc_clear_wrapper
  , .fread = gghlite_enc_fread_raw_wrapper
  , .fwrite = gghlite_enc_fprint_raw_wrapper
  , .size = sizeof(mmap_enc)
  , .set = gghlite_enc_set_wrapper
  , .add = gghlite_enc_add_wrapper
  , .mul = gghlite_enc_mul_wrapper
  , .is_zero = gghlite_enc_is_zero_wrapper
  , .encode = gghlite_enc_set_gghlite_clr_wrapper
  };

const mmap_vtable gghlite_vtable =
  { .pp  = &gghlite_pp_vtable
  , .sk  = &gghlite_sk_vtable
  , .enc = &gghlite_enc_vtable
  };

void gghlite_params_clear_read_wrapper(mmap_pp *pp)
  { gghlite_params_clear_read(pp->gghlite_self); }
void fread_gghlite_params_wrapper(mmap_pp *const pp, FILE *const fp)
  { fread_gghlite_params(fp, pp->gghlite_self); }
void fwrite_gghlite_params_wrapper(const mmap_pp *const pp, FILE *const fp)
  { fwrite_gghlite_params(fp, pp->gghlite_self); }

void gghlite_jigsaw_init_gamma_wrapper(mmap_sk *const sk, size_t lambda, size_t kappa, size_t gamma, aes_randstate_t randstate) {
  const gghlite_flag_t flags = GGHLITE_FLAGS_GOOD_G_INV | GGHLITE_FLAGS_QUIET;
  gghlite_jigsaw_init_gamma(sk->gghlite_self, lambda, kappa, gamma, flags, randstate);
}

void gghlite_sk_clear_wrapper(mmap_sk *const sk)
  { gghlite_sk_clear(sk->gghlite_self, 1); }
void fread_gghlite_sk_wrapper(mmap_sk *const sk, FILE *const fp)
  { fread_gghlite_sk(fp, sk->gghlite_self); }
void fwrite_gghlite_sk_wrapper(const mmap_sk *const sk, FILE *const fp)
  { fwrite_gghlite_sk(fp, sk->gghlite_self); }
const mmap_pp *const gghlite_sk_to_pp(const mmap_sk *const sk)
  /* N.B. This cast is strictly speaking probably not okay from a "portable C"
   * standpoint. However it's almost certainly going to be fine with all the
   * compilers we care about... */
  { return (mmap_pp *)sk->gghlite_self->params; }
void fmpz_poly_oz_ideal_norm_wrapper(const mmap_sk *const sk, fmpz_t p_out)
  { fmpz_poly_oz_ideal_norm(p_out, sk->gghlite_self->g, sk->gghlite_self->params->n, 0); }

void gghlite_enc_init_wrapper(mmap_enc *const enc, const mmap_pp *const pp)
  { gghlite_enc_init(enc->gghlite_self, pp->gghlite_self); }
void gghlite_enc_clear_wrapper(mmap_enc *const enc)
  { gghlite_enc_clear(enc->gghlite_self); }

void gghlite_enc_fread_raw_wrapper(mmap_enc *const enc, FILE *const fp) {
  const int tmp = gghlite_enc_fread_raw(fp, enc->gghlite_self);
  assert(tmp > 0);
}

void gghlite_enc_fprint_raw_wrapper(const mmap_enc *const enc, FILE *const fp)
  { gghlite_enc_fprint_raw(fp, enc->gghlite_self); }
void gghlite_enc_set_wrapper(mmap_enc *const dest, const mmap_enc *const src)
  { gghlite_enc_set(dest->gghlite_self, src->gghlite_self); }
void gghlite_enc_add_wrapper(mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b)
  { gghlite_enc_add(dest->gghlite_self, pp->gghlite_self, a->gghlite_self, b->gghlite_self); }
void gghlite_enc_mul_wrapper(mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b)
  { gghlite_enc_mul(dest->gghlite_self, pp->gghlite_self, a->gghlite_self, b->gghlite_self); }
bool gghlite_enc_is_zero_wrapper(const mmap_enc *const enc, const mmap_pp *const pp)
  { return gghlite_enc_is_zero(pp->gghlite_self, enc->gghlite_self); }

void gghlite_enc_set_gghlite_clr_wrapper(mmap_enc *const enc, const mmap_sk *const sk, const fmpz_t plaintext, int *group, aes_randstate_t randstate)
{
  gghlite_clr_t e;
  gghlite_clr_init(e);
  fmpz_poly_set_coeff_fmpz(e, 0, plaintext);
  gghlite_enc_set_gghlite_clr(enc->gghlite_self, sk->gghlite_self, e, 1, group, 0, randstate);
  gghlite_clr_clear(e);
}

void gghlite_params_clear_read(gghlite_params_t gghlite_self) {
  for(int i = 0; i < gghlite_self->gamma; i++) {
    for(int j = 0; j < gghlite_self->kappa; j++) {
      free(gghlite_self->x[i][j]);
    }
    free(gghlite_self->x[i]);
  }
  free(gghlite_self->x);
  free(gghlite_self->y);

  dgsl_rot_mp_clear(gghlite_self->D_sigma_s);
  dgsl_rot_mp_clear(gghlite_self->D_sigma_p);

  fmpz_mod_poly_clear(gghlite_self->pzt);
  mpfr_clear(gghlite_self->xi);
  mpfr_clear(gghlite_self->sigma_s);
  mpfr_clear(gghlite_self->ell_b);
  mpfr_clear(gghlite_self->sigma_p);
  mpfr_clear(gghlite_self->ell_g);
  mpfr_clear(gghlite_self->sigma);
  fmpz_mod_poly_oz_ntt_precomp_clear(gghlite_self->ntt);
  fmpz_clear(gghlite_self->q);
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
  ), 7);

  gghlite_params_initzero(params, lambda, kappa, gamma);
  params->n = n;
  params->ell = ell;
  params->rerand_mask = rerand_mask;
  params->flags = gghlite_flag_int;

  fmpz_inp_raw(params->q, fp);
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

  gghlite_enc_fread_raw(fp, params->pzt);
  CHECK(fscanf(fp, "\n"), 0);
  CHECK(fscanf(fp, "%zd\n", &params->ntt->n), 1);
  gghlite_enc_fread_raw(fp, params->ntt->w);
  CHECK(fscanf(fp, "\n"), 0);
  gghlite_enc_fread_raw(fp, params->ntt->w_inv);
  CHECK(fscanf(fp, "\n"), 0);
  gghlite_enc_fread_raw(fp, params->ntt->phi);
  CHECK(fscanf(fp, "\n"), 0);
  gghlite_enc_fread_raw(fp, params->ntt->phi_inv);

  gghlite_params_set_D_sigmas(params);
}

void fwrite_gghlite_params(FILE *fp, const gghlite_params_t params) {
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
  fmpz_out_raw(fp, params->q);
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
  gghlite_enc_fprint_raw(fp, params->pzt);
  fprintf(fp, "\n");
  fprintf(fp, "%zd\n", params->ntt->n);
  fmpz_mod_poly_fprint_raw(fp, params->ntt->w);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint_raw(fp, params->ntt->w_inv);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint_raw(fp, params->ntt->phi);
  fprintf(fp, "\n");
  fmpz_mod_poly_fprint_raw(fp, params->ntt->phi_inv);
}

void fread_gghlite_sk(FILE *fp, gghlite_sk_t gghlite_self) {
  uint64_t t = ggh_walltime(0);
  timer_printf("Starting reading gghlite params...\n");
  fread_gghlite_params(fp, gghlite_self->params);
  CHECK(fscanf(fp, "\n"), 0);
  timer_printf("Finished reading gghlite params %8.2fs\n",
    ggh_seconds(ggh_walltime(t)));

  t = ggh_walltime(0);
  timer_printf("Starting reading g, g_inv, h...\n");
  fmpz_poly_init(gghlite_self->g);
  fmpz_poly_fread_raw(fp, gghlite_self->g);
  CHECK(fscanf(fp, "\n"), 0);
  fmpq_poly_init(gghlite_self->g_inv);
  fmpq_poly_fread(fp, gghlite_self->g_inv);
  CHECK(fscanf(fp, "\n"), 0);
  fmpz_poly_init(gghlite_self->h);
  fmpz_poly_fread_raw(fp, gghlite_self->h);
  CHECK(fscanf(fp, "\n"), 0);
    timer_printf("Finished reading g, g_inv, h %8.2fs\n",
    ggh_seconds(ggh_walltime(t)));

  t = ggh_walltime(0);
  timer_printf("Starting reading z, z_inv, a, b...\n");
  gghlite_self->z = malloc(gghlite_self->params->gamma * sizeof(gghlite_enc_t));
  gghlite_self->z_inv = malloc(gghlite_self->params->gamma * sizeof(gghlite_enc_t));
  gghlite_self->a = malloc(gghlite_self->params->gamma * sizeof(gghlite_clr_t));
  gghlite_self->b = malloc(gghlite_self->params->gamma * sizeof(gghlite_clr_t **));
  for(int i = 0; i < gghlite_self->params->gamma; i++) {
    fmpz_t p1, p2;
    fmpz_init(p1);
    fmpz_init(p2);
    fmpz_inp_raw(p1, fp);
    CHECK(fscanf(fp, "\n"), 0);
    fmpz_mod_poly_fread_raw(fp, gghlite_self->z[i]);
    CHECK(fscanf(fp, "\n"), 0);
    fmpz_inp_raw(p2, fp);
    CHECK(fscanf(fp, "\n"), 0);
    fmpz_mod_poly_fread_raw(fp, gghlite_self->z_inv[i]);
    CHECK(fscanf(fp, "\n"), 0);
    fmpz_poly_init(gghlite_self->a[i]);
    fmpz_poly_fread_raw(fp, gghlite_self->a[i]);
    CHECK(fscanf(fp, "\n"), 0);
    gghlite_self->b[i] = malloc(gghlite_self->params->kappa * sizeof(gghlite_clr_t *));
    for(int j = 0; j < gghlite_self->params->kappa; j++) {
      gghlite_self->b[i][j] = malloc(2 * sizeof(gghlite_clr_t));
      fmpz_poly_init(gghlite_self->b[i][j][0]);
      fmpz_poly_fread_raw(fp, gghlite_self->b[i][j][0]);
      CHECK(fscanf(fp, "\n"), 0);
      fmpz_poly_init(gghlite_self->b[i][j][1]);
      fmpz_poly_fread_raw(fp, gghlite_self->b[i][j][1]);
      CHECK(fscanf(fp, "\n"), 0);
    }
    fmpz_clear(p1);
    fmpz_clear(p2);
  timer_printf("\r    Progress: [%lu / %lu] %8.2fs",
    i, gghlite_self->params->gamma, ggh_seconds(ggh_walltime(t)));
  }
  timer_printf("\n");
  timer_printf("Finished reading z, z_inv, a, b %8.2fs\n",
    ggh_seconds(ggh_walltime(t)));

  t = ggh_walltime(0);
  timer_printf("Starting setting D_g...\n");
  gghlite_sk_set_D_g(gghlite_self);
  timer_printf("Finished setting D_g %8.2fs\n",
    ggh_seconds(ggh_walltime(t)));
}

void fwrite_gghlite_sk(FILE *fp, const gghlite_sk_t gghlite_self) {
  uint64_t t = ggh_walltime(0);
  timer_printf("Starting writing gghlite params...\n");
  fwrite_gghlite_params(fp, gghlite_self->params);
  fprintf(fp, "\n");
  timer_printf("Finished writing gghlite params %8.2fs\n",
    ggh_seconds(ggh_walltime(t)));

  t = ggh_walltime(0);
  timer_printf("Starting writing g, g_inv, h...\n");
  fmpz_poly_fprint_raw(fp, gghlite_self->g);
  fprintf(fp, "\n");
  fmpq_poly_fprint(fp, gghlite_self->g_inv);
  fprintf(fp, "\n");
  fmpz_poly_fprint_raw(fp, gghlite_self->h);
  fprintf(fp, "\n");
  timer_printf("Finished writing g, g_inv, h %8.2fs\n",
    ggh_seconds(ggh_walltime(t)));

  t = ggh_walltime(0);
  timer_printf("Starting writing z, z_inv, a, b...\n");
  for(int i = 0; i < gghlite_self->params->gamma; i++) {
    fmpz_out_raw(fp, fmpz_mod_poly_modulus(gghlite_self->z[i]));
    fprintf(fp, "\n");
    fmpz_mod_poly_fprint_raw(fp, gghlite_self->z[i]);
    fprintf(fp, "\n");
    fmpz_out_raw(fp, fmpz_mod_poly_modulus(gghlite_self->z_inv[i]));
    fprintf(fp, "\n");
    fmpz_mod_poly_fprint_raw(fp, gghlite_self->z_inv[i]);
    fprintf(fp, "\n");
    fmpz_poly_fprint_raw(fp, gghlite_self->a[i]);
    fprintf(fp, "\n");
    for(int j = 0; j < gghlite_self->params->kappa; j++) {
      fmpz_poly_fprint_raw(fp, gghlite_self->b[i][j][0]);
      fprintf(fp, "\n");
      fmpz_poly_fprint_raw(fp, gghlite_self->b[i][j][1]);
      fprintf(fp, "\n");
    }
  timer_printf("\r    Progress: [%lu / %lu] %8.2fs",
    i, gghlite_self->params->gamma, ggh_seconds(ggh_walltime(t)));
  }
  timer_printf("\n");
  timer_printf("Finished writing z, z_inv, a, b %8.2fs\n",
  ggh_seconds(ggh_walltime(t)));
}

int fmpz_poly_fread_raw(FILE * file, fmpz_poly_t poly)
{
    int r;
    slong i, len;
    mpz_t t;

    mpz_init(t);
    r = mpz_inp_str(t, file, 10);
    if (r == 0)
    {
        mpz_clear(t);
        return 0;
    }
    if (!mpz_fits_slong_p(t))
    {
        flint_printf("Exception (fmpz_poly_fread). Length does not fit into a slong.\n");
        abort();
    }
    len = flint_mpz_get_si(t);
    mpz_clear(t);

    fmpz_poly_fit_length(poly, len);

    for (i = 0; i < len; i++)
    {
        r = fmpz_inp_raw(poly->coeffs + i, file);
        if (r <= 0)
            return r;
    }

    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);

    return 1;
}


int _fmpz_poly_fprint_raw(FILE * file, const fmpz * vec, slong len)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%li", len);
    if ((len > 0) && (r > 0))
    {
        for (i = 0; (i < len) && (r > 0); i++)
        {
            r = fmpz_out_raw(file, vec + i);
        }
    }

    return r;
}

int fmpz_poly_fprint_raw(FILE * file, const fmpz_poly_t poly)
{
    return _fmpz_poly_fprint_raw(file, poly->coeffs, poly->length);
}



int _fmpz_mod_poly_fprint_raw(FILE * file, const fmpz *poly, slong len,
                          const fmpz_t p)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%wd ", len);
    if (r <= 0)
        return r;

    r = fmpz_out_raw(file, p);
    if (r <= 0)
        return r;

    if (len == 0)
        return r;

    r = flint_fprintf(file, " ");
    if (r <= 0)
        return r;

    for (i = 0; (r > 0) && (i < len); i++)
    {
        r = flint_fprintf(file, " ");
        if (r <= 0)
            return r;
        r = fmpz_out_raw(file, poly + i);
        if (r <= 0)
            return r;
    }

    return r;
}

int fmpz_mod_poly_fprint_raw(FILE * file, const fmpz_mod_poly_t poly)
{
    return _fmpz_mod_poly_fprint_raw(file, poly->coeffs, poly->length,
        &(poly->p));
}

/* copied from fmpz_mod_poly_fread, mostly (but fixed) */
int fmpz_mod_poly_fread_raw(FILE * f, fmpz_mod_poly_t poly)
{
    slong i, length;
    fmpz_t coeff;
    ulong res;

    fmpz_init(coeff);
    if (flint_fscanf(f, "%wd ", &length) != 1) {
        fmpz_clear(coeff);
        return 0;
    }

    fmpz_inp_raw(coeff,f);
    fmpz_mod_poly_init(poly, coeff);
    fmpz_mod_poly_fit_length(poly, length);

    poly->length = length;
    flint_fscanf(f, " ");

    for (i = 0; i < length; i++)
    {
        flint_fscanf(f, " ");
        res = fmpz_inp_raw(coeff, f);

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
