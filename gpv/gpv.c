#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint-addons/flint-addons.h>
#include "gso.h"

int dgs_disc_gauss_lattice_mp_call_identity(fmpz *rop,  const dgs_disc_gauss_lattice_mp_t *self, gmp_randstate_t state) {
  assert(rop); assert(self);

  const long n = fmpz_mat_ncols(self->B);
  mpz_t  tmp_g;  mpz_init(tmp_g);

  _fmpz_vec_zero(rop, n);

  for(long i=0; i<n; i++) {
    self->D[0]->call(tmp_g, self->D[0], state);
    fmpz_set_mpz(rop+i, tmp_g);
  }

  mpz_clear(tmp_g);

  return 0;
}

int dgs_disc_gauss_lattice_mp_call_inlattice(fmpz *rop,  const dgs_disc_gauss_lattice_mp_t *self, gmp_randstate_t state) {
  assert(rop); assert(self);

  const long m = fmpz_mat_nrows(self->B);
  const long n = fmpz_mat_ncols(self->B);
  mpz_t  tmp_g; mpz_init(tmp_g);
  fmpz_t tmp_f; fmpz_init(tmp_f);

  _fmpz_vec_zero(rop, n);

  if (fmpz_mat_nrows(self->B) == 1) {
    dgs_die("not implemented");
  } else {
    for(long i=0; i<m; i++) {
      self->D[i]->call(tmp_g, self->D[i], state);
      fmpz_set_mpz(tmp_f, tmp_g);
      _fmpz_vec_scalar_addmul_fmpz(rop, self->B->rows[i], n, tmp_f);
    }
    if (self->c)
      _fmpz_vec_add(rop, rop, self->c, n);
  }

  mpz_clear(tmp_g);
  fmpz_clear(tmp_f);

  return 0;
}

dgs_disc_gauss_lattice_mp_t *dgs_disc_gauss_lattice_mp_init(const fmpz_mat_t B, const mpfr_t sigma,
                                                  const fmpz *c, const dgs_gauss_lattice_alg_t algorithm) {
  assert(mpfr_cmp_ui(sigma, 0) > 0);

  dgs_disc_gauss_lattice_mp_t *self = (dgs_disc_gauss_lattice_mp_t*)calloc(1, sizeof(dgs_disc_gauss_lattice_mp_t));
  if(!self) dgs_die("out of memory");

  dgs_gauss_lattice_alg_t alg = algorithm;

  const long m = fmpz_mat_nrows(B);
  const long n = fmpz_mat_ncols(B);
  const mpfr_prec_t prec = mpfr_get_prec(sigma);

  fmpz_mat_init_set(self->B, B);
  self->c = _fmpz_vec_init(fmpz_mat_ncols(self->B));
  if (c)
    _fmpz_vec_set(self->c, c, fmpz_mat_ncols(self->B));
  mpfr_init2(self->sigma, prec);
  mpfr_set(self->sigma, sigma, MPFR_RNDN);

  if (alg == DGS_LATTICE_DETECT) {
    if (fmpz_mat_is_one(self->B) || fmpz_mat_rot_is_one(self->B)) {
      alg = DGS_LATTICE_IDENTITY;
    }    else if (_fmpz_vec_is_zero(self->c, n))
      alg = DGS_LATTICE_INLATTICE;
    else
      alg = DGS_LATTICE_COSET; //TODO: we could test for lattice membership here
  }

  mpfr_t c_;
  mpfr_init2(c_, prec);

  switch(alg) {
  case DGS_LATTICE_IDENTITY:
    self->D = (dgs_disc_gauss_mp_t**)calloc(1, sizeof(dgs_disc_gauss_mp_t*));
    mpfr_set_d(c_, 0.0, MPFR_RNDN);
    self->D[0] = dgs_disc_gauss_mp_init(self->sigma, c_, 6, DGS_DISC_GAUSS_DEFAULT);
    self->call = dgs_disc_gauss_lattice_mp_call_identity;
    break;

  case DGS_LATTICE_INLATTICE:
    self->D = (dgs_disc_gauss_mp_t**)calloc(m, sizeof(dgs_disc_gauss_mp_t*));

    mpfr_mat_t G;
    mpfr_mat_init(G, m, n, prec);
    mpfr_mat_set_fmpz_mat(G, B);

    mpfr_t sigma_;
    mpfr_init2(sigma_, prec);

    mpfr_t norm;
    mpfr_init2(norm, prec);

    mpfr_set_d(c_, 0.0, MPFR_RNDN);

    for(long i=0; i<m; i++) {
      _mpfr_vec_2norm(norm, (const mpfr_t*)G->rows[i], n);
      mpfr_div(sigma_, self->sigma, norm, MPFR_RNDN);
      self->D[i] = dgs_disc_gauss_mp_init(sigma_, c_, 6, DGS_DISC_GAUSS_DEFAULT);
    }

    self->call = dgs_disc_gauss_lattice_mp_call_inlattice;

    mpfr_clear(sigma_);
    mpfr_clear(norm);
    mpfr_mat_clear(G);
    break;
  default:
    dgs_die("not implemented");
  }

  mpfr_clear(c_);

  return self;
}


void dgs_disc_gauss_lattice_mp_clear(dgs_disc_gauss_lattice_mp_t *self) {
  if (self->c)
    _fmpz_vec_clear(self->c, fmpz_mat_ncols(self->B));
  if (!fmpz_mat_is_empty(self->B)) {
    fmpz_mat_clear(self->B);
  }

  if(self->call == dgs_disc_gauss_lattice_mp_call_identity) {
    if(self->D) {
      dgs_disc_gauss_mp_clear(self->D[0]);
      free(self->D);
    }
  }
  mpfr_clear(self->sigma);
  free(self);
}
