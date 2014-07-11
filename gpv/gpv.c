#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>

typedef enum {
  DGS_LATTICE_DETECT               = 0x0,
  DGS_LATTICE_SIMPLE               = 0x2,
  DGS_LATTICE_INLATTICE            = 0x3,
  DGS_LATTICE_COSET                = 0x4,
} dgs_gauss_lattice_alg_t;

struct _dgs_disc_gauss_lattice_mp_t;

typedef struct _dgs_disc_gauss_lattice_mp_t{
  fmpz_mat_t B;
  fmpz *c;
  mpfr_t sigma;
  dgs_disc_gauss_mp_t **D;
  int (*call)(fmpz *rop,  const struct _dgs_disc_gauss_lattice_mp_t *self, gmp_randstate_t state);
  
} dgs_disc_gauss_lattice_mp_t;

int fmpz_mat_is_one(const fmpz_mat_t mat) {
  if (mat->r != mat->c)
    return 0;
  for(long i=0; i<mat->r; i++) {
    for(long j=0; j<i; j++)
      if(!mat->rows[i][j])
        return 0;
    if(!mat->rows[i][i])
      return 0;
    for(long j=i+1; j<mat->c; j++)
      if(!mat->rows[i][j])
        return 0;
  }
  return 1;
}

int dgs_disc_gauss_lattice_mp_call_simple(fmpz *rop,  const dgs_disc_gauss_lattice_mp_t *self, gmp_randstate_t state) {
  assert(rop); assert(self);
  
  const long n = fmpz_mat_ncols(self->B);
  mpz_t  tmp_g;  mpz_init(tmp_g);
  fmpz_t tmp_f;  fmpz_init(tmp_f);

  _fmpz_vec_zero(rop, n);
  
  for(long i=0; i<fmpz_mat_nrows(self->B); i++) {
    self->D[0]->call(tmp_g, self->D[0], state); fmpz_set_mpz(tmp_f, tmp_g);
    _fmpz_vec_scalar_addmul_fmpz(rop, self->B->rows[i], n, tmp_f);
  }

  fmpz_clear(tmp_f);
  mpz_clear(tmp_g);
  
  return 0;
}

dgs_disc_gauss_lattice_mp_t *dgs_disc_gauss_lattice_mp_init(const fmpz_mat_t B, const mpfr_t sigma,
                                                  const fmpz *c, const dgs_gauss_lattice_alg_t algorithm) {
  dgs_disc_gauss_lattice_mp_t *self = (dgs_disc_gauss_lattice_mp_t*)calloc(1, sizeof(dgs_disc_gauss_lattice_mp_t));
  if(!self) dgs_die("out of memory");

  dgs_gauss_lattice_alg_t alg = algorithm;
  
  assert(fmpz_mat_rank(B));
  assert(mpfr_cmp_ui(sigma, 0) > 0);

  
  fmpz_mat_init_set(self->B, B);
  self->c = _fmpz_vec_init(fmpz_mat_ncols(self->B));
  _fmpz_vec_set(self->c, c, fmpz_mat_ncols(self->B));
  mpfr_init_set(self->sigma, sigma, MPFR_RNDN);


  if (alg == DGS_LATTICE_DETECT) {
    if (fmpz_mat_is_one(self->B))
      alg = DGS_LATTICE_SIMPLE;
    else if (_fmpz_vec_is_zero(c, fmpz_mat_ncols(self->B)))
      alg = DGS_LATTICE_INLATTICE;
    else
      alg = DGS_LATTICE_COSET; //TODO: we could test for lattice membership here
  }

  switch(alg) {
  case DGS_LATTICE_SIMPLE:
    self->D = (dgs_disc_gauss_mp_t**)calloc(2, sizeof(dgs_disc_gauss_mp_t*));
    mpfr_t c;
    mpfr_init(c);
    self->D[0] = dgs_disc_gauss_mp_init(self->sigma, c, 6, DGS_DISC_GAUSS_DEFAULT);
    mpfr_clear(c);
    self->call = dgs_disc_gauss_lattice_mp_call_simple;
    break;
  default:
    dgs_die("not implemented");
  }
}


void dgs_disc_gauss_lattice_mp_clear(dgs_disc_gauss_lattice_mp_t *self) {
  if (!fmpz_mat_is_empty(self->B)) {
    _fmpz_vec_clear(self->c, fmpz_mat_ncols(self->B));
    fmpz_mat_clear(self->B);
  }
  if(self->D) {
    dgs_disc_gauss_mp_clear(self->D[0]);
    free(self->D);
  }
  mpfr_clear(self->sigma);
  free(self);
}
