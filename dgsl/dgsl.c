#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint-addons/flint-addons.h>
#include <flint-addons/cyclotomic-2power.h>
#include "dgsl.h"
#include "gso.h"

int dgsl_mp_call_identity(fmpz *rop,  const dgsl_mp_t *self, gmp_randstate_t state) {
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

int dgsl_mp_call_inlattice(fmpz *rop,  const dgsl_mp_t *self, gmp_randstate_t state) {
  assert(rop); assert(self);

  const long m = fmpz_mat_nrows(self->B);
  const long n = fmpz_mat_ncols(self->B);
  mpz_t  tmp_g; mpz_init(tmp_g);
  fmpz_t tmp_f; fmpz_init(tmp_f);

  _fmpz_vec_zero(rop, n);

  for(long i=0; i<m; i++) {
    self->D[i]->call(tmp_g, self->D[i], state);
    fmpz_set_mpz(tmp_f, tmp_g);
    _fmpz_vec_scalar_addmul_fmpz(rop, self->B->rows[i], n, tmp_f);
  }
  _fmpz_vec_add(rop, rop, self->c_z, n);

  mpz_clear(tmp_g);
  fmpz_clear(tmp_f);

  return 0;
}

int dgsl_mp_call_coset(fmpz *rop, const dgsl_mp_t *self, gmp_randstate_t state) {
  assert(rop); assert(self);

  const long n = fmpz_mat_ncols(self->B);
  _fmpz_vec_zero(rop, n);

  mpfr_t *c = _mpfr_vec_init(n, mpfr_get_prec(self->sigma));
  _mpfr_vec_set(c, self->c, n, MPFR_RNDN);

  mpfr_t c_prime;
  mpfr_init2(c_prime, mpfr_get_prec(self->sigma));

  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(self->sigma));

  mpfr_t sigma_prime;
  mpfr_init2(sigma_prime, mpfr_get_prec(self->sigma));

  mpz_t z;
  mpz_init(z);
  
  mpfr_t z_mpfr;
  mpfr_init2(z_mpfr, mpfr_get_prec(self->sigma));

  fmpz_t z_fmpz;
  fmpz_init(z_fmpz);
  
  size_t tau = 3;
  if (ceil(sqrt(log2((double)n))) > tau)
    tau = ceil(sqrt(log2((double)n)));

  mpfr_t *b = _mpfr_vec_init(n, mpfr_get_prec(self->sigma));

  if (fmpz_mat_nrows(self->B)>1) {
    const long m = fmpz_mat_nrows(self->B);
    for(long j=0; j<m; j++) {
      long i = m-i-1;
      _mpfr_vec_dot_product(c_prime, c,                self->G->rows[i], n, MPFR_RNDN);
      _mpfr_vec_dot_product(tmp,     self->G->rows[i], self->G->rows[i], n, MPFR_RNDN);
      mpfr_div(c_prime, c_prime, tmp, MPFR_RNDN);

      mpfr_sqrt(tmp, tmp, MPFR_RNDN);
      mpfr_div(sigma_prime, self->sigma, tmp, MPFR_RNDN);

      assert(mpfr_cmp_d(sigma_prime, 0.0) > 0);
      dgs_disc_gauss_mp_t *D = dgs_disc_gauss_mp_init(sigma_prime, c_prime, tau, DGS_DISC_GAUSS_UNIFORM_ONLINE);
      D->call(z, D, state);
      dgs_disc_gauss_mp_clear(D);

      mpfr_set_z(z_mpfr, z, MPFR_RNDN);
      mpfr_neg(z_mpfr, z_mpfr, MPFR_RNDN);
      _mpfr_vec_set_fmpz_vec(b, self->B->rows[i], n, MPFR_RNDN);
      _mpfr_vec_scalar_addmul_mpfr(c, b, n, z_mpfr, MPFR_RNDN);

      fmpz_set_mpz(z_fmpz, z);     
      _fmpz_vec_scalar_addmul_fmpz(rop, self->B->rows[i], n, z_fmpz);
    }
  } else {
    dgs_die("not implemented");
  }

  fmpz_clear(z_fmpz);
  mpfr_clear(z_mpfr);
  mpfr_clear(sigma_prime);
  mpfr_clear(tmp);
  mpfr_clear(c_prime);
  _mpfr_vec_clear(c, n);
  _mpfr_vec_clear(b, n);
  return 0;
}

dgsl_mp_t *dgsl_mp_init(const fmpz_mat_t B, mpfr_t sigma,
                      mpfr_t *c, const dgsl_alg_t algorithm) {
  assert(mpfr_cmp_ui(sigma, 0) > 0);

  dgsl_mp_t *self = (dgsl_mp_t*)calloc(1, sizeof(dgsl_mp_t));
  if(!self) dgs_die("out of memory");

  dgsl_alg_t alg = algorithm;

  long m = fmpz_mat_nrows(B);
  long n = fmpz_mat_ncols(B);

  const mpfr_prec_t prec = mpfr_get_prec(sigma);

  fmpz_mat_init_set(self->B, B);
  self->c_z = _fmpz_vec_init(n);
  self->c   = _mpfr_vec_init(n, prec);
  mpfr_init2(self->sigma, prec);
  mpfr_set(self->sigma, sigma, MPFR_RNDN);

  if (alg == DGSL_DETECT) {
    if (fmpz_mat_is_one(self->B)) {
      alg = DGSL_IDENTITY;
    } else if (_mpfr_vec_is_zero(c, n))
      alg = DGSL_INLATTICE;
    else
      alg = DGSL_COSET; //TODO: we could test for lattice membership here
  }

  mpfr_t c_;
  mpfr_init2(c_, prec);

  size_t tau = 3;
  if (2*ceil(sqrt(log2((double)n))) > tau)
    tau = 2*ceil(sqrt(log2((double)n)));

  switch(alg) {
  case DGSL_IDENTITY:
    self->D = (dgs_disc_gauss_mp_t**)calloc(1, sizeof(dgs_disc_gauss_mp_t*));
    mpfr_set_d(c_, 0.0, MPFR_RNDN);
    self->D[0] = dgs_disc_gauss_mp_init(self->sigma, c_, tau, DGS_DISC_GAUSS_DEFAULT);

    self->call = dgsl_mp_call_identity;
    break;
    
  case DGSL_INLATTICE:
    self->D = (dgs_disc_gauss_mp_t**)calloc(m, sizeof(dgs_disc_gauss_mp_t*));

    if (c)
      _fmpz_vec_set_mpfr_vec(self->c_z, c, n);
    
    mpfr_mat_t G;
    mpfr_mat_init(G, m, n, prec);
    mpfr_mat_set_fmpz_mat(G, B);

    mpfr_mat_gso(G, MPFR_RNDN);
    
    mpfr_t sigma_;
    mpfr_init2(sigma_, prec);

    mpfr_t norm;
    mpfr_init2(norm, prec);

    mpfr_set_d(c_, 0.0, MPFR_RNDN);

    for(long i=0; i<m; i++) {
      _mpfr_vec_2norm(norm, G->rows[i], n, MPFR_RNDN);
      assert(mpfr_cmp_d(norm, 0.0) > 0);
      mpfr_div(sigma_, self->sigma, norm, MPFR_RNDN);
      assert(mpfr_cmp_d(sigma_, 0.0) > 0);
      self->D[i] = dgs_disc_gauss_mp_init(sigma_, c_, tau, DGS_DISC_GAUSS_DEFAULT);
    }

    mpfr_clear(sigma_);
    mpfr_clear(norm);
    mpfr_mat_clear(G);

    self->call = dgsl_mp_call_inlattice;
    break;
    
  case DGSL_COSET:
    mpfr_mat_init(self->G, m, n, prec);
    mpfr_mat_set_fmpz_mat(self->G, B);
    mpfr_mat_gso(self->G, MPFR_RNDN);

    self->call = dgsl_mp_call_coset;
    break;
  default:
    dgs_die("not implemented");
  }

  mpfr_clear(c_);

  return self;
}


void dgsl_mp_clear(dgsl_mp_t *self) {
  if (!self)
    return;
  _fmpz_vec_clear(self->c_z, fmpz_mat_ncols(self->B));
  _mpfr_vec_clear(self->c,   fmpz_mat_ncols(self->B));
  
  if (!fmpz_mat_is_empty(self->B)) {
    fmpz_mat_clear(self->B);
  }

  if (!mpfr_mat_is_empty(self->G)) {
    mpfr_mat_clear(self->G);
  }
  if(self->call == dgsl_mp_call_identity) {
    if(self->D) {
      dgs_disc_gauss_mp_clear(self->D[0]);
      free(self->D);
    }
  }
  mpfr_clear(self->sigma);
  free(self);
}

/** Rotational Basis **/

dgsl_rot_mp_t *dgsl_rot_mp_init(const fmpz_mat_t B, mpfr_t sigma,
                            mpfr_t *c, const dgsl_alg_t algorithm) {
  assert(mpfr_cmp_ui(sigma, 0) > 0);

  dgsl_rot_mp_t *self = (dgsl_rot_mp_t*)calloc(1, sizeof(dgsl_mp_t));
  if(!self) dgs_die("out of memory");

  dgsl_alg_t alg = algorithm;

  long n = fmpz_mat_ncols(B);

  if (fmpz_mat_nrows(B) != 1)
    dgs_die("only rotational bases are supported, call dgsl_mp_init instead.");

  const mpfr_prec_t prec = mpfr_get_prec(sigma);

  fmpz_mat_init_set(self->B, B);
  self->c_z = _fmpz_vec_init(n);
  self->c   = _mpfr_vec_init(n, prec);
  mpfr_init2(self->sigma, prec);
  mpfr_set(self->sigma, sigma, MPFR_RNDN);

  if (alg == DGSL_DETECT) {
    if (fmpz_mat_rot_is_one(self->B)) {
      alg = DGSL_IDENTITY;
    } else if (_mpfr_vec_is_zero(c, n))
      alg = DGSL_INLATTICE;
    else
      alg = DGSL_COSET; //TODO: we could test for lattice membership here
  }

  mpfr_t c_;
  mpfr_init2(c_, prec);

  size_t tau = 3;
  if (2*ceil(sqrt(log2((double)n))) > tau)
    tau = 2*ceil(sqrt(log2((double)n)));

  switch(alg) {
  case DGSL_IDENTITY:
    self->D = (dgs_disc_gauss_mp_t**)calloc(1, sizeof(dgs_disc_gauss_mp_t*));
    mpfr_set_d(c_, 0.0, MPFR_RNDN);
    self->D[0] = dgs_disc_gauss_mp_init(self->sigma, c_, tau, DGS_DISC_GAUSS_DEFAULT);

    self->call = dgsl_rot_mp_call_identity;
    break;
    
  case DGSL_INLATTICE:
    self->D = (dgs_disc_gauss_mp_t**)calloc(n, sizeof(dgs_disc_gauss_mp_t*));

    if (c)
      _fmpz_vec_set_mpfr_vec(self->c_z, c, n);
    
    mpfr_mat_t G;
    mpfr_mat_init(G, n, n, prec);
    mpfr_mat_set_fmpz_mat_rot(G, B);
    mpfr_mat_gso(G, MPFR_RNDN);
    
    mpfr_t sigma_;
    mpfr_init2(sigma_, prec);

    mpfr_t norm;
    mpfr_init2(norm, prec);

    mpfr_set_d(c_, 0.0, MPFR_RNDN);

    for(long i=0; i<n; i++) {
      _mpfr_vec_2norm(norm, G->rows[i], n, MPFR_RNDN);
      assert(mpfr_cmp_d(norm, 0.0) > 0);
      mpfr_div(sigma_, self->sigma, norm, MPFR_RNDN);
      assert(mpfr_cmp_d(sigma_, 0.0) > 0);
      self->D[i] = dgs_disc_gauss_mp_init(sigma_, c_, tau, DGS_DISC_GAUSS_DEFAULT);
    }

    mpfr_clear(sigma_);
    mpfr_clear(norm);
    mpfr_mat_clear(G);

    self->call = dgsl_rot_mp_call_inlattice;
    break;
    
  case DGSL_COSET:
    dgs_die("not implemented");
    break;
  default:
    dgs_die("not implemented");
  }

  mpfr_clear(c_);

  return self;
}



/**
   sqrt(Σ_2) with Σ_2 = Σ - Σ_1 = σ^2·I - r·g·g^T
*/

void _dgsl_rot_mp_sqrt_sigma_2(fmpq_poly_t rop, const fmpq_poly_t g, const mpfr_t sigma, const int r, const long n, const mpfr_prec_t prec) { 
  fmpq_t tmp;
  fmpq_init(tmp);

  mpf_t tmp_f;
  mpq_t tmp_q;
  mpf_init2(tmp_f, mpfr_get_prec(sigma));
  mpq_init(tmp_q);
  
  // TODO: might as well get mpq out
  mpfr_get_f(tmp_f, sigma, MPFR_RNDN);
  mpq_set_f(tmp_q, tmp_f);
  fmpq_set_mpq(tmp, tmp_q);
  fmpq_pow_si(tmp, tmp, 2);
  mpf_clear(tmp_f);
  mpq_clear(tmp_q);

  fmpq_poly_set_coeff_fmpq(rop, 0, tmp);
  fmpq_clear(tmp);
  
  fmpq_poly_t ggT;
  fmpq_poly_init(ggT);
  fmpq_poly_mul_mod_cnf2pow_transpose(ggT, g, n);
  fmpq_poly_scalar_mul_si(ggT, ggT, r);
  fmpq_poly_scalar_mul_si(ggT, ggT, r);
  fmpq_poly_sub(rop, rop, ggT);
  fmpq_poly_clear(ggT);  

  fmpq_poly_sqrt_mod_cnf2power_approx(rop, rop, n, 2*prec, prec);
}

int dgsl_rot_mp_call_identity(fmpz *rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state) {
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

int dgsl_rot_mp_call_inlattice(fmpz *rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state) {
  assert(rop); assert(self);

  const long n = fmpz_mat_ncols(self->B);
  mpz_t  tmp_g; mpz_init(tmp_g);
  fmpz_t tmp_f; fmpz_init(tmp_f);

  _fmpz_vec_zero(rop, n);

  fmpz* gen = _fmpz_vec_init(n);
  _fmpz_vec_set(gen, self->B->rows[0], n);
  for(long i=0; i<n; i++) {
    self->D[i]->call(tmp_g, self->D[i], state);
    fmpz_set_mpz(tmp_f, tmp_g);
    _fmpz_vec_scalar_addmul_fmpz(rop, gen, n, tmp_f);
    _fmpz_vec_rot_left_neg(gen, gen, n); //TODO: user should be able to config this
  }
  _fmpz_vec_clear(gen, n);
  _fmpz_vec_add(rop, rop, self->c_z, n);

  mpz_clear(tmp_g);
  fmpz_clear(tmp_f);

  return 0;
}

void dgsl_rot_mp_clear(dgsl_rot_mp_t *self) {
  if (!self)
    return;
  _fmpz_vec_clear(self->c_z, fmpz_mat_ncols(self->B));
  _mpfr_vec_clear(self->c,   fmpz_mat_ncols(self->B));
  
  if (!fmpz_mat_is_empty(self->B)) {
    fmpz_mat_clear(self->B);
  }

  if(self->call == dgsl_rot_mp_call_identity) {
    if(self->D) {
      dgs_disc_gauss_mp_clear(self->D[0]);
      free(self->D);
    }
  }
  mpfr_clear(self->sigma);
  free(self);
}
