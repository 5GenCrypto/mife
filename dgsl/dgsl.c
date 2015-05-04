#include <stdio.h>
#include <mpfr.h>
#include <omp.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <oz/flint-addons.h>
#include <oz/oz.h>
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

  const long m = fmpz_mat_nrows(self->B);
  for(long j=0; j<m; j++) {
    long i = m-j-1;
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

/**
   Rotational Basis
**/

dgsl_rot_mp_t *dgsl_rot_mp_init(const long n, const fmpz_poly_t B, mpfr_t sigma, fmpq_poly_t c, const dgsl_alg_t algorithm, const oz_flag_t flags) {
  assert(mpfr_cmp_ui(sigma, 0) > 0);

  dgsl_rot_mp_t *self = (dgsl_rot_mp_t*)calloc(1, sizeof(dgsl_rot_mp_t));
  if(!self) dgs_die("out of memory");

  dgsl_alg_t alg = algorithm;

  self->n = n;

  self->prec = mpfr_get_prec(sigma);

  fmpz_poly_init(self->B);
  fmpz_poly_set(self->B, B);
  if(fmpz_poly_length(self->B) > n)
    dgs_die("polynomial is longer than length n");
  else
    fmpz_poly_realloc(self->B, n);


  fmpz_poly_init(self->c_z);
  fmpq_poly_init(self->c);

  mpfr_init2(self->sigma, self->prec);
  mpfr_set(self->sigma, sigma, MPFR_RNDN);

  if (alg == DGSL_DETECT) {
    if (fmpz_poly_is_one(self->B) && (c && fmpq_poly_is_zero(c))) {
      alg = DGSL_IDENTITY;
    } else if (c && fmpq_poly_is_zero(c))
      alg = DGSL_INLATTICE;
    else
      alg = DGSL_COSET; //TODO: we could test for lattice membership here
  }

  size_t tau = 3;
  if (2*ceil(sqrt(log2((double)n))) > tau)
    tau = 2*ceil(sqrt(log2((double)n)));

  switch(alg) {
  case DGSL_IDENTITY: {
    self->D = (dgs_disc_gauss_mp_t**)calloc(1, sizeof(dgs_disc_gauss_mp_t*));
    mpfr_t c_;
    mpfr_init2(c_, self->prec);
    mpfr_set_d(c_, 0.0, MPFR_RNDN);
    self->D[0] = dgs_disc_gauss_mp_init(self->sigma, c_, tau, DGS_DISC_GAUSS_DEFAULT);
    self->call = dgsl_rot_mp_call_identity;
    mpfr_clear(c_);
    break;
  }
  case DGSL_GPV_INLATTICE: {
    self->D = (dgs_disc_gauss_mp_t**)calloc(n, sizeof(dgs_disc_gauss_mp_t*));

    if (c && !fmpq_poly_is_zero(c)) {
      fmpq_t c_i;
      fmpq_init(c_i);
      for(int i=0; i<n; i++) {
        fmpq_poly_get_coeff_fmpq(c_i, c, i);
        fmpz_poly_set_coeff_fmpz(self->c_z, i, fmpq_numref(c_i));
      }
      fmpq_clear(c_i);
    }
    mpfr_mat_t G;
    mpfr_mat_init(G, n, n, self->prec);
    mpfr_mat_set_fmpz_poly(G, B);
    mpfr_mat_gso(G, MPFR_RNDN);

    mpfr_t sigma_;
    mpfr_init2(sigma_, self->prec);

    mpfr_t norm;
    mpfr_init2(norm, self->prec);

    mpfr_t c_;
    mpfr_init2(c_, self->prec);
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
    mpfr_clear(c_);
    mpfr_mat_clear(G);

    self->call = dgsl_rot_mp_call_gpv_inlattice;
    break;
  }
  case DGSL_INLATTICE: {
    fmpq_poly_init(self->sigma_sqrt);
    long r= 2*ceil(sqrt(log(n)));

    fmpq_poly_t Bq;    fmpq_poly_init(Bq);
    fmpq_poly_set_fmpz_poly(Bq, self->B);
    fmpq_poly_oz_invert_approx(self->B_inv, Bq, n, self->prec, flags);
    fmpq_poly_clear(Bq);

    _dgsl_rot_mp_sqrt_sigma_2(self->sigma_sqrt, self->B, sigma, r, n, self->prec, flags);

    mpfr_init2(self->r_f, self->prec);
    mpfr_set_ui(self->r_f, r, MPFR_RNDN);

    self->call = dgsl_rot_mp_call_inlattice;
    break;
  }
  case DGSL_COSET:
    dgs_die("not implemented");

  default:
    dgs_die("not implemented");
  }


  return self;
}


int dgsl_rot_mp_call_identity(fmpz_poly_t rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state) {
  assert(rop); assert(self);

  const long n = self->n;
  mpz_t  tmp_g;  mpz_init(tmp_g);

  fmpz_poly_zero(rop);
  fmpz_poly_realloc(rop, n);

  for(long i=0; i<n; i++) {
    self->D[0]->call(tmp_g, self->D[0], state);
    fmpz_poly_set_coeff_mpz(rop, i, tmp_g);
  }
  _fmpz_poly_set_length(rop, n);
  _fmpz_poly_normalise(rop);

  mpz_clear(tmp_g);

  return 0;
}

int dgsl_rot_mp_call_gpv_inlattice(fmpz_poly_t rop,  const dgsl_rot_mp_t *self, gmp_randstate_t state) {
  assert(rop); assert(self);

  const long n = self->n;
  mpz_t  tmp_z; mpz_init(tmp_z);
  fmpz_t tmp; fmpz_init(tmp);

  fmpz_poly_zero(rop);

  fmpz_poly_t tmp_poly;
  fmpz_poly_init(tmp_poly);
  fmpz_poly_set(tmp_poly, self->B);
  fmpz_poly_realloc(tmp_poly, n);
  tmp_poly->length = n;

  fmpz_poly_t tmp2;
  fmpz_poly_init(tmp2);

  for(long i=0; i<n; i++) {
    self->D[i]->call(tmp_z, self->D[i], state); fmpz_set_mpz(tmp, tmp_z);
    fmpz_poly_scalar_mul_fmpz(tmp2, tmp_poly, tmp);
    fmpz_poly_add(rop, rop, tmp2);
    _fmpz_vec_rot_left_neg(tmp_poly->coeffs, tmp_poly->coeffs, n);
  }
  fmpz_poly_clear(tmp_poly);
  fmpz_poly_add(rop, rop, self->c_z);

  fmpz_poly_clear(tmp2);
  mpz_clear(tmp_z);
  fmpz_clear(tmp);
  return 0;
}

int dgsl_rot_mp_call_plus1(fmpz_poly_t rop, const dgsl_rot_mp_t *self, gmp_randstate_t state) {
  const long n = self->n;
  fmpq_poly_t x;
  fmpq_poly_init(x);
  fmpq_poly_sample_D1(x, n, self->prec, state);
  fmpq_poly_oz_mul(x, self->sigma_sqrt, x, self->n);
  fmpq_poly_neg(x, x); // (-x2)

  // we sample with centre c = -1/g
  fmpq_poly_sub(x, x, self->B_inv); // (c+x2)

  fmpz_poly_disc_gauss_rounding(rop, x, self->r_f, state);
  fmpz_poly_neg(rop, rop);
  fmpz_poly_oz_mul(rop, self->B, rop, self->n);
  fmpz_add_ui(rop->coeffs, rop->coeffs, 1);

  fmpq_poly_clear(x);
  return 0;
}

int dgsl_rot_mp_call_plus_fmpz_poly(fmpz_poly_t rop, const dgsl_rot_mp_t *self, const fmpz_poly_t c, gmp_randstate_t state) {
  fmpz_poly_t t;  fmpz_poly_init(t);
  fmpz_poly_set(t, c);
  fmpq_poly_t tq; fmpq_poly_init(tq); // == 0
  fmpq_poly_set_fmpz_poly(tq, t);
  fmpq_poly_neg(tq, tq);
  dgsl_rot_mp_call_recenter_fmpq_poly(rop, self, tq, state);
  fmpz_poly_add(rop, rop, t);
  fmpq_poly_clear(tq);
  fmpz_poly_clear(t);
  return 0;
}

int dgsl_rot_mp_call_recenter_fmpq_poly(fmpz_poly_t rop, const dgsl_rot_mp_t *self, const fmpq_poly_t c, gmp_randstate_t state) {
  fmpq_poly_t x;
  fmpq_poly_init(x);
  fmpq_poly_sample_D1(x, self->n, self->prec, state);
  fmpq_poly_oz_mul(x, self->sigma_sqrt, x, self->n);
  fmpq_poly_neg(x, x);

  // we sample with centre c/g
  fmpq_poly_t t;  fmpq_poly_init(t);
  fmpq_poly_oz_mul(t, c, self->B_inv, self->n);
  fmpq_poly_add(x, t, x); // (c+x2)
  fmpq_poly_clear(t);

  fmpz_poly_disc_gauss_rounding(rop, x, self->r_f, state);
  fmpz_poly_oz_mul(rop, self->B, rop, self->n);

  fmpq_poly_clear(x);
  return 0;

}

int _dgsl_rot_mp_call_inlattice_multiplier(fmpz_poly_t rop, const dgsl_rot_mp_t *self, gmp_randstate_t state) {
  const long n = self->n;
  fmpq_poly_t x;
  fmpq_poly_init(x);
  fmpq_poly_sample_D1(x, n, self->prec, state);
  fmpq_poly_oz_mul(x, self->sigma_sqrt, x, self->n);
  fmpq_poly_neg(x, x);;
  fmpz_poly_disc_gauss_rounding(rop, x, self->r_f, state);
  fmpz_poly_neg(rop, rop);;
  fmpq_poly_clear(x);
  return 0;
}

int dgsl_rot_mp_call_inlattice(fmpz_poly_t rop, const dgsl_rot_mp_t *self, gmp_randstate_t state) {
  _dgsl_rot_mp_call_inlattice_multiplier(rop, self, state);
  fmpz_poly_oz_mul(rop, self->B, rop, self->n);
  return 0;
}



void dgsl_rot_mp_clear(dgsl_rot_mp_t *self) {
  if (!self)
    return;
  fmpz_poly_clear(self->c_z);
  fmpq_poly_clear(self->c);

  if (!fmpz_poly_is_zero(self->B))
     fmpz_poly_clear(self->B);
  if (!fmpq_poly_is_zero(self->B_inv))
     fmpq_poly_clear(self->B_inv);

  if(self->call == dgsl_rot_mp_call_identity) {
    if(self->D) {
      dgs_disc_gauss_mp_clear(self->D[0]);
      free(self->D);
    }
  }
  if(self->call == dgsl_rot_mp_call_gpv_inlattice) {

  }
  if(self->call == dgsl_rot_mp_call_inlattice) {
    mpfr_clear(self->r_f);
    fmpq_poly_clear(self->sigma_sqrt);
  }
  mpfr_clear(self->sigma);
  free(self);
}

void fmpz_poly_disc_gauss_rounding(fmpz_poly_t rop, const fmpq_poly_t x, const mpfr_t r_f, gmp_randstate_t randstate) {
  mpfr_t xi;  mpfr_init2(xi, mpfr_get_prec(r_f));
  mpf_t xi_f; mpf_init2(xi_f, mpfr_get_prec(r_f));
  mpq_t xi_q; mpq_init(xi_q);
  mpz_t s_z;  mpz_init(s_z);

  const long n = fmpq_poly_length(x);
  const size_t tau = (ceil(2*sqrt(log2((double)n))) > 3) ? ceil(2*sqrt(log2((double)n))) : 3;

  fmpz_poly_zero(rop);

  for(int i=0; i<n; i++) {
    fmpq_poly_get_coeff_mpq(xi_q, x, i);
    mpf_set_q(xi_f, xi_q);
    mpfr_set_f(xi, xi_f, MPFR_RNDN);

    dgs_disc_gauss_mp_t *D = dgs_disc_gauss_mp_init(r_f, xi, tau, DGS_DISC_GAUSS_UNIFORM_ONLINE);
    D->call(s_z, D, randstate);
    dgs_disc_gauss_mp_clear(D);

    fmpz_poly_set_coeff_mpz(rop, i, s_z);
  }
  mpz_clear(s_z);
  mpq_clear(xi_q);
  mpf_clear(xi_f);
  mpfr_clear(xi);
}

void fmpq_poly_sample_D1(fmpq_poly_t f, int n, mpfr_prec_t prec, gmp_randstate_t state) {
  mpfr_t u1; mpfr_init2(u1, prec);
  mpfr_t u2; mpfr_init2(u2, prec);
  mpfr_t z1; mpfr_init2(z1, prec);
  mpfr_t z2; mpfr_init2(z2, prec);

  mpfr_t pi2; mpfr_init2(pi2, prec);
  mpfr_const_pi(pi2, MPFR_RNDN);
  mpfr_mul_si(pi2, pi2, 2, MPFR_RNDN);

  mpf_t tmp_f;
  mpq_t tmp_q;
  mpf_init(tmp_f);
  mpq_init(tmp_q);

  assert(n%2==0);

  for(long i=0; i<n; i+=2) {
    mpfr_urandomb(u1, state);
    mpfr_urandomb(u2, state);
    mpfr_log(u1, u1, MPFR_RNDN);
    mpfr_mul_si(u1, u1, -2, MPFR_RNDN);
    mpfr_sqrt(u1, u1, MPFR_RNDN);
    mpfr_mul(u2, pi2, u2, MPFR_RNDN);
    mpfr_cos(z1, u2, MPFR_RNDN);
    mpfr_mul(z1, z1, u1, MPFR_RNDN); //z1 = sqrt(-2*log(u1)) * cos(2*pi*u2)
    mpfr_sin(z2, u2, MPFR_RNDN);
    mpfr_mul(z2, z2, u1, MPFR_RNDN); //z1 = sqrt(-2*log(u1)) * sin(2*pi*U2)

    mpfr_get_f(tmp_f, z1, MPFR_RNDN);
    mpq_set_f(tmp_q, tmp_f);
    fmpq_poly_set_coeff_mpq(f, i, tmp_q);

    mpfr_get_f(tmp_f, z2, MPFR_RNDN);
    mpq_set_f(tmp_q, tmp_f);
    fmpq_poly_set_coeff_mpq(f, i+1, tmp_q);
  }
  mpf_clear(tmp_f);
  mpq_clear(tmp_q);
  mpfr_clear(pi2);
  mpfr_clear(u1);
  mpfr_clear(u2);
  mpfr_clear(z1);
  mpfr_clear(z2);
}

/**
   sqrt(Σ_2) with Σ_2 = Σ - Σ_1 = σ^2·g^-T·g^-1 - r^2·I
*/

void _dgsl_rot_mp_sqrt_sigma_2(fmpq_poly_t rop, const fmpz_poly_t g, const mpfr_t sigma,
                              const int r, const long n, const mpfr_prec_t prec, const oz_flag_t flags) {
  fmpq_poly_zero(rop);

  fmpq_t r_q2;
  fmpq_init(r_q2);
  fmpq_set_si(r_q2, r, 1);
  fmpq_mul(r_q2, r_q2, r_q2);
  fmpq_neg(r_q2, r_q2);
  fmpq_poly_set_coeff_fmpq(rop, 0, r_q2);
  fmpq_clear(r_q2);

  fmpq_poly_t g_q; fmpq_poly_init(g_q);
  fmpq_poly_set_fmpz_poly(g_q, g);

  fmpq_poly_t ng; fmpq_poly_init(ng);
  fmpq_poly_oz_invert_approx(ng, g_q, n, prec, flags);

  fmpq_poly_t ngt;
  fmpq_poly_init(ngt);
  fmpq_poly_oz_conjugate(ngt, ng, n);

  fmpq_poly_t nggt;
  fmpq_poly_init(nggt);
  fmpq_poly_oz_mul(nggt, ng, ngt, n);

  /**
     We compute sqrt(g^-T · g^-1) to use it as the starting point for
     convergence on sqrt(σ^2 · g^-T · g^-1 - r^2) below. We can compute the
     former with less precision than the latter
  */

  mpfr_t norm;
  mpfr_init2(norm, prec);
  fmpz_poly_2norm_mpfr(norm, g, MPFR_RNDN);
  double p = mpfr_get_d(norm, MPFR_RNDN);

  /**
     |g^-1| ~= 1/|g|
     |g^-T| ~= |g^-1|
     |g^-1·g^-T| ~= sqrt(n)·|g^-T|·|g^-1|
  */
  fmpq_poly_t sqrt_start; fmpq_poly_init(sqrt_start);
  p = log2(n) + 4*log2(p);
  int fail = -1;
  while (fail) {
    p = 2*p;
    if (fail<0)
      fail = fmpq_poly_oz_sqrt_approx_db(sqrt_start, nggt, n, p, prec/2, flags, NULL);
    else
      fail = fmpq_poly_oz_sqrt_approx_db(sqrt_start, nggt, n, p, prec/2, flags, sqrt_start);
    if(fail)
      fprintf(stderr, "FAILED for precision %7.1f with code (%d), doubling precision.\n", p, fail);
  }

  fmpq_t sigma2;
  fmpq_init(sigma2);
  fmpq_set_mpfr(sigma2, sigma, MPFR_RNDN);
  fmpq_poly_scalar_mul_fmpq(sqrt_start, sqrt_start, sigma2);

  fmpq_mul(sigma2, sigma2, sigma2);
  fmpq_poly_scalar_mul_fmpq(nggt, nggt, sigma2);
  fmpq_clear(sigma2);

  fmpq_poly_add(rop, rop, nggt);

  p = p + 2*log2(mpfr_get_d(sigma, MPFR_RNDN));

  fmpq_poly_oz_sqrt_approx_babylonian(rop, rop, n, p, prec, flags, sqrt_start);

  mpfr_clear(norm);
  fmpq_poly_clear(g_q);
  fmpq_poly_clear(ng);
  fmpq_poly_clear(ngt);
  fmpq_poly_clear(nggt);
  fmpq_poly_clear(sqrt_start);
}
