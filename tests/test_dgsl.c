#include <flint/fmpz_mat.h>
#include <flint-addons/flint-addons.h>
#include <dgsl/dgsl.h>
#include <dgsl/gso.h>
#include <math.h>

int test_dist_coset(long nrows, long ncols, mp_bitcnt_t bits, double sigma, long c, size_t ntrials, flint_rand_t state) {
  assert(nrows>0);
  assert(ncols>0);
  assert(bits>0);
  assert(sigma>0);
  assert(ntrials>0);
  fmpz_mat_t B;
  fmpz_mat_init(B, nrows, ncols);
  fmpz_mat_randtest(B, state, bits);
  
  mpfr_t sigma_;
  mpfr_init_set_d(sigma_, sigma, MPFR_RNDN);

  mpfr_t *c_ = _mpfr_vec_init(ncols, 53);
  mpfr_set_si(c_[0], c, MPFR_RNDN);
  
  dgsl_mp_t *D = dgsl_mp_init(B, sigma_, c_, DGSL_INLATTICE);

  _mpfr_vec_clear(c_, ncols);

  fmpz *v = _fmpz_vec_init(ncols);

  mpfr_t tmp, acc;
  mpfr_init(tmp);
  mpfr_init_set_ui(acc, 0, MPFR_RNDN);
  
  for(size_t i=0; i<ntrials; i++) {
    D->call(v, D, state->gmp_state);
    _fmpz_vec_2norm_mpfr(tmp, v, ncols);
    mpfr_add(acc, acc, tmp, MPFR_RNDN);
  }

  mpfr_div_ui(tmp, acc, ntrials, MPFR_RNDN);

  double left;
  if(nrows == 1)
    left = sqrt((double)ncols)*mpfr_get_d(sigma_, MPFR_RNDN);
  else
    left = sqrt((double)nrows)*mpfr_get_d(sigma_, MPFR_RNDN);
  double rght = mpfr_get_d(tmp, MPFR_RNDN);
  double quality = fabs(log2(left/rght));

  printf("     coset:: m: %4ld, n: %4ld, bits: %3ld, σ: %10.4lf :: want: %8.2lf, have: %8.2lf, |log₂(want/have)|: %8.4f", nrows, ncols, bits, sigma, left, rght, quality);

  dgsl_mp_clear(D);
  _fmpz_vec_clear(v, ncols);

  mpfr_clear(sigma_);
  mpfr_clear(tmp);
  mpfr_clear(acc);
  if (quality < 0.1)
    return 0;
  else
    return 1;
}

double *dist_norms(dgsl_mp_t *D, size_t ntrials, flint_rand_t state) {

  const size_t n = fmpz_mat_ncols(D->B);

  double *norms = (double*)calloc(sizeof(double), n+1);
  fmpz *v = _fmpz_vec_init(n);

  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(D->sigma));
  
  fmpz_t tmp_z;
  fmpz_init(tmp_z);

  for(size_t i=0; i<ntrials; i++) {
    D->call(v, D, state->gmp_state);
    _fmpz_vec_2norm_mpfr(tmp, v, n);
    norms[0] += mpfr_get_d(tmp, MPFR_RNDN);
    for (size_t j=0; j<n; j++) {
      fmpz_set(tmp_z, v + j);
      fmpz_pow_ui(tmp_z, tmp_z, 2);
      norms[j+1] += fmpz_get_d(tmp_z);
    }
  }

  norms[0] /= ntrials;
  for (size_t j=1; j<=n; j++) {
    norms[j] /= ntrials;
    norms[j] = sqrt(norms[j]);
  }
  
  mpfr_clear(tmp);
  fmpz_clear(tmp_z);
  _fmpz_vec_clear(v, n);
  return norms;
}

int test_dist_inlattice(long nrows, long ncols, mp_bitcnt_t bits, double sigma, size_t ntrials, flint_rand_t state) {
  assert(nrows>0);
  assert(ncols>0);
  assert(bits>0);
  assert(sigma>0);
  assert(ntrials>0);

  fmpz_mat_t B;
  fmpz_mat_init(B, nrows, ncols);
  fmpz_mat_randrank(B, state, nrows, bits);
  fmpz_mat_randops(B, state, 100);

  mpfr_t sigma_;
  mpfr_init_set_d(sigma_, sigma, MPFR_RNDN);

  dgsl_mp_t *D = dgsl_mp_init(B, sigma_, NULL, DGSL_INLATTICE);

  double *norms = dist_norms(D, ntrials, state);

  printf(" inlattice:: m: %4ld, n: %4ld, bits: %3ld, σ: %10.4lf :: log(E[||v||]): %8.2lf, log(E[||v_0||]): %8.2f, log(E[||v_{n-1}||]): %8.2lf",
         nrows, ncols, bits, sigma, log2(norms[0]), log2(norms[1]), log2(norms[ncols]));

  free(norms);
  dgsl_mp_clear(D);
  mpfr_clear(sigma_);
  fmpz_mat_clear(B);
  return 0;
}

int test_dist_simple(long nrows, long ncols, double sigma, size_t ntrials, flint_rand_t state) {
  assert(nrows>0);
  assert(ncols>0);
  assert(sigma>0);
  assert(ntrials>0);
  fmpz_mat_t B;
  fmpz_mat_init(B, nrows, ncols);
  fmpz_mat_one(B);

  mpfr_t sigma_;
  mpfr_init_set_d(sigma_, sigma, MPFR_RNDN);

  dgsl_mp_t *D = dgsl_mp_init(B, sigma_, NULL, DGSL_DETECT);

  double *norms = dist_norms(D, ntrials, state);

  printf("    simple:: m: %4ld, n: %4ld, σ: %10.4lf :: log(E[||v||]): %8.2lf, log(E[||v_0||]): %8.2f, log(E[||v_{n-1}||]): %8.2lf",
         nrows, ncols, sigma, log2(norms[0]), log2(norms[1]), log2(norms[ncols]));

  free(norms);
  dgsl_mp_clear(D);
  mpfr_clear(sigma_);
  fmpz_mat_clear(B);
  return 0;
}

double *dist_rot_norms(dgsl_rot_mp_t *D, size_t ntrials, flint_rand_t state) {

  const size_t n = fmpz_mat_ncols(D->B);

  double *norms = (double*)calloc(sizeof(double), n+1);
  fmpz *v = _fmpz_vec_init(n);

  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(D->sigma));
  
  fmpz_t tmp_z;
  fmpz_init(tmp_z);

  for(size_t i=0; i<ntrials; i++) {
    D->call(v, D, state->gmp_state);
    _fmpz_vec_2norm_mpfr(tmp, v, n);
    norms[0] += mpfr_get_d(tmp, MPFR_RNDN);
    for (size_t j=0; j<n; j++) {
      fmpz_set(tmp_z, v + j);
      fmpz_pow_ui(tmp_z, tmp_z, 2);
      norms[j+1] += fmpz_get_d(tmp_z);
    }
  }

  norms[0] /= ntrials;
  for (size_t j=1; j<=n; j++) {
    norms[j] /= ntrials;
    norms[j] = sqrt(norms[j]);
  }
  
  mpfr_clear(tmp);
  fmpz_clear(tmp_z);
  _fmpz_vec_clear(v, n);
  return norms;
}


int test_dist_rot_simple(long ncols, double sigma, size_t ntrials, flint_rand_t state) {
  assert(ncols>0);
  assert(sigma>0);
  assert(ntrials>0);
  fmpz_mat_t B;
  fmpz_mat_init(B, 1, ncols);
  fmpz_mat_one(B);

  mpfr_t sigma_;
  mpfr_init_set_d(sigma_, sigma, MPFR_RNDN);

  dgsl_rot_mp_t *D = dgsl_rot_mp_init(B, sigma_, NULL, DGSL_DETECT);

  double *norms = dist_rot_norms(D, ntrials, state);

  printf("simple rot:: n: %4ld, σ: %10.4lf :: log(E[||v||]): %8.2lf, log(E[||v_0||]): %8.2f, log(E[||v_{n-1}||]): %8.2lf",
         ncols, sigma, log2(norms[0]), log2(norms[1]), log2(norms[ncols]));

  free(norms);
  dgsl_rot_mp_clear(D);
  mpfr_clear(sigma_);
  fmpz_mat_clear(B);
  return 0;
}

int test_dist_rot_inlattice(long ncols, mp_bitcnt_t bits, double sigma, size_t ntrials, flint_rand_t state) {
  assert(ncols>0);
  assert(bits>0);
  assert(sigma>0);
  assert(ntrials>0);

  fmpz_mat_t B;
  fmpz_mat_init(B, 1, ncols);
  fmpz_mat_randrank(B, state, 1, bits);
  fmpz_mat_randops(B, state, 100);

  mpfr_t sigma_;
  mpfr_init_set_d(sigma_, sigma, MPFR_RNDN);

  dgsl_rot_mp_t *D = dgsl_rot_mp_init(B, sigma_, NULL, DGSL_INLATTICE);

  double *norms = dist_rot_norms(D, ntrials, state);

  printf(" inlattice:: n: %4ld, bits: %3ld, σ: %10.4lf :: log(E[||v||]): %8.2lf, log(E[||v_0||]): %8.2f, log(E[||v_{n-1}||]): %8.2lf",
         ncols, bits, sigma, log2(norms[0]), log2(norms[1]), log2(norms[ncols]));

  free(norms);
  dgsl_rot_mp_clear(D);
  mpfr_clear(sigma_);
  fmpz_mat_clear(B);
  return 0;
}

int test_gso(long m, long n, long *M, double *GSO) {
  double quality;
  fmpz_mat_t B;
  fmpz_mat_init(B, m, n);

  long x = 0;
  for(long i=0; i<m; i++) {
    for(long j=0; j<n; j++) {
      fmpz_set_si(B->rows[i] + j, M[x]);
      x++;
    }
  }

  mpfr_mat_t G;
  if (m != 1) {
    mpfr_mat_init(G, m, n, 53);
    mpfr_mat_set_fmpz_mat(G, B);
  } else {
    mpfr_mat_init(G, n, n, 53);
    mpfr_mat_set_fmpz_mat_rot(G, B);
  }

  mpfr_mat_gso(G, MPFR_RNDN);

  x = 0;
  for(long i=0; i<G->r; i++) {
    for(long j=0; j<G->c; j++) {
      quality += fabs(mpfr_get_d(G->rows[i][j], MPFR_RNDN) -  GSO[x]);
      x++;
    }
  }
  printf("      gso:: m: %4ld, n: %4ld, dist: %8.4f", m, n, quality);

  mpfr_mat_clear(G);
  fmpz_mat_clear(B);
  return (int)quality;
}

int test_dgsl_run(int status) {
  if (status)
    printf(" FAIL\n");
  else
    printf(" PASS\n");
  return status;
}

int main(int argc, char *argv[]) {
  int status = 0;

  flint_rand_t randstate;
  flint_randinit_seed(randstate, 0x1337, 1);
  
  status += test_dgsl_run( test_dist_inlattice(30, 30,  2,  1849.0, 1<<14, randstate) );
  status += test_dgsl_run( test_dist_inlattice(10, 10,  3,   255.0, 1<<14, randstate) );
  status += test_dgsl_run( test_dist_inlattice(13, 13,  4, 24145.0, 1<<14, randstate) );
  status += test_dgsl_run( test_dist_inlattice(20, 20, 10, 24145.0, 1<<14, randstate) );
  printf("\n");

  status += test_dgsl_run( test_dist_rot_inlattice( 30,  2,   255.0, 1<<13, randstate) );
  status += test_dgsl_run( test_dist_rot_inlattice( 10,  3,    25.0, 1<<13, randstate) );
  status += test_dgsl_run( test_dist_rot_inlattice( 13,  4,  2445.0, 1<<13, randstate) );
  printf("\n");

  status += test_dgsl_run( test_dist_simple(  10,   10,  3.0, 1<<14, randstate) );
  status += test_dgsl_run( test_dist_simple( 100,  100,  3.0, 1<<14, randstate) );
  status += test_dgsl_run( test_dist_simple( 100,  100, 25.2, 1<<12, randstate) );
  status += test_dgsl_run( test_dist_simple( 512,  512, 30.0, 1<<10, randstate) );
  printf("\n");

  status += test_dgsl_run( test_dist_rot_simple(  10,  3.0, 1<<14, randstate) );
  status += test_dgsl_run( test_dist_rot_simple( 100,  3.0, 1<<14, randstate) );
  status += test_dgsl_run( test_dist_rot_simple( 100, 25.2, 1<<12, randstate) );
  status += test_dgsl_run( test_dist_rot_simple( 512, 30.0, 1<<10, randstate) );
  printf("\n");

  {
    long   M[16] = {-1, 0,-1,-1,
                     5,-7,-1,-8,
                    -1,-2, 5, 2,
                     0, 2,-1,-5};
    double G[16] = {-1.00000000000000,  0.00000000000000,  -1.00000000000000, -1.00000000000000,
                     6.33333333333333, -7.00000000000000,   0.33333333333333, -6.66666666666667,
                    -2.81047381546135, -2.20947630922693,   3.00997506234414, -0.199501246882793,
                     0.27364941873717,  2.34556644631867,   1.83736038294962, -2.11100980168680};
    status += test_dgsl_run( test_gso(4,4, M, G) );
  }

  {
    long   M[9] = {1,2,2, -1,0,2, 0,0,1};
    double G[9] = { 1.00000000000000,  2.000000000000000, 2.00000000000000,
                   -1.33333333333333, -0.666666666666667, 1.33333333333333,
                    0.22222222222222, -0.222222222222222, 0.11111111111111};
    status += test_dgsl_run( test_gso(3,3, M, G) );
  }

  {
    long   M[4] = {-11, 6, 1, 1};
    double G[16] = {-11.0000000000000,  6.00000000000000,  1.00000000000000,  1.00000000000000,
                    -4.32075471698113, -9.18867924528302,  6.30188679245283,  1.30188679245283,
                    -2.43517430473952, -4.05209557383470, -8.90677634155895,  6.43243243243243,
                    -2.98113207547170, -3.11320754716981, -4.62264150943396, -9.49056603773585};

    status += test_dgsl_run( test_gso(1,4, M, G) );
  }


  return status;
}
