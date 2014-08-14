#include <flint/fmpz_mat.h>
#include <flint-addons/flint-addons.h>
#include <gpv/gpv.h>
#include <gpv/gso.h>
#include <math.h>

int test_dist_inlattice(long nrows, long ncols, mp_bitcnt_t bits, double sigma, size_t ntrials, flint_rand_t state) {
  assert(nrows>0);
  assert(ncols>0);
  assert(sigma>0);
  assert(ntrials>0);
  fmpz_mat_t B;
  fmpz_mat_init(B, nrows, ncols);
  fmpz_mat_randbits(B, state, bits);

  mpfr_t sigma_;
  mpfr_init_set_d(sigma_, sigma, MPFR_RNDN);

  dgs_disc_gauss_lattice_mp_t *D = dgs_disc_gauss_lattice_mp_init(B, sigma_, NULL, DGS_LATTICE_INLATTICE);

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

  double left = sqrt((double)nrows)*mpfr_get_d(sigma_, MPFR_RNDN);
  double rght = mpfr_get_d(tmp, MPFR_RNDN);
  double quality = fabs(log2(left/rght));

  printf("inlattice:: m: %4ld, n: %4ld, bits: %3ld, σ: %10.4lf :: want: %8.2lf, have: %8.2lf, |log₂(want/have)|: %8.4f\n", nrows, ncols, bits, sigma, left, rght, quality);

  dgs_disc_gauss_lattice_mp_clear(D);
  _fmpz_vec_clear(v, ncols);

  mpfr_clear(sigma_);
  mpfr_clear(tmp);
  mpfr_clear(acc);
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

  dgs_disc_gauss_lattice_mp_t *D = dgs_disc_gauss_lattice_mp_init(B, sigma_, NULL, DGS_LATTICE_DETECT);

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

  double left = sqrt((double)nrows)*mpfr_get_d(sigma_, MPFR_RNDN);
  double rght = mpfr_get_d(tmp, MPFR_RNDN);
  double quality = fabs(log2(left/rght));

  printf("   simple:: m: %4ld, n: %4ld, σ: %10.4lf :: want: %8.2lf, have: %8.2lf, |log₂(want/have)|: %8.4f\n", nrows, ncols, sigma, left, rght, quality);

  dgs_disc_gauss_lattice_mp_clear(D);
  _fmpz_vec_clear(v, ncols);

  mpfr_clear(sigma_);
  mpfr_clear(tmp);
  mpfr_clear(acc);
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
  mpfr_mat_init(G, m, n, 53);
  mpfr_mat_set_fmpz_mat(G, B);

  mpfr_mat_gso(G);

  x = 0;
  for(long i=0; i<m; i++) {
    for(long j=0; j<n; j++) {
      quality += fabs(mpfr_get_d(G->rows[i][j], MPFR_RNDN) -  GSO[x]);
      x++;
    }
  }
  printf("      gso:: m: %4ld, n: %4ld, dist: %8.4f\n", m, n, quality);

  mpfr_mat_clear(G);
  fmpz_mat_clear(B);
  return (int)quality;
}

int main(int argc, char *argv[]) {

  flint_rand_t state;
  flint_randinit(state);
  _flint_rand_init_gmp(state);

  test_dist_inlattice(30, 30,  2,   255.0, 1<<13, state);
  test_dist_inlattice(10, 10,  3,   255.0, 1<<13, state);
  test_dist_inlattice( 3, 13,  4, 24145.0, 1<<13, state);
  printf("\n");
  
  test_dist_simple(  10,   10,  3.0, 1<<14, state);
  test_dist_simple( 100,  100,  3.0, 1<<14, state);
  test_dist_simple( 100,  100, 25.2, 1<<12, state);
  test_dist_simple( 512,  512, 30.0, 1<<10, state);

  test_dist_simple( 1,   10,  3.0, 1<<14, state);
  test_dist_simple( 1,  100,  3.0, 1<<14, state);
  test_dist_simple( 1,  100, 25.2, 1<<12, state);
  test_dist_simple( 1,  512, 30.0, 1<<10, state);

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
    test_gso(4,4, M, G);
  }

  {
    long   M[9] = {1,2,2, -1,0,2, 0,0,1};
    double G[9] = { 1.00000000000000,  2.000000000000000, 2.00000000000000,
                   -1.33333333333333, -0.666666666666667, 1.33333333333333,
                    0.22222222222222, -0.222222222222222, 0.11111111111111};
    test_gso(3,3, M, G);
  }



  return 0;
}
