#include <flint/fmpz_mat.h>
#include <flint-addons/flint-addons.h>
#include <gpv/gpv.h>
#include <math.h>

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

  double left = sqrt((double)ncols)*mpfr_get_d(sigma_, MPFR_RNDN);
  double rght = mpfr_get_d(tmp, MPFR_RNDN);
  double quality = fabs(log2(left/rght));

  printf("m: %4ld, n: %4ld, σ: %8.4lf :: want: %8.2lf, have: %8.2lf, |log₂(want/have)|: %8.4f\n", nrows, ncols, sigma, left, rght, quality);

  dgs_disc_gauss_lattice_mp_clear(D);
  _fmpz_vec_clear(v, ncols);

  mpfr_clear(sigma_);
  mpfr_clear(tmp);
  mpfr_clear(acc);
  return 0;
}



int main(int argc, char *argv[]) {

  flint_rand_t state;
  flint_randinit(state);
  _flint_rand_init_gmp(state);
  test_dist_simple(  10,   10,  3.0, 1<<14, state);
  test_dist_simple( 100,  100,  3.0, 1<<14, state);
  test_dist_simple( 100,  100, 25.2, 1<<12, state);
  test_dist_simple( 512,  512, 30.0, 1<<10, state);

  test_dist_simple( 1,   10,  3.0, 1<<14, state);
  test_dist_simple( 1,  100,  3.0, 1<<14, state);
  test_dist_simple( 1,  100, 25.2, 1<<12, state);
  test_dist_simple( 1,  512, 30.0, 1<<10, state);


  return 0;
}
