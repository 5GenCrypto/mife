#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <oz/flint-addons.h>
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
  mpfr_mul_d(sigma_, sigma_, 0.398942280401433, MPFR_RNDN);

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
    _fmpz_vec_eucl_norm_mpfr(tmp, v, ncols, MPFR_RNDN);
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

  printf(" B+c:: m: %4ld, n: %4ld, bits: %3ld, log(σ): %6.2lf :: want: %8.2lf, have: %8.2lf, |log₂(want/have)|: %8.4f", nrows, ncols, bits, log2(sigma), left, rght, quality);

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
    _fmpz_vec_eucl_norm_mpfr(tmp, v, n, MPFR_RNDN);
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
  mpfr_mul_d(sigma_, sigma_, 0.398942280401433, MPFR_RNDN);

  dgsl_mp_t *D = dgsl_mp_init(B, sigma_, NULL, DGSL_INLATTICE);

  double *norms = dist_norms(D, ntrials, state);

  printf("  B:: m: %4ld, n: %4ld, bits: %3ld, log(σ): %6.2lf :: log(E[|v|]): %8.2lf, log(E[|v_0|]): %8.2f, log(E[|v_{n-1}|]): %8.2lf",
         nrows, ncols, bits, log2(sigma), log2(norms[0]), log2(norms[1]), log2(norms[ncols]));

  free(norms);
  dgsl_mp_clear(D);
  mpfr_clear(sigma_);
  fmpz_mat_clear(B);
  return 0;
}

int test_dist_identity(long nrows, long ncols, double sigma, size_t ntrials, flint_rand_t state) {
  assert(nrows>0);
  assert(ncols>0);
  assert(sigma>0);
  assert(ntrials>0);
  fmpz_mat_t B;
  fmpz_mat_init(B, nrows, ncols);
  fmpz_mat_one(B);

  mpfr_t sigma_;
  mpfr_init_set_d(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, 0.398942280401433, MPFR_RNDN);

  dgsl_mp_t *D = dgsl_mp_init(B, sigma_, NULL, DGSL_DETECT);

  double *norms = dist_norms(D, ntrials, state);

  printf("ZZ^n:: m: %4ld, n: %4ld, log(σ): %6.2lf :: log(E[|v|]): %6.2lf, log(E[|v_0|]): %6.2f, log(E[|v_{n-1}|]): %6.2lf",
         nrows, ncols, log2(sigma), log2(norms[0]), log2(norms[1]), log2(norms[ncols]));

  free(norms);
  dgsl_mp_clear(D);
  mpfr_clear(sigma_);
  fmpz_mat_clear(B);
  return 0;
}

double *dist_rot_norms(dgsl_rot_mp_t *D, size_t ntrials, const fmpz_poly_t c, flint_rand_t state) {

  const size_t n = D->n;

  double *norms = (double*)calloc(sizeof(double), n+1);
  fmpz_poly_t v;
  fmpz_poly_init(v);

  mpfr_t tmp;
  mpfr_init2(tmp, mpfr_get_prec(D->sigma));

  fmpz_t tmp_z;
  fmpz_init(tmp_z);

  for(size_t i=0; i<ntrials; i++) {
    if (fmpz_poly_is_zero(c))
      D->call(v, D, state->gmp_state);
    else if(fmpz_poly_is_one(c))
      dgsl_rot_mp_call_plus1(v, D, state->gmp_state);
    else
      dgsl_rot_mp_call_plus_fmpz_poly(v, D, c, state->gmp_state);

    fmpz_poly_eucl_norm_mpfr(tmp, v, MPFR_RNDN);
    norms[0] += mpfr_get_d(tmp, MPFR_RNDN);
    for (size_t j=0; j<n; j++) {
      fmpz_set(tmp_z, v->coeffs + j);
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
  fmpz_poly_clear(v);
  return norms;
}

double max_min_norm_ratio(double *norms, long n, double *min, double *max) {
  double min_norm = norms[1];
  double max_norm = norms[1];
  for (size_t j=1; j<=n; j++) {
    if (norms[j] < min_norm)
      min_norm = norms[j];
    if (norms[j] > max_norm)
      max_norm = norms[j];
  }
  *min = min_norm;
  *max = max_norm;
  double ratio = max_norm/min_norm;
  return ratio;
}

int test_dist_rot_identity(long ncols, double sigma, size_t ntrials, flint_rand_t state) {
  assert(ncols>0);
  assert(sigma>0);
  assert(ntrials>0);
  fmpz_poly_t B;
  fmpz_poly_init(B);
  fmpz_poly_one(B);

  mpfr_t sigma_;
  mpfr_init2(sigma_, 80);
  mpfr_set_d(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, 0.398942280401433, MPFR_RNDN);

  dgsl_rot_mp_t *D = dgsl_rot_mp_init(ncols, B, sigma_, NULL, DGSL_IDENTITY, 0);

  fmpz_poly_t c; fmpz_poly_init2(c, ncols);
  double *norms = dist_rot_norms(D, ntrials, c, state);
  fmpz_poly_clear(c);
  double min; double max;
  double ratio = max_min_norm_ratio(norms, ncols, &min, &max);

  printf("ZZ^n::          n: %4ld, log(σ): %6.2lf :: log(E[|v|]): %6.2lf, log(E[v_min]): %6.2f, log(E[|v_max|]): %6.2lf, log(E[v_max]/E[|v_min|]): %6.2lf",
         ncols, log2(sigma), log2(norms[0]), log2(min), log2(max), log2(ratio));

  free(norms);
  dgsl_rot_mp_clear(D);
  mpfr_clear(sigma_);
  fmpz_poly_clear(B);
  if (ratio < 1.2)
    return 0;
  else
    return 1;
}

int test_dist_rot_inlattice(long ncols, double sigma, double sigma_p, size_t ntrials, flint_rand_t state, int gpv) {
  assert(ncols>0);
  assert(sigma>0);
  assert(ntrials>0);

  fmpz_poly_t B;
  fmpz_poly_init(B);

  mpfr_t sigma_;
  mpfr_init2(sigma_, 320);
  mpfr_set_d(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, 0.398942280401433, MPFR_RNDN);

  fmpz_poly_sample_sigma(B, ncols, sigma_, state);

  mpfr_t sigma_p_;
  mpfr_init2(sigma_p_, 160);
  mpfr_init_set_d(sigma_p_, sigma_p, MPFR_RNDN);
  mpfr_mul_d(sigma_p_, sigma_p_, 0.398942280401433, MPFR_RNDN);

  dgsl_rot_mp_t *D;
  if (gpv)
    D = dgsl_rot_mp_init(ncols, B, sigma_p_, NULL, DGSL_GPV_INLATTICE, 0);
  else
    D = dgsl_rot_mp_init(ncols, B, sigma_p_, NULL, DGSL_INLATTICE, 0);

  fmpz_poly_t c; fmpz_poly_init2(c, ncols);
  double *norms = dist_rot_norms(D, ntrials, c, state);
  double min; double max;
  double ratio = max_min_norm_ratio(norms, ncols, &min, &max);

  printf("  <g>:: n: %4ld, log(σ): %6.2lf, log(σ'): %6.2lf", ncols, log2(sigma), log2(sigma_p));
  printf(" p:: log(E[|v|]): %6.2lf ?< σ'·sqrt(n): %6.2lf,",log2(norms[0]), log2(sigma_p*sqrt(ncols)));
  printf(" log(E[|v_min|]): %6.2lf, log(E[|v_max|]): %6.2lf, log(E[|v_max|/E[|v_min|]): %6.2lf, gpv: %d",log2(min), log2(max), log2(ratio), gpv);

  free(norms);
  dgsl_rot_mp_clear(D);
  mpfr_clear(sigma_);
  mpfr_clear(sigma_p_);
  fmpz_poly_clear(B);
  fmpz_poly_clear(c);
  if (ratio < 1.5)
    return 0;
  else
    return 1;
}

int test_dist_rot_plus1(long ncols, double sigma, double sigma_p, size_t ntrials, flint_rand_t state) {
  assert(ncols>0);
  assert(sigma>0);
  assert(ntrials>0);

  fmpz_poly_t B;
  fmpz_poly_init(B);

  mpfr_t sigma_;
  mpfr_init2(sigma_, 160);
  mpfr_set_d(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, 0.398942280401433, MPFR_RNDN);

  fmpz_poly_sample_sigma(B, ncols, sigma_, state);

  mpfr_t sigma_p_;
  mpfr_init2(sigma_p_, 160);
  mpfr_init_set_d(sigma_p_, sigma_p, MPFR_RNDN);
  mpfr_mul_d(sigma_p_, sigma_p_, 0.398942280401433, MPFR_RNDN);

  dgsl_rot_mp_t *D;
  D = dgsl_rot_mp_init(ncols, B, sigma_p_, NULL, DGSL_INLATTICE, 0);

  fmpz_poly_t c; fmpz_poly_init2(c, ncols);
  fmpz_poly_set_coeff_si(c, 0, 1);
  double *norms = dist_rot_norms(D, ntrials, c, state);
  double min; double max;
  double ratio = max_min_norm_ratio(norms, ncols, &min, &max);

  printf("<g>+1:: n: %4ld, log(σ): %6.2lf, log(σ'): %6.2lf", ncols, log2(sigma), log2(sigma_p));
  printf(" p:: log(E[|v|]): %6.2lf ?< σ'·sqrt(n): %6.2lf,",log2(norms[0]), log2(sigma_p*sqrt(ncols)));
  printf(" log(E[|v_min|]): %6.2lf, log(E[|v_max|]): %6.2lf, log(E[|v_max|/E[|v_min|]): %6.2lf",log2(min), log2(max), log2(ratio));

  free(norms);
  dgsl_rot_mp_clear(D);
  mpfr_clear(sigma_);
  mpfr_clear(sigma_p_);
  fmpz_poly_clear(B);
  fmpz_poly_clear(c);
  if (ratio < 1.4)
    return 0;
  else
    return 1;
}


int test_dist_rot_plus2(long ncols, double sigma, double sigma_p, size_t ntrials, flint_rand_t state) {
  assert(ncols>0);
  assert(sigma>0);
  assert(ntrials>0);

  fmpz_poly_t B;
  fmpz_poly_init(B);

  mpfr_t sigma_;
  mpfr_init2(sigma_, 160);
  mpfr_set_d(sigma_, sigma, MPFR_RNDN);
  mpfr_mul_d(sigma_, sigma_, 0.398942280401433, MPFR_RNDN);

  fmpz_poly_sample_sigma(B, ncols, sigma_, state);

  mpfr_t sigma_p_;
  mpfr_init2(sigma_p_, 160);
  mpfr_init_set_d(sigma_p_, sigma_p, MPFR_RNDN);
  mpfr_mul_d(sigma_p_, sigma_p_, 0.398942280401433, MPFR_RNDN);

  dgsl_rot_mp_t *D;
  D = dgsl_rot_mp_init(ncols, B, sigma_p_, NULL, DGSL_INLATTICE, 0);

  fmpz_poly_t c; fmpz_poly_init2(c, ncols);
  for(long i=0; i<ncols; i++)
    fmpz_poly_set_coeff_si(c, i, i);

  double *norms = dist_rot_norms(D, ntrials, c, state);
  double min; double max;
  double ratio = max_min_norm_ratio(norms, ncols, &min, &max);

  printf("<g>+c:: n: %4ld, log(σ): %6.2lf, log(σ'): %6.2lf", ncols, log2(sigma), log2(sigma_p));
  printf(" p:: log(E[|v|]): %6.2lf ?< σ'·sqrt(n): %6.2lf,",log2(norms[0]), log2(sigma_p*sqrt(ncols)));
  printf(" log(E[|v_min|]): %6.2lf, log(E[|v_max|]): %6.2lf, log(E[|v_max|/E[|v_min|]): %6.2lf",log2(min), log2(max), log2(ratio));

  free(norms);
  dgsl_rot_mp_clear(D);
  mpfr_clear(sigma_);
  mpfr_clear(sigma_p_);
  fmpz_poly_clear(B);
  fmpz_poly_clear(c);
  if (ratio < 1.4)
    return 0;
  else
    return 1;
}


int test_gso(long m, long n, long *M, double *GSO) {
  double quality = 0.0;
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

  status += test_dgsl_run( test_dist_identity(  16,  16,    100000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_identity(  32,  32,  10000000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_identity(  64,  64,  10000000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_identity( 128, 128, 100000000.0, 1<<10, randstate) );
  printf("\n");

  status += test_dgsl_run( test_dist_rot_identity(  16,    100000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_identity(  32,  10000000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_identity(  64,  10000000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_identity( 128, 100000000.0, 1<<10, randstate) );
  printf("\n");

  status += test_dgsl_run( test_dist_rot_inlattice(     16,  1000.0,    100000.0, 1<<10, randstate, 0) );
  status += test_dgsl_run( test_dist_rot_inlattice(     32, 10000.0,  10000000.0, 1<<10, randstate, 0) );
  status += test_dgsl_run( test_dist_rot_inlattice(     64,  1000.0,  10000000.0, 1<<10, randstate, 0) );
  status += test_dgsl_run( test_dist_rot_inlattice(    128,  1000.0, 100000000.0, 1<<10, randstate, 0) );
  printf("\n");

  status += test_dgsl_run( test_dist_rot_plus1(     16,  1000.0,    100000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_plus1(     32, 10000.0,  10000000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_plus1(     64,  1000.0,  10000000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_plus1(    128,  1000.0, 100000000.0, 1<<10, randstate) );
  printf("\n");

  status += test_dgsl_run( test_dist_rot_plus2(     16,  1000.0,    100000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_plus2(     32, 10000.0,  10000000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_plus2(     64,  1000.0,  10000000.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_rot_plus2(    128,  1000.0, 100000000.0, 1<<10, randstate) );
  printf("\n");

  status += test_dgsl_run( test_dist_rot_inlattice(     16,  1000.0,    100000.0, 1<<10, randstate, 1) );
  status += test_dgsl_run( test_dist_rot_inlattice(     32, 10000.0,  10000000.0, 1<<10, randstate, 1) );
  status += test_dgsl_run( test_dist_rot_inlattice(     64, 10000.0,  10000000.0, 1<<10, randstate, 1) );
  status += test_dgsl_run( test_dist_rot_inlattice(    128, 10000.0, 100000000.0, 1<<10, randstate, 1) );
  printf("\n");

  status += test_dgsl_run( test_dist_inlattice(30, 30,  2,  1849.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_inlattice(10, 10,  3,   255.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_inlattice(13, 13,  4, 24145.0, 1<<10, randstate) );
  status += test_dgsl_run( test_dist_inlattice(20, 20, 10, 24145.0, 1<<10, randstate) );
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

  flint_cleanup();
  return status;
}
