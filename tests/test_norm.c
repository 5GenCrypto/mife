#include <oz/oz.h>
#include <oz/util.h>
#include <math.h>

int test_nmod_poly_oz_ideal_norm(slong n, flint_rand_t state) {
  nmod_poly_t f;
  nmod_poly_t g;

  mp_limb_t q = _n_next_oz_good_probaprime(2*n+1, 2*n);

  nmod_poly_init2(f, q, n);
  nmod_poly_init(g, q);
  nmod_poly_set_coeff_ui(g, 0, 1);
  nmod_poly_set_coeff_ui(g, n, 1);

  nmod_poly_randtest(f, state, n);
  while(!nmod_poly_get_coeff_ui(f, n-1))
      nmod_poly_randtest(f, state, n);

  uint64_t t0 = oz_walltime(0);
  mp_limb_t r0 = nmod_poly_resultant(f, g);
  t0 = oz_walltime(t0);

  uint64_t t1 = oz_walltime(0);
  mp_limb_t r1 = nmod_poly_oz_resultant(f, n);
  t1 = oz_walltime(t1);

  int r = (r0 == r1);

  printf("n: %4ld, q: %4ld, flint: %7.2fs, fft: %7.2fs, flint/fft: %8.2f ", n, q,
         t0/1000000.0, t1/1000000.0, (double)t0/(double)t1);
  if (r)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  nmod_poly_clear(f);
  nmod_poly_clear(g);
  return !r;
}

int test_fmpz_poly_oz_ideal_norm(slong n, mp_bitcnt_t bits, flint_rand_t state) {
  fmpz_poly_t f;
  fmpz_poly_t g;
  fmpz_poly_init(f);
  fmpz_poly_init_oz_modulus(g, n);

  fmpz_poly_randtest(f, state, n, bits);
  while(!fmpz_poly_get_coeff_ptr(f, n-1))
      fmpz_poly_randtest(f, state, n, bits);

  fmpz_t r0, r1, r2;
  fmpz_init(r0);
  fmpz_init(r1);
  fmpz_init(r2);

  uint64_t t0 = oz_walltime(0);
  fmpz_poly_resultant_modular(r0, f, g);
  t0 = oz_walltime(t0);

  uint64_t t1 = oz_walltime(0);
  fmpz_poly_oz_ideal_norm(r1, f, n, 0);
  t1 = oz_walltime(t1);

  uint64_t t2 = oz_walltime(0);
  fmpz_poly_oz_ideal_norm(r2, f, n, 40);
  t2 = oz_walltime(t2);

  int r = fmpz_equal(r0,r1);
  if (r!=1) {
    fmpz_print(r0); printf("\n");
    fmpz_print(r1); printf("\n");
    exit(0);
  }

  printf("n: %4ld, bits: %4ld, flint: %7.2fs, oz: %7.2fs, approx: %8.2fs, flint/bounded: %8.2f, oz/approx: %8.2f ", n, bits,
         t0/1000000.0, t1/1000000.0, t2/1000000.0, (double)t0/(double)t1, (double)t0/(double)t2);
  if (r)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  fmpz_poly_clear(f);
  fmpz_poly_clear(g);
  fmpz_clear(r0);
  fmpz_clear(r1);
  fmpz_clear(r2);
  return !r;
}

int test_fmpq_poly_oz_ideal_norm(slong n, mp_bitcnt_t bits, flint_rand_t state) {
  fmpq_poly_t f;
  fmpq_poly_t g;
  fmpq_poly_init(f);
  fmpq_poly_init_oz_modulus(g, n);

  fmpq_poly_randtest(f, state, n, bits);
  while(fmpq_poly_degree(f) != n-1)
      fmpq_poly_randtest(f, state, n, bits);

  fmpq_t r0, r1, r2;
  fmpq_init(r0);
  fmpq_init(r1);
  fmpq_init(r2);

  uint64_t t0 = oz_walltime(0);
  fmpq_poly_oz_ideal_norm(r0, f, n, 0);
  t0 = oz_walltime(t0);

  uint64_t t1 = oz_walltime(0);
  fmpq_poly_oz_ideal_norm(r1, f, n, -n);
  t1 = oz_walltime(t1);

  uint64_t t2 = oz_walltime(0);
  fmpq_poly_oz_ideal_norm(r2, f, n, n);
  t2 = oz_walltime(t2);


  fmpq_div(r0, r0, r2);
  mpq_t tmp;
  mpq_init(tmp);
  fmpq_get_mpq(tmp, r0);
  double ratio = mpq_get_d(tmp);
  mpq_clear(tmp);

  int r = (fabs(ratio - 1.0) < 0.1);

  printf("n: %4ld, bits: %4ld, exact: %7.2fs, upper: %7.2fs, truncated: %8.2f, exact/upper: %8.2f, exact/approx: %8.2f ratio: %8.2f", n, bits,
         t0/1000000.0, t1/1000000.0, t2/1000000.0, (double)t0/(double)t1, (double)t0/(double)t2, ratio);
  if (r)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  fmpq_poly_clear(f);
  fmpq_poly_clear(g);
  fmpq_clear(r0);
  fmpq_clear(r1);
  fmpq_clear(r2);
  return !r;
}


int main(int argc, char *argv[]) {

  flint_rand_t state;
  flint_randinit(state);

  int status = 0;

  unsigned long n[4] = {32,64,128,0};

  for(int i=0; n[i]; i++) {
    for(mp_bitcnt_t bits=2; bits<=2*n[i]; bits=2*bits) {
      status += test_fmpz_poly_oz_ideal_norm(n[i], bits, state);
    }
  }
  printf("\n");

  for(int i=0; n[i]; i++) {
    for(mp_bitcnt_t bits=2; bits<=2*n[i]; bits=2*bits) {
      status += test_fmpq_poly_oz_ideal_norm(n[i], bits, state);
    }
  }
  printf("\n");

  for(int i=0; n[i]; i++) {
    status += test_nmod_poly_oz_ideal_norm(n[i],state);
  }

  flint_randclear(state);
  flint_cleanup();
  return status;
}
