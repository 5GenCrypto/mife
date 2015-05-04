#include <oz/oz.h>
#include <oz/util.h>
#include <oz/flint-addons.h>

int test_fmpz_mod_poly_oz_invert(long n, long q_, flint_rand_t state) {
  fmpz_t q;
  fmpz_init(q);
  fmpz_set_si(q, q_);

  fmpz_mod_poly_t f;  fmpz_mod_poly_init(f, q);
  fmpz_mod_poly_t g;  fmpz_mod_poly_init_oz_modulus(g, q, n);

  fmpz_mod_poly_randtest(f, state, n);
  while (fmpz_mod_poly_degree(f) < n-1)
    fmpz_mod_poly_randtest(f, state, n);

  fmpz_mod_poly_t r0, r1;
  fmpz_mod_poly_init(r0, q);
  fmpz_mod_poly_init(r1, q);

  uint64_t t0 = oz_walltime(0);
  fmpz_mod_poly_invert_mod(r0, f, g);
  t0 = oz_walltime(t0);

  uint64_t t1 = oz_walltime(0);
  fmpz_mod_poly_oz_invert(r1, f, n);
  t1 = oz_walltime(t1);

  int r= fmpz_mod_poly_equal(r0, r1);

  if(!r) {
    fmpz_mod_poly_print_pretty(r0, "x"); printf("\n");
    fmpz_mod_poly_print_pretty(r1, "x"); printf("\n");
  }

  printf("n: %4ld,    q: %4ld, xgcd: %7.2fs, oz: %7.2fs, xgcd/oz: %7.2f ", n, q_,
         oz_seconds(t0), oz_seconds(t1), (double)t0/(double)t1);
  if (r)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  fmpz_mod_poly_clear(f);
  fmpz_mod_poly_clear(g);
  fmpz_mod_poly_clear(r0);
  fmpz_mod_poly_clear(r1);
  return !r;
}

int test_fmpq_poly_oz_invert(long n, mp_bitcnt_t bits, flint_rand_t state) {
  fmpq_poly_t f;  fmpq_poly_init(f);
  fmpq_poly_t g;  fmpq_poly_init_oz_modulus(g, n);

  fmpq_poly_randtest(f, state, n, bits);
  while (fmpq_poly_degree(f) < n-1)
    fmpq_poly_randtest(f, state, n, bits);

  fmpq_poly_t r0, r1;
  fmpq_poly_init(r0);
  fmpq_poly_init(r1);

  uint64_t t0 = oz_walltime(0);
  fmpq_poly_invert_mod(r0, f, g);
  t0 = oz_walltime(t0);

  uint64_t t1 = oz_walltime(0);
  fmpq_poly_oz_invert_approx(r1, f, n, 0, 0);
  t1 = oz_walltime(t1);

  int r= fmpq_poly_equal(r0, r1);

  if(!r) {
    fmpq_poly_print_pretty(r0, "x"); printf("\n");
    fmpq_poly_print_pretty(r1, "x"); printf("\n");
    exit(1);
  }

  printf("n: %4ld, bits: %4ld, xgcd: %7.2fs, oz: %7.2fs, xgcd/oz: %7.2f ", n, bits,
         oz_seconds(t0), oz_seconds(t1), (double)t0/(double)t1);
  if (r)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  fmpq_poly_clear(f);
  fmpq_poly_clear(g);
  fmpq_poly_clear(r0);
  fmpq_poly_clear(r1);
  return !r;
}


int main(int argc, char *argv[]) {

  flint_rand_t state;
  flint_randinit(state);

  int status = 0;

  unsigned long n[5] = {16,32,64,128,0};

  for(int i=0; n[i]; i++)
    for(long q=n_nextprime(1,0); q<100; q = n_nextprime(q, 0))
      status += test_fmpz_mod_poly_oz_invert(n[i], q, state);

  printf("\n");
  for(int i=0; n[i]; i++)
    for(mp_bitcnt_t bits=1; bits < n[i]; bits=2*bits)
      status += test_fmpq_poly_oz_invert(n[i], bits, state);

  flint_randclear(state);
  flint_cleanup();
  return status;
}
