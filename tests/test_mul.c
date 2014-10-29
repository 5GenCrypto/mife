#include <oz/oz.h>
#include <oz/util.h>
#include <oz/flint-addons.h>

int test_fmpz_mod_poly_oz_mul(long n, long q_, flint_rand_t state) {
  fmpz_t q;
  fmpz_init(q);
  fmpz_set_si(q, q_);

  fmpz_mod_poly_t f0;  fmpz_mod_poly_init(f0, q);
  fmpz_mod_poly_t f1;  fmpz_mod_poly_init(f1, q);
  fmpz_mod_poly_t g;  fmpz_mod_poly_init_oz_modulus(g, q, n);

  fmpz_mod_poly_randtest(f0, state, n);
  while (fmpz_mod_poly_degree(f0) < n-1)
    fmpz_mod_poly_randtest(f0, state, n);

  fmpz_mod_poly_randtest(f1, state, n);
  while (fmpz_mod_poly_degree(f1) < n-1)
    fmpz_mod_poly_randtest(f1, state, n);

  fmpz_mod_poly_t r0, r1;
  fmpz_mod_poly_init(r0, q);
  fmpz_mod_poly_init(r1, q);

  uint64_t t0 = oz_walltime(0);
  fmpz_mod_poly_mulmod(r0, f0, f1, g);
  t0 = oz_walltime(t0);

  uint64_t t1 = oz_walltime(0);
  fmpz_mod_poly_oz_mul(r1, f0, f1, n);
  t1 = oz_walltime(t1);

  int r= fmpz_mod_poly_equal(r0, r1);

  if(!r) {
    fmpz_mod_poly_print_pretty(r0, "x"); printf("\n");
    fmpz_mod_poly_print_pretty(r1, "x"); printf("\n");
  }

  printf("n: %6ld, q: %6ld, flint: %7.2fs, oz: %7.2fs, flint/oz: %7.2f, passes: %d\n", n, q_,
         t0/1000000.0, t1/1000000.0, (double)t0/(double)t1, r);

  fmpz_mod_poly_clear(f0);
  fmpz_mod_poly_clear(f1);
  fmpz_mod_poly_clear(g);
  fmpz_mod_poly_clear(r0);
  fmpz_mod_poly_clear(r1);
  return !r;
}

int main(int argc, char *argv[]) {

  flint_rand_t state;
  flint_randinit(state);

  int status = 0;

  mp_bitcnt_t bits[5] = {11,12,13,14,0};

  for(int i=0; bits[i]; i++) {
    unsigned long n = ((unsigned long)1)<<bits[i];
    for(unsigned long q=n_nextprime(n,0); q<n+100; q = n_nextprime(q, 0))
      status += test_fmpz_mod_poly_oz_mul(n, q, state);
  }
  flint_randclear(state);

  return status;
}
