#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include <oz/oz.h>

int test_randtest(slong n, mp_bitcnt_t bits, flint_rand_t state) {
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

  uint64_t t0 = ggh_walltime(0);
  fmpz_poly_resultant_modular(r0, f, g);
  t0 = ggh_walltime(t0);

  uint64_t t1 = ggh_walltime(0);
  fmpz_poly_oz_ideal_norm(r1, f, n, 0);
  t1 = ggh_walltime(t1);

  int r= fmpz_equal(r0,r1);

  printf("n: %4ld, bits: %4ld, flint mulmod: %7.2fs, mulmod bound: %7.2fs, ratio mulmod/bound: %7.2f, passes: %d\n", n, bits,
         t0/1000000.0, t1/1000000.0, (double)t0/(double)t1, r);

  fmpz_poly_clear(f);
  fmpz_poly_clear(g);
  fmpz_clear(r0);
  fmpz_clear(r1);
  fmpz_clear(r2);
  return r;
}

int main(int argc, char *argv[]) {

  flint_rand_t state;
  flint_randinit(state);

  int status = 0;

  unsigned long n[4] = {32,64,128,0};

  for(int i=0; n[i]; i++) {
    for(mp_bitcnt_t bits=2; bits<=2*n[i]; bits=2*bits) {
      status += test_randtest(n[i], bits, state);
    }
  }

  flint_randclear(state);

  return status;
}
