#include <oz/oz.h>
#include <oz/util.h>
#include <oz/flint-addons.h>

int test_fmpz_mod_poly_oz_mul_fftnwc(long n, mp_bitcnt_t bits, aes_randstate_t state) {

  fmpz_t q;
  fmpz_init(q);
  fmpz_zero(q);
  fmpz_randbits_aes(q, state, bits);
  fmpz_abs(q, q);
  fmpz_fdiv_q_2exp(q, q, n_flog(n,2)+1);
  fmpz_mul_2exp(q, q, n_flog(n,2)+1);
  fmpz_add_ui(q, q, 1);
  for(int i=0; ; i++) {
    if (fmpz_is_probabprime(q))
      break;
    fmpz_add_ui(q, q, 2*n);
  }
  fmpz_mod_poly_t f0;  fmpz_mod_poly_init(f0, q);
  fmpz_mod_poly_t f1;  fmpz_mod_poly_init(f1, q);
  fmpz_mod_poly_t g;  fmpz_mod_poly_init_oz_modulus(g, q, n);

  fmpz_mod_poly_randtest_aes(f0, state, n);
  while (fmpz_mod_poly_degree(f0) < n-1)
    fmpz_mod_poly_randtest_aes(f0, state, n);

  fmpz_mod_poly_randtest_aes(f1, state, n);
  while (fmpz_mod_poly_degree(f1) < n-1)
    fmpz_mod_poly_randtest_aes(f1, state, n);

  fmpz_mod_poly_t r0, r1;
  fmpz_mod_poly_init(r0, q);
  fmpz_mod_poly_init(r1, q);

  uint64_t t0 = oz_walltime(0);
  fmpz_mod_poly_oz_mul(r0, f0, f1, n);
  t0 = oz_walltime(t0);

  fmpz_mod_poly_oz_ntt_precomp_t precomp;
  fmpz_mod_poly_oz_ntt_precomp_init(precomp, n, q);

  uint64_t t1 = oz_walltime(0);
  _fmpz_mod_poly_oz_mul_nttnwc(r1, f0, f1, precomp);
  t1 = oz_walltime(t1);

  fmpz_mod_poly_oz_ntt_precomp_clear(precomp);

  int r = fmpz_mod_poly_equal(r0, r1);

  if(!r) {
    fmpz_mod_poly_print_pretty(r0, "x"); printf("\n");
    fmpz_mod_poly_print_pretty(r1, "x"); printf("\n");
  }
  printf("n: %6ld, log(q): %6ld, flint: %7.2fs, oz: %7.2fs, flint/oz: %7.2f ", n, fmpz_sizeinbase(q,2),
         oz_seconds(t0), oz_seconds(t1), (double)t0/(double)t1);
  if (r)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  fmpz_mod_poly_clear(f0);
  fmpz_mod_poly_clear(f1);
  fmpz_mod_poly_clear(g);
  fmpz_mod_poly_clear(r0);
  fmpz_mod_poly_clear(r1);
  return !r;
}

int test_fmpz_mod_poly_oz_mul(long n, long q_, aes_randstate_t state) {
  fmpz_t q;
  fmpz_init(q);
  fmpz_set_si(q, q_);

  fmpz_mod_poly_t f0;  fmpz_mod_poly_init(f0, q);
  fmpz_mod_poly_t f1;  fmpz_mod_poly_init(f1, q);
  fmpz_mod_poly_t g;  fmpz_mod_poly_init_oz_modulus(g, q, n);

  fmpz_mod_poly_randtest_aes(f0, state, n);
  while (fmpz_mod_poly_degree(f0) < n-1)
    fmpz_mod_poly_randtest_aes(f0, state, n);

  fmpz_mod_poly_randtest_aes(f1, state, n);
  while (fmpz_mod_poly_degree(f1) < n-1)
    fmpz_mod_poly_randtest_aes(f1, state, n);

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

  printf("n: %6ld,      q: %6ld, flint: %7.2fs, oz: %7.2fs, flint/oz: %7.2f ", n, q_,
         oz_seconds(t0), oz_seconds(t1), (double)t0/(double)t1);
  if (r)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  fmpz_mod_poly_clear(f0);
  fmpz_mod_poly_clear(f1);
  fmpz_mod_poly_clear(g);
  fmpz_mod_poly_clear(r0);
  fmpz_mod_poly_clear(r1);
  return !r;
}

int main(int argc, char *argv[]) {

  aes_randstate_t state;
  aes_randinit(state);

  int status = 0;

  mp_bitcnt_t bits[5] = {10,11,12,13,0};

  for(int i=0; bits[i]; i++) {
    unsigned long n = ((unsigned long)1)<<bits[i];
    status += test_fmpz_mod_poly_oz_mul_fftnwc(n, n/2, state);
  }

  for(int i=0; bits[i]; i++) {
    unsigned long n = ((unsigned long)1)<<bits[i];
    for(unsigned long q=n_nextprime(n,0); q<n+100; q = n_nextprime(q, 0)) {
      status += test_fmpz_mod_poly_oz_mul(n, q, state);
    }
  }
  aes_randclear(state);
  flint_cleanup();
  return status;
}
