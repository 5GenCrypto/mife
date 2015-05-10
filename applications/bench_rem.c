#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include <oz/oz.h>

void bench_rem(const fmpz_poly_t g, const fmpz_t a, const long n, const mp_bitcnt_t prec) {
  fmpq_poly_t gq;
  fmpq_poly_init(gq);
  fmpq_poly_set_fmpz_poly(gq, g);

  fmpq_poly_t g_inv;
  fmpq_poly_init(g_inv);

  fmpq_poly_oz_invert_approx(g_inv, gq, n, prec, 0);

  fmpz_poly_t f; fmpz_poly_init(f);
  fmpz_poly_t r; fmpz_poly_init(r);
  fmpz_poly_set_coeff_fmpz(f, 0, a);

  uint64_t t = ggh_walltime(0);
  _fmpz_poly_oz_rem_small_iter(r, f, g, n, g_inv, prec, OZ_VERBOSE);
  printf("n: %6ld, log σ: %4zu, prec: %5ld, t: %10.4f\n", n, labs(fmpz_poly_max_bits(g)), prec, ggh_seconds(ggh_walltime(t)));

  fmpz_poly_clear(f);
  fmpz_poly_clear(r);
  fmpq_poly_clear(g_inv);
  fmpq_poly_clear(gq);
}

int main(int argc, char *argv[]) {
  assert(argc==2);
  const long n = atol(argv[1]);

  flint_rand_t randstate;
  flint_randinit_seed(randstate, 0x1337+1, 1);

  mpfr_t sigma;
  mpfr_init2(sigma, 80);
  mpfr_set_d(sigma, _gghlite_sigma(n), MPFR_RNDN);


  fmpz_poly_t g;
  fmpz_poly_init(g);
  fmpz_poly_sample_sigma(g, n, sigma, randstate);

  fmpz_t p;  fmpz_init(p);
  fmpz_poly_oz_ideal_norm(p, g, n, 0);

  fmpz_t a;  fmpz_init(a);
  fmpz_randm(a, randstate, p);

  for(size_t i=11; i<=log2(n)+1; i++) {
    bench_rem(g, a, n,  1ULL<<i);
  }

  fmpz_clear(p);
  fmpz_clear(a);
  fmpz_poly_clear(g);
  mpfr_clear(sigma);

  flint_randclear(randstate);
  flint_cleanup();
  return 0;
}
