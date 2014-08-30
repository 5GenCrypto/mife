#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>
#include <flint-addons/flint-addons.h>

int test_randtest(slong n, mp_bitcnt_t bits, flint_rand_t state) {
  fmpz_poly_t f;
  fmpz_poly_t g;
  fmpz_poly_init(f);
  fmpz_poly_init(g);


  fmpz_poly_randtest(f, state, n, bits);
  while(!fmpz_poly_get_coeff_ptr(f, n-1))
      fmpz_poly_randtest(f, state, n, bits);
  fmpz_poly_set_coeff_si(g, n, 1);
  fmpz_poly_set_coeff_si(g, 0, 1);

  fmpz_t r0, r1, r2, r3;
  fmpz_init(r0);
  fmpz_init(r1);
  fmpz_init(r2);
  fmpz_init(r3);

  uint64_t t0 = ggh_walltime(0);
  fmpz_poly_resultant_modular(r0, f, g);
  t0 = ggh_walltime(t0);

  uint64_t t1 = ggh_walltime(0);
  fmpz_poly_ideal_norm(r1, f, g);
  t1 = ggh_walltime(t1);

  uint64_t t2 = ggh_walltime(0);
  //fmpz_poly_resultant_euclidean(r2, f, g);
  t2 = ggh_walltime(t2);

  int r= fmpz_equal(r0,r1);

  printf("n: %4ld, bits: %4ld, flint mulmod: %7.2fs, mulmod bound: %7.2fs, flint euclidean: %7.2fs, ratio mulmod/bound: %7.2f, ratio euclidean/bound: %7.2f, passes: %d\n", n, bits,
         t0/1000000.0, t1/1000000.0, t2/1000000.0, (double)t0/(double)t1, (double)t2/(double)t1, r);


  if (!r) {
    fmpz_print(r0);
    printf("\n");
    fmpz_print(r1);
    printf("\n");

    fmpz_mat_t M;
    fmpz_mat_init(M, n, n);
    fmpz_poly_ideal_rot_basis(M, f);
    fmpz_mat_det(r3, M);

    fmpz_print(r3);
    printf("\n");
    fmpz_clear(r3);
    fmpz_mat_clear(M);
  }

  fmpz_poly_clear(f);
  fmpz_poly_clear(g);
  fmpz_clear(r0);
  fmpz_clear(r1);
  fmpz_clear(r2);
  return r;
}

#define N 3

int main(int argc, char *argv[]) {

  flint_rand_t state;
  flint_randinit(state);

  int status =0;

  int n[N] = {32,64,128};

  for(int i=0; i<N; i++) {
    for(mp_bitcnt_t bits=2; bits<=2*n[i]; bits=2*bits) {
      status += test_randtest(n[i], bits, state);
    }
  }

  flint_randclear(state);

  return 0;
}
