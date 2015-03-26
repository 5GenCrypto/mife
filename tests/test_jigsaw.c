#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>

int test_jigsaw(const size_t lambda, const size_t kappa, flint_rand_t randstate) {

  printf("λ: %4zu, κ: %2zu", lambda, kappa);

  gghlite_t self;
  const gghlite_flag_t flags = GGHLITE_FLAGS_GDDH_HARD | GGHLITE_FLAGS_QUIET;
  uint64_t rerand = 1;
  gghlite_init(self, lambda, kappa, rerand, flags, randstate);

  fmpz_t p; fmpz_init(p);
  fmpz_poly_oz_ideal_norm(p, self->g, self->pk->n, 0);

  fmpz_t        a[kappa];
  gghlite_clr_t e[kappa];
  gghlite_enc_t u[kappa];

  fmpz_t a_acc;  fmpz_init(a_acc);
  fmpz_set_ui(a_acc, 1);

  gghlite_enc_t u_acc;
  gghlite_enc_init(u_acc, self->pk);
  gghlite_enc_set_ui(u_acc, 1, self->pk);

  for(size_t i=0; i<kappa; i++) {
    // sample random elements in ZZ_p
    fmpz_init(a[i]);
    fmpz_randm(a[i], randstate, p);

    fmpz_mul(a_acc, a_acc, a[i]);
    fmpz_mod(a_acc, a_acc, p);

    gghlite_clr_init(e[i]);
    fmpz_poly_set_coeff_fmpz(e[i], 0, a[i]);
    gghlite_enc_init(u[i], self->pk);
    gghlite_enc0(u[i], self, e[i], randstate);
    gghlite_elevate(u[i], self->pk, u[i], 1, 0, 0, 0, randstate);

    gghlite_mul(u_acc, self->pk, u_acc, u[i]);
  }

  gghlite_clr_t rc;  gghlite_clr_init(rc);
  fmpz_poly_set_coeff_fmpz(rc, 0, a_acc);
  gghlite_enc_t re; gghlite_enc_init(re, self->pk);
  gghlite_enc0(re, self, rc, randstate);
  gghlite_elevate(re, self->pk, re, kappa, 0, 0, 0, randstate);

  gghlite_sub(re, self->pk, re, u_acc);

  int status = 1 - gghlite_is_zero(self->pk, re);

  gghlite_extract(rc, self->pk, re);

  fmpz_poly_print_pretty(rc, "x");

  gghlite_enc_clear(re);
  gghlite_clr_clear(rc);

  for(size_t i=0; i<kappa; i++) {
    fmpz_clear(a[i]);
    gghlite_clr_clear(e[i]);
    gghlite_enc_clear(u[i]);
  }

  gghlite_enc_clear(u_acc);
  fmpz_clear(a_acc);
  fmpz_clear(p);
  gghlite_clear(self, 1);

  if (status == 0)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  return status;
}


int main(int argc, char *argv[]) {
  flint_rand_t randstate;
  flint_randinit_seed(randstate, 0x1337, 1);

  int status = 0;

  status += test_jigsaw(20, 2, randstate);

  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
  return status;
}
