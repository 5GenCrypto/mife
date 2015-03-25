#include <gghlite/gghlite.h>

int test_instgen_symm(const size_t lambda, const size_t kappa, const uint64_t rerand, flint_rand_t randstate) {

  printf("symm: 1, λ: %4zu, κ: %2zu, rerand: 0x%016zx", lambda, kappa, rerand);

  gghlite_t self;
  gghlite_pk_init_params(self->pk, lambda, kappa, rerand, GGHLITE_FLAGS_GDDH_HARD | GGHLITE_FLAGS_QUIET);
  gghlite_init_instance(self, randstate);

  int status = 0;

  /* The elements a & y should be zero if there is no rerandomisation, and non-zero otherwise.  */

  if(rerand & 1) {
    if (fmpz_poly_is_zero(self->a[0]))          status++;
    if (fmpz_mod_poly_is_zero(self->pk->y[0]))  status++;
  } else {
    if (!fmpz_poly_is_zero(self->a[0]))         status++;
    if (!fmpz_mod_poly_is_zero(self->pk->y[0])) status++;
  }

  /* We want rerandomisers for any level where gghlite_pk_have_rerand(self->pk, k) is true.  */

  for(size_t k=0; k< kappa; k++) {
    if (gghlite_pk_have_rerand(self->pk, k)) {
      if (fmpz_poly_is_zero(self->b[k][0]))         status++;
      if (fmpz_poly_is_zero(self->b[k][1]))         status++;
      if (fmpz_mod_poly_is_zero(self->pk->x[k][0])) status++;
      if (fmpz_mod_poly_is_zero(self->pk->x[k][1])) status++;
    } else {
      if (!fmpz_poly_is_zero(self->b[k][0]))         status++;
      if (!fmpz_poly_is_zero(self->b[k][1]))         status++;
      if (!fmpz_mod_poly_is_zero(self->pk->x[k][0])) status++;
      if (!fmpz_mod_poly_is_zero(self->pk->x[k][0])) status++;
    }
  }

  /* In the symmetric setting only z[0] and z_inv[0] should be used */

  if (!fmpz_mod_poly_is_zero(self->z[1]))      status++;
  if (!fmpz_mod_poly_is_zero(self->z_inv[1]))  status++;

  if (status == 0)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  gghlite_clear(self, 1);
  return status;
}


int test_instgen_asymm(const size_t lambda, const size_t kappa, const uint64_t rerand, flint_rand_t randstate) {

  printf("symm: 0, λ: %4zu, κ: %2zu, rerand: 0x%016zx", lambda, kappa, rerand);

  gghlite_t self;
  const gghlite_flag_t flags = GGHLITE_FLAGS_GDDH_HARD | GGHLITE_FLAGS_QUIET | GGHLITE_FLAGS_ASYMMETRIC;
  gghlite_pk_init_params(self->pk, lambda, kappa, rerand, flags);
  gghlite_init_instance(self, randstate);

  int status = 0;

  for(size_t k=0; k< kappa; k++) {

    /* We want rerandomisers for any source group where gghlite_pk_have_rerand(self->pk, k) is
       true and not anywhere else */

    if (gghlite_pk_have_rerand(self->pk, k)) {

      if (fmpz_poly_is_zero(self->a[k]))            status++;
      if (fmpz_poly_is_zero(self->b[k][0]))         status++;
      if (fmpz_poly_is_zero(self->b[k][1]))         status++;

      if (fmpz_mod_poly_is_zero(self->pk->y[k]))    status++;
      if (fmpz_mod_poly_is_zero(self->pk->x[k][0])) status++;
      if (fmpz_mod_poly_is_zero(self->pk->x[k][1])) status++;

    } else {

      if (!fmpz_poly_is_zero(self->a[k]))           status++;
      if (!fmpz_poly_is_zero(self->b[k][0]))        status++;
      if (!fmpz_poly_is_zero(self->b[k][1]))        status++;

      if (!fmpz_mod_poly_is_zero(self->pk->y[k])) status++;
      if (!fmpz_mod_poly_is_zero(self->pk->x[k][0])) status++;
      if (!fmpz_mod_poly_is_zero(self->pk->x[k][0])) status++;

    }

    /* we want all z[i] and z_inv[i] to be != 0 */

    if (fmpz_mod_poly_is_zero(self->z[k]))      status++;
    if (fmpz_mod_poly_is_zero(self->z_inv[k]))  status++;
  }


  if (status == 0)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  gghlite_clear(self, 1);
  return status;
}


int main(int argc, char *argv[]) {

  flint_rand_t randstate;
  flint_randinit_seed(randstate, 0x1337, 1);

  int status = 0;

  status += test_instgen_symm(20, 2, 0x0, randstate);
  status += test_instgen_symm(20, 2, 0x1, randstate);
  status += test_instgen_symm(20, 4, 0x0, randstate);
  status += test_instgen_symm(20, 4, 0x1, randstate);

  status += test_instgen_asymm(20, 2, 0x0, randstate);
  status += test_instgen_asymm(20, 2, 0x1, randstate);
  status += test_instgen_asymm(20, 4, 0x0, randstate);
  status += test_instgen_asymm(20, 4, 0xf, randstate);

  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
  return status;
}
