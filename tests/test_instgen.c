#include <gghlite/gghlite.h>

int test_instgen_symm(const size_t lambda, const size_t kappa, const uint64_t rerand, aes_randstate_t randstate) {

  printf("symm: 1, λ: %4zu, κ: %2zu, rerand: 0x%016zx", lambda, kappa, rerand);

  gghlite_sk_t self;
  gghlite_flag_t flags = GGHLITE_FLAGS_GDDH_HARD | GGHLITE_FLAGS_QUIET;
  gghlite_init(self, lambda, kappa, 1 /* gamma */, rerand, flags, randstate);

  int status = 0;

  /* The elements a & y should be zero if there is no rerandomisation, and non-zero otherwise.  */

  if (gghlite_params_have_rerand(self->params, 0)) {
    if (fmpz_poly_is_zero(self->a[0]))              status++;
    if (fmpz_mod_poly_is_zero(self->params->y[0]))  status++;
  } else {
    if (!fmpz_poly_is_zero(self->a[0]))              status++;
    if (!fmpz_mod_poly_is_zero(self->params->y[0]))  status++;
  }

  /* We want rerandomisers for any level where gghlite_pk_have_rerand(self->params, k) is true.  */

  for(size_t k=0; k< kappa; k++) {
    if (gghlite_params_have_rerand(self->params, k)) {
      if (fmpz_poly_is_zero(self->b[0][k][0]))         status++;
      if (fmpz_poly_is_zero(self->b[0][k][1]))         status++;
      if (fmpz_mod_poly_is_zero(self->params->x[0][k][0])) status++;
      if (fmpz_mod_poly_is_zero(self->params->x[0][k][1])) status++;
    } else {
      if (!fmpz_poly_is_zero(self->b[0][k][0]))         status++;
      if (!fmpz_poly_is_zero(self->b[0][k][1]))         status++;
      if (!fmpz_mod_poly_is_zero(self->params->x[0][k][0])) status++;
      if (!fmpz_mod_poly_is_zero(self->params->x[0][k][1])) status++;
    }
  }

  /* In the symmetric setting only z[0] and z_inv[0] should be used */

  if (fmpz_mod_poly_is_zero(self->z[0]))       status++;
  if (fmpz_mod_poly_is_zero(self->z_inv[0]))   status++;
  if (!fmpz_mod_poly_is_zero(self->z[1]))      status++;
  if (!fmpz_mod_poly_is_zero(self->z_inv[1]))  status++;

  if (status == 0)
    printf(" (%d) PASS\n", status);
  else
    printf(" (%d) FAIL\n", status);

  gghlite_sk_clear(self, 1);
  return status;
}


int test_instgen_asymm(const size_t lambda, const size_t kappa, const uint64_t rerand, aes_randstate_t randstate) {

  printf("symm: 0, λ: %4zu, κ: %2zu, rerand: 0x%016zx", lambda, kappa, rerand);

  gghlite_sk_t self;
  const gghlite_flag_t flags = GGHLITE_FLAGS_GDDH_HARD | GGHLITE_FLAGS_QUIET | GGHLITE_FLAGS_ASYMMETRIC;
  gghlite_init(self, lambda, kappa, kappa /* gamma */, rerand, flags, randstate);

  int status = 0;

  for(size_t i=0; i< kappa; i++) {
    /* we want all z[i], z_inv[i] to be != 0 */
    if (fmpz_mod_poly_is_zero(self->z[i]))         status++;
    if (fmpz_mod_poly_is_zero(self->z_inv[i]))     status++;

    if (gghlite_params_have_rerand(self->params, 0)) {
      if (fmpz_poly_is_zero(self->a[i]))             status++;
      if (fmpz_mod_poly_is_zero(self->params->y[i])) status++;
    } else {
      if (!fmpz_poly_is_zero(self->a[i]))             status++;
      if (!fmpz_mod_poly_is_zero(self->params->y[i])) status++;
    }

    for(size_t k=0; k< kappa; k++) {
    /* We want rerandomisers for any source group where gghlite_pk_have_rerand(self->params, k) is true
       and not anywhere else */

      int s = 0;
      if (fmpz_poly_is_zero(self->b[i][k][0]))         s++;
      if (fmpz_poly_is_zero(self->b[i][k][1]))         s++;

      if (fmpz_mod_poly_is_zero(self->params->x[i][k][0])) s++;
      if (fmpz_mod_poly_is_zero(self->params->x[i][k][1])) s++;

      if (gghlite_params_have_rerand(self->params, k)) {
        status += s;
      } else {
        status += (4-s);
      }
    }
  }

  if (status == 0)
    printf(" (%d) PASS\n", status);
  else
    printf(" (%d) FAIL\n", status);

  gghlite_sk_clear(self, 1);
  return status;
}


int main(int argc, char *argv[]) {

  aes_randstate_t randstate;
  aes_randinit(randstate);

  int status = 0;

  status += test_instgen_symm(20, 2, 0x0, randstate);
  status += test_instgen_symm(20, 2, 0x1, randstate);
  status += test_instgen_symm(20, 4, 0x0, randstate);
  status += test_instgen_symm(20, 4, 0x1, randstate);

  status += test_instgen_asymm(20, 2, 0x0, randstate);
  status += test_instgen_asymm(20, 2, 0x1, randstate);
  status += test_instgen_asymm(20, 4, 0x0, randstate);
  status += test_instgen_asymm(20, 4, 0x1, randstate);

  aes_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
  return status;
}
