#include <gghlite/gghlite.h>
#include <gghlite/gghlite-internals.h>

/**
 * Testing the flexibility of large index sets with a smaller kappa
 */
int test_jigsaw_indices(const size_t lambda, const size_t kappa, const size_t gamma, flint_rand_t randstate) {

	printf("lambda: %d, kappa: %d, gamma: %d", (int) lambda, (int) kappa, (int) gamma);

  gghlite_sk_t self;
  gghlite_flag_t flags = GGHLITE_FLAGS_QUIET;
  gghlite_jigsaw_init_gamma(self, lambda, kappa, gamma, flags, randstate);

  fmpz_t p; fmpz_init(p);
  fmpz_poly_oz_ideal_norm(p, self->g, self->params->n, 0);

	/* partitioning of the universe set */
	int partition[GAMMA];
	for (int i = 0; i < gamma; i++) {
		partition[i] = rand() % kappa;
	}

  fmpz_t a[kappa];
  fmpz_t acc;  fmpz_init(acc);
  fmpz_set_ui(acc, 1);

  for(size_t k=0; k<kappa; k++) {
    fmpz_init(a[k]);
    fmpz_randm(a[k], randstate, p);
    fmpz_mul(acc, acc, a[k]);
    fmpz_mod(acc, acc, p);
  }

  gghlite_clr_t e[kappa];
  gghlite_enc_t u[kappa];

  for(size_t k=0; k<kappa; k++) {
    gghlite_clr_init(e[k]);
    gghlite_enc_init(u[k], self->params);
  }

  gghlite_enc_t left;
  gghlite_enc_init(left, self->params);
  gghlite_enc_set_ui0(left, 1, self->params);

  for(size_t k=0; k<kappa; k++) {
    fmpz_poly_set_coeff_fmpz(e[k], 0, a[k]);
    const int i = (gghlite_sk_is_symmetric(self)) ? 0 : k;
		int group[GAMMA];
		memset(group, 0, GAMMA * sizeof(int));
		for (int j = 0; j < gamma; j++) {
			if (partition[j] == k) {
				group[j] = 1;
			}
		}
    gghlite_enc_set_gghlite_clr(u[k], self, e[k], 1, group, 1, randstate);
    gghlite_enc_mul(left, self->params, left, u[k]);
  }

  gghlite_enc_t rght;
  gghlite_enc_init(rght, self->params);
  gghlite_enc_set_ui0(rght, 1, self->params);

  fmpz_poly_t tmp; fmpz_poly_init(tmp);
  fmpz_poly_set_coeff_fmpz(tmp, 0, acc);
  gghlite_enc_set_gghlite_clr0(rght, self, tmp, randstate);

  for(size_t k=0; k<gamma; k++) {
    gghlite_enc_mul(rght, self->params, rght, self->z_inv[k]);
  }
 
  gghlite_enc_sub(rght, self->params, rght, left);
  int status = 1 - gghlite_enc_is_zero(self->params, rght);

  for(size_t i=0; i<kappa; i++) {
    fmpz_clear(a[i]);
    gghlite_clr_clear(e[i]);
    gghlite_enc_clear(u[i]);
  }

  gghlite_enc_clear(left);
  gghlite_enc_clear(rght);
  gghlite_clr_clear(tmp);
  fmpz_clear(p);
  gghlite_sk_clear(self, 1);

  if (status == 0)
    printf(" PASS\n");
  else
    printf(" FAIL\n");

  return status;
}

int test_jigsaw(const size_t lambda, const size_t kappa, int symmetric, flint_rand_t randstate) {

  printf("λ: %4zu, κ: %2zu, symmetric: %d …", lambda, kappa, symmetric);

  gghlite_sk_t self;
  gghlite_flag_t flags = GGHLITE_FLAGS_QUIET;
  if (symmetric)
    gghlite_init(self, lambda, kappa, kappa, 0x0, flags | GGHLITE_FLAGS_GOOD_G_INV, randstate);
  else
    gghlite_jigsaw_init(self, lambda, kappa, flags, randstate);

  fmpz_t p; fmpz_init(p);
  fmpz_poly_oz_ideal_norm(p, self->g, self->params->n, 0);


  fmpz_t a[kappa];
  fmpz_t acc;  fmpz_init(acc);
  fmpz_set_ui(acc, 1);

  for(size_t k=0; k<kappa; k++) {
    fmpz_init(a[k]);
    fmpz_randm(a[k], randstate, p);
    fmpz_mul(acc, acc, a[k]);
    fmpz_mod(acc, acc, p);
  }

  gghlite_clr_t e[kappa];
  gghlite_enc_t u[kappa];

  for(size_t k=0; k<kappa; k++) {
    gghlite_clr_init(e[k]);
    gghlite_enc_init(u[k], self->params);
  }

  gghlite_enc_t left;
  gghlite_enc_init(left, self->params);
  gghlite_enc_set_ui0(left, 1, self->params);

  for(size_t k=0; k<kappa; k++) {
    fmpz_poly_set_coeff_fmpz(e[k], 0, a[k]);
    const int i = (gghlite_sk_is_symmetric(self)) ? 0 : k;
		int group[GAMMA];
		memset(group, 0, GAMMA * sizeof(int));
		group[i] = 1;
    gghlite_enc_set_gghlite_clr(u[k], self, e[k], 1, group, 1, randstate);
    gghlite_enc_mul(left, self->params, left, u[k]);
  }

  gghlite_enc_t rght;
  gghlite_enc_init(rght, self->params);
  gghlite_enc_set_ui0(rght, 1, self->params);

  fmpz_poly_t tmp; fmpz_poly_init(tmp);
  fmpz_poly_set_coeff_fmpz(tmp, 0, acc);
  gghlite_enc_set_gghlite_clr0(rght, self, tmp, randstate);

  for(size_t k=0; k<kappa; k++) {
    const int i = (gghlite_sk_is_symmetric(self)) ? 0 : k;
    gghlite_enc_mul(rght, self->params, rght, self->z_inv[i]);
  }

  gghlite_enc_sub(rght, self->params, rght, left);
  int status = 1 - gghlite_enc_is_zero(self->params, rght);

  for(size_t i=0; i<kappa; i++) {
    fmpz_clear(a[i]);
    gghlite_clr_clear(e[i]);
    gghlite_enc_clear(u[i]);
  }

  gghlite_enc_clear(left);
  gghlite_enc_clear(rght);
  gghlite_clr_clear(tmp);
  fmpz_clear(p);
  gghlite_sk_clear(self, 1);

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

  status += test_jigsaw(20, 2, 1, randstate);
  status += test_jigsaw(20, 3, 1, randstate);
  status += test_jigsaw(20, 4, 1, randstate);

  status += test_jigsaw(20, 2, 0, randstate);
  status += test_jigsaw(20, 3, 0, randstate);
  status += test_jigsaw(20, 4, 0, randstate);

  status += test_jigsaw_indices(20, 4, 90, randstate);
  status += test_jigsaw_indices(20, 20, 60, randstate);



  flint_randclear(randstate);
  flint_cleanup();
  mpfr_free_cache();
  return status;
}
