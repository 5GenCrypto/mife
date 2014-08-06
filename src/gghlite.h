#ifndef GGHLITE__H
#define GGHLITE__H

#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>

#include <stdint.h>
#include <math.h>

/**
   Structs
*/

struct _gghlite_pk_struct {
  size_t lambda;
  size_t kappa;

  size_t n;
  fmpz_t q;
  mpfr_t sigma;   //< σ
  mpfr_t sigma_p; //< σ'
  mpfr_t sigma_s; //< σ^*
  mpfr_t ell_b;  //< l_b

  fmpz_poly_t cyclotomic_polynomial;

};

typedef struct _gghlite_pk_struct gghlite_pk_t[1];

struct _gghlite_struct {
  gghlite_pk_t pk;

  fmpz_poly_t g;
  fmpz_poly_t z;
  fmpz_poly_t h;
};

typedef struct _gghlite_struct gghlite_t[1];

/**
   Auxilery functions
*/

static inline int64_t _gghlite_log_q(const int64_t log_n, const int64_t kappa) {
  int64_t log_q = (8.5*log_n + log2(log2(kappa))/2.0 + log2(log_n)) * (8*kappa);
  return log_q;
}

static inline int _gghlite_check_sec(int64_t log_q, size_t n, size_t lambda) {
  return (3.0/8.0 * log_q < n/(double)lambda);
}


static inline int gghlite_check_sec(const gghlite_pk_t self) {
  return _gghlite_check_sec(fmpz_sizeinbase(self->q, 2), self->n, self->lambda);
}

static inline long _gghlite_prec(const gghlite_pk_t self) {
  if (self->lambda < 53)
    return 53;
  else
    return self->lambda;
}

/**
   Decclarations
*/
void gghlite_init_step1(gghlite_t self, int64_t lambda, int64_t kappa);
void gghlite_init_step2(gghlite_t self, flint_rand_t randstate);
void gghlite_init(gghlite_t self, const int64_t lambda, const int64_t kappa, flint_rand_t randstate);

void gghlite_print(const gghlite_t self);

void gghlite_pk_clear(gghlite_pk_t self);
void gghlite_clear(gghlite_t self, int clear_pk);

#endif //GGHLITE__H
