#ifndef GPV__H
#define GPV__H

#include <mpfr.h>
#include <dgs/dgs.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>

/**
   Sampling algorithms
*/

typedef enum {
  DGS_LATTICE_DETECT               = 0x0, //< detect which algorithm to use
  DGS_LATTICE_IDENTITY             = 0x2, //< identity lattice
  DGS_LATTICE_INLATTICE            = 0x3, //< c is in the lattice
  DGS_LATTICE_COSET                = 0x4, //< c is not in the lattice
} dgs_gauss_lattice_alg_t;

struct _dgs_disc_gauss_lattice_mp_t;

typedef struct _dgs_disc_gauss_lattice_mp_t{
  fmpz_mat_t B; //< basis matrix
  fmpz *c; //< center
  mpfr_t sigma; //< Gaussian parameter
  dgs_disc_gauss_mp_t **D; //< storage for internal samplers
  int (*call)(fmpz *rop,  const struct _dgs_disc_gauss_lattice_mp_t *self, gmp_randstate_t state); //< call this function

} dgs_disc_gauss_lattice_mp_t;

/**
   @param B basis matrix (copied), if matrix is 1 x n it is assumed that it represents a rotational basis mod X^n + 1
   @param sigma Gaussian width parameter (copied)
   @param c center
   @param algorithm
*/

dgs_disc_gauss_lattice_mp_t *dgs_disc_gauss_lattice_mp_init(const fmpz_mat_t B, const mpfr_t sigma,
                                                            const fmpz *c, const dgs_gauss_lattice_alg_t algorithm);

/**
   Simple sampling when B is the identity

   @param rop return value (pre-allocated)
   @param self sampler
   @param state entropy source

*/

int dgs_disc_gauss_lattice_mp_call_simple(fmpz *rop,  const dgs_disc_gauss_lattice_mp_t *self, gmp_randstate_t state);


void dgs_disc_gauss_lattice_mp_clear(dgs_disc_gauss_lattice_mp_t *self);

#endif
