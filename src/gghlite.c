#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <dgs/dgs.h>

int main(int argc, char *argv[])
{
  fmpz_poly_t f;
  dgs_disc_gauss_dp_t *D = dgs_disc_gauss_dp_init(3.0, 0, 6, DGS_DISC_GAUSS_SIGMA2_LOGTABLE);
  D->call(D);
  dgs_disc_gauss_dp_clear(D);
}
