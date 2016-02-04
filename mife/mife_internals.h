#ifndef _MIFE_INTERNALS_H_
#define _MIFE_INTERNALS_H_

#include "mife_defs.h"

/* MIFE internal functions */
void reset_T();
float get_T();
void set_NUM_ENC(int val);
int get_NUM_ENC();
void mife_apply_randomizer(mife_pp_t pp, aes_randstate_t randstate, fmpz_mat_t m);
void mife_apply_randomizers(mife_mat_clr_t met, mife_pp_t pp, mife_sk_t sk,
    aes_randstate_t randstate);
int ***mife_partitions(mife_pp_t pp, fmpz_t index);
void mife_partitions_clear(mife_pp_t pp, int ***partitions);
void mife_apply_kilian(mife_pp_t pp, mife_sk_t sk, fmpz_mat_t m, int global_index);
void mife_set_encodings(mife_ciphertext_t ct, mife_mat_clr_t met, fmpz_t index,
    mife_pp_t pp, mife_sk_t sk, aes_randstate_t randstate);
void mife_clear_pp_read(mife_pp_t pp);
void mife_clear_pp(mife_pp_t pp);
void mife_clear_sk(mife_sk_t sk);
void mife_mat_clr_clear(mife_pp_t pp, mife_mat_clr_t met);
void mife_gen_partitioning(int *partitioning, fmpz_t index, int L, int nu);
void mife_mat_encode(mife_pp_t pp, mife_sk_t sk, gghlite_enc_mat_t enc,
    fmpz_mat_t m, int *group, aes_randstate_t randstate);
void gghlite_enc_mat_zeros_print(mife_pp_t pp, gghlite_enc_mat_t m);
void gghlite_enc_mat_mul(gghlite_params_t params, gghlite_enc_mat_t r,
    gghlite_enc_mat_t m1, gghlite_enc_mat_t m2);
void mife_ciphertext_clear(mife_pp_t pp, mife_ciphertext_t ct);
void message_to_dary(ulong *dary, int bitstring_len, fmpz_t message, int d);

/* functions dealing with fmpz types and matrix multiplications mod fmpz_t */
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);
void fmpz_mat_scalar_mul_modp(fmpz_mat_t m, fmpz_t scalar, fmpz_t modp);
void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n,
    fmpz_t p);
void gghlite_enc_mat_init(gghlite_params_t params, gghlite_enc_mat_t m,
    int nrows, int ncols);
void gghlite_enc_mat_clear(gghlite_enc_mat_t m);
void fmpz_init_exp(fmpz_t exp, int base, int n);

/* functions for computing matrix inverse mod fmpz_t */
void fmpz_mat_modp(fmpz_mat_t m, int dim, fmpz_t p);
void fmpz_mat_det_modp(fmpz_t det, fmpz_mat_t a, int n, fmpz_t p);
void fmpz_mat_cofactor_modp(fmpz_mat_t b, fmpz_mat_t a, int n, fmpz_t p);
void fmpz_modp_matrix_inverse(fmpz_mat_t inv, fmpz_mat_t a, int dim, fmpz_t p);


#endif /* _MIFE_INTERNALS_H_ */
