#ifndef _MIFE_H_
#define _MIFE_H_

#include "mife_internals.h"
#include "mife_io.h"

int g_parallel;

/* front end functions which handle command line parsing, file IO, and all */
int mife_keygen_main (const_mmap_vtable mmap, int argc, char **argv);
int mife_encrypt_main(const_mmap_vtable mmap, int argc, char **argv);
int mife_eval_main   (const_mmap_vtable mmap, int argc, char **argv);

/* MIFE interface */
void mife_init_params(mife_pp_t pp, mife_flag_t flags);
void mife_mbp_set(
    void *mbp_params,
    mife_pp_t pp,
    int num_inputs,
    int (*paramfn)  (mife_pp_t, int),
    void (*kilianfn)(mife_pp_t, int *),
    void (*orderfn) (mife_pp_t, int, int *, int *),
    void (*setfn)   (mife_pp_t, mife_mat_clr_t, void *),
    int (*parsefn)  (mife_pp_t, f2_matrix)
    );
void mife_setup(const_mmap_vtable mmap, mife_pp_t pp, mife_sk_t sk, int L, int lambda,
    aes_randstate_t randstate);
void mife_encrypt(const_mmap_vtable mmap, mife_ciphertext_t ct, void *message, mife_pp_t pp,
    mife_sk_t sk, aes_randstate_t randstate);
int mife_evaluate(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t *cts);

/* memory-efficient MIFE interface */
void mife_encrypt_setup(mife_pp_t pp, fmpz_t uid, void *message,
    mife_mat_clr_t out_clr, int ****out_partitions);
void mife_encrypt_single(const_mmap_vtable mmap, mife_pp_t pp, mife_sk_t sk, aes_randstate_t randstate,
    int global_index, mife_mat_clr_t clr, int ***partitions,
    mmap_enc_mat_t out_ct);
void mife_encrypt_clear(mife_pp_t pp, mife_mat_clr_t clr, int ***out_partitions);
f2_matrix mife_zt_all(const_mmap_vtable mmap, const mife_pp_t pp, mmap_enc_mat_t m);

#endif /* _MIFE_H_ */
