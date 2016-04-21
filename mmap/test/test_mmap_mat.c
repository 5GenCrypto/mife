#include <mmap.h>
#include <mmap_clt.h>
#include <mmap_gghlite.h>
#include <flint/fmpz.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include "utils.h"

ulong nzs = 2;
ulong lambda = 10;
ulong kappa = 2;

static void encode(const mmap_vtable *vtable, mmap_sk *sk,
                   mmap_enc_mat_t out, fmpz_mat_t in, int idx, int nrows,
                   int ncols, aes_randstate_t rand)
{
    int *pows;

    pows = calloc(nzs, sizeof(int));
    /* for (ulong i = 0; i < nzs; i++) pows[i] = 1; */
    pows[idx] = 1;

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            vtable->enc->encode(out->m[i][j], sk, 1, fmpz_mat_entry(in, i, j),
                                pows, rand);
        }
    }
    free(pows);
}

static int test(const mmap_vtable *vtable)
{
    int ok = 1;
    mmap_sk sk;
    const mmap_pp *pp;

    aes_randstate_t rng;
    aes_randinit(rng);

    /* fmpz_t mod; */
    /* fmpz_init(mod); */

    vtable->sk->init(&sk, lambda, kappa, nzs, rng);
    /* vtable->sk->plaintext_field(&sk, mod); */
    pp = vtable->sk->pp(&sk);

    fmpz_mat_t zero_1, one_1, zero_2, one_2, res;
    fmpz_mat_init(zero_1, 1, 2);
    fmpz_mat_init(one_1,  1, 2);
    fmpz_mat_init(zero_2, 2, 2);
    fmpz_mat_init(one_2,  2, 2);
    fmpz_mat_init(res,    1, 2);

    fmpz_set_ui(fmpz_mat_entry(zero_1, 0, 0), 1);
    fmpz_set_ui(fmpz_mat_entry(zero_1, 0, 1), 0);

    fmpz_set_ui(fmpz_mat_entry(one_1, 0, 0), 1);
    fmpz_set_ui(fmpz_mat_entry(one_1, 0, 1), 1);

    fmpz_set_ui(fmpz_mat_entry(zero_2, 0, 0), 1);
    fmpz_set_ui(fmpz_mat_entry(zero_2, 0, 1), 0);
    fmpz_set_ui(fmpz_mat_entry(zero_2, 1, 0), 0);
    fmpz_set_ui(fmpz_mat_entry(zero_2, 1, 1), 0);

    fmpz_set_ui(fmpz_mat_entry(one_2, 0, 0), 1);
    fmpz_set_ui(fmpz_mat_entry(one_2, 0, 1), 0);
    fmpz_set_ui(fmpz_mat_entry(one_2, 1, 0), 0);
    fmpz_set_ui(fmpz_mat_entry(one_2, 1, 1), 1);


    mmap_enc_mat_t zero_enc_1, one_enc_1, zero_enc_2, one_enc_2, result;
    mmap_enc_mat_init(vtable, pp, zero_enc_1, 1, 2);
    mmap_enc_mat_init(vtable, pp, one_enc_1,  1, 2);
    mmap_enc_mat_init(vtable, pp, zero_enc_2, 2, 2);
    mmap_enc_mat_init(vtable, pp, one_enc_2,  2, 2);
    mmap_enc_mat_init(vtable, pp, result,     1, 2);

    encode(vtable, &sk, zero_enc_1, zero_1, 0, 1, 2, rng);
    encode(vtable, &sk, one_enc_1,  one_1,  0, 1, 2, rng);
    encode(vtable, &sk, zero_enc_2, zero_2, 1, 2, 2, rng);
    encode(vtable, &sk, one_enc_2,  one_2,  1, 2, 2, rng);

    mmap_enc_mat_mul(vtable, pp, result, zero_enc_1, zero_enc_2);
    ok &= expect("[1 0] * [1 0][0 0]", 1, vtable->enc->is_zero(result->m[0][1], pp));
    mmap_enc_mat_mul(vtable, pp, result, zero_enc_1, one_enc_2);
    ok &= expect("[1 0] * [1 0][0 1]", 1, vtable->enc->is_zero(result->m[0][1], pp));
    mmap_enc_mat_mul(vtable, pp, result, one_enc_1, zero_enc_2);
    ok &= expect("[1 1] * [1 0][0 0]", 1, vtable->enc->is_zero(result->m[0][1], pp));
    mmap_enc_mat_mul(vtable, pp, result, one_enc_1, one_enc_2);
    ok &= expect("[1 1] * [1 0][0 1]", 0, vtable->enc->is_zero(result->m[0][1], pp));

    return !ok;
}

int main(void)
{
    printf("* CLT13\n");
    if (test(&clt_vtable) == 1)
        return 1;
    printf("* GGHLite\n");
    if (test(&gghlite_vtable) == 1)
        return 1;
    return 0;
}
