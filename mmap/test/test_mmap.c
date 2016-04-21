#include <mmap.h>
#include <mmap_clt.h>
#include <mmap_gghlite.h>
#include <flint/fmpz.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

int expect(char * desc, int expected, int recieved);

static int test(const mmap_vtable *mmap)
{
    srand(time(NULL));

    ulong nzs     = 10;
    ulong lambda  = 10;
    ulong kappa   = 2;

    aes_randstate_t rng;
    aes_randinit(rng);

    int pows [nzs];
    for (ulong i = 0; i < nzs; i++) pows[i] = 1;

    FILE *sk_f = tmpfile();
    if (sk_f == NULL) {
        fprintf(stderr, "Couldn't open test.map!\n");
        exit(1);
    }

    FILE *pp_f = tmpfile();
    if (pp_f == NULL) {
        fprintf(stderr, "Couldn't open test.pp!\n");
        exit(1);
    }

    // test initialization & serialization
    mmap_sk *sk_ = malloc(mmap->sk->size);
    mmap->sk->init(sk_, lambda, kappa, nzs, rng);
    mmap->sk->fwrite(sk_, sk_f);
    rewind(sk_f);
    mmap->sk->clear(sk_);
    free(sk_);
    mmap_sk *sk = malloc(mmap->sk->size);
    mmap->sk->fread(sk, sk_f);

    mmap_pp *pp_ = mmap->sk->pp(sk);
    mmap->pp->fwrite(pp_, pp_f);
    rewind(pp_f);
    mmap_pp *pp = malloc(mmap->pp->size);
    mmap->pp->fread(pp, pp_f);

    fmpz_t x [1];
    fmpz_init_set_ui(x[0], 0);
    while (fmpz_cmp_ui(x[0], 0) <= 0) {
        fmpz_set_ui(x[0], rand());
        fmpz_t mod;
        fmpz_init(mod);
        mmap->sk->plaintext_field(sk, mod);
        fmpz_mod(x[0], x[0], mod);
    }
    printf("x = ");
    fmpz_print(x[0]);
    puts("");

    fmpz_t zero [1];
    fmpz_init_set_ui(zero[0], 0);

    fmpz_t one [1];
    fmpz_init_set_ui(one[0], 1);

    int top_level [nzs];
    for (ulong i = 0; i < nzs; i++) {
        top_level[i] = 1;
    }

    mmap_enc x0, x1, xp;
    mmap->enc->init(&x0, pp);
    mmap->enc->init(&x1, pp);
    mmap->enc->init(&xp, pp);
    mmap->enc->encode(&x0, sk, 1, zero, top_level, rng);
    mmap->enc->encode(&x1, sk, 1, zero, top_level, rng);
    mmap->enc->add(&xp, pp, &x0, &x1);
    int ok = expect("is_zero(0 + 0)", 1, mmap->enc->is_zero(&xp, pp));

    mmap->enc->encode(&x0, sk, 1, zero, top_level, rng);
    mmap->enc->encode(&x1, sk, 1, one,  top_level, rng);
    mmap->enc->add(&xp, pp, &x0, &x1);
    ok &= expect("is_zero(0 + 1)", 0, mmap->enc->is_zero(&xp, pp));

    mmap->enc->encode(&x0, sk, 1, zero, top_level, rng);
    mmap->enc->encode(&x1, sk, 1, x,    top_level, rng);
    mmap->enc->add(&xp, pp, &x0, &x1);
    ok &= expect("is_zero(0 + x)", 0, mmap->enc->is_zero(&xp, pp));

    int ix0 [nzs];
    int ix1 [nzs];
    for (ulong i = 0; i < nzs; i++) {
        if (i < nzs / 2) {
            ix0[i] = 1;
            ix1[i] = 0;
        } else {
            ix0[i] = 0;
            ix1[i] = 1;
        }
    }
    mmap->enc->encode(&x0, sk, 1, x   , ix0, rng);
    mmap->enc->encode(&x1, sk, 1, zero, ix1, rng);
    mmap->enc->mul(&xp, pp, &x0, &x1);
    ok &= expect("is_zero(x * 0)", 1, mmap->enc->is_zero(&xp, pp));

    mmap->enc->encode(&x0, sk, 1, x   , ix0, rng);
    mmap->enc->encode(&x1, sk, 1, one, ix1, rng);
    mmap->enc->mul(&xp, pp, &x0, &x1);
    ok &= expect("is_zero(x * 1)", 0, mmap->enc->is_zero(&xp, pp));

    mmap->enc->encode(&x0, sk, 1, x   , ix0, rng);
    mmap->enc->encode(&x1, sk, 1, x, ix1, rng);
    mmap->enc->mul(&xp, pp, &x0, &x1);
    ok &= expect("is_zero(x * x)", 0, mmap->enc->is_zero(&xp, pp));

    mmap->enc->clear(&x0);
    mmap->enc->clear(&x1);
    mmap->enc->clear(&xp);
    free(sk);
    return !ok;
}

int main(void)
{
    printf("* CLT13\n");
    if (test(&clt13_vtable) == 1)
        return 1;
    printf("* GGHLite\n");
    if (test(&gghlite_vtable) == 1)
        return 1;
    return 0;
}

int expect(char * desc, int expected, int recieved)
{
    if (expected != recieved) {
        printf("\033[1;41m");
    }
    printf("%s = %d", desc, recieved);
    if (expected != recieved) {
        printf("\033[0m");
    }
    puts("");
    return expected == recieved;
}
