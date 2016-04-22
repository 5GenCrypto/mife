#include "mmap.h"
#include "mmap_clt.h"

#include <gmp.h>

static void clt_pp_clear_wrapper (mmap_pp *pp);
static void clt_pp_read_wrapper  (mmap_pp *const pp, FILE *const fp);
static void clt_pp_save_wrapper  (const mmap_pp *const pp, FILE *const fp);

static const mmap_pp_vtable clt_pp_vtable =
  { .clear  = clt_pp_clear_wrapper
  , .fread  = clt_pp_read_wrapper
  , .fwrite = clt_pp_save_wrapper
  , .size   = sizeof(mmap_pp)
  };

static void clt_state_init_wrapper  (mmap_sk *const sk, size_t lambda, size_t kappa, size_t gamma, aes_randstate_t randstate, bool verbose);
static void clt_state_clear_wrapper (mmap_sk *const sk);
static void clt_state_read_wrapper  (mmap_sk *const sk, FILE *const fp);
static void clt_state_save_wrapper  (const mmap_sk *const sk, FILE *const fp);
static void clt_state_get_modulus   (const mmap_sk *const sk, fmpz_t p_out);
static const mmap_pp *const clt_pp_init_wrapper (const mmap_sk *const sk);

static const mmap_sk_vtable clt_sk_vtable =
  { .init   = clt_state_init_wrapper
  , .clear  = clt_state_clear_wrapper
  , .fread  = clt_state_read_wrapper
  , .fwrite = clt_state_save_wrapper
  , .pp     = clt_pp_init_wrapper
  , .size   = sizeof(clt_state)
  , .plaintext_field = clt_state_get_modulus
  };

static void clt_enc_init_wrapper    (mmap_enc *const enc, const mmap_pp *const pp);
static void clt_enc_clear_wrapper   (mmap_enc *const enc);
static void clt_enc_fread_wrapper   (mmap_enc *const enc, FILE *const fp);
static void clt_enc_fwrite_wrapper  (const mmap_enc *const enc, FILE *const fp);
static void clt_enc_set_wrapper     (mmap_enc *const dest, const mmap_enc *const src);
static void clt_enc_add_wrapper     (mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b);
static void clt_enc_mul_wrapper     (mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b);
static bool clt_enc_is_zero_wrapper (const mmap_enc *const enc, const mmap_pp *const pp);
static void clt_encode_wrapper (mmap_enc *const enc, const mmap_sk *const sk, int n, const fmpz_t *plaintext, int *group, aes_randstate_t rng);

static const mmap_enc_vtable clt_enc_vtable =
  { .init    = clt_enc_init_wrapper
  , .clear   = clt_enc_clear_wrapper
  , .fread   = clt_enc_fread_wrapper
  , .fwrite  = clt_enc_fwrite_wrapper
  , .set     = clt_enc_set_wrapper
  , .add     = clt_enc_add_wrapper
  , .mul     = clt_enc_mul_wrapper
  , .is_zero = clt_enc_is_zero_wrapper
  , .encode  = clt_encode_wrapper
  , .size    = sizeof(mmap_enc)
  };

const mmap_vtable clt_vtable =
  { .pp  = &clt_pp_vtable
  , .sk  = &clt_sk_vtable
  , .enc = &clt_enc_vtable
  };

////////////////////////////////////////////////////////////////////////////////
// implementation

static void clt_pp_clear_wrapper (mmap_pp *pp)
{
    clt_pp_clear(&(pp->clt_self));
}

static void clt_pp_read_wrapper (mmap_pp *const pp, FILE *const fp)
{
    clt_pp_fread(fp, &(pp->clt_self));
}

static void clt_pp_save_wrapper (const mmap_pp *const pp, FILE *const fp)
{
    clt_pp_fsave(fp, &(pp->clt_self));
}

static void clt_state_init_wrapper (mmap_sk *const sk, size_t lambda, size_t kappa,
                                    size_t gamma, aes_randstate_t rng, bool verbose)
{
    int flags = CLT_FLAG_DEFAULT | CLT_FLAG_OPT_PARALLEL_ENCODE;
    int *pows = malloc(gamma * sizeof(int));
    for (size_t i = 0; i < gamma; i++) {
        pows[i] = 1;
    }
    if (verbose)
        flags |= CLT_FLAG_VERBOSE;
    clt_state_init(&(sk->clt_self), kappa, lambda, gamma, pows, flags, rng);
    free(pows);
}

static void clt_state_clear_wrapper (mmap_sk *const sk)
{
    clt_state_clear(&(sk->clt_self));
}

static void clt_state_read_wrapper (mmap_sk *const sk, FILE *const fp)
{
    clt_state_fread(fp, &(sk->clt_self));
}

static void clt_state_save_wrapper (const mmap_sk *const sk, FILE *const fp)
{
    clt_state_fsave(fp, &(sk->clt_self));
}

static void clt_state_get_modulus (const mmap_sk *const sk, fmpz_t p_out)
{
    fmpz_set_mpz(p_out, sk->clt_self.gs[0]);
}

static const mmap_pp *const clt_pp_init_wrapper (const mmap_sk *const sk)
{
    mmap_pp *pp = malloc(sizeof(mmap_pp));
    clt_pp_init(&(pp->clt_self), &(sk->clt_self));
    return pp;
}

static void clt_enc_init_wrapper (mmap_enc *const enc, const mmap_pp *const pp)
{
    mpz_init(enc->clt_self);
}

static void clt_enc_clear_wrapper (mmap_enc *const enc)
{
    mpz_clear(enc->clt_self);
}

static void clt_enc_fread_wrapper (mmap_enc *enc, FILE *const fp)
{
    mpz_init(enc->clt_self);
    mpz_inp_raw(enc->clt_self, fp);
}

static void clt_enc_fwrite_wrapper (const mmap_enc *const enc, FILE *const fp)
{
    mpz_out_raw(fp, enc->clt_self);
}

static void clt_enc_set_wrapper (mmap_enc *const dest, const mmap_enc *const src)
{
    mpz_set(dest->clt_self, src->clt_self);
}

static void clt_enc_add_wrapper (mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b)
{
    mpz_add(dest->clt_self, a->clt_self, b->clt_self);
    mpz_mod(dest->clt_self, dest->clt_self, pp->clt_self.x0);
}

static void clt_enc_mul_wrapper (mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b)
{
    mpz_mul(dest->clt_self, a->clt_self, b->clt_self);
    mpz_mod(dest->clt_self, dest->clt_self, pp->clt_self.x0);
}

static bool clt_enc_is_zero_wrapper (const mmap_enc *const enc, const mmap_pp *const pp)
{
    return clt_is_zero(&(pp->clt_self), enc->clt_self);
}

static void
clt_encode_wrapper (mmap_enc *const enc, const mmap_sk *const sk, int n,
                    const fmpz_t *plaintext, int *group, aes_randstate_t rng)
{
    mpz_t *ins;

    ins = calloc(n, sizeof(mpz_t));
    for (int i = 0; i < n; ++i) {
        mpz_init(ins[i]);
        fmpz_get_mpz(ins[i], plaintext[i]);
    }
    clt_encode(enc->clt_self, &(sk->clt_self), n, ins, group, rng);
    for (int i = 0; i < n; ++i) {
        mpz_clear(ins[i]);
    }
    free(ins);
}
