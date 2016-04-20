#include "mife_glue_clt.h"

#include <gmp.h>
#include "mife_glue.h"

void clt_pp_clear_wrapper (mmap_pp *pp);
void clt_pp_read_wrapper  (mmap_pp *const pp, FILE *const fp);
void clt_pp_save_wrapper  (const mmap_pp *const pp, FILE *const fp);

const mmap_pp_vtable clt13_pp_vtable =
  { .clear  = clt_pp_clear_wrapper
  , .fread  = clt_pp_read_wrapper
  , .fwrite = clt_pp_save_wrapper
  , .size   = sizeof(mmap_pp)
  };

void clt_state_init_wrapper  (mmap_sk *const sk, size_t lambda, size_t kappa, size_t gamma, aes_randstate_t randstate);
void clt_state_clear_wrapper (mmap_sk *const sk);
void clt_state_read_wrapper  (mmap_sk *const sk, FILE *const fp);
void clt_state_save_wrapper  (const mmap_sk *const sk, FILE *const fp);
void clt_state_get_modulus   (const mmap_sk *const sk, fmpz_t p_out);
const mmap_pp *const clt_pp_init_wrapper (const mmap_sk *const sk);

const mmap_sk_vtable clt13_sk_vtable =
  { .init   = clt_state_init_wrapper
  , .clear  = clt_state_clear_wrapper
  , .fread  = clt_state_read_wrapper
  , .fwrite = clt_state_save_wrapper
  , .pp     = clt_pp_init_wrapper
  , .size   = sizeof(clt_state)
  , .plaintext_field = clt_state_get_modulus
  };

void clt_enc_init_wrapper    (mmap_enc *const enc, const mmap_pp *const pp);
void clt_enc_clear_wrapper   (mmap_enc *const enc);
void clt_enc_fread_wrapper   (mmap_enc *const enc, FILE *const fp);
void clt_enc_fwrite_wrapper  (const mmap_enc *const enc, FILE *const fp);
void clt_enc_set_wrapper     (mmap_enc *const dest, const mmap_enc *const src);
void clt_enc_add_wrapper     (mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b);
void clt_enc_mul_wrapper     (mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b);
bool clt_enc_is_zero_wrapper (const mmap_enc *const enc, const mmap_pp *const pp);

void clt_encode_wrapper (
    mmap_enc *const enc,
    const mmap_sk *const sk,
    const fmpz_t plaintext,
    int *group,
    aes_randstate_t randstate
);

const mmap_enc_vtable clt13_enc_vtable =
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

const mmap_vtable clt13_vtable =
  { .pp  = &clt13_pp_vtable
  , .sk  = &clt13_sk_vtable
  , .enc = &clt13_enc_vtable
  };

////////////////////////////////////////////////////////////////////////////////
// implementation

void clt_pp_clear_wrapper (mmap_pp *pp)
{
    clt_pp_clear(&(pp->clt13_self));
}

void clt_pp_read_wrapper (mmap_pp *const pp, FILE *const fp)
{
    clt_pp_fread(fp, &(pp->clt13_self));
}

void clt_pp_save_wrapper (const mmap_pp *const pp, FILE *const fp)
{
    clt_pp_fsave(fp, &(pp->clt13_self));
}

void clt_state_init_wrapper (mmap_sk *const sk, size_t lambda, size_t kappa,
                             size_t gamma, aes_randstate_t rng)
{
    int *pows = malloc(gamma * sizeof(int));
    for (int i = 0; i < gamma; i++) {
        pows[i] = 1;
    }
    clt_state_init(&(sk->clt13_self), kappa, lambda, gamma, pows,
                   CLT_FLAG_VERBOSE & CLT_FLAG_DEFAULT & CLT_FLAG_OPT_PARALLEL_ENCODE, rng);
    free(pows);
}

void clt_state_clear_wrapper (mmap_sk *const sk)
{
    clt_state_clear(&(sk->clt13_self));
}

void clt_state_read_wrapper (mmap_sk *const sk, FILE *const fp)
{
    clt_state_fread(fp, &(sk->clt13_self));
}

void clt_state_save_wrapper (const mmap_sk *const sk, FILE *const fp)
{
    clt_state_fsave(fp, &(sk->clt13_self));
}

void clt_state_get_modulus (const mmap_sk *const sk, fmpz_t p_out)
{
    fmpz_set_mpz(p_out, sk->clt13_self.gs[0]);
}

const mmap_pp *const clt_pp_init_wrapper (const mmap_sk *const sk)
{
    mmap_pp *pp = malloc(sizeof(mmap_pp));
    clt_pp_init(&(pp->clt13_self), &(sk->clt13_self));
    return pp;
}

void clt_enc_init_wrapper (mmap_enc *const enc, const mmap_pp *const pp)
{
    mpz_init(enc->clt13_self);
}

void clt_enc_clear_wrapper (mmap_enc *const enc)
{
    mpz_clear(enc->clt13_self);
}

void clt_enc_fread_wrapper (mmap_enc *enc, FILE *const fp)
{
    mpz_inp_raw(enc->clt13_self, fp);
}

void clt_enc_fwrite_wrapper (const mmap_enc *const enc, FILE *const fp)
{
    clt_elem_fsave(fp, enc->clt13_self);
}

void clt_enc_set_wrapper (mmap_enc *const dest, const mmap_enc *const src)
{
    mpz_set(dest->clt13_self, src->clt13_self);
}

void clt_enc_add_wrapper (mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b)
{
    mpz_add(dest->clt13_self, a->clt13_self, b->clt13_self);
    mpz_mod(dest->clt13_self, dest->clt13_self, pp->clt13_self.x0);
}

void clt_enc_mul_wrapper (mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b)
{
    mpz_mul(dest->clt13_self, a->clt13_self, b->clt13_self);
    mpz_mod(dest->clt13_self, dest->clt13_self, pp->clt13_self.x0);
}

bool clt_enc_is_zero_wrapper (const mmap_enc *const enc, const mmap_pp *const pp)
{
    return clt_is_zero(&(pp->clt13_self), enc->clt13_self);
}

void
clt_encode_wrapper (mmap_enc *const enc, const mmap_sk *const sk,
                    const fmpz_t plaintext, int *group, aes_randstate_t rng)
{
    mpz_t ins[1];
    mpz_init(ins[0]);
    fmpz_get_mpz(ins[0], plaintext);
    clt_encode(enc->clt13_self, &(sk->clt13_self), 1, ins, group, rng);
    mpz_clear(ins[0]);
}
