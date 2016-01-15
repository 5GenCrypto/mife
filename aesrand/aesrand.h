#ifndef _AESRAND_H_
#define _AESRAND_H_

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "../flint/fmpz_poly.h"
#include "../flint/fmpz_mod_poly.h"
#include <openssl/evp.h>
#include <openssl/sha.h>

struct _aes_randstate_struct {
  char aes_init;
  unsigned long ctr;
  EVP_CIPHER_CTX *ctx;
  unsigned char key[SHA256_DIGEST_LENGTH];
  unsigned char *iv;
};

typedef struct _aes_randstate_struct aes_randstate_t[1];

void aes_randinit(aes_randstate_t state);
void aes_randinit_seed(aes_randstate_t state, char *seed, char *additional);
void aes_randclear(aes_randstate_t state);

void fmpz_mod_poly_randtest_aes(fmpz_mod_poly_t f, aes_randstate_t state,
    slong len);
void fmpz_randm_aes(fmpz_t f, aes_randstate_t state, const fmpz_t m);
void mpfr_urandomb_aes(mpfr_t rop, aes_randstate_t state);
void mpz_urandomb_aes(mpz_t rop, aes_randstate_t state, mp_bitcnt_t n);
void mpz_urandomm_aes(mpz_t rop, aes_randstate_t state, const mpz_t n);


#endif /* _AESRAND_H_ */
