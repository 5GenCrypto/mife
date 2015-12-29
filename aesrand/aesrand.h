#ifndef _AESRAND_H_
#define _AESRAND_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <openssl/evp.h>

typedef struct {
  unsigned long ctr;
  EVP_CIPHER_CTX *ctx;
  unsigned char *key;
  unsigned char *iv;
} aes_randstate_t[1];

void aesrand_init(aes_randstate_t state);
void aesrand_clear(aes_randstate_t state);

void mpz_urandomb_aes(mpz_t rop, aes_randstate_t state, mp_bitcnt_t n);
void mpz_urandomm_aes(mpz_t rop, aes_randstate_t state, const mpz_t n);


#endif /* _AESRAND_H_ */
