#include "aesrand.h"

void aesrand_init(aes_randstate_t state) {
  state->ctr = 0;
  state->ctx = EVP_CIPHER_CTX_new();
  state->key = "12345678901234567890123456789012";
  state->iv = "000000000000";
  EVP_EncryptInit_ex (state->ctx, EVP_aes_256_gcm(), NULL, state->key, state->iv);
}

void aesrand_clear(aes_randstate_t state) {
  EVP_CIPHER_CTX_cleanup(state->ctx);
}

void mpz_urandomm_aes(mpz_t rop, aes_randstate_t state, const mpz_t n) {
  unsigned long size = mpz_sizeinbase(n, 2);

  while(1) {
    mpz_urandomb_aes(rop, state, size);
    if(mpz_cmp(rop, n) < 0) {
      break;
    }
  }
}

void mpz_urandomb_aes(mpz_t rop, aes_randstate_t state, mp_bitcnt_t n) {
  int nb = n/8+1; // number of bytes

  unsigned char *in = malloc(nb);
  for(int i = 0; i < nb; i++) {
    in[i] = state->ctr++;
  }

  unsigned char output[256];
  int outlen, bytelen;
  EVP_EncryptUpdate(state->ctx, output, &outlen, in, nb);

  if(outlen != nb) {
    printf("ERROR: with randomness being generated.\n");
  }
  
  bytelen = outlen;
  int true_len = bytelen + 4;

  unsigned char *buf = malloc(true_len);
  memset(buf, 0, true_len);
  memcpy(buf+4, output, outlen);
  buf[4] >>= ((nb*8) - (unsigned int) n);
    
  for(int i = 3; i >= 0; i--) {
    buf[i] = (unsigned char) (bytelen % (1 << 8));
    bytelen /= (1 << 8);
  }

  /*
  char *printbuf = malloc(3 * true_len + 1);
  char *buf_ptr = printbuf;
  memset(printbuf, 0, 3 * true_len + 1);
  for(int i = 0; i < true_len; i++) {
    buf_ptr += sprintf(buf_ptr, "%02x ", buf[i]);
  }
  printf("%s\n", printbuf);
  */
  
  // this generates a random n-bit number.
  FILE *fp = fmemopen(buf, outlen+4, "rb");

  mpz_inp_raw(rop, fp);
  //gmp_printf("%Zd\n", p);
   
}
