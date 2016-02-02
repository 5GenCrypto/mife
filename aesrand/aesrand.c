#include "aesrand.h"

#define AES_ALGORITHM EVP_aes_256_ctr()

int X = 0;
int Y = 0;

int VERBOSE = 0;

void aes_randinit(aes_randstate_t state) {
  char *default_seed = "12345678901234567890123456789012";
  aes_randinit_seed(state, default_seed, NULL);
}

void aes_randinit_seed(aes_randstate_t state, char *seed, char *additional) {
  if(additional == NULL) {
    additional = "";
  }
  aes_randinit_seedn(state, seed, strlen(seed), additional, strlen(additional));
}

void aes_randinit_seedn(aes_randstate_t state, char *seed, size_t seed_len, char *additional, size_t additional_len) {
  state->aes_init = 1;
  state->ctr = 0;
  state->iv = malloc(EVP_CIPHER_iv_length(AES_ALGORITHM));
  memset(state->iv, 0, EVP_CIPHER_iv_length(AES_ALGORITHM));

  SHA256_CTX sha256;
  SHA256_Init(&sha256);
  SHA256_Update(&sha256, seed, seed_len);
  SHA256_Update(&sha256, additional, additional_len);
  SHA256_Final(state->key, &sha256);
}

void aes_randclear(aes_randstate_t state) {
  free(state->iv);
}

void fmpz_mod_poly_randtest_aes(fmpz_mod_poly_t f, aes_randstate_t state,
    slong len) {
    slong i;

    fmpz_mod_poly_fit_length(f, len);

    for (i = 0; i < len; i++) {
      fmpz_randm_aes(f->coeffs + i, state, &(f->p));
    }

    _fmpz_mod_poly_set_length(f, len);
    _fmpz_mod_poly_normalise(f);
}

void fmpz_randm_aes(fmpz_t f, aes_randstate_t state, const fmpz_t m) {
  mpz_t x, rop;
  mpz_init(x);
  mpz_init(rop);
  fmpz_get_mpz(x, m);
  mpz_urandomm_aes(rop, state, x);
  fmpz_set_mpz(f, rop);
  mpz_clear(x);
  mpz_clear(rop);
}

void mpfr_urandomb_aes(mpfr_t rop, aes_randstate_t state) {
  //printf("calling the mpfr randomness function which is weird...\n");
  unsigned long size = mpfr_get_prec(rop);

  mpfr_t num, denom;
  mpz_t mpz_num, mpz_denom;
  mpz_init(mpz_num);
  mpz_init(mpz_denom);

  //printf("size: %lu\n", size);
  //VERBOSE = 1;
  mpz_urandomb_aes(mpz_num, state, size);
  //VERBOSE = 0;
  mpfr_init_set_z(num, mpz_num, MPFR_RNDN);

  //mpfr_printf("mpz_print for mpz_num: %Zd\n", mpz_num);
  //mpfr_printf("mpfr_printf for num: %Zd\n", num);
  mpz_ui_pow_ui(mpz_denom, 2, size);
  mpfr_init_set_z(denom, mpz_denom, MPFR_RNDN);

  mpfr_div(rop, num, denom, MPFR_RNDN);
  mpz_clear(mpz_num);
  mpz_clear(mpz_denom);
  mpfr_clear(num);
  mpfr_clear(denom);
  //mpfr_printf("mpfr_printf: %.128Rf\n", rop);
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


void gen_random_bytes(unsigned char *buf, mp_bitcnt_t n) {
  unsigned long p = 0;

  for( ; p < n; p++) {
    buf[p] = (unsigned char) (rand() % 256);
  }
  /*
  while(p < n) {
    unsigned long num_get = (n-p > 256) ? 256 : n-p;
    syscall(SYS_getrandom, buf+p, num_get, 0);
    p += num_get;
  }*/

}

void fmpz_randbits_aes(fmpz_t out, aes_randstate_t state, mp_bitcnt_t bits) {
  mpz_t rop;
  mpz_init(rop);
  mpz_urandomb_aes(rop, state, bits);
  fmpz_set_mpz(out, rop);
  mpz_clear(rop);
}

void mpz_urandomb_aes(mpz_t rop, aes_randstate_t state, mp_bitcnt_t n) {
  mp_bitcnt_t nb = n/8+1; // number of bytes
  int TRY_MAX = 1000000000; // number of times to attempt getting good randomness
  //printf("nb: %lu\n", nb);

  // update the internal counter, works at most 2^64 times
  memcpy(state->iv, &state->ctr, sizeof(state->ctr)); 

  state->ctx = EVP_CIPHER_CTX_new();
  EVP_EncryptInit_ex (state->ctx, AES_ALGORITHM, NULL, state->key, state->iv);

  unsigned char *output = malloc(2 * (nb + EVP_MAX_IV_LENGTH));
  mp_bitcnt_t outlen = 0;
  int halt_count = 0;

  int in_size = nb;
  unsigned char in[in_size];
  memset(in, 0, in_size);

  while(outlen < nb && halt_count < TRY_MAX) {
    int buflen = 0;
    EVP_EncryptUpdate(state->ctx, output+outlen, &buflen, in, in_size);
    state->ctr++;
    outlen += buflen;
    halt_count++;
  }
  int final_len = 0;
  EVP_EncryptFinal(state->ctx, output+outlen, &final_len);
  outlen += final_len;

  if(halt_count == TRY_MAX) {
    printf("ERROR: randomness being generated is not perfect.\n");
  }

  if(outlen > nb) {
    outlen = nb; // we will only use nb bytes
  }

/*  
  outlen = nb;
  gen_random_bytes(output, nb);
*/
  
  mp_bitcnt_t true_len = outlen + 4;
  mp_bitcnt_t bytelen = outlen;

  unsigned char *buf = malloc(true_len);
  memset(buf, 0, true_len);
  memcpy(buf+4, output, outlen);
  buf[4] >>= ((outlen*8) - (unsigned int) n);

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
  FILE *fp = fmemopen(buf, true_len, "rb");
  if(!fp) {
    printf("Error in generating randomness.\n");
  }

  if(mpz_inp_raw(rop, fp) == 0) {
    printf("Error in parsing randomness.\n");
  }

  //mpfr_printf("rop: %Zd\n", rop);

  fclose(fp); 
  free(output);
  free(buf); 

  EVP_CIPHER_CTX_cleanup(state->ctx);
  free(state->ctx);

}

