#include "aesrand.h"

int X = 0;
int Y = 0;

void aes_randinit(aes_randstate_t state) {
  char *default_seed = "12345678901234567890123456789012";
  aes_randinit_seed(state, default_seed);
}

void aes_randinit_seed(aes_randstate_t state, char *seed) {
  state->aes_init = 1;
  state->ctr = 0;
  state->key = seed;
  state->iv = "000000000000";
  state->ctx = EVP_CIPHER_CTX_new();
  EVP_EncryptInit_ex (state->ctx, EVP_aes_128_gcm(), NULL, state->key,
      state->iv);
}

void aes_randclear(aes_randstate_t state) {
  EVP_CIPHER_CTX_cleanup(state->ctx);
  free(state->ctx);
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
  unsigned long size = mpfr_get_prec(rop);

  mpfr_t num, denom;
  mpz_t mpz_num, mpz_denom;
  mpz_init(mpz_num);
  mpz_init(mpz_denom);

  mpz_urandomb_aes(mpz_num, state, size);
  mpfr_init_set_z(num, mpz_num, MPFR_RNDN);

  mpz_ui_pow_ui(mpz_denom, 2, size);
  mpfr_init_set_z(denom, mpz_denom, MPFR_RNDN);

  mpfr_div(rop, num, denom, MPFR_RNDN);
  mpz_clear(mpz_num);
  mpz_clear(mpz_denom);
  mpfr_clear(num);
  mpfr_clear(denom);
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

  unsigned char *output = malloc(nb);
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
  if(!fp) {
    printf("Error in generating randomness.\n");
  }

  mpz_inp_raw(rop, fp);
  //gmp_printf("%Zd\n", p);

  fclose(fp); 
  free(in);
  free(output);
  free(buf); 
}

