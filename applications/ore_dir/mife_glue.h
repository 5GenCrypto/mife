#ifndef _MIFE_GLUE_H
#define _MIFE_GLUE_H

#include <gghlite/gghlite-defs.h>
#include <mife/mife_defs.h>

struct mmap_pp  { gghlite_params_t self; };
struct mmap_sk  { gghlite_sk_t     self; };
struct mmap_enc { gghlite_enc_t    self; };
const mmap_vtable gghlite_vtable;

void gghlite_params_clear_read_wrapper(mmap_pp *pp);
void fread_gghlite_params_wrapper(mmap_pp *const pp, FILE *const fp);
void fwrite_gghlite_params_wrapper(const mmap_pp *const pp, FILE *const fp);

void gghlite_jigsaw_init_gamma_wrapper(mmap_sk *const sk, size_t lambda, size_t kappa, size_t gamma, aes_randstate_t randstate);
void gghlite_sk_clear_wrapper(mmap_sk *const sk);
void fread_gghlite_sk_wrapper(mmap_sk *const sk, FILE *const fp);
void fwrite_gghlite_sk_wrapper(const mmap_sk *const sk, FILE *const fp);
const mmap_pp *const gghlite_sk_to_pp(const mmap_sk *const sk);
void fmpz_poly_oz_ideal_norm_wrapper(const mmap_sk *const sk, fmpz_t p_out);

void gghlite_enc_init_wrapper(mmap_enc *const enc, const mmap_pp *const pp);
void gghlite_enc_clear_wrapper(mmap_enc *const enc);
void gghlite_enc_fread_raw_wrapper(mmap_enc *const enc, FILE *const fp);
void gghlite_enc_fprint_raw_wrapper(const mmap_enc *const enc, FILE *const fp);
void gghlite_enc_set_wrapper(mmap_enc *const dest, const mmap_enc *const src);
void gghlite_enc_add_wrapper(mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b);
void gghlite_enc_mul_wrapper(mmap_enc *const dest, const mmap_pp *const pp, const mmap_enc *const a, const mmap_enc *const b);
bool gghlite_enc_is_zero_wrapper(const mmap_enc *const enc, const mmap_pp *const pp);
void gghlite_enc_set_gghlite_clr_wrapper(mmap_enc *const enc, const mmap_sk *const sk, const fmpz_t plaintext, int *group, aes_randstate_t randstate);

void gghlite_params_clear_read(gghlite_params_t self);

void fread_gghlite_params(FILE *fp, gghlite_params_t params);
void fwrite_gghlite_params(FILE *fp, const gghlite_params_t params);
void fread_gghlite_sk(FILE *fp, gghlite_sk_t self);
void fwrite_gghlite_sk(FILE *fp, const gghlite_sk_t self);

/* functions dealing with file reading and writing for encodings */
#define gghlite_enc_fprint fmpz_mod_poly_fprint
#define gghlite_enc_fread fmpz_mod_poly_fread
#define gghlite_enc_fprint_raw fmpz_mod_poly_fprint_raw
#define gghlite_enc_fread_raw fmpz_mod_poly_fread_raw
int fmpz_mod_poly_fprint_raw(FILE * file, const fmpz_mod_poly_t poly);
int gghlite_enc_fread_raw(FILE * f, gghlite_enc_t poly);
int fmpz_poly_fprint_raw(FILE * file, const fmpz_poly_t poly);
int fmpz_poly_fread_raw(FILE * file, fmpz_poly_t poly);

#endif /* ifndef _MIFE_GLUE_H */
