#ifndef _MIFE_IO_H_
#define _MIFE_IO_H_

#include "mife_defs.h"
#include "flint_raw_io.h"

/* functions dealing with file reading and writing for encodings */
#define gghlite_enc_fprint fmpz_mod_poly_fprint 
#define gghlite_enc_fread fmpz_mod_poly_fread
void fread_gghlite_params(FILE *fp, gghlite_params_t params);
void fwrite_gghlite_params(FILE *fp, gghlite_params_t params);
void fread_gghlite_sk(FILE *fp, gghlite_sk_t self);
void fwrite_gghlite_sk(FILE *fp, gghlite_sk_t self);
void fwrite_mife_pp(const_mmap_vtable mmap, mife_pp_t pp, char *filepath);
void fread_mife_pp(const_mmap_vtable mmap, mife_pp_t pp, char *filepath);
void fwrite_mife_sk(const_mmap_vtable mmap, mife_sk_t sk, char *filepath);
void fread_mife_sk(const_mmap_vtable mmap, mife_sk_t sk, char *filepath);
void fwrite_mife_ciphertext(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t ct, char *filepath);
void fwrite_gghlite_enc_mat(const_mmap_vtable mmap, mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp);
void fread_mife_ciphertext(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t ct, char *filepath);
void fread_gghlite_enc_mat(const_mmap_vtable mmap, const mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp);


#endif /* _MIFE_IO_H_ */
