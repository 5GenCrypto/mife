#ifndef _MIFE_IO_H_
#define _MIFE_IO_H_

#include "mife_defs.h"

/* functions dealing with file reading and writing for encodings */
#define gghlite_enc_fprint fmpz_mod_poly_fprint 
void fread_gghlite_params(FILE *fp, gghlite_params_t params);
void fwrite_gghlite_params(FILE *fp, gghlite_params_t params);
int gghlite_enc_fread(FILE * f, gghlite_enc_t poly);
void fwrite_mife_pp(mife_pp_t pp, char *filepath);
void fread_mife_pp(mife_pp_t pp, char *filepath);
void fwrite_mife_sk(mife_sk_t sk, char *filepath);
void fread_mife_sk(mife_sk_t sk, char *filepath);
void fwrite_mife_ciphertext(mife_pp_t pp, mife_ciphertext_t ct, char *filepath);
void fwrite_gghlite_enc_mat(mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp);
void fread_mife_ciphertext(mife_pp_t pp, mife_ciphertext_t ct, char *filepath);
void fread_gghlite_enc_mat(const mife_pp_t pp, gghlite_enc_mat_t m, FILE *fp);


#endif /* _MIFE_IO_H_ */

