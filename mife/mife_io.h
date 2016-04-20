#ifndef _MIFE_IO_H_
#define _MIFE_IO_H_

#include "mmap/mmap.h"
#include "mife_defs.h"
#include "flint_raw_io.h"

void fwrite_mife_pp(const_mmap_vtable mmap, mife_pp_t pp, char *filepath);
void fread_mife_pp(const_mmap_vtable mmap, mife_pp_t pp, char *filepath);
void fwrite_mife_sk(const_mmap_vtable mmap, mife_sk_t sk, char *filepath);
void fread_mife_sk(const_mmap_vtable mmap, mife_sk_t sk, char *filepath);
void fwrite_mife_ciphertext(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t ct, char *filepath);
void fwrite_mmap_enc_mat(const_mmap_vtable mmap, mmap_enc_mat_t m, FILE *fp);
void fread_mife_ciphertext(const_mmap_vtable mmap, mife_pp_t pp, mife_ciphertext_t ct, char *filepath);
void fread_mmap_enc_mat(const_mmap_vtable mmap, mmap_enc_mat_t m, FILE *fp);


#endif /* _MIFE_IO_H_ */
