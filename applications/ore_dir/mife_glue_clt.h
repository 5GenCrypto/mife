#ifndef _MIFE_GLUE_CLT_H
#define _MIFE_GLUE_CLT_H

#include <clt13.h>
#include <mife/mife_defs.h>

struct mmap_pp  { clt_pp    self; };
struct mmap_sk  { clt_state self; };
struct mmap_enc { mpz_t     self; };

extern const mmap_vtable clt13_vtable;

#endif /* ifndef _MIFE_GLUE_CLT_H */
