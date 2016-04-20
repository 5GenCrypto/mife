#ifndef _LIBMMAP_MMAP_H_
#define _LIBMMAP_MMAP_H_

#include <gghlite/gghlite-defs.h>
#include <clt13.h>

struct mmap_pp {
  union {
    gghlite_params_t gghlite_self;
    clt_pp clt13_self;
  };
};

struct mmap_sk {
  union {
    gghlite_sk_t gghlite_self;
    clt_state clt13_self;
  };
};

struct mmap_enc {
  union {
    gghlite_enc_t gghlite_self;
    clt_elem_t clt13_self;
  };
};

#endif
