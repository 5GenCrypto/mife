#ifndef _ORE_H_
#define _ORE_H_

#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <gghlite/gghlite-defs.h>

typedef enum {
  ORE_DEFAULT    = 0x00, //!< default behaviour
  ORE_NO_KILIAN    = 0x01, //!< do not multiply kilian randomizers into the encodings
  ORE_NO_RANDOMIZERS    = 0x02, //!< do not multiply the scalar randomizers into the encodings
  ORE_SIMPLE_PARTITIONS  = 0x04, //!< pick a simple partitioning (x[0] is encoded ta the universe, all others are encoded at the empty set.)
} ore_flag_t;



#endif /* _ORE_H_ */
