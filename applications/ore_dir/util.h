#ifndef _ORE_UTILS_H
#define _ORE_UTILS_H

#ifdef ALLOC_FAILS
#error "trying to redefine ALLOC_FAILS"
#else  /* ifdef ALLOC_FAILS */
#define ALLOC_FAILS(path, len) (NULL == ((path) = malloc((len) * sizeof(*(path)))))
#endif /* ifdef ALLOC_FAILS */

#ifdef SEED_SIZE
#error "trying to redefine SEED_SIZE"
#else /* ifdef SEED_SIZE */
#define SEED_SIZE 32
#endif /* ifdef SEED_SIZE */

#endif /* ifndef _ORE_UTILS_H */
