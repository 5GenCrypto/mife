#ifndef _ORE_UTILS_H

#ifdef ALLOC_FAILS
#error "trying to redefine ALLOC_FAILS"
#else  /* ifdef ALLOC_FAILS */
#define ALLOC_FAILS(path, len) (NULL == ((path) = malloc((len) * sizeof(*(path)))))
#endif /* ifdef ALLOC_FAILS */

#endif /* ifndef _ORE_UTILS_H */
