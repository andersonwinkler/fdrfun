#ifndef EX_BKY_H
#define EX_BKY_H

#ifdef  __cplusplus
//extern "C" {
#endif

#include <stdbool.h>

#ifdef DT32
#define flote float
#else
#define flote double
#endif

flote bky(flote *pvals, flote *padj, size_t nvox, double qval);

#ifdef  __cplusplus
//}
#endif

#endif // EX_BKY_H
