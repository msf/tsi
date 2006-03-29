#ifndef _TSI_DEBUG_H
#define _TSI_DEBUG_H

#include "memdebug.h"

#ifdef TSI_DEBUG
#define printf_dbg printf
#define tsi_malloc debug_malloc
#define tsi_free   debug_free
#else
#define printf_dbg //
#define tsi_malloc malloc
#define tsi_free   free
#endif

#ifdef TSI_DEBUG2
#define printf_dbg2 printf
#undef  printf_dbg
#undef  tsi_malloc
#undef  tsi_free
#define printf_dbg printf
#define tsi_malloc debug_malloc
#define tsi_free   debug_free
#else
#define printf_dbg2 //
#endif

#endif /* _TSI_DEBUG_H */

