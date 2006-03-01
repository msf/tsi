#ifndef _TSI_DEBUG_H
#define _TSI_DEBUG_H

#ifdef TSI_DEBUG
#define printf_dbg printf
#else
#define printf_dbg //
#endif

#ifdef TSI_DEBUG2
#define printf_dbg2 printf
#undef  printf_dbg
#define printf_dbg printf
#else
#define printf_dbg2 //
#endif

#endif //_TSI_DEBUG_H

