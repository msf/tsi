#ifndef _MEMDEBUG_H
#define _MEMDEBUG_H


#include <stdlib.h>

void * debug_malloc(unsigned int size);
void debug_free(void *);
void debug_check();

#endif
