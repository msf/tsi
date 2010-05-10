#include <stdio.h>
#include <stdlib.h>
#include "memdebug.h"
#include "list.h"

/* this code was adapted from
 * v9fs project
 * it is therefore GPL
 */

#define MAGIC_HEAD	0x5EAD5EAD
#define MAGIC_TAIL	0xDEA5DEA5

#define MAGIC_FREE	0xDEADDEAD

struct tsi_mem {
	unsigned int		magichead;
	unsigned int 		size;
	void*			callerpc;
	struct list_head	mem_list;
};

LIST_HEAD(tsi_mem_list);

/* memdebug prototypes */
void debug_check(void);

void *debug_malloc(unsigned int size)
{
	struct tsi_mem *m;

	m = malloc(sizeof(struct tsi_mem) + size + sizeof(unsigned int));
	if (!m)
		return NULL;

	m->magichead = MAGIC_HEAD;
	m->size = size;
	m->callerpc = __builtin_return_address(0);
	INIT_LIST_HEAD(&m->mem_list);
	*(unsigned int *)((char *)m + sizeof(struct tsi_mem) + size) = MAGIC_TAIL;
	list_add(&m->mem_list, &tsi_mem_list);

	return (char *)m + sizeof(struct tsi_mem);
}


void debug_free(void *p)
{
	struct tsi_mem *m;
	unsigned int *magictail;

	if (!p)
		return;
	// printf( "debug_free(): %p callerpc %p\n", p, __builtin_return_address(0));

	m = (struct tsi_mem*) ((char *)p - sizeof(struct tsi_mem));
	magictail = (unsigned int *)((char *)m + sizeof(struct tsi_mem) + m->size);

	if (m->magichead == MAGIC_FREE && *magictail == MAGIC_FREE) {
		printf( "debug_free(): block %p already freed from function %p\n",
			p, m->callerpc);
		return;
	}

	if (m->magichead != MAGIC_HEAD) {
		printf( "debug_free(): freeing memory that wasn't allocated (or buffer has been underrun) %p function %p\n",
			p, __builtin_return_address(0));
		return;
	}

	if (*magictail != MAGIC_TAIL) {
		printf( "debug_free(): buffer %p has been overriden allocated from %p\n",
			p, m->callerpc);
		return;
	}

	m->magichead = MAGIC_FREE;
	*magictail = MAGIC_FREE;
	m->callerpc = __builtin_return_address(0);

	list_del(&m->mem_list);
	free(m);
}

void debug_check(void)
{
	struct tsi_mem *m, *mtmp;

	list_for_each_entry_safe(m, mtmp, &tsi_mem_list, mem_list) {
		printf("debug_check(): leak %p (size: %lubytes) allocated from %p\n", 
			(char *)m + sizeof(struct tsi_mem), m->size*sizeof(char), m->callerpc);
	}
}
