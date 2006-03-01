#ifndef _TSI_PARALLEL_H
#define _TSI_PARALLEL_H

#include "registry.h"

typedef struct tsi_parallel_type {
    int rank;
    int comsize;
} tsi_parallel;

tsi_parallel *new_tsi_parallel(registry *r);

void delete_tsi_parallel(tsi_parallel *);


#endif /* _TSI_PARALLEL_H */
