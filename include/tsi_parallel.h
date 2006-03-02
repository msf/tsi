#ifndef TSI_PARALLEL_H
#define TSI_PARALLEL_H

#include "registry.h"

typedef struct tsi_parallel_type {
    int rank;
    int comsize;
} tsi;

tsi *new_tsi(registry *r);

void delete_tsi(tsi_parallel *);


#endif /* TSI_PARALLEL_H */
