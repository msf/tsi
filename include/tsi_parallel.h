#ifndef _TSI_PARALLEL_H
#define _TSI_PARALLEL_H

#include "tsi.h"

int new_tsi_parallel(int *n_procs, int *proc_id);

int delete_tsi_parallel();

/* TSI functions */
int tsi_is_best_parallel(tsi *t);

int tsi_compare_parallel(tsi *t);

int tsi_generate_layers_parallel(tsi *t, int *nlayers, int **layer_size);


#endif /* _TSI_PARALLEL_H */
