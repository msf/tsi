#ifdef TSI_MPI
#include <mpi.h>
#endif

#include "debug.h"
#include "log.h"
#include "tsi.h"
#include "tsi_parallel.h"
#include "tsi_math.h"



int new_tsi_parallel(int *n_procs, int *proc_id) {
#ifdef TSI_MPI
    if (MPI_Init(NULL, NULL) != MPI_SUCCESS) {      /* MPI-2 */
        printf_dbg("new_tsi_parallel(): MPI_Init failed!\n");
        return -1;
    } else if (MPI_Comm_size(MPI_COMM_WORLD, n_procs) != MPI_SUCCESS) {
        printf_dbg("new_tsi_parallel(): MPI_Comm_size failed!\n");
        return -1;
    } else if (MPI_Comm_rank(MPI_COMM_WORLD, proc_id) != MPI_SUCCESS) {
        printf_dbg("new_tsi_parallel(): MPI_Comm_rank failed!\n");    
        return -1;
    } else {
        printf_dbg("new_tsi_parallel(): MPI started...\n");
        return 1;
    }
#endif /* TSI_MPI */
    *proc_id = 0;
    *n_procs = 1;
    return 0;
} /* new_tsi_parallel */



int delete_tsi_parallel(tsi *t) {
#ifdef TSI_MPI
    int mpi_flag = 0;

    if (MPI_Initialized(&mpi_flag) != MPI_SUCCESS) {
        printf_dbg("delete_tsi_parallel(): failed to detect MPI status!\n");
        fflush(stdout);
        return -1;
    } else if (mpi_flag) {
        printf_dbg("delete_tsi_parallel(): entering MPI_Barrier\n");    
        if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg("delete_tsi_parallel(): MPI_Barrier failed!\n");    
            fflush(stdout);
            return -1;
        }
		printf("delete_tsi_parallel() calling MPI_Finalize()\n");
        if (MPI_Finalize() != MPI_SUCCESS) {
            printf_dbg("delete_tsi_parallel(): MPI_Finalize failed!\n");
            fflush(stdout);
            return -1;
        } 
		printf_dbg("delete_tsi_parallel(): MPI stoped\n");
		printf("delete_tsi_parallel() FINISHED OK\n");
		fflush(stdout);
		return 1;
    } else
		return 1;
#endif /* TSI_MPI */
    return 0;
}



int tsi_set_layers_parallel(tsi *t, cm_grid *g) {
#ifdef TSI_MPI
    int nlayers, *layers_size;

    nlayers = get_nlayers(g);
    layers_size = get_layers(g);

    printf_dbg("tsi_set_layers_parallel(%d): starting broadcast of the number of layers\n", t->proc_id);
    if (MPI_Bcast(&nlayers, 1, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
        printf_dbg2("tsi_set_layers_parallel: failed to broadcast nlayers\n");
            fflush(stdout);
        return 1;
    }


    if (t->proc_id == t->root_id) {
        printf_dbg("tsi_set_layers_parallel(%d): broadcast layers array\n", t->proc_id);
        if (MPI_Bcast(layers_size, nlayers, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg2("tsi_set_layers_parallel: failed to broadcast layers array\n");
            fflush(stdout);
            return 1;
        }
        printf_dbg("tsi_set_layers_parallel(%d): finished sending layers\n", t->proc_id);
    } else {
        printf_dbg("tsi_set_layers_parallel(%d): received %d layers\n", t->proc_id, nlayers);
        if ((layers_size = (int *) tsi_malloc(sizeof(int) * nlayers)) == NULL) {
            printf_dbg2("tsi_set_layers_parallel: failed to allocate layers array\n");
            fflush(stdout);
            return 1;
        }
        printf_dbg("tsi_set_layers_parallel(%d): waiting for layers array\n", t->proc_id);
        if (MPI_Bcast(layers_size, nlayers, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg2("tsi_set_layers_parallel: failed to receive layers array\n");
            fflush(stdout);
            return 1;
        }
        build_cmgrid(g, nlayers, layers_size);
        printf_dbg("tsi_set_layers_parallel(%d): finished receiving layers\n", t->proc_id);
            fflush(stdout);
    }
#endif /* TSI_MPI */
    return 0;
} /* tsi_set_layers_parallel */



int tsi_is_best_parallel(tsi *t) {
#ifdef TSI_MPI
    corr result;

    printf_dbg("tsi_is_best_parallel(%d): start reducing best correlation (%f)\n", t->proc_id, t->global_best.value);
    if (MPI_Allreduce(&t->global_best, &result, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD) != MPI_SUCCESS) {
        printf_dbg("tsi_is_best_parallel: Failed to execute is_best reduce\n");
        return 1;
    }
    t->global_best.value = result.value;
    t->global_best.proc_id = result.proc_id;
    printf_dbg("tsi_is_best_parallel(%d): finished reducing best correlation\n", t->proc_id);
    printf_dbg("tsi_is_best_parallel(%d): result: max %f, loc %d\n", t->proc_id, result.value, result.proc_id);
#endif /* TSI_MPI */
    return 0;
} /* tsi_is_best_parallel */



int tsi_compare_parallel(tsi *t) {
#ifdef TSI_MPI
    corr corr_data, result;
    int z0, z1, *layer_size, nlayers, last_layer, layer, nxy, ai_p;
    unsigned int i, j, grid_size;
    float ai_val;
    
    t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
    //t->nextBCM = load_grid(t->heap, t->nextBCM_idx);

    t->nextBCM = t->nextBCM_c->cg;
    layer_size = t->nextBCM_c->layer_size;
    nlayers = t->nextBCM_c->nlayers;
    nxy = t->nextBCM_c->nxy;
    grid_size = nlayers * nxy;  /* compressed grid size */

    z0 = 0;
    z1 = layer_size[0];
    layer = last_layer = 0;
    log_message(t->l,0,"tsi_compare_parallel()");
    printf_dbg("parallel_compare(%d): layer: %d, z0: %d, z1:%d\n", t->proc_id, layer, z0, z1);
    for (i = 0; i < grid_size; i++) {
        corr_data.value = t->nextBCM[i];
        corr_data.proc_id = t->proc_id;
        if (MPI_Allreduce(&corr_data, &result, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg("tsi_compare_parallel(): Failed to execute CM reduce\n");
	    return 1;
        }
        t->nextBCM[i] = result.value;   /* compressed grid... */

        /* broadcast new best AI values */
        layer = i / nxy;
        if (layer > last_layer) {   /* ugly hack, I know... */
            last_layer = layer;
            z0 = z1;
            z1 += layer_size[layer];
            printf_dbg("parallel_compare(%d): layer: %d, z0: %d, z1:%d\n", t->proc_id, layer, z0, z1);            
        }

        for (j = z0; j < z1; j++) {
            ai_p = (j * nxy) + (i - layer*nxy);
            //printf_dbg("parallel_compare(%d): size: %d, nlayers: %d, layer: %d, i: %d, [i]: %d\n", t->proc_id, grid_size, nlayers, layer, i, (i-layer*nxy));
            //printf_dbg("parallel_compare(%d): i: %d, ai_p: %d\n", t->proc_id, i, ai_p);
            fflush(stdout);
            ai_val = t->nextBAI[ai_p];
            if (MPI_Bcast(&ai_val, 1, MPI_FLOAT, result.proc_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
                printf_dbg("tsi_compare_parallel: failed to broadcast new AI value\n");
                fflush(stdout);
                return 1;
            }
            t->nextBAI[ai_p] = ai_val;
        }
    }
    log_message(t->l,0,"tsi_compare_parallel() finished");

    dirty_grid(t->heap, t->nextBAI_idx);
    //dirty_grid(t->heap, t->nextBCM_idx);
#endif /* TSI_MPI */
    return 0;
} /* tsi_compare_parallel */

/* end of tsi_parallel.c */
