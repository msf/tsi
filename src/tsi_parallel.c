#include "debug.h"
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



int delete_tsi_parallel() {
#ifdef TSI_MPI
    int mpi_flag = 0;

    if (MPI_Initialized(&mpi_flag) != MPI_SUCCESS) {
        printf_dbg("delete_tsi_parallel(): failed to detect MPI status!\n");
        return -1;
    } else if (mpi_flag) {
        if (MPI_Finalize() != MPI_SUCCESS) {
            printf_dbg("delete_tsi_parallel(): MPI_Finalize failed!\n");
            return -1;
        } else {
            printf_dbg("delete_tsi_parallel(): MPI stoped\n");
            return 1;
        }
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

    if (MPI_Bcast(&nlayers, 1, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
        printf_dbg2("tsi_set_layers_parallel: failed to broadcast nlayers\n");
        return 1;
    }

    if (t->proc_id == t->root_id) {
        if (MPI_Bcast(layers_size, nlayers, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg2("tsi_set_layers_parallel: failed to broadcast layers array\n");
            return 1;
        }
    } else {
        if ((layers_size = (int *) tsi_malloc(sizeof(int) * nlayers)) == NULL) {
            printf_dbg2("tsi_set_layers_parallel: failed to allocate layers array\n");
            return 1;
        }
        if (MPI_Bcast(layers_size, nlayers, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg2("tsi_set_layers_parallel: failed to receive layers array\n");
            return 1;
        }
        build_cmgrid(g, nlayers, layers_size);
    }
#endif /* TSI_MPI */
    return 0;
} /* tsi_set_layers_parallel */



int tsi_is_best_parallel(tsi *t) {
#ifdef TSI_MPI
    corr result;

    if (MPI_Reduce(&t->global_best, &result, 1, MPI_FLOAT_INT, MPI_MAXLOC, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
        printf_dbg("tsi_is_best_parallel: Failed to execute is_best reduce\n");
        return 1;
    }
    t->global_best.value = result.value;
    t->global_best.proc_id = result.proc_id;
#endif /* TSI_MPI */
    return 0;
} /* tsi_is_best_parallel */



int tsi_compare_parallel(tsi *t) {
#ifdef TSI_MPI
    corr corr_data, result;
    int i;
    unsigned int j;
    float ai_val;
    
    t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
    t->nextBCM = load_grid(t->heap, t->nextBAI_idx);
    t->ai = load_grid(t->heap, t->nextBAI_idx);
    t->cm = load_grid(t->heap, t->nextBAI_idx);
    j = 0;
    for (i = 0; i < t->grid_size; i++) {
        corr_data.value = t->cm[i];
        corr_data.proc_id = t->proc_id;
        if (MPI_Reduce(&corr_data, &result, 1, MPI_FLOAT_INT, MPI_MAXLOC, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg("tsi_compare_parallel(): Failed to execute CM reduce\n");
	    return 1;
        }
        if (result.value > corr_data.value) {
            t->cm[i] = result.value;
            ai_val = t->ai[i];
            if (MPI_Bcast(&ai_val, 1, MPI_FLOAT, result.proc_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
                printf_dbg("tsi_compare_parallel: failed to broadcast new AI value\n");
                return 1;
            }
            t->nextBAI[i] = ai_val;
            j++;
        }
    }
#endif /* TSI_MPI */
    return 0;
} /* tsi_compare_parallel */

/* end of tsi_parallel.c */
