#ifdef TSI_MPI
#include <mpi.h>
#endif

#include "debug.h"
#include "log.h"
#include "tsi.h"
#include "tsi_parallel.h"
#include "tsi_math.h"

int tsi_compare_parallel_v1(tsi *t);
int tsi_compare_parallel_v2(tsi *t);
int tsi_compare_parallel_v3(tsi *t);

void local_compare_update(float *lBCM, float *lBAI, float *rBCM, float *rBAI, int size);

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


int tsi_compare_parallel(tsi *t) 
{
	//return tsi_compare_parallel_v1(t);
	//return tsi_compare_parallel_v2(t);
	return tsi_compare_parallel_v3(t);

}


int tsi_compare_parallel_v1(tsi *t) {
#ifdef TSI_MPI
    corr corr_data, result;
    int i;
    unsigned int j;
    float ai_val;
    
    t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
	t->nextBCM_idx = new_grid(t->heap);
    t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
	expand_correlations_grid(t->nextBCM_c, t->nextBCM);
	delete_cmgrid(t->nextBCM_c);

    j = 0;
	log_message(t->l,0,"tsi_compare_parallel()");
    for (i = 0; i < t->grid_size; i++) {
        corr_data.value = t->nextBCM[i];
        corr_data.proc_id = t->proc_id;


		if (MPI_Allreduce(&corr_data, &result, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD) != MPI_SUCCESS) {
            log_string(t->l,"tsi_compare_parallel(): Failed to execute CM reduce\n");
	    	return 1;
        }
		/*
        if (MPI_Reduce(&corr_data, &result, 1, MPI_FLOAT_INT, MPI_MAXLOC, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg("tsi_compare_parallel(): Failed to execute CM reduce\n");
	    return 1;
        }

        // broadcast result of reduce
        if (MPI_Bcast(&result, 1, MPI_FLOAT_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg("tsi_compare_parallel: Failed to execute is_best reduce\n");
            fflush(stdout);
            return 1;
        }
		*/
		
        t->nextBCM[i] = result.value;

        /* broadcast new best AI value */
        ai_val = t->nextBAI[i];
        if (MPI_Bcast(&ai_val, 1, MPI_FLOAT, result.proc_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg("tsi_compare_parallel: failed to broadcast new AI value\n");
            fflush(stdout);
            return 1;
        }
        t->nextBAI[i] = ai_val;
    }
	log_message(t->l,0,"tsi_compare_parallel() finished");
    dirty_grid(t->heap, t->nextBAI_idx);
    dirty_grid(t->heap, t->nextBCM_idx);
#endif /* TSI_MPI */
    return 0;
} /* tsi_compare_parallel_v1 */

int tsi_compare_parallel_v2(tsi *t) 
{
#ifdef TSI_MPI
    corr corr_data, result;
    int z0, z1, *layer_size, nlayers, last_layer, layer, nxy, ai_p;
    unsigned int i, j, compressed_grid_size;
    float ai_val;
	int ret, x, y;
    
    t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
    //t->nextBCM = load_grid(t->heap, t->nextBCM_idx);

	float *compressedBCM = t->nextBCM_c->cg;
    layer_size = t->nextBCM_c->layer_size;
    nlayers = t->nextBCM_c->nlayers;
    nxy = t->nextBCM_c->nxy;
    compressed_grid_size = nlayers * nxy;  /* compressed grid size */

    z0 = 0;
    z1 = layer_size[0];
    layer = last_layer = 0;
    log_message(t->l,0,"tsi_compare_parallel()");
    printf_dbg("parallel_compare(%d): layer: %d, z0: %d, z1:%d\n", t->proc_id, layer, z0, z1);
    for (i = 0; i < compressed_grid_size; i++) {
        corr_data.value = compressedBCM[i];
        corr_data.proc_id = t->proc_id;
		ret = MPI_Allreduce(&corr_data, &result, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        if ( ret != MPI_SUCCESS) {
            printf_dbg("tsi_compare_parallel(): Failed to execute CM reduce\n");
		    return 1;
        }
		/*
        if (MPI_Reduce(&corr_data, &result, 1, MPI_FLOAT_INT, MPI_MAXLOC, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg("tsi_compare_parallel(): Failed to execute CM reduce\n");
	    return 1;
        }

        // broadcast result of reduce
        if (MPI_Bcast(&result, 1, MPI_FLOAT_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            printf_dbg("tsi_compare_parallel: Failed to execute is_best reduce\n");
            fflush(stdout);
            return 1;
        }
		*/

        compressedBCM[i] = result.value;   /* compressed grid... */

        /* broadcast new best AI values */
        layer = i / nxy;
        if (layer > last_layer) { 
            last_layer = layer;
            z0 = z1;
            z1 += layer_size[layer];
            printf_dbg("parallel_compare(%d): layer: %d, z0: %d, z1:%d\n", t->proc_id, layer, z0, z1);            
        }

		y = (i / t->xsize) - (t->ysize * layer);
		x = i % t->xsize;
//		printf_dbg("parallel_compare(%d): corr coord: (%d,%d,%d)\n", t->proc_id, x,y,layer);

        for (j = z0; j < z1; j++) {
            ai_p = (j * nxy) + (i - layer*nxy);
            //printf_dbg("parallel_compare(%d): size: %d, nlayers: %d, layer: %d, i: %d, [i]: %d\n", t->proc_id, grid_size, nlayers, layer, i, (i-layer*nxy));
            //printf_dbg("ai_p: %d\t", ai_p);
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

	/* expand BCM */
    t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
	t->nextBCM_idx = new_grid(t->heap);
    t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
	expand_correlations_grid(t->nextBCM_c, t->nextBCM);
	delete_cmgrid(t->nextBCM_c);

    dirty_grid(t->heap, t->nextBAI_idx);
    dirty_grid(t->heap, t->nextBCM_idx);
#endif /* TSI_MPI */
    return 0;
} /* tsi_compare_parallel_v2 */


/* totally sincronous distributed compate & update */
int tsi_compare_parallel_v3(tsi *t)
{
#ifdef TSI_MPI
	int ret;
	int i;
	int to, from, start;
	float *smallBCM, *smallBAI;
	float *BAI, *BCM;

	// MPI stuff
	int rank, clustersize;
	int size_per_rank;
	int rest;
	MPI_Status  stat;


	//log_message(t->l, 0, "compare_parallel_v3: Called!");

	t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
	t->nextBCM_idx = new_grid(t->heap);
	t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
	expand_correlations_grid(t->nextBCM_c, t->nextBCM);
	delete_cmgrid(t->nextBCM_c);

	BAI = t->nextBAI;
	BCM = t->nextBCM;

	MPI_Comm_size(MPI_COMM_WORLD, &clustersize);

	size_per_rank = t->grid_size / clustersize;
	rest = t->grid_size % clustersize;
	smallBCM = (float *) malloc( size_per_rank * sizeof(float));
	smallBAI = (float *) malloc( size_per_rank * sizeof(float));
	if(rest > 0)
		log_message(t->l, 0, "compare_parallel_v3: WE HAVE REST!");


	/* 
	   circular send & recieve of parcels of BCM + BAI 
	   total DATA through network: grid_size * 2 * sizeof(float) 
	   the last node in group also works the "rest" 
	 */

	to = t->proc_id;
	from = t->proc_id;

	to = (to + 1) % clustersize;
	from = (from  + clustersize - 1) % clustersize;
	start = to * size_per_rank;

	while(to != t->proc_id) {

		ret = MPI_Sendrecv(BCM + start, size_per_rank, MPI_FLOAT, to, to,
				smallBCM, size_per_rank, MPI_FLOAT, from, t->proc_id,
				MPI_COMM_WORLD, &stat);
		if( ret != MPI_SUCCESS) {
			log_string(t->l, "tsi_compare_parallel(): mpi_sendrecv error!\n");
			return 1;
		}
		ret = MPI_Sendrecv(BAI + start, size_per_rank, MPI_FLOAT, to, to+1,
				smallBAI, size_per_rank, MPI_FLOAT, from, t->proc_id+1,
				MPI_COMM_WORLD, &stat);
		if( ret != MPI_SUCCESS) {
			log_string(t->l, "tsi_compare_parallel(): mpi_sendrecv error!\n");
			return 1;
		}

		i = t->proc_id * size_per_rank;
		local_compare_update( BCM + i, BAI + i, smallBCM, smallBAI, size_per_rank);

		if( to == (clustersize - 1) && rest != 0) { /* last gets the rest */
			i = start +  size_per_rank;
			ret = MPI_Send(BCM + i, rest, MPI_FLOAT, to, to, MPI_COMM_WORLD);
			if( ret != MPI_SUCCESS) {
				log_string(t->l, "tsi_compare_parallel(): mpi_send error!\n");
				return 1;
			}
			ret = MPI_Send(BAI + i, rest, MPI_FLOAT, to, to+1, MPI_COMM_WORLD);
			if( ret != MPI_SUCCESS) {
				log_string(t->l, "tsi_compare_parallel(): mpi_send error!\n");
				return 1;
			}
		}
		if( t->proc_id == (clustersize -1)  && rest != 0) {/* last gets the rest */
			i = (t->proc_id *size_per_rank ) + size_per_rank;
			ret = MPI_Recv(smallBCM, rest, MPI_FLOAT, from, t->proc_id, MPI_COMM_WORLD, &stat);
			if( ret != MPI_SUCCESS) {
				log_string(t->l, "tsi_compare_parallel(): mpi_recv error!\n");
				return 1;
			}
			ret = MPI_Recv(smallBAI, rest, MPI_FLOAT, from, t->proc_id+1, MPI_COMM_WORLD, &stat);
			if( ret != MPI_SUCCESS) {
				log_string(t->l, "tsi_compare_parallel(): mpi_recv error!\n");
				return 1;
			}
			local_compare_update( BCM + i, BAI + i, smallBCM, smallBAI, rest);
		}

		/* step to next round */
		to = (to + 1) % clustersize;
		from = (from  + clustersize - 1) % clustersize;
		start = to * size_per_rank;
	} 

	/* broadcast new  BCM & BAI to all
	 * total DATA throught network : grid_size * 2 * sizeof(float) */
	for(i = 0; i < clustersize; i++) {
		start = i * size_per_rank;
		ret = MPI_Bcast(BCM+start, size_per_rank, MPI_FLOAT, i, MPI_COMM_WORLD);
		if( ret != MPI_SUCCESS) {
			log_string(t->l, "tsi_compare_parallel(): mpi_bcast error!\n");
			return 1;
		}
		ret = MPI_Bcast(BAI+start, size_per_rank, MPI_FLOAT, i, MPI_COMM_WORLD);
		if( ret != MPI_SUCCESS) {
			log_string(t->l, "tsi_compare_parallel(): mpi_bcast error!\n");
			return 1;
		}
	}
	/* last node in group */
	if(rest != 0) {
		start += size_per_rank;
		ret = MPI_Bcast(BCM+start, rest, MPI_FLOAT, clustersize-1, MPI_COMM_WORLD);
		if( ret != MPI_SUCCESS) {
			log_string(t->l, "tsi_compare_parallel(): mpi_bcast error!\n");
			return 1;
		}
		ret = MPI_Bcast(BAI+start, rest, MPI_FLOAT, clustersize-1, MPI_COMM_WORLD);
		if( ret != MPI_SUCCESS) {
			log_string(t->l, "tsi_compare_parallel(): mpi_bcast error!\n");
			return 1;
		}
	}

	free(smallBCM);
	free(smallBAI);
	//	log_message(t->l,0,"tsi_compare_parallel_v3() finished");
	dirty_grid(t->heap, t->nextBAI_idx);
	dirty_grid(t->heap, t->nextBCM_idx);

#endif /* TSI_MPI */
	return 0;
} /* end of tsi_compare_parallel_v3 */


void local_compare_update(float *lBCM, float *lBAI, float *rBCM, float *rBAI, int size)
{
	int i;
	for(i = 0; i < size; i++) {
		if(lBCM[i] < rBCM[i]) { // new best correlation 
			lBCM[i] = rBCM[i];
			lBAI[i] = rBAI[i];
		}
	}
}
/* end of tsi_parallel.c */
