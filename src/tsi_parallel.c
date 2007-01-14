#ifdef TSI_MPI
#include <mpi.h>
#endif

#include "debug.h"
#include "log.h"
#include "tsi.h"
#include "tsi_parallel.h"
#include "tsi_math.h"
#include "timer.h"

/* local prototypes */
int expand_correlations_grid(cm_grid *cmg, float *CM);
void local_compare_update(float *lBCM, float *lBAI, float *rBCM, float *rBAI, int size);
int tsi_compare_parallel_direct(tsi *t);
int tsi_compare_parallel_collective(tsi *t);
int tsi_compare_parallel_v3(tsi *t);
int tsi_compare_parallel_v4(tsi *t);
int tsi_compare_parallel_direct_v5(tsi *t);



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
        log_message(t->l, 0, "delete_tsi_parallel(): failed to detect MPI status!\n");
        return -1;
    } else if (mpi_flag) {
        printf_dbg("delete_tsi_parallel(): entering MPI_Barrier\n");    
        if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
            log_message(t->l, 0, "delete_tsi_parallel(): MPI_Barrier failed!\n");    
            return -1;
        }
        printf("delete_tsi_parallel() calling MPI_Finalize()\n");
        if (MPI_Finalize() != MPI_SUCCESS) {
            log_message(t->l, 0, "delete_tsi_parallel(): MPI_Finalize failed!\n");
            return -1;
        } 
        printf_dbg("delete_tsi_parallel(): MPI FINISHED OK\n");
        fflush(stdout);
        return 1;
    } else
        return 1;
#endif /* TSI_MPI */
    return 0;
} /* delete_tsi_parallel */



int tsi_set_layers_parallel(tsi *t, cm_grid *g) {
#ifdef TSI_MPI
    unsigned int nlayers, *layers_size;

    nlayers = get_nlayers(g);
    layers_size = get_layers(g);

    printf_dbg("tsi_set_layers_parallel(%d): starting broadcast of the number of layers\n", t->proc_id);
    if (MPI_Bcast(&nlayers, 1, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
        log_message(t->l, 0, "tsi_set_layers_parallel: failed to broadcast nlayers\n");
        return 1;
    }


    if (t->proc_id == t->root_id) {
        printf_dbg("tsi_set_layers_parallel(%d): broadcast layers array\n", t->proc_id);
        if (MPI_Bcast(layers_size, nlayers, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            log_message(t->l, 0, "tsi_set_layers_parallel: failed to broadcast layers array\n");
            return 1;
        }
        printf_dbg("tsi_set_layers_parallel(%d): finished sending layers\n", t->proc_id);
    } else {
        printf_dbg("tsi_set_layers_parallel(%d): received %d layers\n", t->proc_id, nlayers);
        if ((layers_size = (unsigned int *) tsi_malloc(sizeof(unsigned int) * nlayers)) == NULL) {
            log_message(t->l, 0, "tsi_set_layers_parallel: failed to allocate layers array\n");
            return 1;
        }
        printf_dbg("tsi_set_layers_parallel(%d): waiting for layers array\n", t->proc_id);
        if (MPI_Bcast(layers_size, nlayers, MPI_INT, t->root_id, MPI_COMM_WORLD) != MPI_SUCCESS) {
            log_message(t->l, 0, "tsi_set_layers_parallel: failed to receive layers array\n");
            return 1;
        }
        /* need to build cmg grid */
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

#ifdef TSI_MPI
    int ret;
    struct timeval t1, t2;
    double par_time;

    if (t->compare) {
        ret = tsi_compare_parallel_collective(t);
        if (!ret) {
            t->nextBCM_idx = new_grid(t->heap);
            t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
            load_cmgrid(t->nextBCM_c);
            expand_correlations_grid(t->nextBCM_c, t->nextBCM);
            delete_cmgrid(t->nextBCM_c);
            dirty_grid(t->heap, t->nextBCM_idx);
        }
    } else {
        ret = tsi_compare_parallel_direct(t);
    }
    return ret;
#else
    t->nextBCM_idx = new_grid(t->heap);
    t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
    load_cmgrid(t->nextBCM_c);
    expand_correlations_grid(t->nextBCM_c, t->nextBCM);
    delete_cmgrid(t->nextBCM_c);
    return 0;
#endif /* TSI_MPI */
} /* tsi_compare_parallel */



int tsi_compare_parallel_collective(tsi *t) {
#ifdef TSI_MPI
    cm_grid *bcm;
    float   *bai, *ai_z;
    int     ai_z_idx;

    float *cc,          /* orginal BCM values */
          *rv;          /* results array */
    void  *recv_buf,    /* receive buffer */
          *send_buf;    /* send buffer */
    int   *nv;          /* locations array */

    unsigned int fragment_size,  /* size of each fragment of a grid to be shared */
                 cc_size;        /* size of compressed correlations grid */

    unsigned int n, g, h, i, z0, z1;        /* aux variables */
    unsigned int layer, last_layer;    
    struct timeval t1, t2, t3, t4, t5, t6, t7;
    double par_time, mm_time, run_time;

    getCurrTime(&t1);

    /* distribute all values */
    bcm = t->nextBCM_c;
    load_cmgrid(bcm);
    cc = bcm->cg;

    cc_size = bcm->nxy * bcm->nlayers;
    fragment_size = (cc_size / t->n_procs) + ((cc_size % t->n_procs > 0) ? 1 : 0);
    cc_size = fragment_size * t->n_procs;  /* new "clustered" size */

    //--------------------------- TODO ----------------------------------------------------------------
    if (t->heap->grid_size < (2*cc_size + 2*fragment_size)) {
        printf_dbg("tsi_compare_parallel_collective(): not enough space available on BCM cm_grid...");
        return 1;
    }
    //--------------------------- END OF TODO ---------------------------------------------------------

    send_buf = cc;
    recv_buf = cc + cc_size + fragment_size;
    printf_dbg("tsi_compare_parallel_collective(): performing all to all\n");
    getCurrTime(&t2);
    if (MPI_Alltoall(send_buf, fragment_size, MPI_FLOAT, recv_buf, fragment_size, MPI_FLOAT, MPI_COMM_WORLD) != MPI_SUCCESS) {
        log_string(t->l,"tsi_compare_parallel_collective(): Failed to execute all to all communication\n");
        return 1;
    }
    getCurrTime(&t3);

    /* compare and select best values */
    printf_dbg("tsi_compare_parallel_collective(): begining local compare\n");
    nv = (int *) cc + cc_size;          /* set pointer for nodes array */
    rv = recv_buf;                      /* set pointer for results array (same as recv_buf) */
    for (n = 1; n < t->n_procs; n++) {
        for (g = 0; g < fragment_size; g++) {
            h = n * fragment_size + g;
            if (rv[g] < rv[h]) {
                rv[g] = rv[h];
                nv[g] = n;
            }
        }
    }

    /* distribute results lists */
    getCurrTime(&t4);
    send_buf = rv;    /* get new BCM */
    recv_buf = cc;
    printf_dbg("tsi_compare_parallel_collective(): performing gather float\n");
    if (MPI_Allgather(send_buf, fragment_size, MPI_FLOAT, recv_buf, fragment_size, MPI_FLOAT, MPI_COMM_WORLD) != MPI_SUCCESS) {
        log_string(t->l,"tsi_compare_parallel_collective(): Failed to execute gather all communication\n");
        return 1;
    }

    send_buf = nv;
    recv_buf = nv + fragment_size;
    printf_dbg("tsi_compare_parallel_collective(): performing gather int\n");
    if (MPI_Allgather(send_buf, fragment_size, MPI_INT, recv_buf, fragment_size, MPI_FLOAT, MPI_COMM_WORLD) != MPI_SUCCESS) {
        log_string(t->l,"tsi_compare_parallel_collective(): Failed to execute gather all communication\n");
        return 1;
    }
    getCurrTime(&t5);
    nv += fragment_size;

    /* distribute best AI values */
    bai = load_grid(t->heap, t->nextBAI_idx);
    ai_z_idx = new_grid(t->heap);
    if (ai_z_idx < 0) {
        printf_dbg("tsi_compare_parallel_collective(): failed to allocate ai_z grid\n");
        return 1;
    }
    ai_z = load_grid(t->heap, ai_z_idx);
    if (ai_z == NULL) {
        printf_dbg("tsi_compare_parallel_collective(): failed to load ai_z grid\n");
        return 1;
    }
    cc_size = bcm->nxy * bcm->nlayers;
    getCurrTime(&t6);
    for (n = 0; n < t->n_procs; n++) {
        h = 0;
        if (n == t->proc_id) {      /* broadcast local best AI values*/

            z0 = 0;
            z1 = bcm->layer_size[0];
            layer = last_layer = 0;
            for (g = 0; g < cc_size; g++) {        /* build z vectors AI array */
                layer = g / bcm->nxy;
                if (layer > last_layer) { 
                    last_layer = layer;
                    z0 = z1;
                    z1 += bcm->layer_size[layer];
                    printf_dbg("parallel_compare(%d): layer: %d, z0: %d, z1:%d\n", t->proc_id, layer, z0, z1);            
                }
                if (nv[g] == n) {
                    for (i = z0; i < z1; i++) {
                        ai_z[h] = bai[(i * bcm->nxy) + (g - layer*bcm->nxy)];
                        h++;
                    } /* for */
                }
            } /* for */

            if (h > 0) {
                if (MPI_Bcast(ai_z, h, MPI_FLOAT, n, MPI_COMM_WORLD) != MPI_SUCCESS) {
                    printf_dbg("tsi_compare_parallel: failed to broadcast new AI value\n");
                    return 1;
                }
            }
            printf_dbg("ID: %d - sending size: %d\n",n,h);

        } else {                    /* receive new AI values */

            for (g = 0; g < cc_size; g++) {   /* calc array size of broadcast */
                if (nv[g] == n) {
                    h += bcm->layer_size[g/bcm->nxy];
                }
            } /* for */

            printf_dbg("for ID: %d - recv size: %d\n",n,h);
            if (h > 0) {
                if (MPI_Bcast(ai_z, h, MPI_FLOAT, n, MPI_COMM_WORLD) != MPI_SUCCESS) {
                    printf_dbg("tsi_compare_parallel: failed to broadcast new AI value\n");
                    return 1;
                }

                z0 = 0;
                z1 = bcm->layer_size[0];
                layer = last_layer = 0;
                h = 0;
                for (g = 0; g < cc_size; g++) {        /* unpack z vectors AI array */
                    layer = g / bcm->nxy;
                    if (layer > last_layer) { 
                        last_layer = layer;
                        z0 = z1;
                        z1 += bcm->layer_size[layer];
                        printf_dbg("parallel_compare(%d): layer: %d, z0: %d, z1:%d\n", t->proc_id, layer, z0, z1);            
                    }
                    if (nv[g] == n) {
                        for (i = z0; i < z1; i++) {
                            bai[(i * bcm->nxy) + (g - layer*bcm->nxy)] = ai_z[h];
                            h++;
                        }
                    }
                } /* for */
            }

        }
    } /* for */
    getCurrTime(&t7);

    /* clear grids */
    dirty_cmgrid(bcm);
    dirty_grid(t->heap, t->nextBAI_idx);
    delete_grid(t->heap, ai_z_idx);

    mm_time = getElapsedTime(&t1, &t2);
    run_time = getElapsedTime(&t3, &t4);
    par_time += getElapsedTime(&t4, &t5);
    mm_time += getElapsedTime(&t5, &t6);
    log_action_time(t->l, 0, "parallel_compare_collective(): TIME FOR PARALLEL DISTRIBUTION OF BEST RESULTS", par_time);
    log_action_time(t->l, 0, "parallel_compare_collective(): TIME FOR PARALLEL UPDATE OF AI VALUES", getElapsedTime(&t6, &t7));
    par_time += getElapsedTime(&t6, &t7);
    log_action_time(t->l, 0, "parallel_compare_collective(): TIME SPENT IN MESSAGE PASSING", par_time);
    log_action_time(t->l, 0, "parallel_compare_collective(): time spent in memory management", mm_time);
    log_action_time(t->l, 0, "parallel_compare_collective(): time spent comparing and selecting best results", run_time);
    log_action_time(t->l, 0, "parallel_compare_collective(): execution time", getElapsedTime(&t1,&t7));
    t->par_time += par_time;
    t->mm_time += mm_time;
    t->corr_time += run_time;
#endif /* TSI_MPI */
    return 0;
} /* tsi_compare_parallel_collective */



int tsi_compare_parallel_direct(tsi *t) 
{
    //return tsi_compare_parallel_v3(t);
    return tsi_compare_parallel_v4(t);
    //return tsi_compare_parallel_direct_v5(t);
} /* tsi_compare_parallel_direct */


/* totally sincronous distributed compate & update */
int tsi_compare_parallel_v3(tsi *t)
{
#ifdef TSI_MPI
    int ret;
    int i;
    int to, from, start;
    float *smallBCM, *smallBAI;
    float *BAI, *BCM;

    double par_time, run_time;
    struct timeval t1, t2, t3;

    // MPI stuff
    int clustersize;
    int size_per_rank;
    int rest;
    MPI_Status  stat;


    getCurrTime(&t1);
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
    smallBCM = (float *) tsi_malloc( size_per_rank * sizeof(float));
    smallBAI = (float *) tsi_malloc( size_per_rank * sizeof(float));
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
    getCurrTime(&t2);
    run_time = getElapsedTime(&t1,&t2);
    par_time = 0;

    while(to != t->proc_id) {

        getCurrTime(&t3);
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

        getCurrTime(&t2);
        par_time += getElapsedTime(&t3, &t2);

        i = t->proc_id * size_per_rank;
        local_compare_update( BCM + i, BAI + i, smallBCM, smallBAI, size_per_rank);
        getCurrTime(&t3);
        run_time += getElapsedTime(&t2, &t3);

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
            getCurrTime(&t2);
            par_time += getElapsedTime(&t3, &t2);
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
            getCurrTime(&t2);
            par_time += getElapsedTime(&t3, &t2);
            local_compare_update( BCM + i, BAI + i, smallBCM, smallBAI, rest);
            getCurrTime(&t3);
            run_time += getElapsedTime(&t2, &t3);
        }

        /* step to next round */
        to = (to + 1) % clustersize;
        from = (from  + clustersize - 1) % clustersize;
        start = to * size_per_rank;
    } 

    getCurrTime(&t3);
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

    getCurrTime(&t2);
    par_time += getElapsedTime(&t3, &t2);
    tsi_free(smallBCM);
    tsi_free(smallBAI);
    log_action_time(t->l, 0, "parallel_compare_direct(): TIME SPENT IN MESSAGE PASSING", par_time);
    log_action_time(t->l, 0, "parallel_compare_direct(): time spent comparing and selecting best results", run_time);
    log_action_time(t->l, 0, "parallel_compare_direct(): execution time", getElapsedTime(&t1,&t2));
    t->par_time += par_time;
    t->corr_time += run_time;
    //	log_message(t->l,0,"tsi_compare_parallel_v3() finished");
    dirty_grid(t->heap, t->nextBAI_idx);
    dirty_grid(t->heap, t->nextBCM_idx);

#endif /* TSI_MPI */
    return 0;
} /* tsi_compare_parallel_v3 */



/* totally sincronous distributed compate & update,
 * based on _v3 but:
 * totally point-to-point, no collective operations.
 */
int tsi_compare_parallel_v4(tsi *t)
{
#ifdef TSI_MPI
    int ret;
    int i;
    int to, from, start;
    float *smallBCM, *smallBAI;
    float *BAI, *BCM;

    double par_time, run_time = 0;
    struct timeval t1, t2, t3;
    // MPI stuff
    int clustersize;
    int size_per_rank;
    int rest;
    MPI_Status  stat;


    //log_message(t->l, 0, "compare_parallel_v4: Called!");

    t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
    t->nextBCM_idx = new_grid(t->heap);
    t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
    expand_correlations_grid(t->nextBCM_c, t->nextBCM);
    delete_cmgrid(t->nextBCM_c);

    getCurrTime(&t1);

    BAI = t->nextBAI;
    BCM = t->nextBCM;

    MPI_Comm_size(MPI_COMM_WORLD, &clustersize);

    size_per_rank = t->grid_size / clustersize;
    rest = t->grid_size % clustersize;
    smallBCM = (float *) tsi_malloc( size_per_rank * sizeof(float));
    smallBAI = (float *) tsi_malloc( size_per_rank * sizeof(float));
    if(rest > 0)
        log_message(t->l, 0, "compare_parallel_v4: WE HAVE REST!");


    /* 
       circular send & recieve of parcels of BCM + BAI 
       we send the destination's parcel
       we recive _our_ parcel
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
        getCurrTime(&t2);
        local_compare_update( BCM + i, BAI + i, smallBCM, smallBAI, size_per_rank);
        getCurrTime(&t3);
        run_time += getElapsedTime(&t3, &t2);

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
            getCurrTime(&t2);
            local_compare_update( BCM + i, BAI + i, smallBCM, smallBAI, rest);
            getCurrTime(&t3);
            run_time += getElapsedTime(&t3, &t2);
        }

        /* step to next round */
        to = (to + 1) % clustersize;
        from = (from  + clustersize - 1) % clustersize;
        start = to * size_per_rank;
    } 

    tsi_free(smallBCM);
    tsi_free(smallBAI);


    /* now, update BAI & BCM in the same fashion.. */
    /* 
       circular send & recieve of parcels of BCM + BAI 
       we send _our_ parcel
       we recive the origin's parcel
       total DATA through network: grid_size * 2 * sizeof(float) 
       the last node in group also works the "rest" 
       */
    to = t->proc_id;
    from = t->proc_id;

    to = (to + 1) % clustersize;
    from = (from  + clustersize - 1) % clustersize;
    start = to * size_per_rank;

    while(to != t->proc_id) {

        i = t->proc_id * size_per_rank;
        ret = MPI_Sendrecv(BCM + i, size_per_rank, MPI_FLOAT, to, to,
                BCM + start, size_per_rank, MPI_FLOAT, from, t->proc_id,
                MPI_COMM_WORLD, &stat);
        if( ret != MPI_SUCCESS) {
            log_string(t->l, "tsi_compare_parallel(): mpi_sendrecv error!\n");
            return 1;
        }
        ret = MPI_Sendrecv(BAI + i, size_per_rank, MPI_FLOAT, to, to+1,
                BAI + start, size_per_rank, MPI_FLOAT, from, t->proc_id+1,
                MPI_COMM_WORLD, &stat);
        if( ret != MPI_SUCCESS) {
            log_string(t->l, "tsi_compare_parallel(): mpi_sendrecv error!\n");
            return 1;
        }


        if( to == (clustersize - 1) && rest != 0) { /* last gets the rest */
            i += size_per_rank;
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
            i = start + size_per_rank;
            ret = MPI_Recv(BCM + i, rest, MPI_FLOAT, from, t->proc_id, MPI_COMM_WORLD, &stat);
            if( ret != MPI_SUCCESS) {
                log_string(t->l, "tsi_compare_parallel(): mpi_recv error!\n");
                return 1;
            }
            ret = MPI_Recv(BAI + i, rest, MPI_FLOAT, from, t->proc_id+1, MPI_COMM_WORLD, &stat);
            if( ret != MPI_SUCCESS) {
                log_string(t->l, "tsi_compare_parallel(): mpi_recv error!\n");
                return 1;
            }
        }

        /* step to next round */
        to = (to + 1) % clustersize;
        from = (from  + clustersize - 1) % clustersize;
        start = to * size_per_rank;
    } 


    getCurrTime(&t2);
    par_time = getElapsedTime(&t1, &t2) - run_time;
    tsi_free(smallBCM);
    tsi_free(smallBAI);
    log_action_time(t->l, 0, "parallel_compare_direct(): TIME SPENT IN MESSAGE PASSING", par_time);
    log_action_time(t->l, 0, "parallel_compare_direct(): time spent comparing and selecting best results", run_time);
    log_action_time(t->l, 0, "parallel_compare_direct(): execution time", getElapsedTime(&t1,&t2));
    t->par_time += par_time;
    t->corr_time += run_time;
    //	log_message(t->l,0,"tsi_compare_parallel_v4() finished");
    dirty_grid(t->heap, t->nextBAI_idx);
    dirty_grid(t->heap, t->nextBCM_idx);

#endif /* TSI_MPI */
    return 0;
} /* tsi_compare_parallel_v4 */


/* based on tsi_compare_parallel_collective() but:
 * totally point-to-point, no collective operations.
 */
int tsi_compare_parallel_direct_v5(tsi *t) {
#ifdef TSI_MPI
    cm_grid *bcm;
    float   *bai, *ai_z;
    int     ai_z_idx;

    float *cc,          /* orginal BCM values */
          *rv;          /* results array */
    void  *recv_buf,    /* receive buffer */
          *send_buf;    /* send buffer */
    int   *nv;          /* locations array */

    unsigned int fragment_size,  /* size of each fragment of a grid to be shared */
                 cc_size;        /* size of compressed correlations grid */

    unsigned int n, g, h, i, z0, z1;        /* aux variables */
    unsigned int layer, last_layer;    
    struct timeval t1, t2, t3, t4, t5;
    double par_time, run_time;

    /* mpi info */
    int clustersize;
    MPI_Status  stat;

    getCurrTime(&t1);
    MPI_Comm_size(MPI_COMM_WORLD, &clustersize);

    /* distribute all values */
    bcm = t->nextBCM_c;
    load_cmgrid(bcm);
    cc = bcm->cg;

    cc_size = bcm->nxy * bcm->nlayers;
    fragment_size = (cc_size / t->n_procs) + ((cc_size % t->n_procs > 0) ? 1 : 0);
    cc_size = fragment_size * t->n_procs;  /* new "clustered" size */

    //--------------------------- TODO ----------------------------------------------------------------
    if (t->heap->grid_size < (2*cc_size + 2*fragment_size)) {
        printf_dbg("tsi_compare_parallel_direct(): not enough space available on BCM cm_grid...");
        return 1;
    }
    //--------------------------- END OF TODO ---------------------------------------------------------

    send_buf = cc;
    recv_buf = cc + cc_size + fragment_size;
    nv = (int *) cc + cc_size;          /* set pointer for nodes array */
    rv = recv_buf;                      /* set pointer for results array (same as recv_buf) */

    printf_dbg("tsi_compare_parallel_direct(): performing all to all\n");
    getCurrTime(&t2);
    /* this code implements ALLTOALL using point to point paired sendrecvs 
     * this is specially nice on clusters with fullduplex point to point links and no (or bad) multicast support
     */

    int to, from, start, ret;
    to = t->proc_id;
    from = t->proc_id;

    to = (to + 1) % clustersize;
    from = (from  + clustersize - 1) % clustersize;
    start = to * fragment_size;

    while(to != t->proc_id) {

        i = from * fragment_size;
        ret = MPI_Sendrecv(send_buf + start, fragment_size, MPI_FLOAT, to, to,
                recv_buf + i, fragment_size, MPI_FLOAT, from, t->proc_id,
                MPI_COMM_WORLD, &stat);
        if( ret != MPI_SUCCESS) {
            log_string(t->l, "tsi_compare_parallel(): mpi_sendrecv error!\n");
            return 1;
        }

        /* step to next round */
        to = (to + 1) % clustersize;
        from = (from  + clustersize - 1) % clustersize;
        start = to * fragment_size;
    } 
    /* we must copy our part of the cc grid to our block in "rv".. */
    for(i = 0; i < fragment_size; i++) {
        g = (t->proc_id * fragment_size) + i;
        rv[i] = cc[g];
    }


    getCurrTime(&t3);
    par_time = getElapsedTime(&t2, &t3);

    /* compare and select best values */
    printf_dbg("tsi_compare_parallel_direct(): begining local compare\n");

    for (n = 1; n < t->n_procs; n++) {
        for (g = 0; g < fragment_size; g++) {
            h = n * fragment_size + g;
            if (rv[g] < rv[h]) {
                rv[g] = rv[h];
                nv[g] = n;
            }
        }
    }

    /* distribute results lists */
    getCurrTime(&t2);
    run_time = getElapsedTime(&t3, &t2);

    /* this code implements ALLGATHER using point to point paired sendrecvs 
     * this is specially nice on clusters with fullduplex point to point links and no (or bad) multicast support
     */
    to = t->proc_id;
    from = t->proc_id;

    to = (to + 1) % clustersize;
    from = (from  + clustersize - 1) % clustersize;

    printf_dbg("tsi_compare_parallel_direct(): performing gather of new BCM\n");
    send_buf = rv;    /* get new BCM */
    recv_buf = cc;
    start = from * fragment_size;
    while(to != t->proc_id) {

        ret = MPI_Sendrecv(send_buf, fragment_size, MPI_FLOAT, to, to,
                recv_buf + start, fragment_size, MPI_FLOAT, from, t->proc_id,
                MPI_COMM_WORLD, &stat);
        if( ret != MPI_SUCCESS) {
            log_string(t->l, "tsi_compare_parallel(): mpi_sendrecv error!\n");
            return 1;
        }

        /* step to next round */
        to = (to + 1) % clustersize;
        from = (from  + clustersize - 1) % clustersize;
        start = from * fragment_size;
    } 
    /* we must copy our part of the cc grid back to cc grid .. */
    for(i = 0; i < fragment_size; i++) {
        g = (t->proc_id * fragment_size) + i;
        cc[g] = rv[i];
    }

    /* this code implements ALLGATHER using point to point paired sendrecvs 
     * this is specially nice on clusters with fullduplex point to point links and no (or bad) multicast support
     */

    to = t->proc_id;
    from = t->proc_id;

    to = (to + 1) % clustersize;
    from = (from  + clustersize - 1) % clustersize;

    printf_dbg("tsi_compare_parallel_direct(): performing gather of new best AI locations\n");
    send_buf = nv;
    recv_buf = rv;
    start = from * fragment_size;
    while(to != t->proc_id) {

        ret = MPI_Sendrecv(send_buf, fragment_size, MPI_FLOAT, to, to,
                recv_buf + start, fragment_size, MPI_FLOAT, from, t->proc_id,
                MPI_COMM_WORLD, &stat);
        if( ret != MPI_SUCCESS) {
            log_string(t->l, "tsi_compare_parallel(): mpi_sendrecv error!\n");
            return 1;
        }

        /* step to next round */
        to = (to + 1) % clustersize;
        from = (from  + clustersize - 1) % clustersize;
        start = from * fragment_size;
    } 
    /* copy our part */
    for(i = 0; i < fragment_size; i++) {
        rv[i] = nv[i];
    }

    getCurrTime(&t3);
    par_time += getElapsedTime(&t2, &t3);
    nv += fragment_size;

    /* distribute best AI values */
    bai = load_grid(t->heap, t->nextBAI_idx);
    ai_z_idx = new_grid(t->heap);
    if (ai_z_idx < 0) {
        printf_dbg("tsi_compare_parallel_direct(): failed to allocate ai_z grid\n");
        return 1;
    }
    ai_z = load_grid(t->heap, ai_z_idx);
    if (ai_z == NULL) {
        printf_dbg("tsi_compare_parallel_direct(): failed to load ai_z grid\n");
        return 1;
    }
    cc_size = bcm->nxy * bcm->nlayers;


    /* lets calculate its size and create the array of our best AI values */
    int send_size = 0;
    z0 = 0;
    z1 = bcm->layer_size[0];
    layer = last_layer = 0;
    for (g = 0; g < cc_size; g++) {        /* build z vectors AI array */
        layer = g / bcm->nxy;
        if (layer > last_layer) { 
            last_layer = layer;
            z0 = z1;
            z1 += bcm->layer_size[layer];
            printf_dbg("parallel_compare(%d): layer: %d, z0: %d, z1:%d\n", t->proc_id, layer, z0, z1);            
        }
        if (nv[g] == t->proc_id) {
            for (i = z0; i < z1; i++) {
                ai_z[send_size] = bai[(i * bcm->nxy) + (g - layer*bcm->nxy)];
                send_size++;
            } /* for */
        }
    } /* for */



    getCurrTime(&t2);
    run_time += getElapsedTime(&t2,&t3);
    /* this code implements ALLGATHER using point to point paired sendrecvs 
     * this is specially nice on clusters with fullduplex point to point links and no (or bad) multicast support
     */

    to = t->proc_id;
    from = t->proc_id;

    to = (to + 1) % clustersize;
    from = (from  + clustersize - 1) % clustersize;

    send_buf = ai_z;
    recv_buf = ai_z + send_size;
    int from_size, total = send_size;
    double r_time;
    r_time = 0;
    printf_dbg("parallel_compare(%d): send_size: %d\n", t->proc_id, send_size);            
    while(to != t->proc_id) {

        from_size = 0;
        for (g = 0; g < cc_size; g++) {   /* calc array size to be received */
            if (nv[g] == from) {
                from_size += bcm->layer_size[g/bcm->nxy];
            }
        } /* for */

        printf_dbg("parallel_compare(%d): sendrecv, send_size: %d\t got so far: %d\t next recv_size: %d from %d\n", t->proc_id, send_size, total, from_size, from);
        ret = MPI_Sendrecv(send_buf , send_size, MPI_FLOAT, to, to,
                recv_buf, from_size, MPI_FLOAT, from, t->proc_id,
                MPI_COMM_WORLD, &stat);

        getCurrTime(&t4);
        if(from_size > 0) {
            h = 0;
            z0 = 0;
            z1 = bcm->layer_size[0];
            layer = last_layer = 0;
            for (g = 0; g < cc_size; g++) {        /* unpack z vectors AI array */
                layer = g / bcm->nxy;
                if (layer > last_layer) { 
                    last_layer = layer;
                    z0 = z1;
                    z1 += bcm->layer_size[layer];
                    printf_dbg("parallel_compare(%d): layer: %d, z0: %d, z1:%d\n", t->proc_id, layer, z0, z1);            
                }
                if (nv[g] == from) {
                    for (i = z0; i < z1; i++) {
                        bai[(i * bcm->nxy) + (g - layer*bcm->nxy)] = ai_z[send_size + h];
                        h++;
                    }
                }
            } /* for */
        }
        total += from_size;

        getCurrTime(&t5);
        r_time += getElapsedTime(&t5, &t4);

        /* step to next round */
        to = (to + 1) % clustersize;
        from = (from  + clustersize - 1) % clustersize;
    } 
    run_time += r_time;

    getCurrTime(&t3);

    /* clear grids */
    dirty_cmgrid(bcm);
    dirty_grid(t->heap, t->nextBAI_idx);
    delete_grid(t->heap, ai_z_idx);

    par_time += getElapsedTime(&t2, &t3) - r_time;
    log_action_time(t->l, 0, "parallel_compare_direct(): TIME SPENT IN MESSAGE PASSING", par_time);
    log_action_time(t->l, 0, "parallel_compare_direct(): execution time", run_time);
    t->par_time += par_time;
    t->corr_time += run_time;

    /* this is done here because the other compare_parallel_direct() don't use compressed CorrMaps */
    t->nextBCM_idx = new_grid(t->heap);
    t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
    load_cmgrid(t->nextBCM_c);
    expand_correlations_grid(t->nextBCM_c, t->nextBCM);
    delete_cmgrid(t->nextBCM_c);
    dirty_grid(t->heap, t->nextBCM_idx);

#endif /* TSI_MPI */
    return 0;
} /* tsi_compare_parallel_v5*/

void local_compare_update(float *lBCM, float *lBAI, float *rBCM, float *rBAI, int size)
{
    int i;
    for(i = 0; i < size; i++) {
        if(lBCM[i] < rBCM[i]) { // new best correlation 
            lBCM[i] = rBCM[i];
            lBAI[i] = rBAI[i];
        }
    }
} /* local_compare_update */

/* end of tsi_parallel.c */
