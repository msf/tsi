#ifndef TSI_MPI
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
#include "dss_wrapper.h"
#include "si_wrapper.h"
#include "tsi.h"

#define HEAP_SIZE 15

tsi *new_tsi(registry *reg) {
    ////////////////////////////////////////////////// FINISH
    tsi *t;
    reg_key *k;
    int n_grids;
    char *filename;

    t = (tsi *) malloc(sizeof(tsi));
    t->reg = reg;
    t->heap = NULL;

    /* get initial data from registry */
    k = get_key(reg, "GLOBAL", "ITERATIONS");
    if (k) {
       t->iterations = get_int(k);
       if (t->iterations < 1) {
           printf("Incoeherent number of iterations!");
           delete_tsi(t);
           return NULL;
       }
    } else {
       printf_dbg("new_tsi(): failed to get number of iterations from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GLOBAL", "SIMULATIONS");
    if (k) {
       t->simulations = get_int(k);
       if (t->simulations < 1) {
           printf("Incoeherent number of simulations!");
           delete_tsi(t);
           return NULL;
       }
    } else {
       printf_dbg("new_tsi(): failed to get number of simulations from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GRID", "XNUMBER");
    if (k)
       t->xsize = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get size of X from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GRID", "YNUMBER");
    if (k)
       t->ysize = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get size of Y from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GRID", "ZNUMBER");
    if (k)
       t->zsize = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get size of Z from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    t->grid_size = (unsigned int)t->zsize * (unsigned int)t->ysize * (unsigned int)t->xsize;

    /* get heap data */
    k = get_key(reg, "GLOBAL", "USEFS");
    if (k)
       t->usefs = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get usefs flag from the registry!\n");
       delete_tsi(t);
       return NULL;
    }

    t->heap = new_heap(HEAP_SIZE, 1, 0, t->usefs, t->xsize, t->ysize, t->zsize);
    if (!t->heap) {
       printf_dbg("new_tsi(): failed to start heap\n");
       delete_tsi(t);
       return NULL;
    }

    /* init grids */
    if ((t->seismic_idx = new_grid(t->heap)) < 0) return 0;
    if ((t->bestAI_idx  = new_grid(t->heap)) < 0) return 0;
    if ((t->currBAI_idx = new_grid(t->heap)) < 0) return 0;
    if ((t->currBCM_idx = new_grid(t->heap)) < 0) return 0;
    if ((t->nextBAI_idx = new_grid(t->heap)) < 0) return 0;
    if ((t->nextBCM_idx = new_grid(t->heap)) < 0) return 0;
    if ((t->ai_idx = new_grid(t->heap)) < 0) return 0;
    if ((t->cm_idx = new_grid(t->heap)) < 0) return 0;
    if ((t->sy_idx = new_grid(t->heap)) < 0) return 0;

    t->dss_eng = new_dss(t->reg, t->heap);
    if (!t->dss_eng) {
       printf_dbg("new_tsi(): failed to start dss engine\n");
       delete_tsi(t);
       return NULL;
    }
    t->si_eng = new_si(t->reg, t->heap);
    if (!t->si_eng) {
       printf_dbg("new_tsi(): failed to start si engine\n");
       delete_tsi(t);
       return NULL;
    }

/*
    k = get_key(reg, "SEISMIC", "FILENAME");
    if (k)
       filename = get_string(k);
    else {
    }
    t->seismic_idx = new_grid(t->heap);
    if (t->seismic_idx) {
        t->seismic = load_grid(t->heap, t->seismic_idx);
        if (load_ascii_grid_from_file(t->seismic, filename)) {
            dirty_grid(t->heap, t->seismic_idx);
        } else {
            printf_dbg("Failed to load seismic file!\n");
            delete_tsi(t);
            return NULL;
        }
    } else {
       printf_dbg("Failed to start the seismic grid!\n");
       delete_tsi(t);
       return NULL;
    }
*/

    return t;
} /* new_tsi */



void delete_tsi(tsi *t) {
    ////////////////////////////////////////////////// TODO
} /* delete_tsi */



float grid_correlation (float *a, float *b) {
    ////////////////////////////////////////////////// TODO
    return 0;
} /* grid_correlation */



int tsi_compare (tsi *t) {
    ////////////////////////////////////////////////// TODO
    t->ai = load_grid(t->heap, t->ai_idx);
    t->cm = load_grid(t->heap, t->cm_idx);
    t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
    t->nextBCM = load_grid(t->heap, t->nextBCM_idx);

    /* execute Compare */    

    dirty_grid(t->heap, t->nextBAI_idx);
    dirty_grid(t->heap, t->nextBCM_idx);
    clear_grid(t->heap, t->ai_idx);
    clear_grid(t->heap, t->cm_idx);
    return 1;
} /* tsi_compare */


void grid_copy(float *a, float *b, unsigned int size)
{
    memcpy(a, b, size*sizeof(float));
}


int run_tsi(tsi *t) {
    int result;
    int i, j;
    float best_corr, ai_corr;
    
    best_corr = ai_corr = 0;

    /* FIRST ITERATION */
    printf_dbg("run_tsi(): first iteration\n");
    for (i = 0; i < t->simulations; i++) {

        /* DSS */
        printf_dbg("run_tsi(): loading AI grid\n");
        t->ai = load_grid(t->heap, t->ai_idx);
        if (t->ai) {
            printf_dbg("run_tsi(): running DSS\n");
            result = run_dss(t->dss_eng, t->ai);         /* 6 grids */
        } else {
            printf_dbg("run_tsi(): failed to load AI for DSS\n");
            return 0;
        }

        /* SI */
        if (result) {
            printf_dbg("run_tsi(): loading CM\n");
            t->cm = load_grid(t->heap, t->cm_idx);
            printf_dbg("run_tsi(): loading seismic\n");
            t->seismic = load_grid(t->heap, t->seismic_idx);
            printf_dbg("run_tsi(): loading SY\n");
            t->sy = load_grid(t->heap, t->sy_idx);
            if (t->cm && t->seismic && t->sy) {
                printf_dbg("run_tsi(): running SI \n");
                result = si_simulation(t->si_eng, t->ai, t->seismic, t->cm, t->sy);  /* 4 grids */
            } else {
                printf_dbg("run_tsi(): failed to load CM, Seismic or SY for SI!\n");
                return 0;
            }
        } else {
            printf_dbg("run_tsi(): DSS failed!\n");
            return 0;
        }


        /* Is best? */                                 /* 4 grids */
        if (result) {
            printf_dbg("run_tsi(): checking best global correlation\n");
            ai_corr = grid_correlation(t->seismic, t->sy);
            if (ai_corr > best_corr) {
                best_corr = ai_corr;
                delete_grid(t->heap, t->bestAI_idx);   /* new bestAI, delete old */
                t->bestAI_idx = t->ai_idx;
            }
            clear_grid(t->heap, t->seismic_idx);
            clear_grid(t->heap, t->sy_idx);
        } else {
            printf_dbg("run_tsi(): SI failed!\n");
            return 0;
        }
        
        /* Compare */
        printf_dbg("run_tsi(): running Compare\n");
        result = tsi_compare(t);   /* builds nextBAI and nextBCM */     /* 4 grids */
        if (result) {
            printf_dbg("run_tsi(): final clean-up\n");
            if (t->ai_idx == t->bestAI_idx) {
                dirty_grid(t->heap, t->bestAI_idx);
                t->ai_idx = new_grid(t->heap);
            }
            clear_grid(t->heap, t->cm_idx);
        } else {
            printf_dbg("run_tsi(): Compare failed!\n");
            return 0;
        }
    } /* for */
    clear_grid(t->heap, t->ai_idx);

    /* NEXT ITERATIONS */
    for (j = 1; j < t->iterations; j++) {

        printf_dbg("run_tsi(): interation %d\n", j+1);
        t->currBAI = load_grid(t->heap, t->currBAI_idx);
        t->currBCM = load_grid(t->heap, t->currBCM_idx);
        t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
        t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
        grid_copy(t->currBAI, t->nextBAI, t->grid_size);
        grid_copy(t->currBCM, t->nextBCM, t->grid_size);
        clear_grid(t->heap, t->nextBAI_idx);
        clear_grid(t->heap, t->nextBCM_idx);

        for (i = 0; i < t->simulations; i++) {

            /* CoDSS */
            printf_dbg("run_tsi(): loading AI grid\n");
            t->ai = load_grid(t->heap, t->ai_idx);
            if (t->ai && t->currBAI && t->currBCM) {
                result = codss_simulation(t->dss_eng, t->currBAI, t->currBCM, t->ai);
                clear_grid(t->heap, t->currBAI_idx);
                clear_grid(t->heap, t->currBCM_idx);
            } else {
                printf_dbg("run_tsi(): failed to load AI, currBAI or currBCM for CoDSS!\n");
                return 0;
            }

            /* SI */
            if (result) {
                t->cm = load_grid(t->heap, t->cm_idx);
                t->seismic = load_grid(t->heap, t->seismic_idx);
                if (t->cm && t->seismic) {
                    printf_dbg("run_tsi(): running CoDSS\n");
                    result = si_simulation(t->si_eng, t->ai, t->seismic, t->cm, t->sy);
                    clear_grid(t->heap, t->seismic_idx);
                } else {
                    printf_dbg("run_tsi(): failed to load CM, Seismic or SY for SI!\n");
                    return 0;
                }
            } else {
                printf_dbg("run_tsi(): CoDSS simulation failed!\n");
                return 0;
            }

            /* Is best? */
            if (result) {
                ai_corr = grid_correlation(t->seismic, t->sy);
                if (ai_corr > best_corr) {
                    delete_grid(t->heap, t->bestAI_idx);
                    clear_grid(t->heap, t->sy_idx);
                    t->bestAI_idx = t->ai_idx;
                }
            } else {
                printf_dbg("run_tsi(): si simulation failed!\n");
                return 0;
            }

            /* Compare */
            result = tsi_compare(t);   /* builds nextBAI and nextBCM */
            if (result) {
                if (t->ai_idx == t->bestAI_idx) {
                    dirty_grid(t->heap, t->bestAI_idx);
                    t->ai_idx = new_grid(t->heap);
                    t->ai = load_grid(t->heap, t->ai_idx);
                }
            } else {
                printf_dbg("run_tsi(): Compare failed!\n");
                return 0;
            }
        } /* for */
        
    } /* for */
    
    printf_dbg("run_tsi(): heap performance: R=%d W=%d G=%d\n", t->heap->reads, t->heap->writes, t->heap->curr_grids);
    return 1;
} /* run_tsi */



#endif /* TSI_MPI */

/* end of file */
