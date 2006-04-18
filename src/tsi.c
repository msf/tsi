#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "debug.h"
#include "memdebug.h"
#include "registry.h"
#include "grid_heap.h"
#include "dss.h"
#include "si.h"
#include "tsi.h"
#include "tsi_math.h"
#include "tsi_parallel.h"
#include "timer.h"


/* additional prototypes to TSI */
int   tsi_compare(unsigned int size, float *AI, float *CM, float *nextBAI, float *BCM);
void  grid_copy(float *a, float *b, unsigned int grid_size);



tsi *new_tsi(registry *reg) {
    tsi *t;
    reg_key *k;
    int usefs, heap_size, swap_thr, new_sims, over_sims, over_procs, i;
    TSI_FILE *fp;

    t = (tsi *) tsi_malloc(sizeof(tsi));
    if (!t) return NULL;
    t->reg = reg;
    t->heap = NULL;
    t->l = new_log(NULL);    /************* FINISH **************/

    /* get machine parameters */
    if (new_tsi_parallel(&t->n_procs, &t->proc_id) < 0) {
        printf("Machine failed to start!\n");
        return NULL;
    }
    t->root_id = t->proc_id;

    t->global_best.value = -999;
    t->global_best.proc_id = t->proc_id;
    
    /* get initial data from registry */
    if ((k = get_key(reg, "GLOBAL", "ITERATIONS")) == NULL) {
       delete_tsi(t);
       return NULL;
    }
    t->iterations = get_int(k);
    if (t->iterations < 1) {
        printf("Incoeherent number of iterations!\n");
        delete_tsi(t);
        return NULL;
    }
    if ((k = get_key(reg, "GLOBAL", "SIMULATIONS")) == NULL) {
       delete_tsi(t);
       return NULL;
    }
    t->simulations = get_int(k);
    if (t->simulations < 1) {
        printf("Incoeherent number of simulations!\n");
        delete_tsi(t);
        return NULL;
    }

	printf_dbg("TSI: %d Simulations x %d Iterations\n",t->simulations, t->iterations);
    k = get_key(reg, "GLOBAL", "OPTIMIZE");
    if (k) {
       t->optimize = get_int(k);
    } else {
       printf_dbg("new_tsi(%d): failed to get optimize flag from the registry! Using defaults...\n", t->proc_id);
       t->optimize = 1;
    }

    /* change simulation data to optimize resource use during execution */
    printf_dbg("new_tsi(%d): number of processes=%d\n", t->proc_id, t->n_procs);
    if (t->n_procs > 1) {
        new_sims = 0;
        if (t->simulations > t->n_procs) new_sims++;
        new_sims += t->simulations / t->n_procs;  /* initial set of sims for each process */
        over_sims = new_sims * t->n_procs - t->simulations;   /* number of exceeding simulations */

        /* calculate the optimal number of CPUs based on the number of simulations */
        if (t->optimize) {
            over_procs = over_sims / new_sims;
            t->n_procs -= over_procs;
        } else {
            if (new_sims == 1) t->n_procs = t->simulations;
        }
        over_sims = new_sims * t->n_procs - t->simulations;
        printf_dbg("new_tsi(%d): new number of processes=%d\n", t->proc_id, t->n_procs);

        /* terminate processes that aren't needed */
        if (t->n_procs <= t->proc_id) {
            exit(0);
            //return NULL;
        }
        
        /* correct the number of simulations to be executed */
        t->simulations = new_sims;
        if (t->proc_id < over_sims) t->simulations--;

        /* set root node */
        k = get_key(reg, "GLOBAL", "MIDDLEZERO");
        if (k) {
            t->root_id = get_int(k);
        } else {
            printf_dbg("new_tsi(%d): failed to get optimize flag from the registry! Using defaults...\n", t->proc_id);
            t->root_id = 1;
        }
        if (t->root_id) t->root_id = t->n_procs / 2;
    }
    printf_dbg("new_tsi(%d): number of simulations=%d\n", t->proc_id, t->simulations);


    /* get heap data */
    k = get_key(reg, "HEAP", "USEFS");
    if (k)
       usefs = get_int(k);
    else {
       printf_dbg("new_tsi(%d): failed to get usefs flag from the registry! Using defaults...\n", t->proc_id);
       usefs = 0;
    }
    k = get_key(reg, "HEAP", "SIZE");
    if (k)
       heap_size = get_int(k);
    else {
       printf_dbg("new_tsi(%d): failed to get heap size flag from the registry! Using defaults...\n", t->proc_id);
       heap_size = 14;
    }
    k = get_key(reg, "HEAP", "THRESHOLD");
    if (k)
       swap_thr = get_int(k);
    else {
       printf_dbg("new_tsi(%d): failed to get swap threshold from the registry! Using defaults...\n", t->proc_id);
       swap_thr = 8;
    }

    /* get correlation data */
    k = get_key(reg, "CORR", "BCM_ROOT");
    if (k)
       t->root = get_int(k);
    else {
       printf_dbg("new_tsi(%d): failed to get BCM root from the registry! Using defaults...\n", t->proc_id);
       t->root = 1;
    }

    /* get grid parameters */
    k = get_key(reg, "GRID", "XNUMBER");
    if (k)
       t->xsize = get_int(k);
    else {
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GRID", "YNUMBER");
    if (k)
       t->ysize = get_int(k);
    else {
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GRID", "ZNUMBER");
    if (k)
       t->zsize = get_int(k);
    else {
       delete_tsi(t);
       return NULL;
    }
    t->even_size = t->grid_size = (unsigned int)t->zsize * (unsigned int)t->ysize * (unsigned int)t->xsize;
    if (t->grid_size % 2) t->even_size++;

	printf_dbg("TSI: grid size is %d\n",t->grid_size);

    /* start grid heap */
    printf_dbg("new_tsi(%d): starting heap\n", t->proc_id);
    t->heap = new_heap(t->n_procs, t->proc_id, heap_size, swap_thr, usefs, t->even_size);
    if (!t->heap) {
       printf_dbg("new_tsi(%d): failed to start heap\n", t->proc_id);
       delete_tsi(t);
       return NULL;
    }

    /* start DSS engine */
    t->dss_eng = new_dss(t->reg, t->heap);
    if (!t->dss_eng) {
       printf_dbg("new_tsi(%d): failed to start dss engine\n", t->proc_id);
       delete_tsi(t);
       return NULL;
    }
	printf_dbg("new_tsi(): DSS engine loaded\n");
    
    /* start SI engine */
    t->si_eng = new_si(t->reg, t->heap);
    if (!t->si_eng) {
       printf_dbg("new_tsi(%d): failed to start si engine\n", t->proc_id);
       delete_tsi(t);
       return NULL;
    }
	printf_dbg("new_tsi(): SI engine loaded\n");


    /* load seismic grid */
    if ((t->seismic_idx = new_grid(t->heap)) < 0) {
        printf_dbg("new_tsi(%d): failed to allocate seismic grid\n", t->proc_id);
        delete_tsi(t);
        return NULL;
    }
    t->seismic = load_grid(t->heap, t->seismic_idx);
    

    if ((k = get_key(reg, "SEISMIC", "FILENAME")) == NULL) {
        printf_dbg("new_tsi(%d): failed to allocate seismic grid\n", t->proc_id);
        delete_tsi(t);
        return NULL;
    }
    if ((fp = open_file(get_string(k))) == NULL) {
       printf_dbg("Failed to open the seismic grid file!\n");
       delete_tsi(t);
       return NULL;
    }
    if (!read_ascii_grid_file(fp, t->seismic, t->grid_size)) {
        printf_dbg("Failed to load seismic file!\n");
        delete_tsi(t);
        return NULL;
    }
	printf_dbg("new_tsi(): Seismic Data loaded\n");
    dirty_grid(t->heap, t->seismic_idx);

    return t;
} /* new_tsi */



void delete_tsi(tsi *t) {
    if (t) {
        if (t->si_eng) delete_si(t->si_eng);
        if (t->dss_eng) delete_dss(t->dss_eng);
        if (t->heap) delete_heap(t->heap);
        if (t->reg) delete_registry(t->reg);
        tsi_free(t);
    }
    delete_tsi_parallel();
} /* delete_tsi */



int run_tsi(tsi *t) {
    int result;
    int i, j;
    float ai_corr;
    struct timeval t1,t2;
    TSI_FILE *fp;
    cm_grid *cmg;

    ai_corr = -1;

    printf_dbg("run_tsi(%d): first iteration, first_simulation\n", t->proc_id);

    /* DSS 0/0 */
    setup_dss(t->dss_eng, NULL);
    printf_dbg2("run_tsi(%d): loading AI grid 0/0\n", t->proc_id);
    if ((t->ai_idx = new_grid(t->heap)) < 0) {
        printf_dbg("new_tsi(%d): failed to allocate AI grid 0/0\n", t->proc_id);
        delete_tsi(t);
        return 0;
    }
    t->ai = load_grid(t->heap, t->ai_idx);
	printf_dbg("\nrun_tsi(%d): (0/0) DSS is starting \n",t->proc_id);
    if (t->ai) {
        printf_dbg("run_tsi(%d): running DSS\n", t->proc_id);
        getCurrTime(&t1);
        result = run_dss(t->dss_eng, t->ai);
        getCurrTime(&t2);
    } else {
        printf_dbg("run_tsi(%d): failed to load AI for DSS 0/0\n", t->proc_id);
        return 0;
    }

    printf_dbg("run_tsi(%d): (0/0) DSS terminated (%f secs)\n\n",t->proc_id,getElapsedTime(&t1,&t2));
    

    /* SI 0/0 */
    if (result) {
        /* random layers */
        getCurrTime(&t1);
        if (t->proc_id == t->root_id) { // MPI
            cmg = new_cmgrid(t->si_eng, 0);    /* generate new set of layers */
            printf("New set of layers:\n");
        } else {
            cmg = new_cmgrid(t->si_eng, 1);    /* return an empty cm_grid */
        }
		
        if (tsi_set_layers_parallel(t, cmg)) {
            printf_dbg("run_tsi(%d): failed to set layers\n", t->proc_id);
            delete_tsi(t);
            return 0;
		}
		print_layers(cmg);
        store_cmgrid(t->si_eng, cmg);  /* adds CM/C grid to SI engine */
        getCurrTime(&t2);
        printf_dbg("run_tsi(%d): New set of layers took %f secs\n\n",t->proc_id,getElapsedTime(&t1,&t2));

        /* load data and grids for SI */
        printf_dbg2("run_tsi(%d): loading seismic 0/0\n", t->proc_id);
        t->seismic = load_grid(t->heap, t->seismic_idx);
        if ((t->cm_idx = new_grid(t->heap)) < 0) {
            printf_dbg("run_tsi(%d): failed to allocate CM grid 0/0\n", t->proc_id);
            delete_tsi(t);
            return 0;
        }
        if ((t->sy_idx = new_grid(t->heap)) < 0) {
            printf_dbg("run_tsi(%d): failed to allocate SY grid 0/0\n", t->proc_id);
            delete_tsi(t);
            return 0;
        }
        printf_dbg2("run_tsi(%d): loading CM 0/0\n", t->proc_id);
        t->cm = load_grid(t->heap, t->cm_idx);
        printf_dbg2("run_tsi(%d): loading SY 0/0\n", t->proc_id);
        t->sy = load_grid(t->heap, t->sy_idx);

        /* run Seismic Inversion if all went well */
        if (t->cm && t->seismic && t->sy) {
            printf_dbg2("run_tsi(%d): (0/0) running SI\n", t->proc_id);
            getCurrTime(&t1);
            result = run_si(t->si_eng, t->ai, t->seismic, t->cm, t->sy);  /* 4 grids */
            getCurrTime(&t2);
        } else {
            printf_dbg("run_tsi(%d): failed to load CM, Seismic or SY for SI! 0/0\n", t->proc_id);
            return 0;
        }
    } else {
        printf_dbg("run_tsi(%d): (0/0) DSS failed!\n", t->proc_id);
        return 0;
    }
    printf_dbg("run_tsi(%d): (0/0) SI execution took %f secs\n",t->proc_id,getElapsedTime(&t1,&t2));

    /* first isBest/Compare (optimal execution) 1/1 */
    getCurrTime(&t1);
    t->global_best.value = grid_correlation(t->sy, t->seismic, t->grid_size);    /* first best correlation */
    printf("\n\t\t (0/0) Current best correlation = %f\n\n", t->global_best.value);
    clear_grid(t->heap, t->seismic_idx);

    grid_copy(t->sy, t->ai, t->grid_size);     /* nextBAI = bestAI = AI at this stage */
    t->bestAI_idx = t->ai_idx;
    t->nextBAI_idx = t->sy_idx;
    dirty_grid(t->heap, t->bestAI_idx);
    dirty_grid(t->heap, t->nextBAI_idx);

    t->nextBCM_idx = t->cm_idx;  /* nextBCM = CM at this stage */
    dirty_grid(t->heap, t->nextBCM_idx);
    t->ai_idx = t->cm_idx = t->sy_idx = -1;  /* free aux grids */
    getCurrTime(&t2);
    printf_dbg("run_tsi(%d): Compare execution took %f secs\n",t->proc_id,getElapsedTime(&t1,&t2));

    

    printf_dbg("run_tsi(%d): first iteration, next simulations\n", t->proc_id);
    /* FIRST ITERATION, NEXT SIMULATIONS */
    for (i = 1; i < t->simulations; i++) {

        /* DSS 1/n */
        printf_dbg2("run_tsi(%d): loading AI grid 1/n\n", t->proc_id);
        if ((t->ai_idx = new_grid(t->heap)) < 0) {
            printf_dbg("run_tsi(%d): failed to allocate AI grid 1/n\n", t->proc_id);
            delete_tsi(t);
            return 0;
        }
        printf_dbg2("run_tsi(%d): loading AI grid 1/n\n", t->proc_id);
        t->ai = load_grid(t->heap, t->ai_idx);
        if (t->ai) {
            printf_dbg("run_tsi(%d): (0/%d) running DSS!\n", t->proc_id,i);
            getCurrTime(&t1);
            result = run_dss(t->dss_eng, t->ai);
            getCurrTime(&t2);
        } else {
            printf_dbg("run_tsi(%d): failed to load AI for DSS 1/n\n", t->proc_id);
            return 0;
        }
        printf_dbg("run_tsi(%d): (0/%d) DSS terminated (%f secs)\n",t->proc_id,i,getElapsedTime(&t1,&t2));

        /* SI 1/n */
        if (result) {
			/* load data and grids */
            if ((t->cm_idx = new_grid(t->heap)) < 0) {
                printf_dbg("run_tsi(%d): failed to allocate CM grid 1/n\n", t->proc_id);
                delete_tsi(t);
                return 0;
            }
            if ((t->sy_idx = new_grid(t->heap)) < 0) {
                printf_dbg("run_tsi(%d): failed to allocate SY grid 1/n\n", t->proc_id);
                delete_tsi(t);
                return 0;
            }
            printf_dbg2("run_tsi(): loading CM\n");
            t->cm = load_grid(t->heap, t->cm_idx);
            printf_dbg2("run_tsi(): loading seismic\n");
            t->seismic = load_grid(t->heap, t->seismic_idx);
            printf_dbg2("run_tsi(): loading SY\n");
            t->sy = load_grid(t->heap, t->sy_idx);

			/* run Seismic Inversion if all went well */
            if (t->cm && t->seismic && t->sy) {
                printf_dbg("run_tsi(): (0/%d) running SI \n",i);
                getCurrTime(&t1);
                result = run_si(t->si_eng, t->ai, t->seismic, t->cm, t->sy);
                getCurrTime(&t2);
            } else {
                printf_dbg("run_tsi(): failed to load CM, Seismic or SY for SI!\n");
                return 0;
            }
        } else {
            printf_dbg("run_tsi(): (0/%d) DSS failed!\n",i);
            return 0;
        }
        printf_dbg("run_tsi(%d): (0/%d) SI execution took %f secs\n",t->proc_id,i,getElapsedTime(&t1,&t2));

        /* Is best? 1/n */
        getCurrTime(&t1);
        if (result) {
            printf_dbg("run_tsi(): checking best global correlation\n");
            ai_corr = grid_correlation(t->seismic, t->sy, t->grid_size);   /* between real and synthetic seismic data */
			printf_dbg("isBest(0/%d): grid correlation: %f\n",i,ai_corr);
            if (ai_corr > t->global_best.value) {
                t->global_best.value = ai_corr;
                delete_grid(t->heap, t->bestAI_idx);   /* new bestAI, delete old */
                t->bestAI_idx = t->ai_idx;
                printf("\n\t\t\t(t/%d) NEW BEST CORRELATION = %f\n\n",i, t->global_best.value);
            }
            clear_grid(t->heap, t->seismic_idx);
        } else {
            printf_dbg("run_tsi(): (0/%d) SI failed!\n",i);
            return 0;
        }
        
        /* Compare 1/n */
        printf_dbg("run_tsi(): (0/%d) running Compare\n",i);
        t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
        t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
        if (t->nextBAI && t->nextBCM)
            result = tsi_compare(t->grid_size, t->ai, t->cm, t->nextBAI, t->nextBCM);   /* builds nextBAI and nextBCM */
        else {
            printf_dbg("run_tsi(%d): failed to load nextBAI or nextBCM! 1/n", t->proc_id);
        }
        dirty_grid(t->heap, t->nextBAI_idx);       /* possible optimization -> delete dirty_grid */
        dirty_grid(t->heap, t->nextBCM_idx);
        if (result) {
            printf_dbg("run_tsi(): final clean-up\n");
            if (t->ai_idx == t->bestAI_idx) {
                dirty_grid(t->heap, t->bestAI_idx);
            } else {
                delete_grid(t->heap, t->ai_idx);
            }
        } else {
            printf_dbg("run_tsi(): Compare failed!\n");
            return 0;
        }
        getCurrTime(&t2);
        printf_dbg("run_tsi(%d): Compare execution took %f secs\n",t->proc_id,getElapsedTime(&t1,&t2));
        
		// TODO: REMOVE THIS
		/*
		printf_dbg("run_tsi(): DUMPING BAI & BCM!\n");
		fp = create_file("BAI.out");
		write_ascii_grid_file(fp, t->nextBAI, t->grid_size);
		close_file(fp);
		fp = create_file("BCM.out");
		write_ascii_grid_file(fp, t->nextBCM, t->grid_size);
		close_file(fp);
		*/
		// --end--
		
        delete_grid(t->heap, t->cm_idx);
        delete_grid(t->heap, t->sy_idx);
        t->ai_idx = t->cm_idx = t->sy_idx = -1;  /* free aux grids */
    } /* for 1/n */
    
    if (t->n_procs > 1) {
	/* find which node got best correlation */
        if (tsi_is_best_parallel(t)) return 0;
	
		if (t->global_best.proc_id != t->proc_id) {
		    /* delete local best AI grid */
			delete_grid(t->heap, t->bestAI_idx);
			t->bestAI_idx = -1;
		}

        if (tsi_compare_parallel(t)) return 0;
    }

    /* NEXT ITERATIONS */
    for (j = 1; j < t->iterations; j++) {
        printf_dbg("run_tsi(%d): iteration %d\n", t->proc_id, j);

        /* random layers */
        getCurrTime(&t1);
        delete_cmgrid(cmg);   // conflict with store???
        if (t->proc_id == t->root_id) {
            cmg = new_cmgrid(t->si_eng, 0);    /* generate new set of layers */
            printf("New set of layers:\n");
        } else {
            cmg = new_cmgrid(t->si_eng, 1);    /* return an empty cm_grid */
		}
		print_layers(cmg);
        if (tsi_set_layers_parallel(t, cmg)) {
            printf_dbg("run_tsi(%d): failed to set layers\n", t->proc_id);
            delete_tsi(t);
            return 0;
        }
        store_cmgrid(t->si_eng, cmg);  /* adds CM/C grid to SI engine */
        getCurrTime(&t2);
        printf_dbg("run_tsi(%d): New set of layers took %f secs\n",t->proc_id,getElapsedTime(&t1,&t2));

        t->currBAI_idx = t->nextBAI_idx;
        t->currBCM_idx = t->nextBCM_idx;
        t->nextBAI_idx = -1;
        t->nextBCM_idx = -1;
        t->currBAI = load_grid(t->heap, t->currBAI_idx);
        setup_dss(t->dss_eng, t->currBAI);

        /* CODSS n/1 */
        if ((t->ai_idx = new_grid(t->heap)) < 0) {
            printf_dbg("run_tsi(%d): failed to allocate AI grid n/1\n", t->proc_id);
            delete_tsi(t);
            return 0;
        }
        printf_dbg2("run_tsi(%d): loading AI grid n/1\n", t->proc_id);
        t->ai = load_grid(t->heap, t->ai_idx);
        printf_dbg2("run_tsi(%d): loading currBCM grid n/1\n", t->proc_id);
        t->currBCM = load_grid(t->heap, t->currBCM_idx);
        printf_dbg2("run_tsi(%d): loading currBAI grid n/n\n", t->proc_id);
        t->currBAI = load_grid(t->heap, t->currBAI_idx);
		if (t->ai && t->currBAI && t->currBCM) {
			printf_dbg("run_tsi(%d): (%d/0) running CoDSS \n", t->proc_id, j);
            getCurrTime(&t1);
            result = run_codss(t->dss_eng, t->currBAI, t->currBCM, t->ai);    /* 8 GRIDS -> fuck-up */
            getCurrTime(&t2);
        } else {
            printf_dbg("run_tsi(%d): failed to load AI, currBAI or currBCM for CoDSS! n/1\n", t->proc_id);
            return 0;
        }
        clear_grid(t->heap, t->currBAI_idx);
        clear_grid(t->heap, t->currBCM_idx);
        printf_dbg("run_tsi(%d): (%d/0) CoDSS terminated (%f secs)\n",t->proc_id,j,getElapsedTime(&t1,&t2));


        /* SI */
        if (result) {
            t->seismic = load_grid(t->heap, t->seismic_idx);
            printf_dbg2("run_tsi(%d): loading SY n/n\n", t->proc_id);
            if ((t->cm_idx = new_grid(t->heap)) < 0) {
                printf_dbg("run_tsi(%d): failed to allocate CM grid\n", t->proc_id);
                delete_tsi(t);
                return 0;
            }
            if ((t->sy_idx = new_grid(t->heap)) < 0) {
                printf_dbg("run_tsi(%d): failed to allocate SY grid\n", t->proc_id);
                delete_tsi(t);
                return 0;
            }
            printf_dbg2("run_tsi(%d): loading CM n/n\n", t->proc_id);
            t->cm = load_grid(t->heap, t->cm_idx);
            printf_dbg2("run_tsi(%d): loading seismic n/n\n", t->proc_id);
            t->sy = load_grid(t->heap, t->sy_idx);
            if (t->cm && t->seismic && t->sy) {
                printf_dbg("run_tsi(%d): (%d/0) running SI\n", t->proc_id,j);
                getCurrTime(&t1);
                result = run_si(t->si_eng, t->ai, t->seismic, t->cm, t->sy);
                getCurrTime(&t2);
            } else {
                printf_dbg("run_tsi(%d): failed to load CM, Seismic or SY for SI! n/1\n", t->proc_id);
                return 0;
            }
        } else {
            printf_dbg("run_tsi(%d): SI simulation failed!\n", t->proc_id);
            return 0;
        }
        printf_dbg("run_tsi(%d): (%d/0) SI terminated (%f secs)\n",t->proc_id,j,getElapsedTime(&t1,&t2));

        /* first isBest/Compare (optimal execution) 1/1 */
        getCurrTime(&t1);
        ai_corr = grid_correlation(t->seismic, t->sy, t->grid_size);   /* between real and synthetic seismic data */
		printf_dbg("isBest(%d/0): grid correlation: %f\n",j,ai_corr);
        clear_grid(t->heap, t->seismic_idx);

        if (ai_corr > t->global_best.value) {
            t->global_best.value = ai_corr;
            delete_grid(t->heap, t->bestAI_idx);   /* new bestAI, delete old */
            t->bestAI_idx = t->ai_idx;
            printf("\n\t\t\ti (%d/0) NEW BEST CORRELATION = %f\n\n", j,t->global_best.value);
            grid_copy(t->sy, t->ai, t->grid_size);     /* nextBAI = bestAI = AI at this stage (t->sy as aux grid) */
            t->nextBAI_idx = t->sy_idx;
            dirty_grid(t->heap, t->bestAI_idx);
        } else {
            delete_grid(t->heap, t->sy_idx);
            t->nextBAI_idx = t->ai_idx;
        }
        dirty_grid(t->heap, t->nextBAI_idx);

        t->nextBCM_idx = t->cm_idx;  /* nextBCM = CM at this stage */
        dirty_grid(t->heap, t->nextBCM_idx);
        t->ai_idx = t->cm_idx = t->sy_idx = -1;  /* free aux grids */
        getCurrTime(&t2);
        printf_dbg("run_tsi(%d): Compare execution took %f secs\n",t->proc_id,getElapsedTime(&t1,&t2));


        for (i = 1; i < t->simulations; i++) {
            printf_dbg("run_tsi(%d): iteration %d, simulation %d\n", t->proc_id, j, i);

            /* CODSS n/n */
            if ((t->ai_idx = new_grid(t->heap)) < 0) {
                printf_dbg("run_tsi(%d): failed to allocate AI grid n/1\n", t->proc_id);
                delete_tsi(t);
                return 0;
            }
            printf_dbg2("run_tsi(%d): loading AI grid n/1\n", t->proc_id);
            t->ai = load_grid(t->heap, t->ai_idx);
            printf_dbg2("run_tsi(%d): loading currBCM grid n/1\n", t->proc_id);
            t->currBCM = load_grid(t->heap, t->currBCM_idx);
            printf_dbg2("run_tsi(%d): loading currBAI grid n/n\n", t->proc_id);
            t->currBAI = load_grid(t->heap, t->currBAI_idx);
			if (t->ai && t->currBAI && t->currBCM) {
				printf_dbg("run_tsi(%d): (%d/%d) running CoDSS \n", t->proc_id, j,i);
                getCurrTime(&t1);
                result = run_codss(t->dss_eng, t->currBAI, t->currBCM, t->ai);    /* 8 GRIDS -> fuck-up */
                getCurrTime(&t2);
            } else {
                printf_dbg("run_tsi(%d): failed to load AI, currBAI or currBCM for CoDSS! n/1\n", t->proc_id);
                return 0;
            }
            clear_grid(t->heap, t->currBAI_idx);
            clear_grid(t->heap, t->currBCM_idx);
            printf_dbg("run_tsi(%d): (%d/%d) CoDSS terminated (%f secs)\n",t->proc_id,j,i,getElapsedTime(&t1,&t2));


            /* SI */
            if (result) {
                t->seismic = load_grid(t->heap, t->seismic_idx);
                printf_dbg2("run_tsi(%d): loading SY n/n\n", t->proc_id);
                if ((t->cm_idx = new_grid(t->heap)) < 0) {
                    printf_dbg("run_tsi(%d): failed to allocate CM grid\n", t->proc_id);
                    delete_tsi(t);
                    return 0;
                }
                if ((t->sy_idx = new_grid(t->heap)) < 0) {
                    printf_dbg("run_tsi(%d): failed to allocate SY grid\n", t->proc_id);
                    delete_tsi(t);
                    return 0;
                }
                printf_dbg2("run_tsi(%d): loading CM n/n\n", t->proc_id);
                t->cm = load_grid(t->heap, t->cm_idx);
                printf_dbg2("run_tsi(%d): loading seismic n/n\n", t->proc_id);
                t->sy = load_grid(t->heap, t->sy_idx);
                if (t->cm && t->seismic && t->sy) {
                    printf_dbg("run_tsi(%d): (%d/%d) running SI n/1\n", t->proc_id,j,i);
                    getCurrTime(&t1);
                    result = run_si(t->si_eng, t->ai, t->seismic, t->cm, t->sy);
                    getCurrTime(&t2);
                } else {
                    printf_dbg("run_tsi(%d): failed to load CM, Seismic or SY for SI! n/1\n", t->proc_id);
                    return 0;
                }
            } else {
                printf_dbg("run_tsi(%d): (%d/%d) SI simulation failed!\n", t->proc_id,j,i);
                return 0;
            }
            printf_dbg("run_tsi(%d): (%d/%d) SI terminated (%f secs)\n",t->proc_id,j,i,getElapsedTime(&t1,&t2));

            /* Is best? */
            getCurrTime(&t1);
            if (result) {
                ai_corr = grid_correlation(t->seismic, t->sy, t->grid_size);
				printf_dbg("isBest(%d/%d): grid correlation: %f\n",j,i,ai_corr);
                if (ai_corr > t->global_best.value) {
					if (t->bestAI_idx > -1) 
						delete_grid(t->heap, t->bestAI_idx);
					t->global_best.value = ai_corr;
					t->global_best.proc_id = t->proc_id;
					t->bestAI_idx = t->ai_idx;
                    printf("\n\tNEW BEST CORRELATION = %f\n\n", t->global_best.value);
                }
                delete_grid(t->heap, t->sy_idx);
                clear_grid(t->heap, t->seismic_idx);
            } else {
                printf_dbg("run_tsi(): si simulation failed!\n");
                return 0;
            }

            /* Compare */
            t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
            t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
            result = tsi_compare(t->grid_size, t->ai, t->cm, t->nextBAI, t->nextBCM);   /* builds nextBAI and nextBCM */
            if (result) {
                if (t->ai_idx == t->bestAI_idx) {
                    dirty_grid(t->heap, t->bestAI_idx);
                } else {
                    delete_grid(t->heap, t->ai_idx);
                }
            } else {
                printf_dbg("run_tsi(%d): Compare failed!\n", t->proc_id);
                return 0;
            }
            getCurrTime(&t2);
            delete_grid(t->heap, t->cm_idx);
            dirty_grid(t->heap, t->nextBAI_idx);
            dirty_grid(t->heap, t->nextBCM_idx);
            t->ai_idx = t->cm_idx = t->sy_idx = -1;  /* free aux grids */
            printf_dbg("run_tsi(%d): (%d/%d) Compare terminated (%f secs)\n",t->proc_id,j,i,getElapsedTime(&t1,&t2));

			// TODO: REMOVE THIS
			/*
			printf_dbg("run_tsi(): (%d/%d) DUMPING BAI & BCM!\n",j,i);
			fp = create_file("BAI.out");
			write_ascii_grid_file(fp, t->nextBAI, t->grid_size);
			close_file(fp);
			fp = create_file("BCM.out");
			write_ascii_grid_file(fp, t->nextBCM, t->grid_size);
			close_file(fp);
			*/
			// --end--

        } /* for (simulations) */

        if (t->n_procs > 1) {
            if (tsi_is_best_parallel(t)) return 0;
			
			if (t->global_best.proc_id != t->proc_id) {
				/* delete local best AI grid */
				delete_grid(t->heap, t->bestAI_idx);
				t->bestAI_idx = -1;
			}

            if (tsi_compare_parallel(t)) return 0;
        }

        delete_grid(t->heap, t->currBAI_idx);
        delete_grid(t->heap, t->currBCM_idx);
    } /* for (iterations) */
    
	/* save best AI grid */
	if (t->global_best.proc_id == t->proc_id) {
		fp = create_file("bestAI.out");
		write_ascii_grid_file(fp, t->ai, t->grid_size);
		close_file(fp);
		printf_dbg("run_tsi(%d): bestAI grid dumped\n",t->proc_id);
	}
	
    printf_dbg("run_tsi(%d): heap performance: R=%d ", t->proc_id, t->heap->reads);
    printf_dbg("W=%d G=%d X=%d\n", t->heap->writes, t->heap->curr_grids, t->heap->max_grids);
    printf("\t\t\tFinal correlation = %f\n", t->global_best.value);
    return 1;
} /* run_tsi */



int tsi_compare(unsigned int size, float *AI, float *CM, float *BAI, float *BCM) {
    int i;
    unsigned int x;

    printf_dbg("tsi_compare(): called\n");
    /* execute Compare */
    x = 0;
    for (i = 0; i < size; i++) {
        if (BCM[i] < CM[i]) {
            /* new best correlation */ 
            BCM[i] = CM[i];
            BAI[i] = AI[i];
            x++;
        }
    }
    printf_dbg("tsi_compare(): changed %d points\n", x);
    return 1;
} /* tsi_compare */



void grid_copy(float *a, float *b, unsigned int grid_size)
{
    printf_dbg2("tsi_compare(): called\n");
    memcpy(a, b, grid_size*sizeof(float));
}


/* end of file tsi.c */

