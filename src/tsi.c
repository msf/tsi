#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

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
int   tsi_recurse_iterations(tsi *t, int i, int s);
int   tsi_recurse_simulations(tsi *t, int i, int s);
int   tsi_seismic_inversion(tsi *t, int iteration, int simulation);
int   tsi_direct_sequential_simulation(tsi *t, int iteration, int simulation);
int   tsi_evaluate_best_correlations(tsi *t, int iteration, int simulation);
int   tsi_setup_iteration(tsi *t, int iteration);
int   tsi_finish_iteration(tsi *t, int iteration, int simulation);
int   tsi_save_results(tsi *t);
int   tsi_compare(unsigned int size, float *AI, cm_grid *CM, float *nextBAI, cm_grid *BCM);
void  grid_copy(float *a, float *b, unsigned int grid_size);
int   tsi_write_grid(tsi *t, TSI_FILE *fp, float *grid, int type, char *desc);
int   tsi_read_grid(tsi *t, TSI_FILE *fp, float *grid, int type);



tsi *new_tsi(registry *reg) {
    tsi *t;
    reg_key *k, *kpath;
    int usefs, heap_size, swap_thr, new_sims, over_sims, over_procs;
    TSI_FILE *fp;
    char buf[512];
    unsigned int timeSeed;
    long lTime;

    t = (tsi *) tsi_malloc(sizeof(tsi));
    if (!t) return NULL;
    t->reg = reg;
    t->heap = NULL;

    /* set a starting point for the rand() */
    lTime = time(NULL);
    timeSeed = (unsigned int) lTime / 2;
    srandom(timeSeed);

    t->empty_path = 0; 

    t->input_path = &t->empty_path;
    if ((k = get_key(reg, "GLOBAL", "INPUT_PATH")) != NULL) 
		t->input_path = get_string(k);

    t->output_path = &t->empty_path;
    if ((k = get_key(reg, "GLOBAL", "OUTPUT_PATH")) != NULL) 
		t->output_path = get_string(k);

    t->seismic_path = t->input_path;
    if ((k = get_key(reg, "SEISMIC", "PATH")) != NULL) 
		t->seismic_path = get_string(k);

    t->log_path = t->output_path;
    if ((k = get_key(reg, "GLOBAL", "LOG_PATH")) != NULL) 
		t->log_path = get_string(k);

    t->dump_path = t->output_path;
    if ((k = get_key(reg, "DUMP", "PATH")) != NULL) 
		t->dump_path = get_string(k);


	/* get machine parameters */
	if (new_tsi_parallel(&t->n_procs, &t->proc_id) < 0) {
		printf("Machine failed to start!\n");
		delete_tsi(t);
		return NULL;
	}
	t->root_id = t->proc_id;

	/* create logger */
	t->l = new_log(t->reg, t->proc_id);
	if(t->l == NULL) {
		printf("log creation failed, aborting\n");
		delete_tsi(t);
		return NULL;
	}
	sprintf(buf,"TSI %s - started!\n",TSI_VERSION);
	log_string(t->l, buf); 

	t->global_best.value = -999;
	t->global_best.proc_id = t->proc_id;

	/* get initial data from registry */
	if ((k = get_key(reg, "GLOBAL", "ITERATIONS")) == NULL) {
		log_string(t->l, "ERROR, missing ITERATIONS parameter on config file\n Aborting\n");
		delete_tsi(t);
		return NULL;
	}
    t->iterations = get_int(k);
    if (t->iterations < 1) {
        log_string(t->l,"Incoeherent number of iterations!\n");
        delete_tsi(t);
        return NULL;
    }
    if ((k = get_key(reg, "GLOBAL", "SIMULATIONS")) == NULL) {
		log_string(t->l, "ERROR, missing SIMULATIONS parameter on config file\n Aborting\n");
       delete_tsi(t);
       return NULL;
    }
    t->simulations = get_int(k);
    if (t->simulations < 1) {
        log_string(t->l,"Incoeherent number of simulations!\n");
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

    /* get dump parameters */
    t->dump_ai = 0;
    t->dump_cm = 0;
    t->dump_bai = 0;
    t->dump_bcm = 0;
    if ((k = get_key(reg, "DUMP", "AI"))   != NULL) 
		t->dump_ai = get_int(k);
    if ((k = get_key(reg, "DUMP", "CORR")) != NULL) 
		t->dump_cm = get_int(k);
    if ((k = get_key(reg, "DUMP", "BAI"))  != NULL) 
		t->dump_bai = get_int(k);
    if ((k = get_key(reg, "DUMP", "BCM"))  != NULL) 
		t->dump_bcm = get_int(k);
    
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
       heap_size = 12;
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


	/* log some info */
	sprintf(buf," TSI - %d Simulations x %d Iterations\n",t->simulations, t->iterations);
	log_string(t->l, buf);

	sprintf(buf," TSI - grid size is %u\n",t->grid_size);
	log_string(t->l, buf);

	int size = t->grid_size / (1024*1024) * sizeof(float);
	if(usefs)
		sprintf(buf," TSI - predicted memory use is: <%u MegaBytes\n",(size*swap_thr) + (size/2) );
	else 
		sprintf(buf," TSI - predicted memory use is: <%u MegaBytes\n",(size*heap_size) + (size/2) ); 
	log_string(t->l, buf);


    /* start grid heap */
    printf_dbg2("new_tsi(%d): starting heap\n", t->proc_id);
    kpath = get_key(reg, "GLOBAL", "TMP_PATH");
    if (kpath)
        t->heap = new_heap(t->n_procs, t->proc_id, heap_size, swap_thr, usefs, get_string(kpath), t->even_size);
    else
        t->heap = new_heap(t->n_procs, t->proc_id, heap_size, swap_thr, usefs, NULL, t->even_size);
    if (!t->heap) {
       printf_dbg("new_tsi(%d): failed to start heap\n", t->proc_id);
       delete_tsi(t);
       return NULL;
    }
    t->seismic_idx =
    t->bestAI_idx =
    t->currBAI_idx =
    t->currBCM_idx =
    t->nextBAI_idx =
    t->nextBCM_idx =
    t->ai_idx =
    t->sy_idx =
    t->cm_idx = -1;

    /* start DSS engine */
    t->dss_eng = new_dss(t->reg, t->heap, t->l);
    if (!t->dss_eng) {
       printf_dbg("new_tsi(%d): failed to start dss engine\n", t->proc_id);
       delete_tsi(t);
       return NULL;
    }
	printf_dbg("new_tsi(): DSS engine loaded\n");
    
    /* start SI engine */
    t->si_eng = new_si(t->reg, t->heap, t->l);
    if (!t->si_eng) {
       printf_dbg("new_tsi(%d): failed to start si engine\n", t->proc_id);
       delete_tsi(t);
       return NULL;
    }
	printf_dbg("new_tsi(): SI engine loaded\n");


    /* get file types */
    t->seismic_file = TSI_ASCII_FILE;
    t->dump_file = TSI_BIN_FILE;     /* increases performance for dump/resume features */
    t->result_file = TSI_ASCII_FILE;
    if ((k = get_key(reg, "DUMP", "FILE_TYPE")) != NULL) {
        if (!strcmp(get_string(k), "cart-grid")) t->dump_file = CARTESIAN_FILE;
        else if (!strcmp(get_string(k), "tsi-ascii")) t->dump_file = TSI_ASCII_FILE;
        else if (!strcmp(get_string(k), "tsi-bin")) t->dump_file = TSI_BIN_FILE;
        //else if (!strcmp(get_string(k), "sgems")) t->dump_file = SGEMS_FILE;
        else if (!strcmp(get_string(k), "gslib")) t->dump_file = GSLIB_FILE;
    }
    if ((k = get_key(reg, "GLOBAL", "RESULT_TYPE")) != NULL) {
        if (!strcmp(get_string(k), "cart-grid")) t->result_file = CARTESIAN_FILE;
        else if (!strcmp(get_string(k), "tsi-ascii")) t->result_file = TSI_ASCII_FILE;
        else if (!strcmp(get_string(k), "tsi-bin")) t->result_file = TSI_BIN_FILE;
        //else if (!strcmp(get_string(k), "sgems")) t->result_file = SGEMS_FILE;
        else if (!strcmp(get_string(k), "gslib")) t->result_file = GSLIB_FILE;
    }
    if ((k = get_key(reg, "SEISMIC", "FILE_TYPE")) != NULL) {
        if (!strcmp(get_string(k), "cart-grid")) t->seismic_file = CARTESIAN_FILE;
        else if (!strcmp(get_string(k), "tsi-ascii")) t->seismic_file = TSI_ASCII_FILE;
        else if (!strcmp(get_string(k), "tsi-bin")) t->seismic_file = TSI_BIN_FILE;
        //else if (!strcmp(get_string(k), "sgems")) t->seismic_file = SGEMS_FILE;
        else if (!strcmp(get_string(k), "gslib")) t->seismic_file = GSLIB_FILE;
    }

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
    sprintf(buf, "%s%s", t->seismic_path, get_string(k));
    if ((fp = open_file(buf)) == NULL) {
        printf("ERROR: Failed to open the seismic grid file!\n");
        delete_tsi(t);
        return NULL;
    }
    if (!tsi_read_grid(t, fp, t->seismic, t->seismic_file)) {
        printf("ERROR: Failed to load seismic file!\n");
        delete_tsi(t);
        return NULL;
    }
    printf_dbg("new_tsi(): Seismic Data loaded\n");
    dirty_grid(t->heap, t->seismic_idx);

    /* reset time counters */
    t->mm_time = 0;
    t->dss_time = 0;
    t->si_time = 0;
    t->corr_time = 0;
    t->par_time = 0;
    
    t->nextBCM_c = t->cm_c = NULL;

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

	debug_check();
    delete_tsi_parallel();
} /* delete_tsi */



int run_tsi(tsi *t) {
	struct timeval t1,t2;
	
    printf_dbg("run_tsi(%d,0,0): begin\n", t->proc_id);
	log_separator(t->l);

	getCurrTime(&t1);

    /* Expand inverted execution tree */
    if (tsi_recurse_iterations(t, t->iterations, t->simulations) == 0) return 0;

	getCurrTime(&t2);
    tsi_save_results(t);

	log_result(t->l, 0, "run_tsi() - HEAP READ grids ",t->heap->reads);
	log_result(t->l, 0, "run_tsi() - HEAP WRITE grids",t->heap->writes);
	log_result(t->l, 0, "run_tsi() - HEAP CURR grids ",t->heap->curr_grids);
	log_result(t->l, 0, "run_tsi() - HEAP MAXIMUM grids",t->heap->max_grids);
	log_string(t->l,"\n");
	log_action_time(t->l, 0, "run_tsi() - TOTAL cpuTime for Memory Management",t->mm_time);
	log_action_time(t->l, 0, "run_tsi() - TOTAL cpuTime for Direct Sequential Simulation",t->dss_time);
	log_action_time(t->l, 0, "run_tsi() - TOTAL cpuTime for Seismic Inversion",t->si_time);
	log_action_time(t->l, 0, "run_tsi() - TOTAL cpuTime for Eval Best Correlations",t->corr_time);
	log_action_time(t->l, 0, "run_tsi() - TOTAL cpuTime",getElapsedTime(&t1, &t2));
	log_string(t->l,"\n");
	log_result(t->l, 0, "run_tsi() - FINAL CORRELATION ",t->global_best.value); 

    return 1;
} /* run_tsi */



int tsi_recurse_iterations(tsi *t, int i, int s) {
    if (--i) {
        if (tsi_recurse_iterations(t, i, s))
            if (tsi_recurse_simulations(t, i, s))
                return tsi_finish_iteration(t, i, s);
    } else {
        if (tsi_recurse_simulations(t, i, s))
            return tsi_finish_iteration(t, i, s);
    }
    return 0;
} /* tsi_recurse_iterations */



int tsi_recurse_simulations(tsi *t, int i, int s) {
    if (--s) {
        if (tsi_recurse_simulations(t, i, s))
			return tsi_simulation(t, i, s);
    } else {
        if (tsi_setup_iteration(t, i))
			return tsi_simulation(t, i, s);
    }
    return 0;
} /* tsi_recurse_simulations */



int tsi_simulation(tsi *t, int i, int s)
{
	int r;
	struct timeval t1, t2;

	log_separator(t->l);
	log_simulation_number(t->l, s);
	
	getCurrTime(&t1);
	r = tsi_direct_sequential_simulation(t, i, s);
	if (r) {
		r = tsi_seismic_inversion(t, i, s);
		if (r) 
			r = tsi_evaluate_best_correlations(t, i ,s);
	}
	getCurrTime(&t2);
	log_action_time( t->l, 0, "tsi_simulation() cpuTime", getElapsedTime(&t1,&t2));
	
	return r;
}



int tsi_setup_iteration(tsi *t, int iteration) 
{
    struct timeval t1, t2, t3, t4;
    double mm_time, run_time, par_time;
    cm_grid *cmg;

    /* prepare simulations */
	log_separator(t->l);
	log_iteration_number(t->l, iteration);
	log_message(t->l, 0, "tsi_setup_iteration() setup [Co]-Direct Sequential Simulation");

    getCurrTime(&t1);
    if (iteration == 0) { /* first iteration */
        getCurrTime(&t2);
        t->currBAI_idx = -1; 
        t->currBCM_idx = -1;
        //setup_dss(t->dss_eng, NULL);
    } else {
        t->currBAI_idx = t->nextBAI_idx;
        t->currBCM_idx = t->nextBCM_idx;
        t->nextBAI_idx = -1;
        t->nextBCM_idx = -1;
        //t->currBAI = load_grid(t->heap, t->currBAI_idx);
        
        /* expand compressed correlations grid from previous iteration */
        t->currBCM_idx = new_grid(t->heap);
        t->currBCM = load_grid(t->heap, t->currBCM_idx);
        expand_correlations_grid(t->nextBCM_c, t->currBCM);
        dirty_grid(t->heap, t->currBCM_idx);
        delete_cmgrid(t->nextBCM_c);

        getCurrTime(&t2);
        //setup_dss(t->dss_eng, t->currBAI);
        //dirty_grid(t->heap, t->currBAI_idx);
    }
    getCurrTime(&t3);

    /* prepare seismic inversion (new set of layers) */
    log_message(t->l, 0, "tsi_setup_iteration() setup Seismic Inversion");
    cmg = load_cmgrid(t->si_eng);
    delete_cmgrid(cmg);   /* conflict with store??? */
    if (t->proc_id == t->root_id) {
        cmg = new_cmgrid(t->si_eng, 0);    /* generate new set of layers */
    } else {
        cmg = new_cmgrid(t->si_eng, 1);    /* return an empty cm_grid */
    }
    if (!cmg) return 0;

    if (tsi_set_layers_parallel(t, cmg)) {
        printf_dbg("tsi_setup_iteration(%d,%d,0): failed to set layers\n", t->proc_id, iteration);
        delete_tsi(t);
        return 0;
    }
    store_cmgrid(t->si_eng, cmg);  /* adds CM/C grid to SI engine */
    getCurrTime(&t4);

   // printf_dbg("tsi_setup_iteration: \tNew set of layers:\n");
   // print_layers(cmg);
   // printf("\n");
    
    mm_time  = getElapsedTime(&t1, &t2);
    run_time = getElapsedTime(&t2, &t3);
    par_time = getElapsedTime(&t3, &t4);
    t->mm_time  += mm_time;
    t->dss_time += run_time;
    t->par_time += par_time;

	log_action_time(t->l, 0, "tsi_setup_iteration() cpuTime in DSS setup",run_time - mm_time - par_time); 
	log_action_time(t->l, 0, "tsi_setup_iteration() cpuTime in SI setup",par_time); 
	log_action_time(t->l, 0, "tsi_setup_iteration() cpuTime in memory management",mm_time); 
	log_action_time(t->l, 0, "tsi_setup_iteration() cpuTime",run_time); 

	/*
    printf_dbg("tsi_setup_iteration(%d,%d,0): terminated", t->proc_id, iteration);
    printf_dbg("\t %f secs (%f secs for memory management)\t", run_time, mm_time);
    printf_dbg("\tlayers setup took %f secs\n", t->proc_id, iteration, par_time);
	*/

    return 1;
} /* tsi_setup_iteration */



int tsi_direct_sequential_simulation(tsi *t, int iteration, int simulation) 
{
    struct timeval t1, t2, t3;
    int result;
    double run_time, mm_time;
    int aux_idx;
    float *aux;

    /* prepare simulations */
    if (iteration == 0) { /* first iteration */
	//t->currBAI_idx = -1; 
	//t->currBCM_idx = -1;
        setup_dss(t->dss_eng, NULL);
    } else {
        //t->currBAI_idx = t->nextBAI_idx;
        //t->currBCM_idx = t->nextBCM_idx;
        //t->nextBAI_idx = -1;
        //t->nextBCM_idx = -1;
        t->currBAI = load_grid(t->heap, t->currBAI_idx);
        
        /* ugly hack to reset BAI for each simulation */
        aux_idx = new_grid(t->heap);
        aux = load_grid(t->heap, aux_idx);
        grid_copy(aux, t->currBAI, t->grid_size);
        clear_grid(t->heap, t->currBAI_idx);
        t->currBAI = aux;

        //getCurrTime(&t2);
        setup_dss(t->dss_eng, aux);
        //dirty_grid(t->heap, t->currBAI_idx);
    }

    /* load new AI for simulation result */
    getCurrTime(&t1);
    printf_dbg2("tsi_dss(%d,%d,%d): new AI grid\n", t->proc_id, iteration, simulation);
	if ((t->ai_idx = new_grid(t->heap)) < 0) {
		log_message(t->l, 0, "ERROR tsi_DSSimulation() failed to create new AI grid");
        delete_tsi(t);
        return 0;
    }
    printf_dbg2("tsi_dss(%d,%d,%d): loading AI grid\n", t->proc_id, iteration, simulation);
	if ((t->ai = load_grid(t->heap, t->ai_idx)) == NULL) {
		log_message(t->l, 0, "ERROR tsi_DSSimulation() failed to load new AI grid");
        return 0;
    }

    /* run simulation */
    if (iteration == 0) { /* first iteration */
		log_message(t->l, 1, "tsi_DSSimulation() running Direct Sequential Simulation"); 
        getCurrTime(&t2);
        result = run_dss(t->dss_eng, t->ai);
        getCurrTime(&t3);
        dirty_grid(t->heap, t->ai_idx);
    } else {              /* co-simulation */
        printf_dbg2("tsi_dss(%d,%d,%d): loading currBCM grid\n", t->proc_id, iteration, simulation);
        t->currBCM = load_grid(t->heap, t->currBCM_idx);
        //printf_dbg2("tsi_dss(%d,%d,%d): loading currBAI grid\n", t->proc_id, iteration, simulation);
        //t->currBAI = load_grid(t->heap, t->currBAI_idx);
        if (!t->currBAI || !t->currBCM) {
			log_message(t->l, 0, "ERROR tsi_DSSimulation() failed to load currBAI or currBCM");
            delete_tsi(t);
            return 0;
        }
		log_message(t->l, 1, "tsi_DSSimulation() running Direct Sequential Co-Simulation"); 
        getCurrTime(&t2);
        result = run_codss(t->dss_eng, t->currBAI, t->currBCM, t->ai);    /* 8 GRIDS -> too much  :-( */
        getCurrTime(&t3);
		printf_dbg("tsi_dss() - coDSS ended\n");
        //clear_grid(t->heap, t->currBAI_idx);
        delete_grid(t->heap, aux_idx);
        clear_grid(t->heap, t->currBCM_idx);
        dirty_grid(t->heap, t->ai_idx);
		printf_dbg("tsi_dss() - coDSS ended, grids curr* grids cleared\n");
    }

	TSI_FILE *fp;

	/*
	printf_dbg("DSS(%d,%d,%d) DUMPING AI!\n",t->proc_id,iteration,simulation);
	fp = create_file("AI.out");
	write_ascii_grid_file(fp, t->ai, t->grid_size);
	close_file(fp);
	*/

    /* evaluate simulation result */
    run_time = getElapsedTime(&t2,&t3);
    mm_time = getElapsedTime(&t1,&t2);
    t->dss_time += run_time;
    t->mm_time += mm_time;
    if (!result) {
        printf_dbg("tsi_dss(%d,%d,%d): [Co]DSS run failed!", t->proc_id, iteration, simulation);
        delete_tsi(t);
        return 0;
    }

	if(mm_time > 0.01)
		log_action_time(t->l, 1, "tsi_DSSimulation() cpuTime in memory management", mm_time);
	log_action_time(t->l, 1, "tsi_DSSimulation() cpuTime",run_time); 
	log_string(t->l,"\n");

    return 1;
} /* tsi_direct_sequential_simulation */



int tsi_seismic_inversion(tsi *t, int iteration, int simulation) 
{
    struct timeval t1, t2, t3;
    double mm_time, run_time;
    int result;

    getCurrTime(&t1);
    /* allocate new grids for CM and SY */
    if ((t->cm_idx = new_grid(t->heap)) < 0) {
		log_message(t->l, 0, "ERROR tsi_seismic_inversion() FAILED TO ALLOCATED CM grid");
        delete_tsi(t);
        return 0;
    }
    if ((t->sy_idx = new_grid(t->heap)) < 0) {
		log_message(t->l, 0, "ERROR tsi_seismic_inversion() FAILED TO ALLOCATED SY grid");
        delete_tsi(t);
        return 0;
    }

    /* load all grids needed for seismic inversion */
    printf_dbg2("tsi_seismic_inversion(%d,%d,%d): loading CM\n", t->proc_id, iteration, simulation);
    t->cm = load_grid(t->heap, t->cm_idx);
    printf_dbg2("tsi_seismic_inversion(%d,%d,%d): loading SY\n", t->proc_id, iteration, simulation);
    t->sy = load_grid(t->heap, t->sy_idx);
    printf_dbg2("tsi_seismic_inversion(%d,%d,%d): loading seismic\n", t->proc_id, iteration, simulation);
    t->seismic = load_grid(t->heap, t->seismic_idx);
    printf_dbg2("tsi_seismic_inversion(%d,%d,%d): loading AI\n", t->proc_id, iteration, simulation);
    t->ai = load_grid(t->heap, t->ai_idx);

    /* run seismic inversion */
    if (t->cm && t->seismic && t->sy && t->ai) {
		log_message(t->l, 1, "tsi_seismic_inversion() running Seismic Inversion");
        getCurrTime(&t2);
        result = run_si(t->si_eng, t->ai, t->seismic, t->cm, t->sy);  /* 4 grids */
        getCurrTime(&t3);
    } else {
        printf_dbg("tsi_seismic_inversion(%d,%d,%d): failed to load CM, Seismic, AI or SY for SI!\n", t->proc_id, iteration, simulation);
        return 0;
    }

	/*
	printf_dbg("Seismic_Inversion(%d,%d,%d) DUMPING SY!\n",t->proc_id,iteration,simulation);
	fp = create_file("SY.out");
	write_ascii_grid_file(fp, t->sy, t->grid_size);
	close_file(fp);
	*/

    clear_grid(t->heap, t->seismic_idx);
    clear_grid(t->heap, t->ai_idx);
    //dirty_grid(t->heap, t->cm_idx); /* SI output - needed for BCM update */
    delete_grid(t->heap, t->cm_idx); /* used only as aux grid */
    dirty_grid(t->heap, t->sy_idx); /* SI output - needed for grid_correlation() */

    /* evaluate result */
    run_time = getElapsedTime(&t2,&t3);
    mm_time = getElapsedTime(&t1,&t2);
    t->si_time += run_time;
    t->mm_time += mm_time;
    if (!result) {
        printf_dbg("tsi_seismic_inversion(%d,%d,%d): SI run failed!", t->proc_id, iteration, simulation);
        delete_tsi(t);
        return 0;
    }
	if(mm_time > 0.01)
		log_action_time(t->l, 1, "tsi_seismic_inversion() cpuTime in memory management", mm_time);
	log_action_time(t->l, 1, "tsi_seismic_inversion() cpuTime",run_time);
	/*
    printf_dbg("tsi_seismic_inversion(%d,%d,%d): \t\t\t SI terminated", t->proc_id, iteration, simulation);
    printf_dbg("\t %f secs (%f secs for memory management)\n", run_time, mm_time);
	*/
	log_string(t->l,"\n");

    return 1;
} /* tsi_seismic_inversion */



/**
 * evaluate best correlations is done in each simulation and does the following:
 * - calculate global correlation for AI grid from this simulation
 * - update Best Correlations Map and Best Acoustic Impedance Map with this AI & Corr Map
 */
int tsi_evaluate_best_correlations(tsi *t, int iteration, int simulation) 
{
    struct timeval t1, t2, t3, t4, t5, t6;
    double mm_time, run_time;
    float corr;
    int result;
    cm_grid *cmg;

    /* load grids needed */
    getCurrTime(&t1);

    printf_dbg2("tsi_eval_best_corr(%d,%d,%d): loading SY\n", t->proc_id, iteration, simulation);
    t->sy = load_grid(t->heap, t->sy_idx);
    printf_dbg2("tsi_eval_best_corr(%d,%d,%d): loading seismic\n", t->proc_id, iteration, simulation);
    t->seismic = load_grid(t->heap, t->seismic_idx);
    printf_dbg2("tsi_eval_best_corr(%d,%d,%d): loading AI\n", t->proc_id, iteration, simulation);
    t->ai = load_grid(t->heap, t->ai_idx);
    //printf_dbg2("tsi_eval_best_corr(%d,%d,%d): loading CM\n", t->proc_id, iteration, simulation);
    //t->cm = load_grid(t->heap, t->cm_idx);

    if (!t->seismic || !t->sy || !t->ai/* || !t->cm*/) {
        printf_dbg("tsi_eval_best_corr(%d,%d,%d):", t->proc_id, iteration, simulation);
        printf_dbg(" failed to load Seismic, SY, CM or AI for correlation evaluation!\n");
        return 0;
    }
    
	log_message(t->l, 1, "tsi_eval_best_correlations() computing grid_correlation");
    /* evaluate best global correlation */
    getCurrTime(&t2);    /* BENCHMARKING... */
    corr = grid_correlation(t->seismic, t->sy, t->grid_size);   /* between real and synthetic seismic data */
    getCurrTime(&t3);
    clear_grid(t->heap, t->seismic_idx);
	delete_grid(t->heap, t->sy_idx); /* not needed anymore, delete it */
	t->sy_idx = -1; 

	log_result(t->l, 1, "tsi_eval_best_correlations() grid global correlation",corr);

    /* is best? */
    if (corr > t->global_best.value) {
        t->global_best.value = corr;
        t->global_best.proc_id = t->proc_id;
        if (t->bestAI_idx > -1) 
			delete_grid(t->heap, t->bestAI_idx);   /* new bestAI, delete old */

		t->bestAI_idx = new_grid(t->heap);
		t->bestAI = load_grid(t->heap, t->bestAI_idx);
		grid_copy(t->bestAI, t->ai, t->grid_size);
		dirty_grid(t->heap, t->bestAI_idx);
		log_result(t->l, 1,"tsi_eval_best_correlations() \tNEW BEST AI\t CORRELATION",t->global_best.value);
    }

    /* update best values found for AI and CM */
	log_message(t->l, 1, "tsi_eval_best_correlations() running Compare & Update");
    if (simulation == 0) { /* first simulation */

        getCurrTime(&t4);
        getCurrTime(&t5);

        t->nextBAI_idx = t->ai_idx;    /* nextBxx = xx at this stage */
        //t->nextBCM_idx = t->cm_idx;
        t->ai_idx = t->cm_idx = -1;  /* free aux grids */

        cmg = load_cmgrid(t->si_eng);
        t->nextBCM_c = clone_cmgrid(cmg);

    } else {               /* compare all values */

        printf_dbg2("tsi_eval_best_corr(%d,%d,%d): loading nextBAI + nextBCM\n", t->proc_id, iteration, simulation);
        t->nextBAI = load_grid(t->heap, t->nextBAI_idx);
        //t->nextBCM = load_grid(t->heap, t->nextBCM_idx);
        if (!t->nextBAI /*|| !t->nextBCM*/) {
            printf_dbg("tsi_eval_best_corr(%d,%d,%d): failed to load nextBAI or nextBCM!", t->proc_id);
            return 0;
        }

        getCurrTime(&t4);
		/* update  nextBAI and nextBCM */
        result = tsi_compare(t->grid_size, t->ai, load_cmgrid(t->si_eng), t->nextBAI, t->nextBCM_c);
        getCurrTime(&t5);
		//delete_grid(t->heap, t->cm_idx); /* not needed anymore, delete it */
		delete_grid(t->heap, t->ai_idx);
		t->ai_idx = t->cm_idx = -1;

		log_result(t->l, 2, "tsi_compare() - Compare & Update changed points",result);
        if (result > 0) {
            printf_dbg2("tsi_eval_best_corr(%d,%d,%d): final clean-up\n");
        } else {
            printf_dbg("tsi_eval_best_corr(%d,%d,%d): Compare&Update failed!\n");
            return 0;
        }
    } /* if first simulation */



    dirty_grid(t->heap, t->nextBAI_idx);
    //dirty_grid(t->heap, t->nextBCM_idx);

    getCurrTime(&t6);

    mm_time = getElapsedTime(&t1, &t2);
    run_time = getElapsedTime(&t2, &t3);
    mm_time += getElapsedTime(&t3, &t4);
    run_time += getElapsedTime(&t4, &t5);
    mm_time += getElapsedTime(&t5, &t6);
    t->mm_time += mm_time;
    t->corr_time += run_time;
	if(mm_time > 0.01)
		log_action_time(t->l, 1, "tsi_eval_best_correlation() cpuTime for memory management",mm_time);
	log_action_time(t->l, 1, "tsi_eval_best_correlation() cpuTime",run_time);

	log_string(t->l,"\n");
    return 1;    
} /* tsi_evaluate_best_correlations */



int tsi_finish_iteration(tsi *t, int iteration, int simulation) 
{
    struct timeval t1, t2;
    char fname[32];
    TSI_FILE *fp;

    getCurrTime(&t1);

    if (t->n_procs > 1) {
		log_message(t->l, 0, "tsi_finish_iteration running distributed Compare & Update");
        if (tsi_is_best_parallel(t))
            return 0;

        if (t->global_best.proc_id != t->proc_id) {
	    /* delete local best AI grid */
	    delete_grid(t->heap, t->bestAI_idx);
	    t->bestAI_idx = -1;
        }
    }

    if (tsi_compare_parallel(t))
        return 0;

	log_message(t->l, 0, "tsi_finish_iteration deleting currBAI & currBCM");
    if (iteration > 0) {
        delete_grid(t->heap, t->currBAI_idx);
        delete_grid(t->heap, t->currBCM_idx);
        t->currBAI_idx = t->currBCM_idx = -1;
    }

    getCurrTime(&t2);
	log_action_time(t->l, 0,"tsi_finish_iteration cpuTime", getElapsedTime(&t1, &t2));

    return 1;
} /* tsi_eval_best_correlations */



int tsi_save_results(tsi *t) {
    reg_key *k;
    TSI_FILE *fp;
    char filename[128];

    /* save best AI grid */
    if (t->global_best.proc_id == t->proc_id) {
        t->bestAI = load_grid(t->heap, t->bestAI_idx);

        if ((k = get_key(t->reg, "GLOBAL", "RESULT_FILE")) == NULL)
            sprintf(filename, "%sBestAI.tsi", t->output_path);
        else
            sprintf(filename, "%s%s", t->output_path, get_string(k));

        if ((fp = create_file(filename)) == NULL) {
            printf("ERROR: Failed to create the result grid file!\n");
            delete_tsi(t);
            return 0;
        }

        sprintf(filename, "bestAI");
        if (!tsi_write_grid(t, fp, t->bestAI, t->result_file, filename)) {
            printf("ERROR: Failed to write result file!\n");
            delete_tsi(t);
            return 0;
        }

        close_file(fp);
		log_message(t->l, 0, "save_results wrote bestAI do disk");
    }
	
/*
        t->seismic = load_grid(t->heap, t->seismic_idx);
        sprintf(filename, "%stest_end-19x19x200.tsi", t->output_path);

        if ((fp = create_file(filename)) == NULL) {
            printf("ERROR: Failed to create the result grid file!\n");
            delete_tsi(t);
            return 0;
        }
        sprintf(filename, "seismic");
        if (!tsi_write_grid(t, fp, t->seismic, t->result_file, filename)) {
            printf("ERROR: Failed to write result file!\n");
            delete_tsi(t);
            return 0;
        }
        close_file(fp);
*/

    return 1;
} /* tsi_save_results */



int tsi_compare(unsigned int size, float *AI, cm_grid *CM_C, float *BAI, cm_grid *BCM_C) {
    unsigned int i, j, k, l, z0, z1, nxy;
    int nlayers, *layer;
    float *CM, *BCM;

    nlayers = get_nlayers(CM_C);
    layer = get_layers(CM_C);
    CM = CM_C->cg;
    BCM = BCM_C->cg;
    nxy = CM_C->nxy;

    /* execute Compare */
    l = 0;
    z0 = z1 = 0;
    for (i = 0; i < nlayers; i++) {
        z1 += layer[i];
        for (j = 0; j < nxy; j++) {
            if (BCM[i*nxy+j] < CM[i*nxy+j]) {
                /* new best correlation */ 
                BCM[i*nxy+j] = CM[i*nxy+j];
                for (k = z0; k < z1; k++) {
                    BAI[k*nxy + j] = AI[k*nxy + j];
                    l++;
                }
            } /* if */
        } /* for(j) */
        z0 = z1;
    } /* for(i) */
    return l;
} /* tsi_compare */



void grid_copy(float *a, float *b, unsigned int grid_size)
{
    printf_dbg2("\tgrid_copy(): called\n");
    memcpy(a, b, grid_size*sizeof(float));
}



int tsi_read_grid(tsi *t, TSI_FILE *fp, float *grid, int type) {
    switch(type) {
        case CARTESIAN_FILE:
            return (read_cartesian_grid(fp, grid, t->grid_size));
        case TSI_ASCII_FILE:
        case TSI_BIN_FILE:
            return (read_tsi_grid(fp, grid, t->xsize, t->ysize, t->zsize));
        case GSLIB_FILE:
            return (read_gslib_grid(fp, grid, t->grid_size));
/*
        case SGEMS_FILE:
            return read_sgems_grid(fp, grid, t->grid_size));
*/
        default:
            fprintf(stderr, "ERROR: Unknown grid file type!\n");
            return 0;
    } /* switch */ 
}



int tsi_write_grid(tsi *t, TSI_FILE *fp, float *grid, int type, char *desc) {
    switch(type) {
        case CARTESIAN_FILE:
            return (write_cartesian_grid(fp, grid, t->grid_size));
        case TSI_ASCII_FILE:
        case TSI_BIN_FILE:
            return (write_tsi_grid(fp, type, grid, t->xsize, t->ysize, t->zsize));
        case GSLIB_FILE:
            return (write_gslib_grid(fp, grid, t->xsize, t->ysize, t->zsize, desc));
/*
        case SGEMS_FILE:
            return write_sgems_grid(fp, grid, t->grid_size));
*/
        default:
            fprintf(stderr, "ERROR: Unknown grid file type!\n");
            return 0;
    } /* switch */ 
}

/* end of file tsi.c */
