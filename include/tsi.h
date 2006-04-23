#ifndef _TSI_H
#define _TSI_H

#include "registry.h"
#include "grid_heap.h"
#include "dss.h"
#include "si.h"
#include "log.h"


#define CARTESIAN_FILE  0
#define TSI_ASCII_FILE  1
#define TSI_BIN_FILE    2
#define GSLIB_FILE      3
#define SGEMS_FILE      4


typedef struct corr_type {
    float value;         /* current best global correlation value */
    int   proc_id;       /* process that holds best AI data */
} corr;


typedef struct tsi_type {
    /* auxiliar objects */
    registry *reg;            /* reference to the registry */
    grid_heap *heap;          /* reference to the grid heap */
    log_t *l;                   /* reference to the log */

    /* geostats engines */
    dss *dss_eng;             /* reference to the DSS and CoDSS engines */
    si  *si_eng;              /* reference to the SI engine */

    /* TSI execution parameters*/
    int iterations,           /* number of iterations */
        simulations;          /* number of simulations for each iteration */
    int root;
    
    /* execution time counters */
    double mm_time;           /* time spent in memory management functions */
    double dss_time;          /* time for simulations */
    double si_time;           /* time for seismic inversion functions */
    double corr_time;         /* time for correlation results evaluation */
    double par_time;          /* time spent executing parallel functions */

    /* parallel execution related data */
    unsigned int grid_segsize;   /* grid segment size = (int)grid_size/n_procs */
    corr global_best;            /* data to be exchanged during parallel isBestAI? */
    int n_procs,                 /* number of processes running */
        proc_id,                 /* process ID */
        root_id,                 /* root process ID */
        optimize;                /* runtime optimization flag */

    /* grid size */
    int xsize,
        ysize,
        zsize;
    unsigned int grid_size;   /* = xsize * ysize * zsize < 2^32 */
    unsigned int even_size;

    /* grid heap references */
    int seismic_idx,
        bestAI_idx,
        currBAI_idx,
        currBCM_idx,
        nextBAI_idx,
        nextBCM_idx,
        ai_idx,
        cm_idx,
        sy_idx;
        
    /* data pointers */
    float *seismic,
          *bestAI,
          *currBAI,
          *currBCM,
          *nextBAI,
          *nextBCM,
          *ai,
          *cm,
          *sy;

    /* execution log */
    int dump_ai,
        dump_cm,
        dump_bai,
        dump_bcm;
    char *dump_path,
         *input_path,
         *output_path,
         *log_path,
         *seismic_path,
         *empty_path;
    
    /* file formats */
    int result_file,
        seismic_file,
        dump_file;
} tsi;

tsi *new_tsi(registry *reg);

void delete_tsi(tsi *t);

int run_tsi(tsi *t);

#endif /* _TSI_H */
