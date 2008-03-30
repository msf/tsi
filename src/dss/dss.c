#include <stdlib.h>
#include <string.h>

#include "debug.h"
#include "memdebug.h"
#include "log.h"
#include "grid_heap.h"
#include "registry.h"
#include "dss.h"
#include "dss_legacy.h"
#include "tsi_io.h"


dss *new_dss(registry *r, grid_heap *h, log_t *l) {
    dss *d;
    reg_key *k, *kpath;
    /* auxiliar variables for harddata */
    char filename[512];

    printf_dbg2("new_dss(): called\n");
    /* object space allocation */
    d = (dss *) tsi_malloc(sizeof(dss));
    d->general    = (general_vars_t *) tsi_malloc(sizeof(general_vars_t));
    d->search     = (search_vars_t *) tsi_malloc(sizeof(search_vars_t));
    d->simulation = (simulation_vars_t *) tsi_malloc(sizeof(simulation_vars_t));
    d->covariance = (covariance_vars_t *) tsi_malloc(sizeof(covariance_vars_t));
    d->clookup    = (covtable_lookup_vars_t *) tsi_malloc(sizeof(covtable_lookup_vars_t));
    d->krige      = (krige_vars_t *) tsi_malloc(sizeof(krige_vars_t));
	d->harddata	  = (harddata_t *) tsi_malloc(sizeof(harddata_t));
    if (!d ||
        !d->general ||
        !d->search ||
        !d->simulation ||
        !d->covariance ||
        !d->clookup ||
        !d->krige ||
		!d->harddata )
        return NULL;
    d->reg = r;
    d->heap = h;
	d->l = l;
	
	/* 8 and 16 are what it need to work.. :-/
	 * 12 is the precalculated maximum of: na = search->nclose + covtable_lookup->ncnode;
	*/
    d->krige->rr = NULL;
    d->krige->r = NULL;
    d->krige->s = NULL;
    d->krige->a = NULL;
    d->krige->vra = NULL;
    d->krige->vrea = NULL;
    d->krige->last_na = 0;    

    printf_dbg2("new_dss(): Starting new DSS engine.\n Loading dss config settings\n");
	if(dss_parameters(d, r)){
		ERROR(l, "new_dss()", "dss_parameters()");
		return NULL;
	}
	
    printf_dbg2("new_dss(): Reading harddata file.\n");
    /* read harddata file */
    kpath = NULL;
    if ((kpath = get_key(r, "HARDDATA", "PATH")) == NULL)
        kpath = get_key(r, "GLOBAL", "INPUT_PATH");
    if ((k = get_key(r, "HARDDATA", "FILENAME")) == NULL) 
		return NULL;
    if (kpath)
        sprintf(filename, "%s%s", get_string(kpath), get_string(k));
    else
		sprintf(filename, "%s", get_string(k));

	
    if ( load_harddata_file(l, filename, d->harddata) ) {
        printf("new_dss(): ERROR - failed to load harddata file!\n");
        return NULL;
    }
        
    /* parse and prepare harddata */    
    readdata(l, d->harddata, d->general);

    printf_dbg2("new_dss(): DSS engine started sucessfully.\n");
    return d;
} /* new_dss */



int setup_dss(dss *d, int ktype) {
    printf_dbg2("setup_dss(): called\n");

    /* restore initial values */
    d->search->nclose =   0;

    return  0;
} /* setup_dss */



int run_dss(dss *d, float *AI) {
    int *order;
    int ktype = 1;

    printf_dbg2("run_dss(): called\n");
    d->covtab_idx = new_grid(d->heap);
    d->order_idx =  new_grid(d->heap);
    d->clookup->covtab = load_grid(d->heap, d->covtab_idx);
    order = (int *) load_grid(d->heap, d->order_idx);

	d->clookup->ixnode = (short *) tsi_malloc( sizeof(short) * d->general->nxyz);
	d->clookup->iynode = (short *) tsi_malloc( sizeof(short) * d->general->nxyz);
	d->clookup->iznode = (short *) tsi_malloc( sizeof(short) * d->general->nxyz);
    /* restore initial values */
    d->search->nclose =   0;


    /* SIMULATION */
	printf_dbg2("run_dss(): grids allocated, calling dssim()\n");
    dssim(AI, NULL, NULL, order, ktype,
          d->general,
		  d->harddata,
          d->search,
          d->simulation,
          d->covariance,
          d->clookup,
          d->krige);                                                                        

	tsi_free(d->clookup->ixnode);
	tsi_free(d->clookup->iynode);
	tsi_free(d->clookup->iznode);

    delete_grid(d->heap, d->covtab_idx);
    delete_grid(d->heap, d->order_idx);

    
    return 1;
} /* run_dss */



int run_codss(dss *d, float *currBAI, float *currBCM, float *AI) {
    int *order;
    int ktype = 5;

    printf_dbg2("run_codss(): called\n");
    d->covtab_idx = new_grid(d->heap);
    d->order_idx =  new_grid(d->heap);
    d->clookup->covtab = load_grid(d->heap, d->covtab_idx);
    order = (int *) load_grid(d->heap, d->order_idx);

	d->clookup->ixnode = (short *) tsi_malloc( sizeof(short) * d->general->nxyz);
	d->clookup->iynode = (short *) tsi_malloc( sizeof(short) * d->general->nxyz);
	d->clookup->iznode = (short *) tsi_malloc( sizeof(short) * d->general->nxyz);

    /* CO-SIMULATION */
    dssim(AI, currBAI, currBCM, order, ktype,
          d->general,
		  d->harddata,
          d->search,
          d->simulation,
          d->covariance,
          d->clookup,
          d->krige);
                                                    
	tsi_free(d->clookup->ixnode);
	tsi_free(d->clookup->iynode);
	tsi_free(d->clookup->iznode);

    delete_grid(d->heap, d->covtab_idx);
    delete_grid(d->heap, d->order_idx);


    return 1;
} /* run_codss */



void delete_dss(dss *d) {
    printf_dbg2("delete_dss(): called\n");
    if (d) {
        if (d->general) {
            tsi_free(d->general);
        }
        if (d->search) tsi_free(d->search);
        if (d->simulation) tsi_free(d->simulation);
        if (d->covariance) {
            if (d->covariance->variogram) tsi_free(d->covariance->variogram);
        	tsi_free(d->covariance);
        }
        if (d->clookup) tsi_free(d->clookup);
        if (d->krige) {
            if (d->krige->rr) tsi_free(d->krige->rr);
            if (d->krige->r) tsi_free(d->krige->r);
            if (d->krige->s) tsi_free(d->krige->s);
            if (d->krige->a) tsi_free(d->krige->a);
            if (d->krige->vra) tsi_free(d->krige->vra);
            if (d->krige->vrea) tsi_free(d->krige->vrea);            
            tsi_free(d->krige);
        }
		if (d->harddata) {
		   	if(d->harddata->in_grid) tsi_free(d->harddata->in_grid);
			if(d->harddata->point) tsi_free(d->harddata->point);
			tsi_free(d->harddata);
		}
        tsi_free(d);
    }
} /* delete_dss */



