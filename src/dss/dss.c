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

double * load_harddata_file(TSI_FILE *, char *, unsigned int *);

dss *new_dss(registry *r, grid_heap *h, log_t *l) {
    dss *d;
    reg_key *k, *kpath;
    /* auxiliar variables for harddata */
    TSI_FILE *fp;
    char filename[512];
    char hdbuf[64];

    printf_dbg2("new_dss(): called\n");
    /* object space allocation */
    d = (dss *) tsi_malloc(sizeof(dss));
    d->general    = (general_vars_t *) tsi_malloc(sizeof(general_vars_t));
    d->search     = (search_vars_t *) tsi_malloc(sizeof(search_vars_t));
    d->simulation = (simulation_vars_t *) tsi_malloc(sizeof(simulation_vars_t));
    d->covariance = (covariance_vars_t *) tsi_malloc(sizeof(covariance_vars_t));
    d->clookup    = (covtable_lookup_vars_t *) tsi_malloc(sizeof(covtable_lookup_vars_t));
    d->krige      = (krige_vars_t *) tsi_malloc(sizeof(krige_vars_t));
    if (!d ||
        !d->general ||
        !d->search ||
        !d->simulation ||
        !d->covariance ||
        !d->clookup ||
        !d->krige)
        return NULL;
    d->reg = r;
    d->heap = h;
	d->l = l;

    printf_dbg2(d->l,"new_dss(): Starting new DSS engine.\n Loading dss config settings\n");
	if(dss_parameters(d, r)){
		printf("new_dss(): ERROR loaging dss configs\n");
		return NULL;
	}
	
    printf_dbg2("new_dss(): Reading harddata file.\n");
    /* read harddata file */
    d->harddata_size = 0;
    kpath = NULL;
    if ((kpath = get_key(r, "HARDDATA", "PATH")) == NULL)
        kpath = get_key(r, "GLOBAL", "INPUT_PATH");
    if ((k = get_key(r, "HARDDATA", "FILENAME")) == NULL) return NULL;
    if (kpath)
        sprintf(filename, "%s%s", get_string(kpath), get_string(k));
    else
        sprintf(filename, "%s", get_string(k));
    if ((fp = open_file(filename)) == NULL) {
        printf("new_dss(): ERROR - Can't open Hard Data file: %s\n", get_string(k));
        return NULL;
    }
    if ((d->harddata = load_harddata_file(fp, hdbuf, &d->harddata_size)) == NULL) {
        printf("new_dss(): ERROR - failed to load harddata file!\n");
        return NULL;
    }
    close_file(fp);

    /* ugly hack for "readdata" */
    d->general->maxdat = d->harddata_size / 4;
    d->general->x = (float *) tsi_malloc(d->general->maxdat * sizeof(float));
    d->general->y = (float *) tsi_malloc(d->general->maxdat * sizeof(float));
    d->general->z = (float *) tsi_malloc(d->general->maxdat * sizeof(float));
    d->general->vr = (float *) tsi_malloc(d->general->maxdat * sizeof(float));
    d->general->wt = (float *) tsi_malloc(d->general->maxdat * sizeof(float));
    d->general->vrtr = (float *) tsi_malloc(d->general->maxdat * sizeof(float));
    d->general->vrgtr = (float *) tsi_malloc(d->general->maxdat * sizeof(float));
    d->general->sec = (float *) tsi_malloc(d->general->maxdat * sizeof(float));
    if (!d->general->x ||
        !d->general->y ||
        !d->general->z ||
        !d->general->vr ||
        !d->general->wt ||
        !d->general->vrtr ||
        !d->general->vrgtr ||
        !d->general->sec)
        return NULL;
        
    /* read data & readWellsData */    
    readdata(d->harddata, d->harddata_size, d->general, d->search, d->simulation);
	readWellsData(d->general, d->harddata, d->harddata_size);

    printf_dbg2("new_dss(): DSS engine started sucessfully.\n");
    return d;
} /* new_dss */



int setup_dss(dss *d, int ktype) {
    printf_dbg2("setup_dss(): called\n");
    d->general->ktype = 1;
    if (ktype != 1) 
		d->general->ktype = 5;

    /* restore initial values */
    d->search->nclose =   0;
    d->clookup->nodmax =  d->clookup->nodmax_bk;
    d->clookup->ntry =    d->clookup->ntry_bk;
    d->clookup->icmean =  d->clookup->icmean_bk;
    d->clookup->icvar =   d->clookup->icvar_bk;
    d->simulation->nsim = d->simulation->nsim_bk;

    return  0;
} /* setup_dss */



int run_dss(dss *d, float *AI) {
    int *order, *mask;
    struct dss_check test;
    simulation_vars_t *x;

    printf_dbg2("run_dss(): called\n");
    d->covtab_idx = new_grid(d->heap);
    d->ixnode_idx = new_grid(d->heap);
    d->iynode_idx = new_grid(d->heap);
    d->iznode_idx = new_grid(d->heap);
    d->order_idx =  new_grid(d->heap);
    d->clookup->covtab = load_grid(d->heap, d->covtab_idx);
    d->clookup->ixnode = (int *) load_grid(d->heap, d->ixnode_idx);
    d->clookup->iynode = (int *) load_grid(d->heap, d->iynode_idx);
    d->clookup->iznode = (int *) load_grid(d->heap, d->iznode_idx);
    order = (int *) load_grid(d->heap, d->order_idx);
    mask = NULL;

    /* restore initial values */
    d->search->nclose =   0;
    d->clookup->nodmax =  d->clookup->nodmax_bk;
    d->clookup->ntry =    d->clookup->ntry_bk;
    d->clookup->icmean =  d->clookup->icmean_bk;
    d->clookup->icvar =   d->clookup->icvar_bk;
    d->simulation->nsim = d->simulation->nsim_bk;

/*
    test.general = tsi_malloc(sizeof(general_vars_t));
    memcpy(test.general, d->general, sizeof(general_vars_t));
    test.search = tsi_malloc(sizeof(search_vars_t));
    memcpy(test.search, d->search, sizeof(search_vars_t));
    printf("nclose=%d\n", d->search->nclose);
    test.simulation = tsi_malloc(sizeof(simulation_vars_t));
    memcpy(test.simulation, d->simulation, sizeof(simulation_vars_t));
*/

    /* SIMULATION */
	printf_dbg2("run_dss(): grids allocated, calling dssim()\n");
    dssim(AI, NULL, NULL, order, mask,
          d->general,
          d->search,
          d->simulation,
          d->covariance,
          d->clookup,
          d->krige);                                                                        

    delete_grid(d->heap, d->covtab_idx);
    delete_grid(d->heap, d->ixnode_idx);
    delete_grid(d->heap, d->iynode_idx);
    delete_grid(d->heap, d->iznode_idx);
    delete_grid(d->heap, d->order_idx);

/*
    if (memcmp(test.general, d->general, sizeof(general_vars_t))) {
        printf_dbg("failed general\n");
    }
    tsi_free(test.general);
    if (memcmp(test.simulation, d->simulation, sizeof(simulation_vars_t))) {
        printf_dbg("failed simulation\n");
        x = (simulation_vars_t *) test.simulation;
        if (x->nsim == d->simulation->nsim)
            printf_dbg("simulation OK\n");
        else
            printf_dbg("%d\n", d->simulation->nsim);
    }
    tsi_free(test.simulation);
    printf("nclose=%d\n", d->search->nclose);
    if (memcmp(test.search, d->search, sizeof(search_vars_t))) {
        printf_dbg("failed search\n");
    }
    tsi_free(test.search);
*/
    
    return 1;
} /* run_dss */



int run_codss(dss *d, float *currBAI, float *currBCM, float *AI) {
    int *order, *mask, i;
    struct dss_check test;
    general_vars_t *x;

    printf_dbg2("run_codss(): called\n");
    d->covtab_idx = new_grid(d->heap);
    d->ixnode_idx = new_grid(d->heap);
    d->iynode_idx = new_grid(d->heap);
    d->iznode_idx = new_grid(d->heap);
    d->order_idx =  new_grid(d->heap);
    d->clookup->covtab = load_grid(d->heap, d->covtab_idx);
    d->clookup->ixnode = (int *) load_grid(d->heap, d->ixnode_idx);
    d->clookup->iynode = (int *) load_grid(d->heap, d->iynode_idx);
    d->clookup->iznode = (int *) load_grid(d->heap, d->iznode_idx);
    order = (int *) load_grid(d->heap, d->order_idx);
	mask = NULL;
	d->general->ktype = 5;

    /* restore initial values */
    d->search->nclose =   0;
    d->clookup->nodmax =  d->clookup->nodmax_bk;
    d->clookup->ntry =    d->clookup->ntry_bk;
    d->clookup->icmean =  d->clookup->icmean_bk;
    d->clookup->icvar =   d->clookup->icvar_bk;
    d->simulation->nsim = d->simulation->nsim_bk;


/*
    test.general = tsi_malloc(sizeof(general_vars_t));
    memcpy(test.general, d->general, sizeof(general_vars_t));
    test.search = tsi_malloc(sizeof(search_vars_t));
    memcpy(test.search, d->search, sizeof(search_vars_t));
    printf("nclose=%d\n", d->search->nclose);
    test.simulation = tsi_malloc(sizeof(simulation_vars_t));
    memcpy(test.simulation, d->simulation, sizeof(simulation_vars_t));
*/

    /* CO-SIMULATION */
    dssim(AI, currBAI, currBCM, order, mask,
          d->general,
          d->search,
          d->simulation,
          d->covariance,
          d->clookup,
          d->krige);
                                                    
    delete_grid(d->heap, d->covtab_idx);
    delete_grid(d->heap, d->ixnode_idx);
    delete_grid(d->heap, d->iynode_idx);
    delete_grid(d->heap, d->iznode_idx);
    delete_grid(d->heap, d->order_idx);

/*
    if (memcmp(test.general, d->general, sizeof(general_vars_t))) {
        printf_dbg("failed general\n");
        x = (general_vars_t *) test.general;
        if (x->nd   != d->general->nd)   printf_dbg("nd NOK\n");
        if (x->ntr  != d->general->ntr)  printf_dbg("ntr NOK\n");
        if (x->idbg != d->general->idbg) printf_dbg("idbg NOK\n");
        if (x->lin  != d->general->lin)  printf_dbg("lin NOK\n");
        if (x->lout != d->general->lout) printf_dbg("lout NOK\n");
        if (x->ldbg != d->general->ldbg) printf_dbg("ldbg NOK\n");
        if (x->llvm != d->general->llvm) printf_dbg("llvm NOK\n");
        for (i = 0; i < 7000; i++) x->close[i] = d->general->close[i];
        if (memcmp(test.general, d->general, sizeof(general_vars_t))) {
            printf_dbg("general OK\n");
        }
    }
    tsi_free(test.general);
    if (memcmp(test.simulation, d->simulation, sizeof(simulation_vars_t))) {
        printf_dbg("failed simulation\n");
    }
    printf("nclose=%d\n", d->search->nclose);
    tsi_free(test.simulation);
    if (memcmp(test.search, d->search, sizeof(search_vars_t))) {
        printf_dbg("failed search\n");
    }
    tsi_free(test.search);
*/

    return 1;
} /* run_codss */



void delete_dss(dss *d) {
    printf_dbg2("delete_dss(): called\n");
    if (d) {
        if (d->general) {
            if (d->general->x) tsi_free(d->general->x);
            if (d->general->y) tsi_free(d->general->y);
            if (d->general->z) tsi_free(d->general->z);
            if (d->general->vr) tsi_free(d->general->vr);
            if (d->general->wt) tsi_free(d->general->wt);
            if (d->general->vrtr) tsi_free(d->general->vrtr);
            if (d->general->vrgtr) tsi_free(d->general->vrgtr);
            if (d->general->sec) tsi_free(d->general->sec);
			if (d->general->wellsDataVal) tsi_free(d->general->wellsDataVal);
			if (d->general->wellsDataPos) tsi_free(d->general->wellsDataPos);
            tsi_free(d->general);
        }
        if (d->search) tsi_free(d->search);
        if (d->simulation) tsi_free(d->simulation);
        if (d->covariance) {
            if (d->covariance->it) tsi_free(d->covariance->it);
            if (d->covariance->cc) tsi_free(d->covariance->cc);
            if (d->covariance->aa) tsi_free(d->covariance->aa);
            if (d->covariance->ang1) tsi_free(d->covariance->ang1);
            if (d->covariance->ang2) tsi_free(d->covariance->ang2);
            if (d->covariance->ang3) tsi_free(d->covariance->ang3);
            if (d->covariance->anis1) tsi_free(d->covariance->anis1);
            if (d->covariance->anis2) tsi_free(d->covariance->anis2);
            tsi_free(d->covariance);
        }
        if (d->clookup) tsi_free(d->clookup);
        if (d->krige) tsi_free(d->krige);
        tsi_free(d);
    }
} /* delete_dss */



double *load_harddata_file(TSI_FILE *fp, char *buf, unsigned int *size) {
    double x, y, z, val, *ret;
    int i;

	i = *size;
	if(fscanf(fp,"%lf %lf %lf %lf", &x, &y, &z, &val) != EOF) {
        (*size)++;
		ret = load_harddata_file(fp, buf, size);
        if ( ret ) {
            ret[i]   = x;
            ret[i+1] = y;
            ret[i+2] = z;
            ret[i+3] = val;
            return ret;
        } else
            return NULL;
    }
    *size *= 4;
    return (double *) tsi_malloc(*size * sizeof(double));
} /* load_harddata_file */

