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

float * load_harddata_file(log_t *, char *, unsigned int *);

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
	
	/* 8 and 16 are what it need to work.. :-/
	 * 12 is the precalculated maximum of:
	na = search->nclose + covtable_lookup->ncnode;
	*/
    d->krige->rr = NULL;
    d->krige->r__ = NULL;
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
    d->harddata_size = 0;
    kpath = NULL;
    if ((kpath = get_key(r, "HARDDATA", "PATH")) == NULL)
        kpath = get_key(r, "GLOBAL", "INPUT_PATH");
    if ((k = get_key(r, "HARDDATA", "FILENAME")) == NULL) 
		return NULL;
    if (kpath)
        sprintf(filename, "%s%s", get_string(kpath), get_string(k));
    else
        sprintf(filename, "%s", get_string(k));

    if ((d->harddata = load_harddata_file(l, filename, &d->harddata_size)) == NULL) {
        printf("new_dss(): ERROR - failed to load harddata file!\n");
        return NULL;
    }

    /* ugly hack for "readdata" */
    d->general->maxdat = d->harddata_size / d->general->nvari;
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
    struct dss_check test;
    simulation_vars_t *x;
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
    int *order, *mask;

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
 	int i;
    struct dss_check test;
    general_vars_t *x;
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
			if (d->general->wellsDataPos) tsi_free(d->general->wellsDataPos);
			if (d->general->wellsDataVal) tsi_free(d->general->wellsDataVal);
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
        if (d->krige) {
            if (d->krige->rr) tsi_free(d->krige->rr);
            if (d->krige->r__) tsi_free(d->krige->r__);
            if (d->krige->s) tsi_free(d->krige->s);
            if (d->krige->a) tsi_free(d->krige->a);
            if (d->krige->vra) tsi_free(d->krige->vra);
            if (d->krige->vrea) tsi_free(d->krige->vrea);            
            tsi_free(d->krige);
        }
		if (d->harddata) tsi_free(d->harddata);
        tsi_free(d);
    }
} /* delete_dss */



float *load_harddata_file(log_t *l, char *filename,  unsigned int *size) {
    float x, y, z, val;
    int i, m;
	char line[256];
	float *buf;
	TSI_FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) {
		ERROR(l, "fopen()", filename);
        return NULL;
    }
	
    /* ignore gslib header */
    if (fgets(line, 255, fp) == NULL) {
		ERROR(l, "load_harddata_file()", "reading 1st line of harddata file");
		return NULL;
	}
	if (fgets(line, 255, fp) == NULL) {
		ERROR(l, "load_harddata_file()", "reading 2nd line of harddata file");
		return NULL;
	}
	if (fgets(line, 255, fp) == NULL) {
		ERROR(l, "load_harddata_file()", "reading 3rd line of harddata file");
		return NULL;
	}

	/* first find out how many values */
	i = 0;
	while(fgets(line, 255, fp) != NULL)
		i++;
	*size = i * 4; // 4 values per line
	printf_dbg("read_harddata_file(): %d lines\n",i);
	buf = (float *) tsi_malloc(sizeof(float) * *size);
	fseek(fp, 0, SEEK_SET); // back to start of file
	
	i = 0;
	m = 0;
	while( (m = fscanf(fp,"%f %f %f %f", &x, &y, &z, &val)) != EOF ) {
		if(m != 4)
			log_print(l,"load_harddata_file(): ERROR -  %d can't parse data values, line: %d\n",m, 1+(i/4));
		buf[i++] = x;
		buf[i++] = y;
		buf[i++] = z;
		buf[i++] = val;
	}

    fclose(fp);
	return buf;
} /* load_harddata_file */

