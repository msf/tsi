#include <stdlib.h>
#include "debug.h"
#include "grid_heap.h"
#include "registry.h"
#include "dss.h"
#include "dss_legacy.h"
#include "tsi_io.h"


double *load_harddata_file(TSI_FILE *fp, char *buf, int *size);

dss *new_dss(registry *r, grid_heap *h) {
    dss *d;
    reg_key *k;
    /* auxiliar variables for harddata */
    TSI_FILE *fp;
    char hdbuf[64];
    /* auxiliar variables for variograms */
    int i, varnum;
    char variogram[16];
    float aa1, aa2, aa, sill;


    /* object space allocation */
    d = (dss *) malloc(sizeof(dss));
    d->general    = (general_vars_t *) malloc(sizeof(general_vars_t));
    d->search     = (search_vars_t *) malloc(sizeof(search_vars_t));
    d->simulation = (simulation_vars_t *) malloc(sizeof(simulation_vars_t));
    d->covariance = (covariance_vars_t *) malloc(sizeof(covariance_vars_t));
    d->clookup    = (covtable_lookup_vars_t *) malloc(sizeof(covtable_lookup_vars_t));
    d->krige      = (krige_vars_t *) malloc(sizeof(krige_vars_t));
    d->files      = (file_vars_t *) malloc(sizeof(file_vars_t));
    if (!d ||
        !d->general ||
        !d->search ||
        !d->simulation ||
        !d->covariance ||
        !d->clookup ||
        !d->krige ||
        !d->files)
        return NULL;
    d->reg = r;
    d->heap = h;
    printf_dbg("new_dss(): Starting new DSS engine.\n");

    /* harddata parameters */
    if ((k = get_key(r, "HARDDATA", "NVARI")) == NULL) return NULL;
    d->general->nvari = get_int(k);
    
    if ((k = get_key(r, "HARDDATA", "IXL")) == NULL) return NULL;
    d->general->ixl = get_int(k);

    if ((k = get_key(r, "HARDDATA", "IYL")) == NULL) return NULL;
    d->general->iyl = get_int(k);
    
    if ((k = get_key(r, "HARDDATA", "IZL")) == NULL) return NULL;
    d->general->izl = get_int(k);
    
    if ((k = get_key(r, "HARDDATA", "IVRL")) == NULL) return NULL;

    d->general->ivrl = get_int(k);
    if ((k = get_key(r, "HARDDATA", "IWT")) == NULL) return NULL;
    d->general->iwt = get_int(k);

    if ((k = get_key(r, "HARDDATA", "ISECVR")) == NULL) return NULL;
    d->general->isecvr = get_int(k);

    if ((k = get_key(r, "HARDDATA", "TMIN")) == NULL) return NULL;
    d->general->tmin = get_float(k);

    if ((k = get_key(r, "HARDDATA", "TMAX")) == NULL) return NULL;
    d->general->tmax = get_float(k);


    /* harddata transformations */
    if ((k = get_key(r, "HDTRANS", "ITRANS")) == NULL) return NULL;
    d->general->itrans = get_int(k);

    if ((k = get_key(r, "HDTRANS", "ISMOOTH")) == NULL) return NULL;
    d->general->ismooth = get_int(k);

    if ((k = get_key(r, "HDTRANS", "ISVR")) == NULL) return NULL;
    d->general->isvr = get_int(k);

    if ((k = get_key(r, "HDTRANS", "ISWT")) == NULL) return NULL;
    d->general->iswt = get_int(k);

    if ((k = get_key(r, "HDTRANS", "ZMIN")) == NULL) return NULL;
    d->general->zmin = get_float(k);

    if ((k = get_key(r, "HDTRANS", "ZMAX")) == NULL) return NULL;
    d->general->zmax = get_float(k);

    if ((k = get_key(r, "HDTRANS", "LTAIL")) == NULL) return NULL;
    d->general->ltail = get_int(k);

    if ((k = get_key(r, "HDTRANS", "LTPAR")) == NULL) return NULL;
    d->general->ltpar = get_float(k);

    if ((k = get_key(r, "HDTRANS", "UTAIL")) == NULL) return NULL;
    d->general->utail = get_int(k);

    if ((k = get_key(r, "HDTRANS", "UTPAR")) == NULL) return NULL;
    d->general->utpar = get_float(k);

    if ((d->general->ltail != 1) && (d->general->ltail != 2)) {
        printf_dbg("ERROR invalid lower tail option: %d\n", (int) d->general->ltail);
        printf_dbg("\t\tonly allowed 1 or 2 - see manual.\n");
        return NULL;
    }
    if ((d->general->utail != 1) && (d->general->utail != 2) && (d->general->utail != 4)) {
        printf_dbg("ERROR invalid upper tail option: %d\n", (int) d->general->utail);
        printf_dbg("\t\tonly allowed 1, 2 or 4 - see manual.\n");
        return NULL;
    }
    if ((d->general->utail == 4) && (d->general->utpar < 1)) {
        printf_dbg("ERROR invalid power for hyperbolic tail: %f\n", d->general->utpar);
        printf_dbg("\t\tmust be greater than 1.0.\n");
        return NULL;
    }
    if ((d->general->ltail == 2) && (d->general->ltpar < 0)) {
        printf_dbg("ERROR invalid power for power model: %f\n",d->general->ltpar);
        printf_dbg("\t\tmust be greater than 0.0.\n");
        return NULL;
    }
    if ((d->general->utail == 2) && (d->general->utpar < 0)) {
        printf_dbg("ERROR invalid power for power model: %f\n",d->general->utpar);
        printf_dbg("\t\tmust be greater than 0.0.\n");
        return NULL;
    }

    
    /* simulations quality */
    if ((k = get_key(r, "QUALITY", "NTRY")) == NULL) return NULL;
    d->clookup->ntry = get_int(k);

    if ((k = get_key(r, "QUALITY", "ICVAR")) == NULL) return NULL;
    d->clookup->icvar = get_int(k);

    if ((k = get_key(r, "QUALITY", "ICMEAN")) == NULL) return NULL;
    d->clookup->icmean = get_int(k);

    
    /* grid parameters */
    if ((k = get_key(r, "GRID", "XNUMBER")) == NULL) return NULL;
    d->general->nx = get_int(k);

    if ((k = get_key(r, "GRID", "YNUMBER")) == NULL) return NULL;
    d->general->ny = get_int(k);

    if ((k = get_key(r, "GRID", "ZNUMBER")) == NULL) return NULL;
    d->general->nz = get_int(k);

    if ((k = get_key(r, "GRID", "XCOORD")) == NULL) return NULL;
    d->general->xmn = get_float(k);

    if ((k = get_key(r, "GRID", "YCOORD")) == NULL) return NULL;
    d->general->ymn = get_float(k);

    if ((k = get_key(r, "GRID", "ZCOORD")) == NULL) return NULL;
    d->general->zmn = get_float(k);

    if ((k = get_key(r, "GRID", "XSIZE")) == NULL) return NULL;
    d->general->xsiz = get_float(k);

    if ((k = get_key(r, "GRID", "YSIZE")) == NULL) return NULL;
    d->general->ysiz = get_float(k);
    d->general->nxy = d->general->nx * d->general->ny;
    d->general->nxyz = d->general->nxy * d->general->nz;

    
    /* mask */
    if ((k = get_key(r, "MASK", "NULL_VALUE")) == NULL) return NULL;
    d->general->nosvalue = get_float(k);

    if ((k = get_key(r, "MASK", "USE_MASK")) == NULL) return NULL;
    d->general->imask = get_int(k);


    /* search parameters */
    if ((k = get_key(r, "SEARCH", "NDMIN")) == NULL) return NULL;
    d->search->ndmin = get_int(k);

    if ((k = get_key(r, "SEARCH", "NDMAX")) == NULL) return NULL;
    d->search->ndmax = get_int(k);

    if ((k = get_key(r, "SEARCH", "NODMAX")) == NULL) return NULL;
    d->clookup->nodmax = get_int(k);

    if ((k = get_key(r, "SEARCH", "SSTRAT")) == NULL) return NULL;
    d->search->sstrat = get_int(k);
    if (d->search->sstrat == 1) d->search->ndmax = 0;

    if ((k = get_key(r, "SEARCH", "MULTS")) == NULL) return NULL;
    d->search->mults = get_int(k);

    if ((k = get_key(r, "SEARCH", "NMULTS")) == NULL) return NULL;
    d->search->nmult = get_float(k);

    if ((k = get_key(r, "SEARCH", "NOCT")) == NULL) return NULL;
    d->search->noct = get_int(k);

    if ((k = get_key(r, "SEARCH", "RADIUS")) == NULL) return NULL;
    d->search->radius = get_float(k);
    if (d->search->radius < 1e-20) {
        printf_dbg("new_dss(): radius (%f) must be greater than EPSLON (1E-20)\n", d->search->radius);
        return NULL;
    }
    d->search->radsqd = d->search->radius * d->search->radius;

    if ((k = get_key(r, "SEARCH", "RADIUS1")) == NULL) return NULL;
    d->search->sanis1 = get_float(k) / d->search->radius;

    if ((k = get_key(r, "SEARCH", "RADIUS2")) == NULL) return NULL;
    d->search->sanis2 = get_float(k) / d->search->radius;

    if ((k = get_key(r, "SEARCH", "SANG1")) == NULL) return NULL;
    d->search->sang1 = get_float(k);

    if ((k = get_key(r, "SEARCH", "SANG2")) == NULL) return NULL;
    d->search->sang2 = get_float(k);

    if ((k = get_key(r, "SEARCH", "SANG3")) == NULL) return NULL;
    d->search->sang3 = get_float(k);


    /* kriging parameters */
    if ((k = get_key(r, "KRIG", "TYPE")) == NULL) return NULL;
    d->general->ktype = get_int(k);

    if ((k = get_key(r, "KRIG", "COLOCORR")) == NULL) return NULL;
    d->general->colocorr = get_float(k);

    if ((k = get_key(r, "KRIG", "VARRED")) == NULL) return NULL;
    d->general->varred = get_float(k);


    /* softdata parameters */
    if ((k = get_key(r, "SOFT", "NVARIL")) == NULL) return NULL;
    d->general->nvaril = get_int(k);

    if ((k = get_key(r, "SOFT", "ICOLLVM")) == NULL) return NULL;
    d->general->icollvm = get_int(k);


    /* variogram models */
    if ((k = get_key(r, "VARIOGRAM", "NUMBER")) == NULL) return NULL;
    varnum = d->covariance->nst[0] = get_int(k);
    printf_dbg2("new_dss(): varnum = %d\n", varnum);

    if ((k = get_key(r, "VARIOGRAM", "NUGGET")) == NULL) return NULL;
    sill = d->covariance->c0[0] = get_float(k);

    d->covariance->it = (int *) malloc(varnum * sizeof(int));
    d->covariance->cc = (float *) malloc(varnum * sizeof(float));
    d->covariance->aa = (float *) malloc(varnum * sizeof(float));
    d->covariance->ang1 = (float *) malloc(varnum * sizeof(float));
    d->covariance->ang2 = (float *) malloc(varnum * sizeof(float));
    d->covariance->ang3 = (float *) malloc(varnum * sizeof(float));
    d->covariance->anis1 = (float *) malloc(varnum * sizeof(float));
    d->covariance->anis2 = (float *) malloc(varnum * sizeof(float));
    d->covariance->aa = (float *) malloc(varnum * sizeof(float));
    if (!d->covariance->it ||
        !d->covariance->cc ||
        !d->covariance->aa ||
        !d->covariance->ang1 ||
        !d->covariance->ang2 ||
        !d->covariance->ang3 ||
        !d->covariance->anis1 ||
        !d->covariance->anis2 ||
        !d->covariance->aa) {
        printf_dbg("new_dss(): Failed to allocate space for variograms!\n");
        return NULL;
    }

    for (i = 0; i < varnum; i++) {
        sprintf(variogram, "VARIOGRAM%d", i+1);
        if ((k = get_key(r, variogram, "TYPE")) == NULL) return NULL;
        d->covariance->it[i] = get_int(k);
        if (d->covariance->it[i] == 4) {
            printf_dbg("new_dss(): A power model is not allowed! Choose a different model and re-start.\n");
            return NULL;
        }

        if ((k = get_key(r, variogram, "COV")) == NULL) return NULL;
        d->covariance->cc[i] = get_float(k);
        sill += d->covariance->cc[i];

        if ((k = get_key(r, variogram, "ANG1")) == NULL) return NULL;
        d->covariance->ang1[i] = get_float(k);

        if ((k = get_key(r, variogram, "ANG2")) == NULL) return NULL;
        d->covariance->ang2[i] = get_float(k);

        if ((k = get_key(r, variogram, "ANG3")) == NULL) return NULL;
        d->covariance->ang3[i] = get_float(k);

        if ((k = get_key(r, variogram, "AA")) == NULL) return NULL;
        aa = d->covariance->aa[i] = get_float(k);
        if (aa < 1e-20) aa = 1e-20;

        if ((k = get_key(r, variogram, "AA1")) == NULL) return NULL;
        aa1 = get_float(k);
        d->covariance->anis1[i] = aa1 / aa; 

        if ((k = get_key(r, variogram, "AA2")) == NULL) return NULL;
        aa2 = get_float(k);
        d->covariance->anis2[i] = aa2 / aa;
    } /* for */

    if ((sill < 1) || (sill > 1)) {
        printf("WARNING: The sill of the variogram is not 1.0! sill = %f\n", sill);
    }

    d->simulation->nsim = 1;   /* number of simulations */

    printf_dbg("new_dss(): Reading harddata file.\n");
    /* read harddata file */
    d->harddata_size = 0;
    if ((k = get_key(r, "HARDDATA", "FILENAME")) == NULL) return NULL;
    if ((fp = open_file(get_string(k))) == NULL) {
        printf_dbg("new_dss(): Can't open Hard Data file: %s\n", get_string(k));
        return NULL;
    }
    if ((d->harddata = load_harddata_file(fp, hdbuf, &d->harddata_size)) == NULL) {
        printf_dbg("new_dss(): failed to load harddata file!\n");
        return NULL;
    }
    close_file(fp);
    
    /* start acorni */
    d->general->ixv = new_acorni(random());

    /* ugly hack for "readdata" */
    d->general->maxdat = d->harddata_size / 4;
    d->general->x = (float *) malloc(d->general->maxdat * sizeof(float));
    d->general->y = (float *) malloc(d->general->maxdat * sizeof(float));
    d->general->z__ = (float *) malloc(d->general->maxdat * sizeof(float));
    d->general->vr = (float *) malloc(d->general->maxdat * sizeof(float));
    d->general->wt = (float *) malloc(d->general->maxdat * sizeof(float));
    d->general->vrtr = (float *) malloc(d->general->maxdat * sizeof(float));
    d->general->vrgtr = (float *) malloc(d->general->maxdat * sizeof(float));
    d->general->sec = (float *) malloc(d->general->maxdat * sizeof(float));
    if (!d->general->x ||
        !d->general->y ||
        !d->general->z__ ||
        !d->general->vr ||
        !d->general->wt ||
        !d->general->vrtr ||
        !d->general->vrgtr ||
        !d->general->sec)
        return NULL;

    printf_dbg("new_dss(): DSS engine started.\n");
    return d;
} /* new_dss */



int setup_dss(dss *d, float *currBAI) {
    d->general->ktype = 1;
    if (currBAI) d->general->ktype = 5;
    return readdata(currBAI, d->harddata, &d->harddata_size, d->general, d->search, d->simulation);
} /* setup_dss */



int run_dss(dss *d, float *AI) {
    int *order, *mask;

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

    /* SIMULATION */
    sdsim(AI, NULL, NULL, order, mask,
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
    return 1;
} /* run_dss */



int run_codss(dss *d, float *currBAI, float *currBCM, float *AI) {
    int *order, *mask;

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

    /* CO-SIMULATION */
    sdsim(AI, currBAI, currBCM, order, mask,
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
    return 1;
} /* run_codss */



void delete_dss(dss *d) {
    if (d) {
        if (d->general) {
            if (d->general->x) free(d->general->x);
            if (d->general->y) free(d->general->y);
            if (d->general->z__) free(d->general->z__);
            if (d->general->vr) free(d->general->vr);
            if (d->general->wt) free(d->general->wt);
            if (d->general->vrtr) free(d->general->vrtr);
            if (d->general->vrgtr) free(d->general->vrgtr);
            if (d->general->sec) free(d->general->sec);
            free(d->general);
        }
        if (d->search) free(d->search);
        if (d->simulation) free(d->simulation);
        if (d->covariance) {
            if (d->covariance->it) free(d->covariance->it);
            if (d->covariance->cc) free(d->covariance->cc);
            if (d->covariance->aa) free(d->covariance->aa);
            if (d->covariance->ang1) free(d->covariance->ang1);
            if (d->covariance->ang2) free(d->covariance->ang2);
            if (d->covariance->ang3) free(d->covariance->ang3);
            if (d->covariance->anis1) free(d->covariance->anis1);
            if (d->covariance->anis2) free(d->covariance->anis2);
            if (d->covariance->aa) free(d->covariance->aa);
            free(d->covariance);
        }
        if (d->clookup) free(d->clookup);
        if (d->krige) free(d->krige);
        if (d->files) free(d->files);
        free(d);
    }
} /* delete_dss */



double *load_harddata_file(TSI_FILE *fp, char *buf, int *size) {
    double x, y, z, val, *ret;
    int i;

    i = *size;
    if (read_line_file(fp, buf) != EOF) {
        sscanf(buf, "%lf %lf %lf %lf", &x, &y, &z, &val);
        (*size)++;
        if (ret = load_harddata_file(fp, buf, size)) {
            ret[i]   = x;
            ret[i+1] = y;
            ret[i+2] = z;
            ret[i+3] = val;
            return ret;
        } else
            return NULL;
    }
    *size *= 4;
    return (double *) malloc(*size * sizeof(double));
} /* load_harddata_file */

