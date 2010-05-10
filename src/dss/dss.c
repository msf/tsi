
#include "debug.h"
#include "memdebug.h"
#include "log.h"
#include "grid_heap.h"
#include "registry.h"
#include "dss.h"
#include "dss_legacy.h"
#include "tsi_io.h"


mask_t * load_mask_file(log_t *l, general_vars_t *general, char *filename)
{
    unsigned int x, y, z;
    int i, m;
    char line[256];
    float *buf;
    TSI_FILE *fp;
    mask_t *mask = NULL;

    if ((fp = fopen(filename, "r")) == NULL) {
        ERROR(l, "fopen()", filename);
        return NULL;
    }

    /* check gslib header */
    if( !read_gslib_header(l, fp, 3) ) {
        ERROR(l, "load_mask_file()", "read_gslib_header()");
        return NULL;
    }

    mask = mask_new( general->nxyz );

    i = 0;
    m = 0;
    while( (m = fscanf(fp,"%u %u %u", &x, &y, &z)) != EOF ) {
        if(m != 3) {
            log_print(l,"load_mask_file(): ERROR -  %d can't parse data values, line: %d\n",m, 5+i);
            continue;
        }
        if(  x >= general->nx || y >= general->ny || z >= general->nz ) {
            log_print(l,"load_mask_file(): ERROR - position outside of grid, line: %d\n", 5+i);
            continue;
        }

        mask_set(mask, getPos(x, y, z, general->nx, general->nxy) );
        i++;
    }

    fclose(fp);
    return mask;
}


dss *new_dss(registry *r, grid_heap *h, log_t *l)
{
    dss *d;
    reg_key *k, *kpath;
    /* auxiliar variables for harddata */
    char filename[512];

    printf_dbg2("new_dss(): called\n");
    /* object space allocation */
    d = (dss *) tsi_malloc(sizeof(dss));
    d->general    = (general_vars_t *) tsi_malloc(sizeof(general_vars_t));
    d->search     = (search_vars_t *) tsi_malloc(sizeof(search_vars_t));
    d->covariance = (covariance_vars_t *) tsi_malloc(sizeof(covariance_vars_t));
    d->clookup    = (covtable_lookup_vars_t *) tsi_malloc(sizeof(covtable_lookup_vars_t));
    d->krige      = (krige_vars_t *) tsi_malloc(sizeof(krige_vars_t));
    d->harddata	  = (harddata_t *) tsi_malloc(sizeof(harddata_t));
    if (!d ||
            !d->general ||
            !d->search ||
            !d->covariance ||
            !d->clookup ||
            !d->krige ||
            !d->harddata )
        return NULL;
    d->reg = r;
    d->heap = h;
    d->l = l;

    printf_dbg2("new_dss(): Starting new DSS engine.\n Loading dss config settings\n");
    if(dss_parameters(d, l, r)){
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
    snprintf(filename, 512, "%s%s", get_string(kpath), get_string(k));

    if ( load_harddata_file(l, filename, d->harddata) ) {
        ERROR( l, "new_dss()", "load_harddata_file()");
        
        return NULL;
    }

    /* parse and prepare harddata */   
    readdata(l, d->harddata, d->general);

    /* read mask data */
    if ((k = get_key(r, "MASK", "FILENAME")) == NULL) {
        d->mask = NULL;
        filename[0] = '\0';
        printf_dbg("new_dss(): not using mask.\n");
    } else {
        snprintf(filename, 512, "%s%s", get_string(kpath), get_string(k));
        log_print(l, "new_dss(): using mask: %s\n", filename);
        d->mask = load_mask_file(l, d->general, filename);
        if( NULL == d->mask ) {
            ERROR( l, "new_dss()", "load_mask_file()");
            log_print(l, "new_dss(): FAILED to load mask, running without mask.\n");
        }
    }
 
    printf_dbg2("new_dss(): DSS engine started sucessfully.\n");
    return d;
} /* new_dss */


int run_dss(dss *d, float *AI) {
    int *order;
    int ktype = 1;

    printf_dbg2("run_dss(): called\n");

    /* SIMULATION */
    printf_dbg2("run_dss(): grids allocated, calling dssim()\n");
    dssim(AI, NULL, NULL, order, ktype,
            d->mask,
            d->general,
            d->harddata,
            d->search,
            d->covariance,
            d->clookup,
            d->krige);                                                                       


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
            d->mask,
            d->general,
            d->harddata,
            d->search,
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
        if (d->covariance) {
            if (d->covariance->variogram) tsi_free(d->covariance->variogram);
            tsi_free(d->covariance);
        }
        if (d->clookup) tsi_free(d->clookup);
        if (d->krige) {
            tsi_free(d->krige);
        }
        if (d->harddata) {
            if(d->harddata->in_grid) tsi_free(d->harddata->in_grid);
            if(d->harddata->point) tsi_free(d->harddata->point);
            tsi_free(d->harddata);
        }
        if( d->mask ) mask_free( d->mask );
        tsi_free(d);
    }
} /* delete_dss */



