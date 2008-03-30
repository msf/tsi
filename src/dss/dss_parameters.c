#include <stdlib.h>
#include "debug.h"
#include "memdebug.h"
#include "grid_heap.h"
#include "registry.h"
#include "dss.h"
#include "dss_legacy.h"
#include "tsi_io.h"

int  dss_parameters(dss *d, registry *r )
{
    int i, varnum;
    float aa1, aa2, aa, sill;
    reg_key *k;

    /* harddata parameters */
    if ((k = get_key(r, "HARDDATA", "TMIN")) == NULL) return 1;
    d->harddata->min_value = get_float(k);

    if ((k = get_key(r, "HARDDATA", "TMAX")) == NULL) return 1;
    d->harddata->max_value = get_float(k);

    if (get_key(r, "HARDDATA", "NVARI"))
    	printf("[HARDDATA] NVARI is NOT used\n");
    
    if (get_key(r, "HARDDATA", "IXL"))
    	printf("[HARDDATA] IXL is NOT used\n");
    
     if (get_key(r, "HARDDATA", "IYL"))
    	printf("[HARDDATA] IYL is NOT used\n");
    
    if (get_key(r, "HARDDATA", "IZL"))
    	printf("[HARDDATA] IZL is NOT used\n");
    
    if (get_key(r, "HARDDATA", "IVRL"))
    	printf("[HARDDATA] IVRL is NOT used\n");
    
    if (get_key(r, "HARDDATA", "IWT"))
    	printf("[HARDDATA] IWT is NOT used\n");
    
    if (get_key(r, "HARDDATA", "ISECVR"))
    	printf("[HARDDATA] ISECVR is NOT used\n");
    

    /* harddata transformations */
    if (get_key(r, "HDTRANS", "ITRANS"))
    	printf("[HDTRANS] ITRANS is NOT used\n");
 
    if (get_key(r, "HDTRANS", "ISMOOTH"))
    	printf("[HDTRANS] ISMOOTH is NOT used\n");

    if (get_key(r, "HDTRANS", "ISVR"))
    	printf("[HDTRANS] ISVR is NOT used\n");
    
    if (get_key(r, "HDTRANS", "ISWT"))
    	printf("[HDTRANS] ISWT is NOT used\n");

    if (get_key(r, "HDTRANS", "LTAIL"))
    	printf("[HDTRANS] LTAIL is NOT used\n");

    if (get_key(r, "HDTRANS", "UTAIL"))
    	printf("[HDTRANS] UTAIL is NOT used\n");
    
    /* simulations quality */
    if ((k = get_key(r, "QUALITY", "NTRY")) == NULL) return 1;
    d->clookup->ntry = get_int(k);

    if (get_key(r, "QUALITY", "ICVAR"))
    	printf("[QUALITY] ICVAR is NOT used\n");

    if (get_key(r, "QUALITY", "ICMEAN"))
    	printf("[QUALITY] ICMEAN is NOT used\n");

    
    /* grid parameters */
    if ((k = get_key(r, "GRID", "XNUMBER")) == NULL) return 1;
    d->general->nx = get_int(k);

    if ((k = get_key(r, "GRID", "YNUMBER")) == NULL) return 1;
    d->general->ny = get_int(k);

    if ((k = get_key(r, "GRID", "ZNUMBER")) == NULL) return 1;
    d->general->nz = get_int(k);

    if ((k = get_key(r, "GRID", "XCOORD")) == NULL) return 1;
    d->general->xmn = get_float(k);

    if ((k = get_key(r, "GRID", "YCOORD")) == NULL) return 1;
    d->general->ymn = get_float(k);

    if ((k = get_key(r, "GRID", "ZCOORD")) == NULL) return 1;
    d->general->zmn = get_float(k);

    if ((k = get_key(r, "GRID", "XSIZE")) == NULL) return 1;
    d->general->xsiz = get_float(k);

    if ((k = get_key(r, "GRID", "YSIZE")) == NULL) return 1;
    d->general->ysiz = get_float(k);

    if ((k = get_key(r, "GRID", "ZSIZE")) == NULL) return 1;
    d->general->zsiz = get_float(k);

    d->general->nxy = d->general->nx * d->general->ny;
    d->general->nxyz = d->general->nxy * d->general->nz;

    
    /* mask */
    if ((k = get_key(r, "MASK", "NULL_VALUE")) == NULL) return 1;
    d->general->nosim_value = get_float(k);


    /* search parameters */
    if ((k = get_key(r, "SEARCH", "NDMIN")) == NULL) return 1;
    d->search->ndmin = get_int(k);

    if (get_key(r, "SEARCH", "NDMAX"))
    	printf("[SEARCH] NDMAX is NOT used\n");
    d->search->ndmax = 0;

    if ((k = get_key(r, "SEARCH", "NODMAX")) == NULL) return 1;
    d->clookup->nodmax = get_int(k);

    if (get_key(r, "SEARCH", "SSTRAT"))
    	printf("[SEARCH] SSTRAT is NOT used\n");

    if (get_key(r, "SEARCH", "NOCT"))
    	printf("[SEARCH] NOCT is NOT used\n");

    if ((k = get_key(r, "SEARCH", "RADIUS")) == NULL) return 1;
    d->search->radius = get_float(k);
    if (d->search->radius < 1e-20) {
        printf("load_dss_configs(): ERROR - radius (%f) must be greater than EPSLON (1E-20)\n", d->search->radius);
        return 1;
    }
    d->search->radsqd = d->search->radius * d->search->radius;

    if ((k = get_key(r, "SEARCH", "RADIUS1")) == NULL) return 1;
    d->search->sanis1 = get_float(k) / d->search->radius;

    if ((k = get_key(r, "SEARCH", "RADIUS2")) == NULL) return 1;
    d->search->sanis2 = get_float(k) / d->search->radius;

    if ((k = get_key(r, "SEARCH", "SANG1")) == NULL) return 1;
    d->search->sang1 = get_float(k);

    if ((k = get_key(r, "SEARCH", "SANG2")) == NULL) return 1;
    d->search->sang2 = get_float(k);

    if ((k = get_key(r, "SEARCH", "SANG3")) == NULL) return 1;
    d->search->sang3 = get_float(k);


    /* kriging parameters */
    if (get_key(r, "KRIG", "TYPE"))
    	printf("[KRIG] TYPE is NOT used\n");

    if (get_key(r, "KRIG", "COLOCORR"))
    	printf("[KRIG] COLOCORR is NOT used\n");

    if (get_key(r, "KRIG", "VARRED"))
    	printf("[KRIG] VARRED is NOT used\n");

    /* softdata parameters */
    if (get_key(r, "SOFT", "NVARIL"))
    	printf("[SOFT] NVARIL is NOT used\n");
    
    if (get_key(r, "SOFT", "ICOLLVM"))
    	printf("[SOFT] ICOLLVM is NOT used\n");

    /* variogram models */
    if ((k = get_key(r, "VARIOGRAM", "NUMBER")) == NULL) return 1;
    varnum = d->covariance->varnum = get_int(k);
    printf_dbg2("load_dss_configs(): varnum = %d\n", varnum);

    if ((k = get_key(r, "VARIOGRAM", "NUGGET")) == NULL) return 1;
    sill = d->covariance->nugget = get_float(k);


	variogram_t *variogram = (variogram_t *) tsi_malloc(varnum * sizeof(variogram_t));

    if (!variogram ) {
        printf("load_dss_configs(): ERROR - Failed to allocate space for variograms!\n");
        return 1;
    }

    char varname[16];

    for (i = 0; i < varnum; i++) {
        sprintf(varname, "VARIOGRAM%d", i+1);
        if ((k = get_key(r, varname, "TYPE")) == NULL) return 1;
        variogram[i].type = get_int(k);
        
		printf_dbg2("%s - type: %d\t %d\n", varname, get_int(k), variogram[i].type);
        if (variogram[i].type == 4) {
            printf("load_dss_configs(): ERROR - A power model is not allowed! Choose a different model and re-start.\n");
            return 1;
        }

        if ((k = get_key(r, varname, "COV")) == NULL) return 1;
        variogram[i].cov = get_float(k);
        sill += variogram[i].cov;

        if ((k = get_key(r, varname, "ANG1")) == NULL) return 1;
        variogram[i].ang1 = get_float(k);

        if ((k = get_key(r, varname, "ANG2")) == NULL) return 1;
        variogram[i].ang2 = get_float(k);

        if ((k = get_key(r, varname, "ANG3")) == NULL) return 1;
        variogram[i].ang3 = get_float(k);

        if ((k = get_key(r, varname, "AA")) == NULL) return 1;
        aa = variogram[i].aa = get_float(k);
        if (aa < 1e-20)
			aa = (float) 1e-20;

        if ((k = get_key(r, varname, "AA1")) == NULL) return 1;
        aa1 = get_float(k);
        variogram[i].anis1 = aa1 / aa;

        if ((k = get_key(r, varname, "AA2")) == NULL) return 1;
        aa2 = get_float(k);
        variogram[i].anis2 = aa2 / aa;
        
		printf_dbg2("load_dss_configs(): variogram %d\n \ttype:\t%d\n \tcov:\t%.2f\n \tang1:\t%.2f\n \tang2:\t%.2f\n \tang3:\t%.2f\n \taa:\t%.2f\n \taa1:\t%.2f\n \taa2:\t%.2f\n\n",
				i,
				variogram[i].type,
				variogram[i].cov,
				variogram[i].ang1,
				variogram[i].ang2,
				variogram[i].ang3,
				aa,aa1,aa2);
    } /* for */
    d->covariance->variogram = variogram;

    if ((sill < 1) || (sill > 1)) {
        printf("load_dss_configs(): WARNING: The sill of the variogram is not 1.0! sill = %f\n", sill);
    }

    //d->simulation->nsim = d->simulation->nsim_bk = 1;   /* number of simulations */

    return 0;
}

