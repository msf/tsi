/* si.c */

#include <stdlib.h>
#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
#include "si.h"
#include "si_math.h"


si *new_si(registry *r, grid_heap *h) {
    si *s;
    reg_key *k;
    TSI_FILE *fp;
    int i;
    char buf[64];

    s = (si *) tsi_malloc(sizeof(si));
    if (!s) {
        return NULL;
    }
    s->reg = r;
    s->heap = h;
    s->cmg = NULL;

    /* get grid parameters */
    if ((k = get_key(r, "GRID", "XNUMBER")) == NULL) {
       delete_si(s);
       return NULL;
    }
    s->xsize = get_int(k);

    if ((k = get_key(r, "GRID", "YNUMBER")) == NULL) {
       delete_si(s);
       return NULL;
    }
    s->ysize = get_int(k);

    if ((k = get_key(r, "GRID", "ZNUMBER")) == NULL) {
       delete_si(s);
       return NULL;
    }
    s->zsize = get_int(k);

    s->grid_size = (unsigned int)s->zsize * (unsigned int)s->ysize * (unsigned int)s->xsize;
		
    /* load wavelet */
    if ((k = get_key(r, "WAVELET", "USED_VALUES")) == NULL) {
        printf_dbg("new_si(): failed to get WAVELET:USED_VALUES from registry!\n");
        return NULL;
    }

    s->wavelet_used_values = get_int(k);
    s->points = (int *) tsi_malloc((1 + s->wavelet_used_values) * sizeof(int));
    s->values = (float *) tsi_malloc((1 + s->wavelet_used_values) * sizeof(float));
    if (!s->points || !s->values) {
        printf_dbg("new_si(): failed to allocate space for wavelet!\n");
        return NULL;
    }
    
    if ((k = get_key(r, "WAVELET", "FILENAME")) == NULL) {
        printf_dbg("new_si(): failed to get WAVELET:FILENAME from registry!\n");
        return NULL;
    }
    fp = open_file(get_string(k));
    if (!fp) {
        printf_dbg("new_si(): failed to open wavelet file!\n");
        return NULL;
    }

    i = 0;
    while (read_line_file(fp, buf) != EOF) {
        sscanf(buf, "%d %f\n", &(s->points[i]), &(s->values[i]));
        printf_dbg2("Wavelet[%d] = Pair(%d, %f)\n", i, s->points[i], s->values[i]);
        i++;
    }
    close_file(fp);

    if (i < s->wavelet_used_values) {
        printf_dbg("new_si(): EOF found. Uncomplete wavelet.\nThis wavelet must have %d values.\n", s->wavelet_used_values);
        return NULL;
    }
    
	/* layers data */
    if ((k = get_key(r, "CORR", "LAYERS_MIN")) == NULL) {
       printf_dbg("new_si: failed to get number of layers setting from the registry!\n");
       delete_si(s);
       return NULL;
    }
       s->min_number = get_int(k);

    if ((k = get_key(r, "CORR", "LAYER_SIZE_MIN")) == NULL) {
       printf_dbg("new_si: failed to get size of layers setting from the registry!\n");
       delete_si(s);
       return NULL;
    }
       s->min_size = get_int(k);

    return s;
} /* new_si */



int setup_si(si *s) {
	printf_dbg("setup_si(): called but not implemented\n");
    /* generate new set of layers */
    return 1;
}



int run_si(si *s, float *AI, float *seismic, float *CM, float *SY) {
    int result;

    /* build reflections grid (use CM as aux grid) */
    result = make_reflections_grid(s, AI, CM);

    /* build synthetic grid */
    result = make_synthetic_grid(s, CM, SY);   /* CM has reflections grid */
	
    /* build correlations cube, use CM has resulting correlations Grid */
    result = make_correlations_grid(s, seismic, SY, CM);
    return 1;
} /* si */



void delete_si(si *s) {
    if (s) {
       if (s->points) tsi_free(s->points);
       if (s->values) tsi_free(s->values);
       tsi_free(s);
    }
} /* delete_si */


/* end of file si.c */
