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

    s = (si *) malloc(sizeof(si));
    if (!s) {
        return NULL;
    }
    s->reg = r;
    s->heap = h;

    /* load wavelet */
    if ((k = get_key(r, "WAVELET", "USED_VALUES")) == NULL) {
        printf_dbg("new_si(): failed to get WAVELET:USED_VALUES from registry!\n");
        return NULL;
    }

    s->wavelet_used_values = get_int(k);
    s->points = (int *) malloc((1 + s->wavelet_used_values) * sizeof(int));
    s->values = (float *) malloc((1 + s->wavelet_used_values) * sizeof(float));
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

    return s;
} /* new_si */



int setup_si(si *s) {
    /* generate new set of layers */
    return 1;
}



int run_si(si *s, float *AI, float *seismic, float *CM, float *SY) {
    int result;

    /* build reflections grid (use CM as aux grid) */
    result = make_reflections_grid(s, AI, CM);

    /* build synthetic grid */
    result = make_synthetic_grid(s, CM, SY);   /* CM has reflections grid */

    /* build correlations cube */
    result = make_correlations_grid(s, seismic, SY);
    return 1;
} /* si */



void delete_si(si *s) {
    if (s) {
       if (s->points) free(s->points);
       if (s->values) free(s->values);
       free(s);
    }
} /* delete_si */


/* end of file si.c */
