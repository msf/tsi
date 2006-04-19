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

    printf_dbg("\tnew_si(): called\n");
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
        printf_dbg("\tnew_si(): failed to get WAVELET:USED_VALUES from registry!\n");
        return NULL;
    }

    s->wavelet_used_values = get_int(k);
    s->points = (int *) tsi_malloc((1 + s->wavelet_used_values) * sizeof(int));
    s->values = (float *) tsi_malloc((1 + s->wavelet_used_values) * sizeof(float));
    if (!s->points || !s->values) {
        printf_dbg("\tnew_si(): failed to allocate space for wavelet!\n");
        return NULL;
    }
    
    if ((k = get_key(r, "WAVELET", "FILENAME")) == NULL) {
        printf_dbg("\tnew_si(): failed to get WAVELET:FILENAME from registry!\n");
        return NULL;
    }
    fp = open_file(get_string(k));
    if (!fp) {
        printf_dbg("\tnew_si(): failed to open wavelet file!\n");
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
        printf_dbg("\tnew_si(): EOF found. Uncomplete wavelet.\nThis wavelet must have %d values.\n", s->wavelet_used_values);
        return NULL;
    }
    
	/* layers data */
    if ((k = get_key(r, "CORR", "RANDOM_LAYERS")) == NULL) {
       printf_dbg("\tnew_si: failed to get random layers flag from the registry! Using defaults.\n");
       s->random = 1;
    } else {
       s->random = get_int(k);
    }

    if ((k = get_key(r, "CORR", "LAYERS_MIN")) == NULL) {
       printf_dbg("\tnew_si: failed to get number of layers setting from the registry!\n");
       delete_si(s);
       return NULL;
    }
    s->min_number = get_int(k);

    if ((k = get_key(r, "CORR", "LAYER_SIZE_MIN")) == NULL) {
       printf_dbg("\tnew_si: failed to get size of layers setting from the registry!\n");
       delete_si(s);
       return NULL;
    }
    s->min_size = get_int(k);

    return s;
} /* new_si */



int setup_si(si *s) {
    printf_dbg("\tsetup_si(): called but not implemented\n");
    /* generate new set of layers */
    return 1;
}



int run_si(si *s, float *AI, float *seismic, float *CM, float *SY) {
    int result;
    float *RG;

    RG = CM;
    printf_dbg("\trun_si(): called\n");

    /* build reflections grid (use CM as aux grid) */
    result = make_reflections_grid(s, AI, RG);

    /* build synthetic grid */
    result = make_synthetic_grid(s, RG, SY);   /* CM has reflections grid */
	
    /* build correlations cube, use CM has resulting correlations Grid */
    result = make_correlations_grid(s, seismic, SY, CM);
    return 1;
} /* si */



void delete_si(si *s) {
    printf_dbg("\tdelete_si(): called\n");
    if (s) {
       if (s->points) tsi_free(s->points);
       if (s->values) tsi_free(s->values);
       if (s->cmg) delete_cmgrid(s->cmg);
       tsi_free(s);
    }
} /* delete_si */



cm_grid *new_cmgrid(si *s, int empty) {
    cm_grid *g;
    int delta, i, x, n, sum, *temp;
    
    printf_dbg("\tnew_cmgrid(): called\n");
    if ((g = tsi_malloc(sizeof(cm_grid))) == NULL) {
        printf_dbg("\tnew_cmgrid: failed to allocate space for cm_grid\n");
        return NULL;
    }
    g->nlayers = 0;
    g->layer_size = NULL;
    g->cg = NULL;
    g->nxy = s->xsize * s->ysize;
    if (empty) return g;

    if (!s->random) {
        /* generate static layers */
        g->nlayers = s->min_number;
        n = s->zsize / s->min_number;
        x = s->zsize % s->min_number;

        if ((g->layer_size = (int *) tsi_malloc(sizeof(int) * s->min_number)) == NULL) {
            printf_dbg2("new_cmgrid: failed to allocate layers array");
            delete_cmgrid(g);
            return NULL;
        }

        for (i = 0; i < s->min_number; i++) g->layer_size[i] = n;
        g->layer_size += x;

        if ((g->cg = (float *) tsi_malloc(g->nlayers * g->nxy * sizeof(float))) == NULL) {
            printf_dbg2("new_cmgrid: failed to allocate cg array\n");
            delete_cmgrid(g);
            return NULL;
        }
        return g;
    }

    /* generate random layers */
    delta = s->zsize - (s->min_size * s->min_number);
    if (delta < 0) {
        /* too many layers, or layers too big */
        fprintf(stderr, "ERROR, too many layers, or layers too big\n");
        delete_cmgrid(g);
        return NULL;
    } else if (delta == 0) {
        /* all layers must be of min_size */
        g->nlayers = s->min_number;
        if ((g->layer_size = (int *) tsi_malloc(sizeof(int) * g->nlayers)) == NULL) {
            printf_dbg2("new_cmgrid: failed to allocate layers array");
            delete_cmgrid(g);
            return NULL;
        }
        for (i = 0; i < g->nlayers; i++) g->layer_size[i] = s->min_size;
        if ((g->cg = (float *) tsi_malloc(g->nlayers * g->nxy * sizeof(float))) == NULL) {
            printf_dbg2("new_cmgrid: failed to allocate cg array\n");
            delete_cmgrid(g);
            return NULL;
        }
        return g;
    }

    /* proper random layers...
     *
     * number of layers: number < N < max_size / min_size
     */	
    n = s->zsize / s->min_size;
    temp = (int *) tsi_malloc(n * sizeof(int));
    for (i = 0; i < n; i++) temp[i] = 0;

    /* first, we warrantee the minimum number of layers */
    do {
        x = random() % (s->min_size * 2);
        x += s->min_size;
        i = (int) random() % n;
        temp[i] = x;

        x = 0;
        for(i = 0; i < n; i++) {
            if(temp[i] != 0) x++;
        }
    } while (x < s->min_number);

    /* now we assure the max_size */
    do {
        x = random() % (s->min_size);
        i = random() % n;
        if (temp[i] == 0) temp[i] += s->min_size;
        temp[i] += x;
        sum = 0;
        for (i = 0; i < n; i++) sum += temp[i];
    } while (sum < s->zsize);

    while (sum > s->zsize) {
        /* lets find the biggest layer */
        sum = 0;
        for (i = 0; i < n; i++) {
            if (temp[i] > sum) {
                sum = temp[i];
                x = i;
            }
        }
        /* now we reduce its size */
        temp[x] -= s->min_size;
        temp[x] /= 2;
        temp[x] += s->min_size;
        /* lets see if its good now */
        sum = 0;
        for (i = 0; i < n; i++) sum += temp[i];
    }
    
    if (sum != s->zsize) {
          /* add the diference to max_size */
          i = s->zsize - sum;
          temp[x] += i;
    }

    /* now lets count the layers != 0 */
    sum = 0;
    for (i = 0; i < n; i++)
        if (temp[i] != 0)
            sum++;

    /* this is the number of layers we are going to return */
    if ((g->layer_size = (int *) tsi_malloc(sizeof(int) * sum)) == NULL) {
		printf_dbg2("\tnew_cmgrid: failed to allocate layers array\n");
		delete_cmgrid(g);
		return NULL;
	}
    g->nlayers = sum;
	
    x = 0;
    for (i = 0; i < n; i++)
        if (temp[i] != 0)
            g->layer_size[x++] = temp[i];
    
    tsi_free(temp);

    if ((g->cg = (float *) tsi_malloc(g->nlayers * g->nxy * sizeof(float))) == NULL) {
        printf_dbg2("\tnew_cmgrid: failed to allocate cg array\n");
        delete_cmgrid(g);
        return NULL;
    }
    return g;
} /* new_cmgrid */



int build_cmgrid(cm_grid *g, int nlayers, int *layers) {
    printf_dbg("\tbuild_cmgrid(): called\n");
    if (g == NULL) {
        printf_dbg2("\tbuild_cmgrid: unexpected NULL for cm_grid\n");
        return 0;
    }
    g->nlayers = nlayers;
    g->layer_size = layers;
    if ((g->cg = (float *) tsi_malloc(nlayers * g->nxy * sizeof(float))) == NULL) {
        printf_dbg2("\tbuild_cmgrid: failed to allocate cg array\n");
        return 0;
    }
    return 1;
} /* build_cmgrid */



cm_grid *load_cmgrid(si *s) {
    printf_dbg("\tload_cmgrid(): called\n");
    return s->cmg;
} /* load_cmgrid */



int get_nlayers(cm_grid *g) {
    printf_dbg("\tget_nlayers(): called\n");
    if (g) return g->nlayers;
    return 0;
} /* load_cmgrid */



int *get_layers(cm_grid *g) {
    printf_dbg("\tget_layers(): called\n");
    if (g) return g->layer_size;
    return NULL;
} /* load_cmgrid */



int store_cmgrid(si *s, cm_grid *g) {
    //if (s->cmg) delete_cmgrid(s->cmg);   // fuck-up... to fix...
    s->cmg = g;
} /* store_cmgrid */



void delete_cmgrid(cm_grid *g) {
    printf_dbg("\tdelete_cmgrid(): called\n");
    if (g) {
        if (g->layer_size) tsi_free(g->layer_size);
        if (g->cg) tsi_free(g->cg);
        tsi_free(g);
    }
} /* delete_cmgrid */



void print_layers(cm_grid *g) {
    int i;
    for (i = 0; i < g->nlayers; i++)
         printf("layer[%d] = %d\n",i,g->layer_size[i]);
} /* print_layers*/

/* end of file si.c */
