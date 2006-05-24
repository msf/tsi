/* si.c */

#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "log.h"
#include "registry.h"
#include "grid_heap.h"
#include "si.h"
#include "si_math.h"


si *new_si(registry *r, grid_heap *h, log_t *l) {
    si *s;
    reg_key *k, *kpath;
    TSI_FILE *fp;
    char filename[512];
    int i;
    char buf[64];

    printf_dbg("\tnew_si(): called\n");
    s = (si *) tsi_malloc(sizeof(si));
    if (!s) {
        return NULL;
    }
    s->reg = r;
    s->heap = h;
    s->l = l;
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
    
    if ((kpath = get_key(r, "WAVELET", "PATH")) == NULL)
        kpath = get_key(r, "GLOBAL", "INPUT_PATH");
    if ((k = get_key(r, "WAVELET", "FILENAME")) == NULL) {
        printf_dbg("\tnew_si(): failed to get WAVELET:FILENAME from registry!\n");
        return NULL;
    }
    if (kpath)
        sprintf(filename, "%s%s", get_string(kpath), get_string(k));
    else
        sprintf(filename, "%s", get_string(k));
    fp = open_file(filename);
    if (!fp) {
        printf_dbg("\tnew_si(): failed to open wavelet file!\n");
        return NULL;
    }

    i = 0;
    while (fscanf(fp,"%d %f\n", &(s->points[i]), &(s->values[i])) != EOF) {
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
    printf_dbg2("\trun_si(): called\n");

    /* build reflections grid (use CM as aux grid) */
    result = make_reflections_grid(s, AI, RG);

    /* build synthetic grid */
    result = make_synthetic_grid(s, RG, SY);   /* CM has reflections grid */
	
    /* build correlations cube, use CM has resulting correlations Grid */
    result = make_correlations_grid(s, seismic, SY, CM);
    return 1;
} /* si */



void delete_si(si *s) {
    printf_dbg2("\tdelete_si(): called\n");
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
    
    printf_dbg2("\tnew_cmgrid(): called\n");
    if ((g = tsi_malloc(sizeof(cm_grid))) == NULL) {
        printf_dbg("\tnew_cmgrid: failed to allocate space for cm_grid\n");
        return NULL;
    }
    g->nlayers = 0;
    g->layer_size = NULL;
    g->cg = NULL;
    g->cg_idx = -1;
    g->heap = s->heap;
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

        //if ((g->cg = (float *) memalign(16, g->nlayers * g->nxy * sizeof(float))) == NULL) {
        //if ((g->cg = (float *) tsi_malloc(g->nlayers * g->nxy * sizeof(float))) == NULL) {
        if ((g->cg_idx = new_grid(g->heap)) == -1) {
            printf_dbg2("new_cmgrid: failed to allocate cg array\n");
            delete_cmgrid(g);
            return NULL;
        }
        set_grid_size(g->heap, g->cg_idx, g->nxy*g->nlayers);
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
        //if ((g->cg = (float *) memalign(16, g->nlayers * g->nxy * sizeof(float))) == NULL) {
        //if ((g->cg = (float *) tsi_malloc(g->nlayers * g->nxy * sizeof(float))) == NULL) {
        if ((g->cg_idx = new_grid(g->heap)) == -1) {
            printf_dbg2("new_cmgrid: failed to allocate cg array\n");
            delete_cmgrid(g);
            return NULL;
        }
        set_grid_size(g->heap, g->cg_idx, g->nxy*g->nlayers);
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
        /* reduce its size */
        temp[x] -= s->min_size;
        temp[x] /= 2;
        temp[x] += s->min_size;
        /* see if its good now */
        sum = 0;
        for (i = 0; i < n; i++) sum += temp[i];
    }
    
    if (sum != s->zsize) {
          /* add the diference to max_size */
          i = s->zsize - sum;
          temp[x] += i;
    }

    /* count the layers != 0 */
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

    //if ((g->cg = (float *) memalign(16,g->nlayers * g->nxy * sizeof(float))) == NULL) {
    //if ((g->cg = (float *) tsi_malloc(g->nlayers * g->nxy * sizeof(float))) == NULL) {
    if ((g->cg_idx = new_grid(g->heap)) == -1) {
        printf_dbg2("\tnew_cmgrid: failed to allocate cg array\n");
        delete_cmgrid(g);
        return NULL;
    }
    set_grid_size(g->heap, g->cg_idx, g->nxy*g->nlayers);
    return g;
} /* new_cmgrid */



int build_cmgrid(cm_grid *g, int nlayers, int *layers) {
    printf_dbg2("\tbuild_cmgrid(): called\n");
    if (g == NULL) {
        printf_dbg2("\tbuild_cmgrid: unexpected NULL for cm_grid\n");
        return 0;
    }
    g->nlayers = nlayers;
    g->layer_size = layers;
    //if ((g->cg = (float *) memalign(16, nlayers * g->nxy * sizeof(float))) == NULL) {
    //if ((g->cg = (float *) tsi_malloc(nlayers * g->nxy * sizeof(float))) == NULL) {
    if ((g->cg_idx = new_grid(g->heap)) == -1) {
        printf_dbg2("\tbuild_cmgrid: failed to allocate cg array\n");
        return 0;
    }
    set_grid_size(g->heap, g->cg_idx, g->nxy*g->nlayers);
    return 1;
} /* build_cmgrid */



cm_grid *get_cmgrid(si *s) {
    printf_dbg2("\tget_cmgrid(): called\n");
    return s->cmg;
} /* get_cmgrid */



int load_cmgrid(cm_grid *g) {
    printf_dbg2("\tload_cmgrid(): called\n");
    if ((g->cg = load_grid(g->heap, g->cg_idx)) == NULL) {
        printf_dbg("\tload_cmgrid: failed to load cg grid\n");
        return 1;
    }
    return 0;
} /* load_cmgrid */



int clear_cmgrid(cm_grid *g) {
    printf_dbg2("\tclear_cmgrid(): called\n");
    clear_grid(g->heap, g->cg_idx);
    return 0;
} /* clear_cmgrid */



int dirty_cmgrid(cm_grid *g) {
    printf_dbg2("\tdirty_cmgrid(): called\n");
    dirty_grid(g->heap, g->cg_idx);
    return 0;
} /* dirty_cmgrid */



cm_grid *clone_cmgrid(cm_grid *g) {
    cm_grid *clone;

    clone = (cm_grid *) tsi_malloc(sizeof(cm_grid));
    clone->nlayers = g->nlayers;
    clone->nxy = g->nxy;
    clone->heap = g->heap;
    clone->layer_size = (int *) tsi_malloc(g->nlayers * sizeof(int));
    //clone->cg = (float *) memalign(16, g->nlayers * g->nxy * sizeof(float));
    clone->cg = (float *) tsi_malloc(g->nlayers * g->nxy * sizeof(float));
    clone->cg_idx = new_grid(clone->heap);
    if (clone->cg_idx < 0) {
        printf_dbg("clone_cmgrid(): failed to allocate grid for the clone\n");
        delete_cmgrid(clone);
        return NULL;
    }
    memcpy(clone->layer_size, g->layer_size, g->nlayers * sizeof(int));
    load_cmgrid(g);
    load_cmgrid(clone);
    memcpy(clone->cg, g->cg, g->nlayers * g->nxy * sizeof(float));
    clear_cmgrid(g);
    dirty_cmgrid(clone);
    return clone;
} /* clone_cmgrid */



int get_nlayers(cm_grid *g) {
    printf_dbg2("\tget_nlayers(): called\n");
    if (g) return g->nlayers;
    return 0;
} /* load_cmgrid */



int *get_layers(cm_grid *g) {
    printf_dbg2("\tget_layers(): called\n");
    if (g) return g->layer_size;
    return NULL;
} /* load_cmgrid */



int store_cmgrid(si *s, cm_grid *g) {
    //if (s->cmg) delete_cmgrid(s->cmg);   // fuck-up... to fix...
    s->cmg = g;
    return 1;
} /* store_cmgrid */



void delete_cmgrid(cm_grid *g) {
    printf_dbg2("\tdelete_cmgrid(): called\n");
    if (g) {
        if (g->layer_size) tsi_free(g->layer_size);
        //if (g->cg) free(g->cg);
        //if (g->cg) tsi_free(g->cg);
        if (g->cg_idx > -1) delete_grid(g->heap, g->cg_idx);
        tsi_free(g);
    }
} /* delete_cmgrid */



void print_layers(cm_grid *g) {
    int i;
    for (i = 0; i < g->nlayers; i++)
         printf_dbg("layer[%d] = %d\n",i,g->layer_size[i]);
} /* print_layers*/

/* end of file si.c */