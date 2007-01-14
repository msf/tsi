#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "grid_heap.h"
#include "debug.h"
#include "tsi_io.h"


grid_heap *new_heap(int nodes, int rank, int heap_size, int swap_thr, int use_fs, char *path, unsigned int grid_size) {
    grid_heap *h;
    int i;
    char filename[128];
    unsigned int fragment_size;
    
    /* initializes the heap record */
    h = (grid_heap *) tsi_malloc(sizeof(grid_heap));
    if (!h) {
        printf_dbg("new_heap(): failed to allocate heap!\n");
        return NULL;
    }

    /* machine data */
    h->nodes = nodes;
    h->rank = rank;

    /* heap constants */
    h->heap_size = heap_size;   /* number of grids in heap */

    /* ajust grid size to the cluster size */
    fragment_size = (grid_size / nodes) + ((grid_size % nodes) > 0 ? 1:0);
    h->grid_size = nodes * fragment_size;

    h->threshold = swap_thr;
    h->use_fs = use_fs;

    /* free grids stack */
    h->next_grid = 0;      /* free grids stack */
    
    /* grid counters*/
    h->alloc_grids = 0;    /* allocated grids counter */
    h->curr_grids = 0;     /* grids in use counter */
    h->reads = 0;
    h->writes = 0;
    h->max_grids = 0;
    
    /* initializes the grid records array */
    h->g = (grid *) tsi_malloc((size_t)heap_size*sizeof(grid));
    memset(h->g, 0, heap_size*sizeof(grid));
    if (!h->g) {
        printf_dbg("new_heap(): failed to allocate grid data array!\n");
        delete_heap(h);
        return NULL;
    } /* if */

    /* create the swap files for each grid */
    for(i = 0; i < heap_size; i++) {
        if (h->use_fs) {
            if (path)
                sprintf(filename, "%stsi_grid.node%d.%d", path, rank, i);
            else /* fall back to /tmp */
                sprintf(filename, "/tmp/tsi_grid.node%d.%d", rank, i);
            printf_dbg2("new_heap(): filename string >%s<\n", filename);
            h->g[i].fp = create_file(filename);
            h->g[i].filename = strdup(filename);
            
            if (!h->g[i].fp) {
                printf_dbg("new_heap(): failed to create grid swap file!\n");
                delete_heap(h);
                return NULL;
            } /* if */
        }
        h->g[i].next_grid = i+1;  /* set free grid stack */
        h->g[i].swappable = 1;
        h->g[i].valid = 0;
    } /* for */

    printf_dbg2("new_heap(%d): heap started\n", rank);
    return h;
} /* new_heap */



/* allocates a new grid */
int new_grid(grid_heap *h) {
    grid *g;
    int idx;

    idx = 0;
    
    if (h->alloc_grids < h->heap_size) {       /* check if there are any grids available */

        h->alloc_grids++;        /* increase the number of allocated grids */
        if (h->alloc_grids > h->max_grids) 
            h->max_grids = h->alloc_grids;
        idx = h->next_grid;      /* get index from free grids stack */

        g = &(h->g[idx]);
        h->next_grid = g->next_grid;   /* update free grids stack */
        g->dirty = 1;    /* mark new grid to prevent being cleared and overwritten */
        g->size = 0;     /* reset size to default */
        g->valid = 1;    /* mark grid has "allocated" */
        g->first_load = 1;

        printf_dbg2("new_grid(): grid %d, next %d\n", idx, h->next_grid);
        printf_dbg2("new_grid(): curr_grids=%d alloc_grids=%d\n", h->curr_grids, h->alloc_grids);
    } else {
        printf_dbg("new_grid(): no more grids available!\n");
        return -1;
    }

    return idx;
} /* new_grid */



void set_grid_size(grid_heap *h, int idx, unsigned int size)
{
    if (h->g[idx].size) {
        printf_dbg("set_grid_size(): re-sizing grid!!!");
    }
    h->g[idx].size = size;
}



/* return a pointer to the grid */
float *load_grid(grid_heap *h, int idx) {
    int i, j;
    
    printf_dbg2("load_grid(): grid req=%d curr_grids=%d alloc_grids=%d\n", idx, h->curr_grids, h->alloc_grids);
    if (idx < h->heap_size && idx >= 0) {

        if (!h->g[idx].valid ) {
            printf_dbg("load_grid(): ERROR - grid %d is not a valid grid - not created with new_grid\n", idx);
            return NULL;
        }
			
        if (h->g[idx].grid) {   /* check if this grid still has its space allocated */
            h->g[idx].swappable = 0;   /* mark grid in use flag */
            if (h->g[idx].first_load) {
                h->g[idx].first_load = 0;
                memset(h->g[idx].grid, 0, h->grid_size*sizeof(float));
            }
            return h->g[idx].grid;     /* return pointer to grid space */
        } else if (h->curr_grids < h->threshold) {  /* if #grids is bellow threshold */
            h->curr_grids++;
            printf_dbg2("Number of grids raised to %d\n", h->curr_grids);
            //h->g[idx].grid = (float *) memalign(16, h->grid_size*sizeof(float));
            h->g[idx].grid = (float *) tsi_malloc(h->grid_size*sizeof(float));
            memset(h->g[idx].grid, 0, h->grid_size*sizeof(float));
            h->g[idx].pointer = h->g[idx].grid;
            h->g[idx].swappable = 0;
            h->g[idx].first_load = 0;
            return h->g[idx].grid;
        }
        /* we don't have grid in memory, and we've exceeded the threshold, lets try swap some unused grid */
        
        /* seek for a cleared grid */
        j = h->heap_size;
        for (i = 0; i < h->heap_size; i++) {
            if (h->g[i].swappable && h->g[i].grid) {
                j = i;
                if (!h->g[i].dirty) /* if we find a clean, unused grid, lets use it! */
                    break;
            }
        }
        printf_dbg2("load_grid(): j = %d, i = %d\n", j, i);

        if ((h->use_fs) && (i != j)) { /* found dirty grid */
            /* dump j grid */
            h->writes++;
            printf_dbg2("load_grid(): swapping dirty grid %d to disk\n", j);
            if (!dump_binary_grid(h->g[j].fp, h->g[j].grid, (h->g[j].size ? h->g[j].size : h->grid_size))) {
                printf_dbg("load_grid(): failed to dump grid %d\n", j);
                return NULL;
            }
            h->g[j].dirty = 0;
            i = j;
        }
		
        if (i == h->heap_size) {
            /* no grid available, allocate more ram */
            h->curr_grids++;
            printf_dbg("Number of grids raised to %d\n", h->curr_grids);
            h->g[idx].grid = (float *) tsi_malloc(h->grid_size*sizeof(float));
            //h->g[idx].grid = (float *) memalign(16, h->grid_size*sizeof(float));
            memset(h->g[idx].grid, 0, h->grid_size*sizeof(float));
            h->g[idx].first_load = 0;
            h->g[idx].pointer = h->g[idx].grid;
            h->g[idx].swappable = 0;
        } else {
            /* we have a clean available grid, use it! */
            h->g[idx].grid = h->g[j].grid;
            h->g[idx].pointer = h->g[j].pointer;
            h->g[idx].swappable = 0;
            h->g[j].grid = h->g[j].pointer = NULL;
            if (h->g[idx].first_load) {
                h->g[idx].first_load = 0;
                memset(h->g[idx].grid, 0, h->grid_size*sizeof(float));
            }
        }
        
        if ((h->use_fs) && (!h->g[idx].dirty)) {
            /* grid was swapped for disk, we need to load grid from disk */
            h->reads++;
            printf_dbg2("load_grid(): loading grid %d to memory\n", idx);
            if (!load_binary_grid(h->g[idx].fp, h->g[idx].grid, (h->g[idx].size ? h->g[idx].size : h->grid_size))) {
                 printf_dbg("load_grid(): failed to load grid %d\n", idx);
                 return NULL;
            }
        }
        return h->g[idx].grid;
    } else {
        printf_dbg("load_grid(): grid %d off range\n", idx);
    }
    return NULL;
} /* load_grid */



int clear_grid(grid_heap *h, int idx) {
    if ((idx < h->heap_size) && (idx >= 0)) {
        if (!h->g[idx].valid) {
            printf("dirty_grid(): ERROR - grid %d is not a valid grid - not created with new_grid\n",idx);
            return 1;
        }
        if (h->g[idx].grid == NULL) {
            printf("dirty_grid(): ERROR - grid %d space is NULL\n",idx);
            return 1;
        }
        h->g[idx].swappable = 1;   /* mark grid as not needed */
        printf_dbg2("clear_grid(): idx=%d curr_grids=%d alloc_grids=%d\n", idx, h->curr_grids, h->alloc_grids);
    } else {
        printf("clear_grid(): requested grid %d is off range or empty\n", idx);
    	return 1;
    }
    return 0;
} /* clear_grid */



int dirty_grid(grid_heap *h, int idx) {
    if ((idx < h->heap_size)  && (idx >= 0) ) {
        if (!h->g[idx].valid) {
            printf_dbg("dirty_grid(): ERROR - grid %d is not a valid grid - not created with new_grid\n",idx);
            return 1;
        }
        if (h->g[idx].grid == NULL) {
            printf_dbg("dirty_grid(): ERROR - grid %d space is NULL\n",idx);
            return 1;
        }		
        h->g[idx].swappable = 1;   /* mark grid as not needed */
        h->g[idx].dirty = 1;       /* mark grid to be saved before being reused */
        printf_dbg2("dirty_grid(): idx=%d curr_grids=%d alloc_grids=%d\n", idx, h->curr_grids, h->alloc_grids);
    } else {
        printf("dirty_grid(): requested grid %d is off range or empty\n", idx);
        return 1;
    }
    return 0;
} /* dirty_grid */



void delete_grid(grid_heap *h, int idx) {
    if (idx < h->heap_size && idx >= 0) {
        h->alloc_grids--;
        h->g[idx].next_grid = h->next_grid;
        h->next_grid = idx;
        h->g[idx].swappable = 1;
        h->g[idx].dirty = 0;
        h->g[idx].valid = 0;    /* mark grid has invalid */
        if (h->use_fs) {
            if (h->g[idx].fp)
                close_file(h->g[idx].fp);
            h->g[idx].fp = create_file(h->g[idx].filename);
        }
        printf_dbg2("delete_grid(): curr_grids=%d alloc_grids=%d\n", h->curr_grids, h->alloc_grids);
    } else {
        printf_dbg("delete_grid(): requested grid %d is off range\n", idx);
    }
} /* delete_grid */



void delete_heap(grid_heap *h) {
    int i;

    if (h) {
        if (h->g) {
            for (i = 0; i < h->heap_size; i++) {
                if (h->g[i].fp) close_file(h->g[i].fp);
                if (h->g[i].filename) free(h->g[i].filename);
                //if (h->g[i].grid) free(h->g[i].grid);
                if (h->g[i].grid) tsi_free(h->g[i].grid);
            }
            tsi_free(h->g);
        }
        tsi_free(h);
    }
} /* delete_heap */


/* end of grid_heap.c */
